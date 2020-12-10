close all
clear
% clc
dbstop if error

rng('default')

plot_enable = true;

%% Parameters of multi-object dynamic model

%choose a scenario
%scenario 1: targets are born from four seperated locations and move
%randomly
%scenario 2: ten targets are born from seperated locations; they first move
%close to each other and then seperate
scenario = 1;

%survival probability
Ps = 0.99;

%Target kinematic state [x-position,y-position,x-velocity,y-velocity]
%Target extent state [orientation,semi-axis length 1,semi-axis length 2]

%Parameters of a nearly constant velocity motion model
%kinematic state dimension
dxr = 4;
%time interval
Ts = 1;
%transition matrix for kinematic state
Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
%process noise
q = 0.09;
%process noise covariance matrix for kinematic state
Cwr = q*[Ts^3/3 0      Ts^2/2 0;
    0      Ts^3/3 0      Ts^2/2;
    Ts^2/2 0      Ts     0;
    0      Ts^2/2 0      Ts];
%measurement rate parameter used for prediction of gamma distribution
eta = 1.2;
%forgetting factor used for prediction of inverse-Wishart distribution
tau = 20;

%struct representation
motionmodel.Ps = Ps;
motionmodel.Ts = Ts;
motionmodel.dxr = dxr;
motionmodel.Ar = Ar;
motionmodel.Cwr = Cwr;
motionmodel.eta = eta;
motionmodel.tau = tau;

if scenario == 1
    %Poisson birth model with Gaussian mixture intensity
    %number of Gaussian components
    b_nG = 4;
    %Poisson intensity
    birthmodel = repmat(struct('w',log(0.25),'xr',[],'Cr',diag([10 10 2 2])^2,...
        'V',[],'v',100,'alpha',1000,'beta',100),b_nG,1);
    %specify kinematic state (x-position,y-position,x-velocity,y_velocity)
    birthmodel(1).xr = [50 50 0 0]';
    birthmodel(2).xr = [-50 50 0 0]';
    birthmodel(3).xr = [-50 -50 0 0]';
    birthmodel(4).xr = [50 -50 0 0]';
    %specify extent state (orientation,two axis lengths)
    birthmodel(1).V = (100+6)*diag([1 1])^2;
    birthmodel(2).V = (100+6)*diag([1 1])^2;
    birthmodel(3).V = (100+6)*diag([1 1])^2;
    birthmodel(4).V = (100+6)*diag([1 1])^2;
    %Poisson birth rate
    b_lambda = sum(exp([birthmodel.w]));
elseif scenario == 2
    birthmodel = struct('w',log(0.25),'xr',[-100,-100,0,0]','Cr',diag([100 100 4 4])^2,...
        'V',(100+6)*diag([1 1]),'v',100,'alpha',1000,'beta',100);
end

%% Parameters of multi-object measurement model

if scenario == 1
    %detection probability
    Pd = 0.9;
    %Poisson clutter rate
    c_lambda = 60;
elseif scenario == 2
    Pd = 0.3;
    c_lambda = 100;
end

%surveillance area
range_x = [-200 200];
range_y = [-200 200];
%uniform distributed clutter density
c_pdf = 1/(range_x(2)-range_x(1))/(range_y(2)-range_y(1));
%Poisson clutter intensity
c_intensity = c_lambda*c_pdf;

%Parameters of a linear Gaussian measurement model
%measurement dimension
dz = 2;
%observation matrix
H = [1 0 0 0;
    0 1 0 0];
%covariance of the multiplicative noise
Ch = diag([1/4 1/4]);
%covariance of the measurement noise
Cv = diag([1/8 1/8]);

%struct representation
measmodel.dz = dz;
measmodel.H = H;
measmodel.Ch = Ch;
measmodel.Cv= Cv;
measmodel.Pd = Pd;
measmodel.c_intensity = c_intensity;

%% Generating ground truth

if scenario == 1
    %total time steps
    T = 100;
    
    %create memory to store ground truth
    %birth time, death time, target state
    gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);
    
    %number of targets
    nt = 0;
    for i = 1:T-1
        %sample the number of newborn targets
        nb = poissrnd(b_lambda);
        %for the first time step, make sure that at least one target is born
        if i == 1 && nb == 0
            nb = 1;
        end
        for j = 1:nb
            %number of targets increases by 1
            nt = nt + 1;
            %sample a Gaussian component in the Poisson birth intensity
            b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/b_lambda),1) - 1;
            %sample an initial target state
            gt(nt).x_bt = i;
            gt(nt).x_dt = i;
            %assume fixed Poisson rate
            gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
            gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
            gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
        end
        %generate target trajectory for all the newborn targets
        for j = 1:nt
            %termintate the trajectory if the target dies or moves out of the
            %surveillance area
            %also assumes that no target dies if there is only one
            if (rand < Ps || nt==1) && gt(j).xr(1,end) >= range_x(1) ...
                    && gt(j).xr(1,end) <= range_x(2) ...
                    && gt(j).xr(2,end) >= range_y(1) ...
                    && gt(j).xr(2,end) <= range_y(2) && gt(j).x_dt == i
                gt(j).x_dt = i+1;
                gt(j).x_lambda = [gt(j).x_lambda gt(j).x_lambda(end)];
                %add motion noise when generating trajectory
                gt(j).xr = [gt(j).xr mvnrnd((Ar*gt(j).xr(:,end))',Cwr)'];
                gt(j).X = cat(3,gt(j).X,gt(j).X(:,:,end));
            end
        end
    end
elseif scenario == 2
    %trajectories are generated by propagating targets forward and backward
    %from six closely spaced targets states at the mid of the scene at half
    %time
    
    %total time steps
    T = 101;
    
    %number of targets
    nt = 26;
    %specify birth time and death time
    birth_time = 1:2:51;%[1,6,11,16,21,26,31,36,41,46];
    death_time = 51:2:101;%[56,61,66,71,76,81,86,91,96,101];
    %construct trajectory
    for i = 1:nt
        gt(i).x_bt = birth_time(i);
        gt(i).x_dt = death_time(i);
        tra_len = death_time(i)-birth_time(i)+1;
        gt(i).x_lambda = gamrnd(birthmodel.alpha,1/birthmodel.beta)*ones(1,tra_len);
        gt(i).X = repmat(iwishrnd(birthmodel.V,birthmodel.v-3),[1,1,tra_len]);
        %specify target state at time 50
        gt(i).xr = zeros(dxr,tra_len);
        mid_time = 52-birth_time(i);
        gt(i).xr(:,mid_time) = diag([6,6,0,0])*randn(dxr,1) + [0;0;3;3];
        %propagate backward
        j_idx = flip(1:mid_time-1);
        for j = 1:mid_time-1
            gt(i).xr(:,j_idx(j)) = mvnrnd((Ar\gt(i).xr(:,j_idx(j)+1))',Cwr)';
        end
        %propagate forward
        for j = mid_time+1:tra_len
            gt(i).xr(:,j) = mvnrnd((Ar*gt(i).xr(:,j-1))',Cwr)';
        end
    end
    
end

%cardinality of multi-target states
card = zeros(T,1);
for i = 1:T
    for j = 1:nt
        if gt(j).x_bt <= i && gt(j).x_dt >= i
            card(i) = card(i) + 1;
        end
    end
end

%% Plot ground truth

if plot_enable
    figure
    grid on
    hold on
    for i = 1:nt
        %plot trajectory
        gt_ks_plot = plot(gt(i).xr(1,:),gt(i).xr(2,:),'b','linewidth',1);
        %plot extent for every 5 time steps
        for j = 1:4:(gt(i).x_dt - gt(i).x_bt + 1)
            gt_es_plot = plot_extent_iw(gt(i).xr(1:2,j),gt(i).X(:,:,j),'-','r',1);
        end
    end
    xlim(range_x)
    ylim(range_y)
    xlabel('x (m)')
    ylabel('y (m)')
    legend([gt_ks_plot, gt_es_plot], {'Target trajectory','Target extent'});

end

%% Generate measurements

Z = cell(T,1);
%target-generated measurements
for i = 1:nt
    for t = gt(i).x_bt:gt(i).x_dt
        %generate measurements only if the target is detected
        %no misdetection at time 1
        if rand < Pd || t == 1
            nz = poissrnd(gt(i).x_lambda(t-gt(i).x_bt+1));
            Z{t} = [Z{t} mvnrnd(gt(i).xr(1:2,t-gt(i).x_bt+1)',...
                Ch*gt(i).X(:,:,t-gt(i).x_bt+1)+Cv,nz)'];
        end
    end
end

%append clutter
for i = 1:T
    nz = poissrnd(c_lambda);
    zx = rand(nz,1)*(range_x(2)-range_x(1)) + range_x(1);
    zy = rand(nz,1)*(range_y(2)-range_y(1)) + range_y(1);
    Z{i} = [Z{i} [zx zy]'];
end

%% Run filter

%paramter setting

%gate size in probability
paras.gating.Pg = 0.999;
paras.gating.size = chi2inv(paras.gating.Pg,dz);
%hyperparameters in DBSCAN
paras.dbscan.max_dist = 5;
paras.dbscan.min_dist = 0.1;
paras.dbscan.grid_dist = 0.1;
%data association algorithm
%method 1: Murty's algorithm
%method 2: Gibbs sampling
paras.data_assoc = 1;
%number of iterations in Gibbs sampling
paras.gibbs = 100;

%whether to use multi-Bernoulli birth model (the default is ppp birth)
paras.mb_birth = false;

%pruning threshold for global hypotheses
paras.pruning.w = log(1e-2);
%pruning threshold for ppp intensity
paras.pruning.ppp = log(1e-3);
%merging threshold for ppp intensity
paras.merging.ppp = 2;
%cap of global hypotheses
paras.cap.w = 100;
%pruning threshold for Bernoulli
paras.pruning.r = 1e-3;
%whether to perform recycling
paras.recycle = true;
if paras.mb_birth
    paras.recycle = false;
end
if paras.recycle
    paras.pruning.r = 1e-1;
end

%whether to perform MB approximation
paras.mb_approx = true;
%M-best assignments in Murty
paras.M = 20;
%MB approximation methods
%method 1: track-oriented merging
%method 2: SJPDA type merging
%method 3: LP merging
paras.mb_merge = 1;

%two different ways to initiate new tracks
%1: initiate a new track for each measurement
%2: initiate a new track for each cluster (this one is only valid for
%track-oriented merging)
paras.new_track_init = 1;

%convergence threshold for variational approximation
paras.vb_threshold = 1e-2;

%estimator to extract multi-target state
%estimator 1: extract state from global hypothesis with the highest weight
%and Bernoulli components with large enough probability of existence
%estimator 2: MAP cardinality estimator
%threshold to extract state estimate from Bernoulli
paras.estimator = 2;
paras.estimate.r = 0.5;

%parameters of GOSPA metric
gospa.p = 1;
gospa.c = 20;
gospa.alpha = 2;

%initialisation parameters
if ~paras.mb_birth
    %global hypothesis weight in logarithm
    mbm.w = 0;
    %global hypothesis look-up table
    mbm.table = zeros(1,0);
    %local hypothesis trees (collections of single target hypotheses)
    mbm.track = cell(0,1);
    %each Bernoulli is parameterised by 1) existence probability r, 2) mean and
    %covariance of the kinematic state xr, Cr, 3) parameters of the
    %extent state V, v, 4) parameters of gamma distribution alpha, beta.
    %PPP for undetected targets, initialised using birth model
    ppp = birthmodel;
else
    %initial setting for multi-Bernoulli birth model
    mbm.w = 0;
    nb = length(birthmodel);
    mbm.table = ones(1,nb);
    mbm.track = cell(nb,1);
    for i = 1:nb
        mbm.track{i}.r = exp(birthmodel(i).w);
        mbm.track{i}.xr = birthmodel(i).xr;
        mbm.track{i}.Cr = birthmodel(i).Cr;
        mbm.track{i}.V = birthmodel(i).V;
        mbm.track{i}.v = birthmodel(i).v;
        mbm.track{i}.alpha = birthmodel(i).alpha;
        mbm.track{i}.beta = birthmodel(i).beta;
    end
end

%memory to store state estimate
est = cell(T,1);
card_est = zeros(T,1);
t_elapsed = zeros(T,1);
%recursive Bayesian estimation
fprintf('Time step: ')
for t = 1:T
    
    fprintf('%d ',t)
    tic
    
    %ellipsoidal gating for detected targets
    %use a boolean vector to store the gating result of each local
    %hypothesis
    gating_matrix_d = cellfun(@(x) cell2mat(arrayfun(@(x) ...
        ellips_gating(x,Z{t},measmodel,paras.gating.size),...
        x,'uniformoutput',false)),mbm.track,'uniformoutput',false);
    
    if ~paras.mb_birth
        %ellipsoidal gating for undetected targets
        %use a boolean vector to store the gating result of each component
        gating_matrix_u = cell2mat(arrayfun(@(x) ...
            ellips_gating(x,Z{t},measmodel,paras.gating.size),...
            ppp,'uniformoutput',false)');
    end
    
    %remove unused measurements according to the gating result
    %gating result of all targets
    %number of tracks
    n_track = length(mbm.track);
    if paras.mb_birth
        gating = logical(sum(cell2mat(gating_matrix_d'),2));
    else
        if n_track > 0
            gating = logical(sum(gating_matrix_u,2) + ...
                sum(cell2mat(gating_matrix_d'),2));
        else
            gating = logical(sum(gating_matrix_u,2));
        end
    end
    
    %used measurements
    W = Z{t}(:,gating);
    %number of measurements after gating
    nm = size(W,2);
    
    %reconstruct gating matrix
    gating_matrix_d = ...
        cellfun(@(x) x(gating,:),gating_matrix_d,'uniformoutput',false);
    if ~paras.mb_birth
        gating_matrix_u = gating_matrix_u(gating,:);
    end
    %use DBSCAN to obtain multiple partitions
    partitions = gen_partitions(W,paras.dbscan,n_track);
    %number of partitions
    np = length(partitions);
    
    %find all the unique clusters in all measurement partitions
    [clusters,IA,IC] = unique(cell2mat(partitions')','rows');
    clusters = clusters';
    n_clusters = length(IA);
    %reconstruct partitions to let it contain indices of clusters
    nc_p = cellfun(@(x) size(x,2),partitions);
    partitions_indices = cell(np,1);
    idx = 0;
    for i = 1:np
        partitions_indices{i} = IC(idx+1:idx+nc_p(i));
        idx = idx + nc_p(i);
    end
    
    bern_new = repmat(struct('r',0,'xr',zeros(dxr,1),'Cr',ones(dxr,dxr),...
        'V',zeros(2,2),'v',0,'alpha',1,'beta',1),1,n_clusters);
    lik_new = zeros(n_clusters,1);
    if paras.mb_birth
        for c = 1:n_clusters
            lik_new(c) = sum(clusters(:,c))*log(measmodel.c_intensity);
        end
    else
        for c = 1:n_clusters
            %check if the cth cluster is in the gate of any ppp components
            ppp_idx = sum(gating_matrix_u-clusters(:,c)<0) == 0;
            if any(ppp_idx)
                [lik_new(c),bern_new(c)] = ...
                    ppp_upd(ppp(ppp_idx),W(:,clusters(:,c)),measmodel);
            else
                %if not, they are all clutter
                lik_new(c) = sum(clusters(:,c))*log(measmodel.c_intensity);
            end
        end
    end
    
    %create updated single target hypotheses for detected targets
    %number of single target hypotheses per track
    n_local_hypo = cellfun('length',mbm.track);
    tracks_upd = cell(n_track,1);
    for i = 1:n_track
        tracks_upd{i} = cell(n_local_hypo(i),1);
        for j = 1:n_local_hypo(i)
            %misdetection for the jth single target hypothesis under the
            %ith local hypothesis tree
            [l_missed,bern_missed] = bern_miss(mbm.track{i}(j),measmodel);
            tracks_upd{i}{j}.c = 0;
            tracks_upd{i}{j}.lik = l_missed;
            tracks_upd{i}{j}.bern = bern_missed;
            %measurement update for the jth single target hypothesis under
            %the ith local hypothesis tree
            for c = 1:n_clusters
                %check if the cth cluster is in the gate of the
                %corresponding single target hypothesis
                if sum(gating_matrix_d{i}(:,j)-clusters(:,c)<0) == 0
                    [l_upd,bern_updated] = ...
                        bern_upd(mbm.track{i}(j),W(:,clusters(:,c)),measmodel);
                    tracks_upd{i}{j}.c(end+1) = c;
                    tracks_upd{i}{j}.lik(end+1) = l_upd;
                    tracks_upd{i}{j}.bern(end+1) = bern_updated;
                end
            end
        end
    end
    
    if (paras.new_track_init == 2) && paras.mb_approx || paras.mb_birth
        tracks_new = cell(n_clusters,1);
        for c = 1:n_clusters
            tracks_new{c}.c = [0 c];
            tracks_new{c}.lik = [0 lik_new(c)];
            tracks_new{c}.bern = struct('r',0,'xr',zeros(dxr,1),'Cr',...
                ones(dxr,dxr),'V',zeros(2,2),'v',0,'alpha',1,'beta',1);
            tracks_new{c}.bern = [tracks_new{c}.bern bern_new(c)];
        end
    else
        %create updated single target hypotheses for the first detection of
        %undetected targets
        tracks_new = cell(nm,1);
        for i = 1:nm
            %create single target hypothesis for non-existent target
            tracks_new{i}.c = 0;
            tracks_new{i}.lik = 0;
            tracks_new{i}.bern = struct('r',0,'xr',zeros(dxr,1),'Cr',...
                ones(dxr,dxr),'V',zeros(2,2),'v',0,'alpha',1,'beta',1);
        end
        for c = 1:n_clusters
            %create new single target hypothesis
            idx = find(clusters(:,c),1,'last');
            tracks_new{idx}.c = [tracks_new{idx}.c c];
            tracks_new{idx}.lik = [tracks_new{idx}.lik lik_new(c)];
            tracks_new{idx}.bern = [tracks_new{idx}.bern bern_new(c)];
        end
    end
    
    %reset m to the number of new tracks
    m = length(tracks_new);
    
    %update local hypothesis trees
    mbm_upd.track = cell(n_track+m,1);
    for i = 1:n_track
        idx = 0;
        for j = 1:length(tracks_upd{i})
            mbm_upd.track{i} = [mbm_upd.track{i} tracks_upd{i}{j}.bern];
            %use an extra variable to record the index of each new single
            %target hypothesis in local hypothesis tree i
            n_ij = length(tracks_upd{i}{j}.c);
            tracks_upd{i}{j}.idx = (1:n_ij) + idx;
            idx = idx + n_ij;
        end
    end
    for i = 1:m
        mbm_upd.track{n_track+i,1} = tracks_new{i}.bern;
    end
    
    %data association for each global association hypothesis
    mbm_upd.w = [];
    mbm_upd.table = zeros(0,n_track+m);
    A = length(mbm.w);
    for a = 1:A
        a_indices = mbm.table(a,:);
        costs_temp = 0;
        for j = 1:n_track
            if a_indices(j) > 0
                single_hypo = tracks_upd{j}{a_indices(j)};
                costs_temp = costs_temp + single_hypo.lik(1);
            end
        end
        %if there is no measurement partition, all the targets are
        %misdetected
        if n_clusters==0
            table_upd = zeros(1,n_track+m);
            w_upd = 0;
            for j = 1:n_track
                if a_indices(j) > 0
                    table_upd(1,j) = tracks_upd{j}{a_indices(j)}.idx(1);
                    w_upd = w_upd + tracks_upd{j}{a_indices(j)}.lik(1);
                end
            end
            for j = n_track+1:n_track+m
                table_upd(1,j) = 1;
            end
            mbm_upd.w = [mbm_upd.w;w_upd + mbm.w(a)];
            mbm_upd.table = [mbm_upd.table;table_upd];
        end
        %go through each measurement partition
        for p = 1:np
            %find all the clusters under this partition
            p_idx = partitions_indices{p};
            %number of clusters under this partition
            p_n = length(p_idx);
            %construct cost matrix
            C = inf(p_n,n_track+m);
            for j = 1:n_track
                if a_indices(j) > 0
                    single_hypo = tracks_upd{j}{a_indices(j)};
                    %set cost for detected targets
                    [LIA,LOCB] = ismember(single_hypo.c,p_idx);
                    C(LOCB(LOCB>0),j) = -single_hypo.lik(LIA)' + single_hypo.lik(1);
                end
            end
            %set cost for undetected targets
            for j = n_track+1:n_track+m
                single_hypo = tracks_new{j-n_track};
                [LIA,LOCB] = ismember(single_hypo.c,p_idx);
                C(LOCB(LOCB>0),j) = -single_hypo.lik(LIA)';
            end
            
            %find columns of C that contain finite entries
            idx = find(sum(isfinite(C),1) > 0);
            C = C(:,idx);
            if paras.data_assoc == 1
                %find M-best assignments using Murty's algorithm
                [assignments,~,costs] = kBest2DAssign(C,ceil(paras.M*exp(mbm.w(a))));
                assignments = assignments';
                costs = costs';
            elseif paras.data_assoc == 2
                %find M-best assignments using Gibbs sampling
                [assignments,costs] = assign2DByGibbs(C,...
                    paras.gibbs,paras.M);
            end
            %number of assignments
            n_a = size(assignments,1);
            %restore track indices
            for i = 1:n_a
                assignments(i,:) = idx(assignments(i,:));
            end
            
            table_upd = zeros(n_a,n_track+m);
            %update the glocal hypothesis look-up table and weight
            for i = 1:n_a
                %go through each association in a given assignment
                for j = 1:p_n
                    %check if detected targets or undetected targets
                    assoc = assignments(i,j);
                    if assoc <= n_track
                        %find the index of the corresponding cluster
                        table_upd(i,assoc) = tracks_upd{assoc}...
                            {a_indices(assoc)}.idx(tracks_upd{assoc}...
                            {a_indices(assoc)}.c == p_idx(j));
                    else
                        %find the index of the corresponding cluster
                        table_upd(i,assoc) = ...
                            find(tracks_new{assoc-n_track}.c == p_idx(j),1);
                    end
                end
                %go through each unassociated track
                unassign = true(n_track+m,1);
                unassign(assignments(i,:)) = false;
                unassign = find(unassign);
                for j = 1:length(unassign)
                    if unassign(j) <= n_track
                        temp = a_indices(unassign(j));
                        if temp > 0
                            table_upd(i,unassign(j)) = ...
                                tracks_upd{unassign(j)}{temp}.idx(1);
                        else
                            table_upd(i,unassign(j)) = 0;
                        end
                    else
                        table_upd(i,unassign(j)) = 1;
                    end
                end
            end
            
            %when computing the global hypothesis weight, assume that each
            %measurement partition is equal likely, i.e., uniform
            %distributed
            mbm_upd.w = [mbm_upd.w;-costs'+costs_temp+mbm.w(a)];
            mbm_upd.table = [mbm_upd.table;table_upd];
        end
    end
    
    %normalise global hypothesis weight
    mbm_upd.w = normalizeLogWeights(mbm_upd.w);
    %prune updated global hypotheses with small weights
    [mbm_upd_w,order] = sort(exp(mbm_upd.w),'descend');
    pos = find(cumsum(mbm_upd_w)>=1-exp(paras.pruning.w),1);
    mbm_upd.w = mbm_upd.w(order(1:pos));
    mbm_upd.table = mbm_upd.table(order(1:pos),:);
    mbm_upd.w = normalizeLogWeights(mbm_upd.w);
    if ~paras.mb_approx
        %cap the number of global hypotheses
        if length(mbm_upd.w) > paras.cap.w
            [~,idx] = sort(mbm_upd.w,'descend');
            mbm_upd.w = mbm_upd.w(idx(1:paras.cap.w));
            mbm_upd.w = normalizeLogWeights(mbm_upd.w);
            mbm_upd.table = mbm_upd.table(idx(1:paras.cap.w),:);
        end
    end
    
    %remove single target hypotheses with small probability of existence
    for i = 1:length(mbm_upd.track)
        if paras.mb_birth
            idx = find(([mbm_upd.track{i}.r] < paras.pruning.r) | ...
                ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ...
                ([mbm_upd.track{i}.v] < 7));
        else
            idx = find(([mbm_upd.track{i}.r] < exp(paras.pruning.ppp)) | ...
                ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ...
                ([mbm_upd.track{i}.v] < 7));
        end
        for j = 1:length(idx)
            mbm_upd.table(mbm_upd.table(:,i) == idx(j),i) = 0;
        end
        %re-index
        idx_0 = mbm_upd.table(:,i) > 0;
        [idx,~,temp] = unique(mbm_upd.table(idx_0,i));
        mbm_upd.table(idx_0,i) = temp;
        mbm_upd.track{i} = mbm_upd.track{i}(idx);
    end
    %remove empty track
    idx = ~cellfun('isempty',mbm_upd.track);
    mbm_upd.track = mbm_upd.track(idx);
    mbm_upd.table = mbm_upd.table(:,idx);
    if isempty(mbm_upd.table)
        mbm_upd.table = zeros(1,0);
        mbm_upd.track = cell(0,1);
        mbm_upd.w = 0;
    end
    
    %merge rows of global hypothesis look-up table that are the same
    if length(mbm_upd.w) > 1
        [mbm_upd.table,~,IC] = unique(mbm_upd.table,'rows');
        n_a = size(mbm_upd.table,1);
        temp = zeros(n_a,1);
        for i = 1:n_a
            [~,temp(i)] = normalizeLogWeights(mbm_upd.w(IC==i));
        end
        mbm_upd.w = temp;
    end
    
    %misdetection update of ppp
    if ~paras.mb_birth
        ppp = ppp_miss(ppp,measmodel);
    end
    
    %number of global hypotheses and tracks
    n_mb = size(mbm_upd.table,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if paras.mb_approx && n_mb > 1
        %find tracks with only one single target hypothesis being included
        %in any of the global hypotheses
        idx = sum(mbm_upd.table - 1,1) ~= 0;
        mbm_upd.track = [mbm_upd.track(idx);mbm_upd.track(~idx)];
        mbm_upd.table = [mbm_upd.table(:,idx) mbm_upd.table(:,~idx)];
        
        %track-oriented merging
        %number of tracks with conflicts
        n_track_t = length(idx);
        for i = 1:n_track_t
            %compute marginal probability that the single target hypothesis
            %in track i is included in the global hypothesis
            c = unique(mbm_upd.table(:,i));
            c = c(c>0);
            lc = length(c);
            w_margin = zeros(lc,1);
            for j = 1:lc
                %find all the global hypotheses that contain this single
                %target hypothesis and compute the log sum of the weights
                [~,log_sum_w] = normalizeLogWeights(mbm_upd.w(mbm_upd.table(:,i)==c(j)));
                %also take into account the probability of existence
                w_margin(j) = log_sum_w + log(mbm_upd.track{i}(c(j)).r);
            end
            %perform moment matching
            [w_margin_n,log_sum_w] = normalizeLogWeights(w_margin);
            mb_approx(i).r = exp(log_sum_w);
            [mb_approx(i).xr,mb_approx(i).Cr] = ...
                kinematic_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).V,mb_approx(i).v] = ...
                extent_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).alpha,mb_approx(i).beta] = ...
                gamma_merge(mbm_upd.track{i},w_margin_n);
        end
        
        if paras.mb_merge == 1 || n_track_t < 2
            %track-oriented merging
            for i = 1:n_track_t
                mbm_upd.track{i} = mb_approx(i);
            end
        elseif paras.mb_merge == 2 && n_track_t > 1
            %initialized using the result of track-oriented merging
            %only need to consider track with more than one single target
            %hypothesis
            bern0 = struct('r',0,'xr',zeros(dxr,1),...
                'Cr',ones(dxr,dxr),'V',zeros(2,2),'v',0,...
                'alpha',1,'beta',1);
            mbm_hypo = repmat(bern0,[n_mb,n_track_t]);
            for i = 1:n_mb
                for j = 1:n_track_t
                    if mbm_upd.table(i,j) > 0
                        mbm_hypo(i,j) = mbm_upd.track{j}(mbm_upd.table(i,j));
                    end
                end
            end
            nh = cellfun('length',mbm_upd.track(1:n_track_t));
            
            %iterative optimization
            C_optimal = inf;
            mb_hat = mb_approx;
            while(1)
                Cmin = zeros(n_mb,1);
                assignment_mb = zeros(n_track_t,n_mb);
                
                Ch = cell(n_track_t,n_track_t);
                for k = 1:n_track_t
                    for i = 1:n_track_t
                        Cht = zeros(nh(i),1);
                        for j = 1:nh(i)
                            Cht(j) = bern_cross_entropy(mbm_upd.track{i}(j),mb_hat(k));
                        end
                        Ch{k,i} = Cht;
                    end
                end
                
                for i = 1:n_mb
                    %construct cost matrix with columns corresponding
                    %to Bernoullis components in the ith MB and rows
                    %corresponding to the approximated Bernoullis
                    C = zeros(n_track_t,n_track_t);
                    for j = 1:n_track_t
                        for k = 1:n_track_t
                            if mbm_upd.table(i,j) == 0
                                C(k,j) = bern_cross_entropy(bern0,mb_hat(k));
                            else
                                C(k,j) = Ch{k,j}(mbm_upd.table(i,j));
                            end
                        end
                    end
                    %find the optimal assignment, assignment_mb(j,i)
                    %means that the jth approximated Bernoulli should
                    %be merged using the assignment_mb(j,i)th Bernoulli
                    %in the ith MB
                    [assignment_mb(:,i),~,Cmin(i)] = assign2D(C);
                end
                %check the columns of assignment_mb, if they are all of
                %the form 1:n_b, then the iterative optimization
                %converges
                if ~any(assignment_mb - (1:n_track_t)')
                    break;
                else
                    C_hat = Cmin'*exp(mbm_upd.w);
                    if C_optimal-C_hat < paras.vb_threshold
                        break;
                    else
                        C_optimal = C_hat;
                    end
                end
                %if not converged, perform merging to construct new
                %mb_hat
                for j = 1:n_track_t
                    w_margin = zeros(n_mb,1);
                    bern = repmat(bern0,[1,n_mb]);
                    for i = 1:n_mb
                        w_margin(i) = mbm_upd.w(i) + ...
                            log(mbm_hypo(i,assignment_mb(j,i)).r);
                        bern(i) = mbm_hypo(i,assignment_mb(j,i));
                    end
                    [w_margin_n,log_sum_w] = normalizeLogWeights(w_margin);
                    %note that each MB may contain Bernoullis with zero
                    %probability of existence with unvalid pdf. The
                    %merging of these Bernoullis will of course be a
                    %Bernoulli with zero probability of existence
                    if isnan(w_margin_n)
                        mb_hat(j) = bern0;
                    else
                        mb_hat(j).r = exp(log_sum_w);
                        [mb_hat(j).xr,mb_hat(j).Cr] = kinematic_merge(bern,w_margin_n);
                        [mb_hat(j).V,mb_hat(j).v] = extent_merge(bern,w_margin_n);
                        [mb_hat(j).alpha,mb_hat(j).beta] = gamma_merge(bern,w_margin_n);
                    end
                end
            end
            for i = 1:n_track_t
                mbm_upd.track{i} = mb_hat(i);
            end
            
        elseif paras.mb_merge == 3 && n_track_t > 1
            %initialized using the result of track-oriented merging
            %only need to consider track with more than one single target
            %hypothesis
            bern0 = struct('r',0,'xr',zeros(dxr,1),...
                'Cr',ones(dxr,dxr),'V',zeros(2,2),'v',0,...
                'alpha',1,'beta',1);
            nh = cellfun('length',mbm_upd.track(1:n_track_t))+1;
            bern_h = [];
            for i = 1:n_track_t
                bern_h = [bern_h bern0];
                bern_h = [bern_h mbm_upd.track{i}];
            end
            
            n_H = sum(nh);
            %construct constraints
            qhj = zeros(n_H,n_track_t);
            idxi = 0;
            for i = 1:n_track_t
                for j = 1:nh(i)
                    qhj(j+idxi,i) = sum(exp(mbm_upd.w((mbm_upd.table(:,i) == j-1))));
                end
                idxi = idxi + nh(i);
            end
            idx = sum(qhj,2)>0;
            bern_h = bern_h(idx);
            n_H = length(bern_h);
            qhj = qhj(idx,:);
            qj = sum(qhj,1)';
            qh = sum(qhj,2);
            
            %iterative optimization
            C_optimal = inf;
            mb_hat = mb_approx;
            while(1)
                %construct cost matrix
                C = zeros(n_H,n_track_t);
                for i = 1:n_H
                    for j = 1:n_track_t
                        C(i,j) = bern_cross_entropy(bern_h(i),mb_hat(j));
                    end
                end
                %solve the transportation problem
                [C_hat,q_hat] = LP_transport(C,qh,qj);
                %if q_hat and qhj are very similar, then the iterative
                %optimization converges
                if sum(sum(abs(q_hat-qhj))) < 1e-6
                    break;
                end
                if C_optimal-C_hat < paras.vb_threshold
                    break;
                else
                    C_optimal = C_hat;
                end
                %if not converged, perform merging to construct new
                %mb_hat
                for i = 1:n_track_t
                    mb_hat(i).r = [bern_h.r]*q_hat(:,i);
                    if mb_hat(i).r > 0
                        w_margin_n = log([bern_h.r]'.*q_hat(:,i)/mb_hat(i).r);
                        [mb_hat(i).xr,mb_hat(i).Cr] = kinematic_merge(bern_h,w_margin_n);
                        [mb_hat(i).V,mb_hat(i).v] = extent_merge(bern_h,w_margin_n);
                        [mb_hat(i).alpha,mb_hat(i).beta] = gamma_merge(bern_h,w_margin_n);
                    else
                        mb_hat(i) = bern0;
                    end
                end
            end
            for i = 1:n_track_t
                mbm_upd.track{i} = mb_hat(i);
            end
        end
        mbm_upd.w = 0;
        mbm_upd.table = ones(1,n_track_t);
        
        %recycle tracks with small probability of existence
        idx_logical = (cellfun(@(x) x.r, mbm_upd.track) >= paras.pruning.r) & ...
            (cellfun(@(x) x.alpha, mbm_upd.track)./...
            cellfun(@(x) x.beta, mbm_upd.track) > 1) & ...
            (cellfun(@(x) x.v, mbm_upd.track) > 7);
        if paras.recycle && (~paras.mb_birth)
            idx = find(~idx_logical);
            for j = 1:length(idx)
                ppp(end+1,1).w = log(mbm_upd.track{idx(j)}.r);
                ppp(end,1).xr = mbm_upd.track{idx(j)}.xr;
                ppp(end,1).Cr = mbm_upd.track{idx(j)}.Cr;
                ppp(end,1).V = mbm_upd.track{idx(j)}.V;
                ppp(end,1).v = mbm_upd.track{idx(j)}.v;
                ppp(end,1).alpha = mbm_upd.track{idx(j)}.alpha;
                ppp(end,1).beta = mbm_upd.track{idx(j)}.beta;
            end
        end
        %delete these single target hypotheses
        mbm_upd.track = mbm_upd.track(idx_logical);
        mbm_upd.table = mbm_upd.table(idx_logical);
        if isempty(mbm_upd.table)
            mbm_upd.table = zeros(1,0);
            mbm_upd.track = cell(0,1);
            mbm_upd.w = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %recycle single target hypotheses with small probability of existence
        for i = 1:length(mbm_upd.track)
            idx = find(([mbm_upd.track{i}.r] < paras.pruning.r) | ...
                ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ...
                ([mbm_upd.track{i}.v] < 7));
            %for each single target hypothesis to be recycled, find global
            %hypotheses that include it
            for j = 1:length(idx)
                temp = mbm_upd.table(:,i)==idx(j);
                mbm_upd.table(temp,i) = 0;
                if paras.recycle && (~paras.mb_birth)
                    [~,log_sum_w] = normalizeLogWeights(mbm_upd.w(temp));
                    ppp(end+1,1).w = log_sum_w+log(mbm_upd.track{i}(idx(j)).r);
                    ppp(end,1).xr = mbm_upd.track{i}(idx(j)).xr;
                    ppp(end,1).Cr = mbm_upd.track{i}(idx(j)).Cr;
                    ppp(end,1).V = mbm_upd.track{i}(idx(j)).V;
                    ppp(end,1).v = mbm_upd.track{i}(idx(j)).v;
                    ppp(end,1).alpha = mbm_upd.track{i}(idx(j)).alpha;
                    ppp(end,1).beta = mbm_upd.track{i}(idx(j)).beta;
                end
            end
            %delete these single target hypotheses
            %re-index
            idx_0 = mbm_upd.table(:,i) > 0;
            [idx,~,temp] = unique(mbm_upd.table(idx_0,i));
            mbm_upd.table(idx_0,i) = temp;
            mbm_upd.track{i} = mbm_upd.track{i}(idx);
        end
        %remove empty track
        idx = ~cellfun('isempty',mbm_upd.track);
        mbm_upd.track = mbm_upd.track(idx);
        mbm_upd.table = mbm_upd.table(:,idx);
        if isempty(mbm_upd.table)
            mbm_upd.table = zeros(1,0);
            mbm_upd.track = cell(0,1);
            mbm_upd.w = 0;
        end
        
        %merge rows of global hypothesis look-up table that are same
        if length(mbm_upd.w) > 1
            [mbm_upd.table,~,IC] = unique(mbm_upd.table,'rows');
            n_a = size(mbm_upd.table,1);
            temp = zeros(n_a,1);
            for i = 1:n_a
                [~,temp(i)] = normalizeLogWeights(mbm_upd.w(IC==i));
            end
            mbm_upd.w = temp;
        end
    end
    
    if ~paras.mb_birth
        %remove ppp components with small weights
        idx = ([ppp.w] > paras.pruning.ppp) & ([ppp.alpha]./[ppp.beta] > 1) & ([ppp.v] > 7);
        ppp = ppp(idx);
        if ~paras.mb_approx
            %merge similar ppp components
            ppp = mixtureReduction(ppp,paras.merging.ppp);
        end
    end
    
    %multi-target state estimation
    if paras.estimator == 1
        %find the global hypothesis with the highest weight
        [~,I] = max(mbm_upd.w);
        est{t}.xr = zeros(dxr,0);
        est{t}.X = zeros(2,2,0);
        if ~isempty(mbm_upd.table)
            hypo_best = mbm_upd.table(I,:);
            for i = 1:length(hypo_best)
                if mbm_upd.table(I,i) > 0
                    %extract state estimate from Bernoulli component with large
                    %enough probability of existence
                    if mbm_upd.track{i}(mbm_upd.table(I,i)).r > paras.estimate.r
                        card_est(t) = card_est(t) + 1;
                        est{t}.xr = [est{t}.xr mbm_upd.track{i}(mbm_upd.table(I,i)).xr];
                        est{t}.X = cat(3,est{t}.X,...
                            mbm_upd.track{i}(mbm_upd.table(I,i)).V...
                            /(mbm_upd.track{i}(mbm_upd.table(I,i)).v-6));
                    end
                end
            end
        end
        est{t}.card = size(est{t}.xr,2);
    elseif paras.estimator == 2
        %compute the cardinality distribution
        [n_a,n_track] = size(mbm_upd.table);
        card_dist = zeros(1,n_track+1);
        pcard = zeros(n_a,n_track+1);
        for i = 1:n_a
            r = [];
            for j = 1:n_track
                if mbm_upd.table(i,j) > 0
                    r = [r mbm_upd.track{j}(mbm_upd.table(i,j)).r];
                end
            end
            if ~isempty(r)
                %avoid numerical underflow
                r(r>1-1e-6) = 1-1e-6;
                pcard(i,1:length(r)+1) = prod(1-r)*poly(-r./(1-r));
            end
            card_dist = card_dist + pcard(i,:)*exp(mbm_upd.w(i));
        end
        %obtain the maximum cardinality
        [~,card_max] = max(card_dist);
        %find the global hypothesis with the highest weight and the same
        %MAP cardinality estimate
        [~,a_best] = max(pcard(:,card_max));
        r = zeros(n_track,1);
        xr = zeros(dxr,n_track);
        X = zeros(2,2,n_track);
        for i = 1:n_track
            if mbm_upd.table(a_best,i) > 0
                r(i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).r;
                xr(:,i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).xr;
                X(:,:,i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).V...
                    /(mbm_upd.track{i}(mbm_upd.table(a_best,i)).v-6);
            end
        end
        [~,I] = sort(r,'descend');
        card_est(t) = card_max-1;
        est{t}.xr = xr(:,I(1:card_max-1));
        est{t}.X = X(:,:,I(1:card_max-1));
        est{t}.card = sum(card_dist.*(0:length(card_dist)-1));
    end
    
    %prediction step
    %prediction for detected targets
    mbm.w = mbm_upd.w;
    mbm.table = mbm_upd.table;
    mbm.track = cellfun(@(x) arrayfun(@(x) ...
        bern_pred(x,motionmodel),x),mbm_upd.track,'uniformoutput',false);
    
    if paras.mb_birth
        %add multi-Bernoulli birth
        na = length(mbm.w);
        nt = length(mbm.track);
        mbm.table = [mbm.table ones(na,nb)];
        for i = 1:nb
            mbm.track{nt+i,1}.r = exp(birthmodel(i).w);
            mbm.track{nt+i,1}.xr = birthmodel(i).xr;
            mbm.track{nt+i,1}.Cr = birthmodel(i).Cr;
            mbm.track{nt+i,1}.V = birthmodel(i).V;
            mbm.track{nt+i,1}.v = birthmodel(i).v;
            mbm.track{nt+i,1}.alpha = birthmodel(i).alpha;
            mbm.track{nt+i,1}.beta = birthmodel(i).beta;
        end
    else
        %prediction for undetected targets
        ppp = ppp_pred(ppp,motionmodel,birthmodel);
    end
    
    t_elapsed(t) = toc;
end
fprintf('\n')

card_est = cellfun(@(x) x.card,est);
card_err = card_est - card;
%evaluate multi-target filtering performance using GOSPAs
d_gospa = zeros(T,1);
decomposed_cost = repmat(struct('localisation',[],'missed',[],'false',[]),T,1);
for t = 1:T
    x_mat.x = zeros(2,0);
    x_mat.X = zeros(2,2,0);
    for i = 1:length(gt)
        if gt(i).x_bt <= t && gt(i).x_dt >= t
            x_mat.x = [x_mat.x gt(i).xr(1:2,t-gt(i).x_bt+1)];
            x_mat.X = cat(3,x_mat.X,gt(i).X(:,:,t-gt(i).x_bt+1));
        end
    end
    y_mat.x = est{t}.xr;
    y_mat.X = est{t}.X;
    [d_gospa(t),~,decomposed_cost(t)] = ...
        GOSPA_extended(x_mat,y_mat,gospa.p,gospa.c,gospa.alpha);
end

if plot_enable
    
    figure
    plot(1:T,card,'linewidth',2)
    grid on
    hold on
    plot(1:T,card_est,'linewidth',2)
    xlabel('Time step')
    ylabel('Number of targets')
    legend('True cardinality','Estimated cardinality')
    
    figure
    subplot(2,2,1)
    plot(1:T,d_gospa,'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('GOSPA')
    subplot(2,2,2)
    plot(1:T,[decomposed_cost.localisation],'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('Localisation error')
    subplot(2,2,3)
    plot(1:T,[decomposed_cost.missed],'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('Misdetection error')
    subplot(2,2,4)
    plot(1:T,[decomposed_cost.false],'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('False detection error')
    
end

fprintf('Mean GOSPA: %.2f\n',mean(d_gospa))
fprintf('Total runtime: %.2f\n',sum(t_elapsed))

