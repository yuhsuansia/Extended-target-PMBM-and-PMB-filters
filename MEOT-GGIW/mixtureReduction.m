function ppp_hat = mixtureReduction(ppp,threshold)

if length(ppp) == 1
    ppp_hat = ppp;
    return;
end

%Index set of components
I = 1:length(ppp);
el = 1;

%weight of ppp components
w = [ppp.w]';

while ~isempty(I)
    Ij = [];
    %Find the component with the highest weight
    [~,j] = max(w);
    bern2 = ppp(j);
    bern2.r = 1;
    for i = I
        %Find other similar components in the sense of small Kullback
        %Leibler divergence
        if i == j
            Ij= [ Ij i ];
        else
            bern1 = ppp(i);
            bern1.r = 1;
            cross_entropy = (bern_kl(bern2,bern1) + ...
                bern_kl(bern1,bern2))/2;
            if cross_entropy < threshold
                Ij= [ Ij i ];
            end
        end
    end
    
    %Merge components by moment matching
    [temp,ppp_hat(el,1).w] = normalizeLogWeights(w(Ij));
    [ppp_hat(el,1).xr,ppp_hat(el,1).Cr] = kinematic_merge(ppp(Ij),temp);
    [ppp_hat(el,1).V,ppp_hat(el,1).v] = extent_merge(ppp(Ij),temp);
    [ppp_hat(el,1).alpha,ppp_hat(el,1).beta] = gamma_merge(ppp(Ij),temp);
    
    %Remove indices of merged components from index set
    I = setdiff(I,Ij);
    %Set a negative to make sure this component won't be
    %selected again
    w(Ij,1) = -inf;
    el = el+1;
end

end