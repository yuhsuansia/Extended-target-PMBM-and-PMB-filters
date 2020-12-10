function partitions_b = gen_partitions(z,model,n_track)
% generate partitions of given measurement set

m = size(z,2);

if isempty(z)
    partitions_b = cell(0,1);
    return
elseif size(z,2) <= 1
    partitions_b{1} = true(m,1);
    return
end

n_track = min(n_track,m);

np = ceil((model.max_dist-model.min_dist)/model.grid_dist);
dist = linspace(model.min_dist,model.max_dist,np);
partitions = zeros(m,np+n_track);

for i = 1:np
    partitions(:,i) = dbscan(z',dist(i),1);
end
for i = np+1:np+n_track
    partitions(:,i) = kmeans(z',i-np);
end
%select unique partitions
partitions = unique(partitions','rows')';

% convert measurement indices to boolean representation
np = size(partitions,2);
partitions_b = cell(np,1);
for i = 1:np
    np_i = length(unique(partitions(:,i)));
    partitions_b{i} = false(m,np_i);
    for j = 1:np_i
        partitions_b{i}(partitions(:,i)==j,j) = true;
    end
end

end