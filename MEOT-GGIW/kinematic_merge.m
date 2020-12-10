function [xr,Cr] = kinematic_merge(states,w)

w = exp(w);
idx = w > 0;
w = w(idx);
states = states(idx);
xr = [states(:).xr]*w;
Cr = zeros(size(states(1).Cr));
for i = 1:length(w)
    %Add spread of means
    x_diff = states(i).xr - xr;
    Cr = Cr + w(i)*(states(i).Cr + x_diff*x_diff');
end

end

