function [lc_upd,bern_updated] = ppp_upd(ppp,z,measmodel)
%PPP_UPD performs the measurement update of a PPP and computes the
%measurement likelihood

Nu = length(ppp);
l_c = zeros(Nu,1);
for i = 1:Nu
    [l_c(i),states(i)] = ggiw(ppp(i),z,measmodel);
end

[w,temp] = normalizeLogWeights([ppp.w]'+l_c);
[temp,lc_upd] = normalizeLogWeights([temp size(z,2)*log(measmodel.c_intensity)]);

%soft assignment    
[bern_updated.xr,bern_updated.Cr] = kinematic_merge(states,w);
[bern_updated.V,bern_updated.v] = extent_merge(states,w);
[bern_updated.alpha,bern_updated.beta] = gamma_merge(states,w);

%hard assignment
% [~,idx] = max(w);
% bern_updated = states(idx);
bern_updated.r = exp(temp(1));

end

