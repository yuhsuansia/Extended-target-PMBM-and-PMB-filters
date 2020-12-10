function cross_entropy = bern_cross_entropy(bern1,bern2)

%we take the logarithm of bern2.r and 1-bern2.r; thus special care should
%be given here. Also, for numerical accuracy, bern1.r and bern2.r can be
%very close to 1, sometimes even slightly larger

bern2.r(bern2.r>1-1e-6) = 1-1e-6;

if bern1.r == 0 && bern2.r == 0
    cross_entropy = 0;
elseif bern1.r == 0 && bern2.r > 0
    cross_entropy = -log(1-bern2.r);
elseif bern1.r > 0 && bern2.r == 0
    cross_entropy = 1/eps;
else
    ce_kinematic = -2*log(2*pi) - log(det(bern2.Cr))/2 - ...
        trace((bern1.Cr + (bern1.xr-bern2.xr)*(bern1.xr-bern2.xr)')/bern2.Cr)/2;
    
    ce_gamma = bern2.alpha*log(bern2.beta) - gammaln(bern2.alpha) + ...
        (bern2.alpha-1)*(psi(bern1.alpha)-log(bern1.beta)) - ...
        bern2.beta*bern1.alpha/bern1.beta;
    
    ce_extent = -3/2*log(det(bern1.V/2)) + bern2.v/2*(psi((bern1.v-4)/2) + ...
        psi((bern1.v-3)/2)) - log_gamma2((bern2.v-3)/2) - ...
        (bern1.v-3)/2*trace(bern1.V\bern2.V) + (bern2.v-3)/2*log(det(bern1.V\bern2.V));
    
    if bern1.r == 1 && bern2.r == 1
        cross_entropy = -(ce_kinematic + ce_extent + ce_gamma);
    else
        cross_entropy = -(1-bern1.r)*log(1-bern2.r) - bern1.r*log(bern2.r) - ...
            bern1.r*(ce_kinematic + ce_extent + ce_gamma);
    end
    
end


end

