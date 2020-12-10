function kl_div = bern_kl(bern1,bern2)

bern1.r(bern1.r>1-1e-10) = 1;
bern2.r(bern2.r>1-1e-10) = 1;

if bern1.r == 0 && bern2.r == 0
    kl_div = 0;
elseif bern1.r == 0 && bern2.r > 0
    kl_div = log(1/(1-bern2.r));
elseif bern1.r > 0 && bern2.r == 0
    kl_div = 1/eps;
else
    gauss_kl = 1/2*(log(det(bern2.Cr)) - log(det(bern1.Cr)) - 4 ...
        + trace(bern2.Cr\bern1.Cr) + (bern1.xr-bern2.xr)'*inv(bern2.Cr)*(bern1.xr-bern2.xr));
    
    iw_kl = log_gamma2((bern2.v-3)/2) - log_gamma2((bern1.v-3)/2) + ...
        (bern1.v-3)/2*trace(bern1.V\bern2.V) - (bern1.v-3) - ...
        (bern2.v-3)/2*log(det(bern1.V\bern2.V)) - (bern2.v-bern1.v)/2*(psi((bern1.v-4)/2) + psi((bern1.v-3)/2));
    
    gamma_kl = bern1.alpha*log(bern1.beta) - bern2.alpha*log(bern2.beta) + ...
        gammaln(bern2.alpha) - gammaln(bern1.alpha) + bern1.alpha*(bern2.beta/bern1.beta - 1) ...
        + (bern1.alpha-bern2.alpha)*(psi(bern1.alpha)-log(bern1.beta));
    
    if bern1.r == 1 && bern2.r == 1
        kl_div = gauss_kl + iw_kl + gamma_kl;
    else
        r_kl = (1-bern1.r)*log((1-bern1.r)/(1-bern2.r)) ...
            + bern1.r*log(bern1.r/bern2.r);
        kl_div = r_kl + bern1.r*(gauss_kl + iw_kl + gamma_kl);
    end
    
end


end

