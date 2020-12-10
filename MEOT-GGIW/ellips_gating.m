function meas_in_gate = ellips_gating(state,z,measmodel,gating_size)
%ELLIPS_GATING performs ellipsoidal gating

if isempty(z)
    meas_in_gate = false(0,1);
    return
end

nz = size(z,2);
meas_in_gate = false(nz,1);

Cy = measmodel.H*state.Cr*measmodel.H' + ...
    measmodel.Ch*state.V/(state.v-6) + measmodel.Cv;
Cy = (Cy+Cy')/2;

nu = z - repmat(measmodel.H*state.xr,[1 nz]);
dist = diag(nu.'*(Cy\nu));

meas_in_gate(dist<gating_size) = true;

end

