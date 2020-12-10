function [l_upd,state_upd] = ggiw(state,z,measmodel)

m = size(z,2);

state_upd.alpha = state.alpha + m;
state_upd.beta = state.beta + 1;

z_hat = mean(z,2);
epsilon = z_hat - measmodel.H*state.xr;
X_hat = state.V/(state.v-6);
X_hat = (X_hat + X_hat')/2;
R_hat = measmodel.Ch*X_hat + measmodel.Cv;
S = measmodel.H*state.Cr*measmodel.H' + R_hat/m;
S = (S + S')/2;
K = state.Cr*measmodel.H'/S;

state_upd.xr = state.xr + K*epsilon;
state_upd.Cr = state.Cr - K*measmodel.H*state.Cr;

N = epsilon*epsilon';
X_hat_chol = chol(X_hat)';
S_chol = chol(S)';
R_chol = chol(R_hat)';
N_hat = X_hat_chol/S_chol*N/S_chol'*X_hat_chol';
Z = (z-z_hat)*(z-z_hat)';
Z_hat = X_hat_chol/R_chol*Z/R_chol'*X_hat_chol';

state_upd.v = state.v + m;
state_upd.V = state.V + N_hat + Z_hat;

l_upd = -m*log(pi) - log(m);
l_upd = l_upd + (state.v-3)/2*log(det(state.V)) - (state_upd.v-3)/2*log(det(state_upd.V));
l_upd = l_upd + log_gamma2((state_upd.v-3)/2) - log_gamma2((state.v-3)/2);
l_upd = l_upd + m/2*log(det(X_hat)) - (m-1)/2*log(det(R_hat)) - 1/2*log(det(S));
l_upd = l_upd + state.alpha*log(state.beta) - gammaln(state.alpha) ...
    - state_upd.alpha*log(state_upd.beta) + gammaln(state_upd.alpha);
l_upd = l_upd + log(measmodel.Pd);

end

