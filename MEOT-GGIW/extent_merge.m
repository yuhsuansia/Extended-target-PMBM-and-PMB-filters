function [V2,nu2] = extent_merge(states,w)

d = 2;
w = exp(w);

idx = w>0;
w = w(idx);
states = states(idx);

wb = sum(w);
nu = [states.v];
N = length(states);
v_k = mean(nu);

C1 = zeros(d,d);
C2 = 0;
C3 = 0;
for i = 1:N
    % inverse Wisharts
    C1 = C1+w(i)*(nu(i)-d-1)*(states(i).V\eye(d));
    C2 = C2+w(i)*sum(psi(0,(nu(i)-d-(1:d))/2));
    C3 = C3 + w(i)*log(det(states(i).V));
end
C = d*wb*log(wb) -wb*log(det(C1)) +C2 -C3;

iter = 1;
while iter<100
    iter=iter+1;

    h_k  = d*wb*log(v_k-d-1)...
        -wb*sum(psi(0,(v_k-d-(1:d))/2))...
        +C;
    hp_k = d*wb/(v_k-d-1)...
        -0.5*wb*sum(psi(1,(v_k-d-(1:d))/2));
    hb_k = -d*wb/((v_k-d-1)^2)...
        -0.25*wb*sum(psi(2,(v_k-d-(1:d))/2));
   
%     v_new = v_k-h_k/hp_k; % Newtons
    v_new = v_k-(2*h_k*hp_k)/(2*hp_k^2-h_k*hb_k); % Halley's
    
    if abs(v_new-v_k)<1e-2
        v_k = v_new;
        break
    else
        v_k = v_new;
    end
    
    v_k=max(v_k,7);
end
nu2=v_k;
V2 = (nu2-d-1)*wb*(C1\eye(d));

end

