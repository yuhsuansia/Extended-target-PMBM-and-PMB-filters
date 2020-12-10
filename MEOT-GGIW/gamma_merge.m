function [alpha,beta] = gamma_merge(states,w)

w = exp(w);

idx = w > 0;
w = w(idx);
states = states(idx);

wb = sum(w);
a = [states.alpha]';
b = [states.beta]';
c1 = 0;
c2 = 0;
for i = 1:length(w)
    c1 = c1 + w(i)*(psi(0,a(i)) - log(b(i)));
    c2 = c2 + w(i)*a(i)/b(i);
end
c = c1/wb - log(c2/wb);

a_k = sum(w(:).*a(:))/wb;
iter = 1;
while iter < 100
    iter = iter+1;
    
    h_k = log(a_k) - psi(0,a_k)+c;
    hp_k = 1/a_k - psi(1,a_k);
    hb_k = -1/a_k^2 - psi(2,a_k);
    
    a_new = a_k - (2*h_k*hp_k)/(2*hp_k^2 - h_k*hb_k); % Halley's
    
    if abs(a_new - a_k) < 1e-2
        a_k = a_new;
        break
    else
        a_k = a_new;
    end
    
    a_k = max(a_k,1);
end
alpha = a_k;
beta = alpha/(sum(w(:).*a(:)./b(:))/wb);

end

