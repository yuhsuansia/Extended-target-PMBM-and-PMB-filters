function b = log_gamma2(a)
%GAMMA2 evaluates the logarithm of bivariate gamma function at a

b = 2*log(pi) + gammaln(a) + gammaln(a-0.5);

end

