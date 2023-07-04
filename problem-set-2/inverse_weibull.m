function[W] = inverse_weibull(U, lambda, k)

W = lambda * (-log(1-U)).^(1/k);

end