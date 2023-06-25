function [iters, sols, errors] = newton(fun, x0, tol, kmax)

k = 0;
x = x0;
error = tol + 1;

iters = zeros(1, kmax);
sols = zeros(1, kmax);
errors = zeros(1, kmax);

while error > tol && k < kmax
    k = k + 1;
    df = (fun(x + tol) - fun(x))/tol;
    x_new = x - fun(x)/df;
    error = abs(x_new - x);
    x = x_new;
    iters(k) = k;
    sols(k) = x;
    errors(k) = error;
end

top = max(iters);

iters = iters(1:top);
sols = sols(1:top);
errors = errors(1:top);

end
