function [x, iters, residuals, errors] = gauss_seidel(A, b, x0, tol, kmax)
% Arguments:
% A: matrix
% b: vector
% x0: initial guess
% tol: tolerance
% kmax: maximum number of iterations

k = 0;
x = x0;
error = tol + 1;

errors = zeros(1, kmax);
iters = zeros(1, kmax);
residuals = zeros(1, kmax);

L = tril(A);
U = triu(A, 1);

while error > tol && k < kmax
    k = k + 1;
    x_new = L \ (b - U * x);
    error = norm(x_new - x);
    euclidean_norm = norm(b - A*x_new);
    x = x_new;
    errors(k) = error;
    iters(k) = k;
    residuals(k) = euclidean_norm;
end

top = max(iters);

errors = errors(1:top);
iters = iters(1:top);
residuals = residuals(1:top);

end