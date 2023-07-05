tol = 1e-6;
kmax = 1e4;

A = [10 0 1; 0.5 7 1; 1 0 6];
b = [21; 9; 8];

x0 = zeros(size(b));

[x, iterations, residuals, errors] = gauss_seidel(A, b, x0, tol, kmax);

disp('Gauss-Seidel Method Solution:');
disp(x);

figure;
semilogy(iterations, residuals);
title('Gauss-Seidel Method');
xlabel('Iterations');
ylabel('Residuals');

ratios = (errors(2:end) ./ errors(1:end-1));

figure;
plot(1:length(errors) - 1, ratios);
title('Gauss-Seidel Method');
xlabel('Iterations');
ylabel('Ratio of Errors');

[x, iterations, residuals, errors] = jacobi(A, b, x0, tol, kmax);

disp('Jacobi Method Solution:');
disp(x);

figure;
semilogy(iterations, residuals);
title('Jacobi Method');
xlabel('Iterations');
ylabel('Residuals');

ratios = (errors(2:end) ./ errors(1:end-1));

figure;
plot(1:length(errors) - 1, ratios);
title('Jacobi Method');
xlabel('Iterations');
ylabel('Ratio of Errors');