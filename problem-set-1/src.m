% Problem Set 1
% Applied Macroeconometrics
% ITAM
% Summer 2023
% Written by Carlos Leddfama

setup;

% nonlinear_ex;

% Note that the function x^2 + 1 has no real roots,
% so none of the algorithms can find a solution.
% It is also clear that the Newton-Raphson method
% is faster than the Chord method.

linear_ex;

% The Gauss-Seidel method is faster than the Jacobi method.
% This is because the Gauss-Seidel method uses the updated values of the variables in the same iteration.
% It happens taking advantage of the lower triangular matrix.
% The Jacobi method uses the old values of the variables in the same iteration.