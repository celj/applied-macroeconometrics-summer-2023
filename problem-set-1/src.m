% Problem Set 1
% Applied Macroeconometrics
% ITAM
% Summer 2023
% Written by Carlos Leddfama

setup;

nonlinear_ex;

% Note that the function x^2 + 1 has no real roots,
% so none of the algorithms can find a solution.
% It is also clear that the Newton-Raphson method
% is faster than the Chord method.

linear_ex;

% The Gauss-Seidel method is faster than the Jacobi method.
% This is because the Gauss-Seidel method uses the updated values of the variables in the same iteration.
% This occurs by taking advantage of the lower triangular matrix.
% The Jacobi method uses the old values of the variables in the same iteration.

labor_ex;

% The variables calculated stand for the following:
% E: employment
% O: out of the labor force
% U: unemployment
% These variables represent the percentage of the population in each state such that E + U + O = 100%.
% Actually, this is set as a constraint in our linear system of equations.
% Additionally, the gross flows measure the movement of people between the states.
% Thus, f^UE stands for the flow of people from employment to unemployment;
% therefore, we consider it to calculate the unemployment rate.