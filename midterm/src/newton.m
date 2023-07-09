function [iterations,solutions,errors] = newton(x0,tol,crit)
%
% NEWTON.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
% Updated by Carlos Lezama

i = 0;

iterations(1) = i;
solutions(1,:) = x0;
errors(1) = norm(func0(x0),'inf');

while crit > tol
    i = i + 1;
    iterations(i + 1) = i;
    [b,A] = func0(x0);
    s = gaussj(A,-b);
    x1 = x0 + s';
    crit = norm(x1 - x0);
    solutions(i + 1,:) = x1;
    errors(i + 1) = norm(func0(x1),'inf') ./ errors(1);
    x0 = x1;
end

end