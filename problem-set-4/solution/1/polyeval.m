function [pxi] = polyeval(coefs,x,xi)
%%
% POLYEVAL.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.22.13
% Modified on 07.18.13
%
% PURPOSE   Quadratic spline interpolation using [P(coefs,x),x] at xi
% USAGE     s = polyeval(coefs,x,xi)
% INPUTS    coefs   : coefficients of the quadratic polynomial
%           x       : nodes
%           xi      : evaluation point
% OUTPUTS   pxi : interpolation approximation at point xi
%
if xi >= x(1) && xi <= x(2);
    pxi = coefs(1)+coefs(2)*xi;
else
    pxi = coefs(3) + coefs(4)*xi + coefs(5)*xi^2;
end

end



