function [c] = coefs(x,fx)
%%
% COEFS.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.24.13
% Modified on 07.18.13
%
% PURPOSE   Assembles coefficients of the quadratic polynomial given nodes x
%           and function values corresponding to those nodes
% USAGE     c = coefs(x,fx)
% INPUTS    fx  : function evaluated at nodes x
%           it  : number of iterations
% OUTPUTS   c   : coefficients
%
c = zeros(5,1);
c(3) = (fx(1)-fx(2))/(x(1)-x(2));
c(1) = fx(2)-c(3)*x(2);
c(4) = ((x(3)+x(2))/(x(3)-x(2)))*c(3) + ( 2*x(2)*(fx(2)-fx(3)) )/( (x(2)-x(3))^2 );
c(5) = (fx(2)-fx(3))/(x(2)^2-x(3)^2) - c(4)/(x(2)+x(3));
c(2) = fx(3)-c(4)*x(3)-c(5)*x(3)^2;

end