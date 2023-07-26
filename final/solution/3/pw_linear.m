function [fxi] = pw_linear(x,fx,xi)
%%
% PW_LINEAR.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.20.13
% Modified on 07.18.13
% 
% PURPOSE   Performs piecewise linear interpolation
% USAGE     s = pw_linear(x,fx,xi) 
% INPUTS    x   : nodes
%           fx  : function evaluated at nodes x
%           xi  : evaluation point
% OUTPUTS   fxi : interpolation approximation at point xi
%    
n = size(fx,1);
a = zeros(n-1,1);
b = zeros(n-1,1);

for i = 1:n-1;
    b(i) = (fx(i+1)-fx(i))/(x(i+1)-x(i));
    a(i) = fx(i)-b(i)*x(i);
end;

if xi >= x(n);
    fxi = a(n-1)+b(n-1)*xi;
elseif xi < x(1);
    fxi = a(1)+b(1)*xi;
else
    i = 1;
    while x(i) <= xi;
        i = i+1;
    end;
    fxi = a(i-1)+b(i-1)*xi;
end;

end