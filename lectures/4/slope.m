function [s] = slope(x,fx,xi)
%%
% SLOPE.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.25.13
% Modified on 06.21.13
%
% PURPOSE   Computes the slope of a function by linear interpolation
% USAGE     s = slope(x,fx,xi)
% INPUTS    x   : nodes
%           fx  : function evaluated at nodes x
%           xi  : evaluation point
% OUTPUTS   fxi : interpolation approximation at point xi
%
n = size(fx,1);
b = zeros(n-1,1);

for i = 1:n-1;
    b(i) = (fx(i+1)-fx(i))/(x(i+1)-x(i));
end;

if xi >= x(n);
    s = b(n-1);
elseif xi < x(1);
    s = b(1);
else
    i = 1;
    while x(i) <= xi;
        i = i+1;
    end;
    
    if i == 2;
        s = b(i-1);
    else
        if xi == x(i-1);
            s = mean( [b(i-2);b(i-1)] );
        else
            s = b(i-1);
        end
    end
end;

end