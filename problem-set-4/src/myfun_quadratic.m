function [v] = myfun_quadratic(theta,k0,t)
%%
% MYFUN_LINEAR.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.03.13
% Modified on 07.16.13
%
% PURPOSE   Computes the value function following a linear policy rule
% USAGE     v = myfun_linear(theta,k0,t)
% INPUTS    theta : parameters of the policy rule
%           k0    : initial capital value
%           t     : horizon
% OUTPUTS   v     : value function
%
global beta sig

c = zeros(t,1);
k1 = zeros(t,1);
for i = 1:t;
    k1(i) = theta(1) + theta(2)*k0 + 0.5 * theta(3) * k0^2;
    c(i) = prodfunc(k0) - k1(i);
    k0 = k1(i);
end;
v = 0;
for i = 1:t;
    if sig == 1;
        v = v + beta^(i-1)*log(c(i));
    else
        v = v + beta^(i-1)*(1/(1-sig))*(c(i))^(1-sig);
    end;
end;
v = -v;

end

