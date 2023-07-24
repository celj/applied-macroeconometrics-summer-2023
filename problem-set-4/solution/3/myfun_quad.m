function [v] = myfun_quad(theta,k0,t)
%%
% MYFUN_QUAD.M
% Numerical methods class
% Summer 2013
% Written by Gustavo Leyva
% Created on 07.04.13
% Modified on 07.16.13
%
% PURPOSE   Computes the value function following a quadratic policy rule
% USAGE     v = myfun_quad(theta,k0,t)
% INPUTS    theta : parameters of the policy rule
%           k0    : initial capital value
%           t     : horizon
% OUTPUTS   v     : value function
%
global beta sig
c = zeros(t,1);
k1 = zeros(t,1);
for i = 1:t;
    k1(i) = theta(1) + theta(2)*k0 + theta(3)*k0^2;
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

