function [iter_store,sol_store,error_store] = bisection(xl,xr,fx)
%%
% BISECTION.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% Created on July, 2013
%

crit = 100;
tol = 10^-5;

fx = inline(fx);
xx = [xl;xr];
xm_0 = mean(xx);

sol_store = [];
iter_store = [];

f0=feval(fx,xm_0);
if f0 ~= 0;
    i = 0;
    while crit > tol
        i = i+1;
        f0=feval(fx,xm_0);
        fl=feval(fx,xl);
        if f0*fl < 0;
            xr = xm_0;
        elseif f0*fl > 0;
            xl = xm_0;
        end;
        xm_1 = mean([xl;xr]);
        crit = abs(xm_1-xm_0);
        sol_store = [sol_store;xm_1];
        iter_store = [iter_store;i];
        xm_0 = xm_1;
    end;
    error_store = abs(sol_store-xm_0);
end;
