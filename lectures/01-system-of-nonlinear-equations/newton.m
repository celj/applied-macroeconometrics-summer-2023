function [iter_store,sol_store,error_store] = newton(x0)
%%
% NEWTON.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% Created in July, 2013
%

crit = 100;
tol = 10^-5;
sol_store = [];
iter_store = [];

i = 0;
while crit > tol;
    i = i+1;
    [f,df] = func0(x0);
    %disp([f df]);
    dirtn = -df\f;
    x1 = x0+dirtn;
    crit = norm(x1-x0);
    sol_store = [sol_store;x1];
    iter_store = [iter_store;i];
    x0 = x1;
end;

error_store = abs(sol_store-x0);
