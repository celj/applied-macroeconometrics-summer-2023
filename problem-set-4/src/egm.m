function [et,it,vf,p] = egm(k_grid)
%%
% EGM.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 07.10.13
%
% PURPOSE   Performs the basic value function iteration method using an
%           endogenous grid method (by solving for the today's capital in the Euler
%           equation). Notice that k_grid is the grid over the tomorrows' market
%           resources (y = F(k)=c+k')
% USAGE     [et,it,vf,p] = egm(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : today's capital
% USES      retrn.m
%           prodfunc.m
%           fsys_k.m
% NOTES     Remains to include the case delta = 1 (07.10.13)
%
global sig beta

ngrid = length(k_grid);
%--------------------------------------------------------------------------
%Starting value function interation
%--------------------------------------------------------------------------
y_grid = prodfunc(k_grid);
v0 = retrn_egm(y_grid,k_grid)/(1-beta);                                     % value function defined on y_grid/k_grid
crit = 1;
tol = 10^-2;

tic;
it = 0;
while crit > tol*(1-beta);
    it = it + 1;
    v1 = zeros(ngrid,1);
    s = zeros(ngrid,1);
    
    for j = 1:ngrid;
        s(j) = slope(k_grid,v0,k_grid(j));
    end
    c = (beta*s).^(-1/sig);
    y = c + k_grid;                                                         % today's market resources
    v1y = retrn_egm(y,k_grid) + beta*v0;                                    % updating the value function on y ('y endog')
    for i = 1:ngrid;
        v1(i) = pw_linear(y,v1y,y_grid(i));                                 % interpolation on y_grid
    end;
    
    crit = max(abs(v1-v0));
    v0 = v1;
end
toc;
et = toc;
%--------------------------------------------------------------------------
% finding the today's k
%--------------------------------------------------------------------------
p = zeros(ngrid,1);
for i = 1:ngrid;
    k = fsolve(@(k) fsys_k(k,y(i)),0.01*y(i),optimset('MaxIter',50,'MaxFunEvals',50,'Display','Off','Algorithm','Levenberg-Marquardt'));
    p(i) = k;
end;
%vfs = 1/(1-beta)*retrn_egm(y,k_grid);
vf = v0;

end