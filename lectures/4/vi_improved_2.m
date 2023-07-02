function [et,it,vf,p] = vi_improved_2(k_grid)
%%
% VI_IMPROVED_2.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.21.13
%
% PURPOSE   Performs an improved value function iteration method by solving for
%           the tomorrow's capital in the Euler equation
% USAGE     [et,it,vf,p] = vi_improved_2(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : policy function for tomorrow's capital
% USES      retrn.m
%           prodfunc.m
%           fsys.m
% NOTES     Remains to include the case delta = 1 (07.10.13)
%	    And sigma\=1? (08.04.13)
%
global beta

ngrid = length(k_grid);
%--------------------------------------------------------------------------
% starting value function interation
%--------------------------------------------------------------------------
v0 = retrn(k_grid,k_grid)/(1-beta);
crit = 1;
tol = 10^-2;

tic;
it = 0;
while crit > tol*(1-beta);
    it = it + 1;
    v1 = zeros(ngrid,1);
    p = zeros(ngrid,1);
    
    for i = 1:ngrid;
        c = fsolve(@(c) fsys(c,k_grid,v0,k_grid(i)),0.01*prodfunc(k_grid(i)),optimset('MaxIter',50,'MaxFunEvals',50,'Display','Off'));
        kiter = prodfunc(k_grid(i))-c;
        p(i) = kiter;
        v1(i) = retrn(k_grid(i),kiter) + beta*pw_linear(k_grid,v0,kiter);
    end
    
    crit = max(abs(v1-v0));
    disp(crit);
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn(k_grid,p);
vf = v0;

end