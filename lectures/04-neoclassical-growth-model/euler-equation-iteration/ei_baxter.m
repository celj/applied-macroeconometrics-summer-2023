function [et,it,k1] = ei_baxter(k_grid)
%%
% EI_BAXTER.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.03.13
% Modified on 07.03.13
%
% PURPOSE   Performs the euler equation iteration method due to Baxter(1990)
% USAGE     [et,it,vf,p] = ei_baxter(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           k1      : policy function for tomorrow's capital
% USES      retrn_mg.m
%           prodfunc_mg.m
%           prodfunc.m
%
global beta sig

ngrid = length(k_grid);
%--------------------------------------------------------------------------
% starting policy function interation
%--------------------------------------------------------------------------
k0 = zeros(ngrid,1);
k1 = zeros(ngrid,1);
crit = 1;
tol = 10^-2;

tic;
it = 0;
while crit > tol*(1-beta);
    it = it + 1;
    for i = 1:ngrid;
        m = beta*prodfunc_mg(k_grid).*retrn_mg(k_grid,k0);
        z = k_grid + m.^(-1/sig);
        [~,j] = min(abs( prodfunc(k_grid(i))-z ));
        k1(i) = k_grid(j);
    end;
    
    crit = max(abs(k1-k0));
    k0 = k1;
end;
toc;
et = toc;

%vfs = 1/(1-beta)*retrn(k_grid,p);
end