function [it,vf,p] = vi(k_grid)
%%
% VI.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.21.13
%
% PURPOSE   Performs the basic value function iteration method
% USAGE     [et,it,vf,p] = vi(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   it      : number of iterations
%           vf      : value function
%           p       : policy function index for tomorrow's capital
% USES      retrn.m
%           prodfunc.m
%

global beta

ngrid = length(k_grid);

%--------------------------------------------------------------------------
% starting value function interation
%--------------------------------------------------------------------------

v0 = retrn(k_grid,k_grid) / (1 - beta);
crit = 1;
tol = 1e-2;

it = 0;

while crit > tol * (1 - beta)
    it = it + 1;
    v1 = zeros(ngrid,1);
    p = zeros(ngrid,1);
    for i = 1:ngrid
        vstore = -999 * ones(ngrid,1);
        for j = 1:ngrid
            c = prodfunc(k_grid(i)) - k_grid(j);
            if c > 0
                vstore(j) = retrn(k_grid(i),k_grid(j)) + beta * v0(j);
            end
        end
        [vjj,jj] = max(vstore);
        v1(i) = vjj;
        p(i) = jj;
    end
    crit = max(abs(v1 - v0));
    v0 = v1;
end

vf = v0;

end