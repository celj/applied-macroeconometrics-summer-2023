function [et,it,vf,p] = vi_improved_1(k_grid)
%%
% VI_IMPROVED_0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.21.13
%
% PURPOSE   Performs an improved value function iteration method by
%           exploiting both the monotonicity of the policy function and the concavity
%           of the value function as in Judd(1998)
% USAGE     [et,it,vf,p] = vi_improved_0(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : policy function index for tomorrow's capital
% USES      retrn_0.m
%           prodfunc.m
%
global alpha beta delta

ngrid = length(k_grid);
%--------------------------------------------------------------------------
% starting value function interation
%--------------------------------------------------------------------------
v0 = retrn_1(k_grid,k_grid)/(1-beta);
crit = 1;
tol = 10^-2;

tic;
it = 0;
while crit > tol*(1-beta);
    it = it + 1;
    tag = 1;
    v1 = zeros(ngrid,1);
    p = zeros(ngrid,1);
    for i = 1:ngrid;
        vstore = -999*ones(ngrid,1);
        for j = tag:ngrid;                                                  % exploiting monotonicity of the policy function
            if k_grid(j) - ((1 - delta) * k_grid(i)) >= 0
                c = prodfunc(k_grid(i)) - k_grid(j);
            else
                c = k_grid(i) .^ alpha;
            end
            if c > 0;
                vstore(j) = retrn_1(k_grid(i),k_grid(j)) + beta*v0(j);
            end;
            if j > 1 && vstore(j) < vstore(j-1);                            % exploiting concavity of the value function
                break;
            end;
        end;
        [vjj,jj] = max(vstore);
        v1(i) = vjj;
        p(i) = jj;
        tag = jj;                                                           % exploiting monotonicity of the policy function
    end
    crit = max(abs(v1-v0));
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn_0(k_grid,k_grid(p));
vf = v0;

end