function [et,it,vf,p] = vi_interpol_2(k_grid)
%%
% VI_INTERPOL_2.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.21.13
%
% PURPOSE   Performs an improved (exploiting the monotonicity of the
%           policy function) value function iteration method using
%           a quadratic spline interpolation to refine the policy function.
%           The refinement in the optimization step uses a golden
%           search method
% USAGE     [et,it,vf,p] = vi_interpol_2(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : policy function for tomorrow's capital
% USES      retrn.m
%           prodfunc.m
%           golden2.m
%
global beta

ngrid = length(k_grid);
%--------------------------------------------------------------------------
%Starting value function interation
%--------------------------------------------------------------------------
v0 = retrn(k_grid,k_grid)/(1-beta);
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
        for j = tag:ngrid;
            c = prodfunc(k_grid(i))-k_grid(j);
            if c > 0;
                vstore(j) = retrn(k_grid(i),k_grid(j)) + beta*v0(j);
            end;
        end;
        [~,jj] = max(vstore);
        tag = jj;
        if jj == 1;
            vss = [v0(jj);v0(jj+1);v0(jj+2)];
            jss = [jj;jj+1;jj+2];
        elseif jj == ngrid;
            vss = [v0(jj-2);v0(jj-1);v0(jj)];
            jss = [jj-2;jj-1;jj];
        else
            vss = [v0(jj-1);v0(jj);v0(jj+1)];
            jss = [jj-1;jj;jj+1];
        end;
        [kiter viter] = golden2(k_grid(i),k_grid(jss),vss);
        v1(i) = viter;
        p(i) = kiter;
    end;
    crit = max(abs(v1-v0));
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn(k_grid,p);
vf = v0;

end