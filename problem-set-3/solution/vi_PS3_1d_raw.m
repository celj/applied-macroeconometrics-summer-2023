function [et,it,vf,p,n,c] = vi_PS3_1d_raw(k_grid,n_grid)
%%
% VI_PS3_1D_RAW.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

global beta

ngrid = length(k_grid);
%--------------------------------------------------------------------------
% starting value function interation
%--------------------------------------------------------------------------
v0 = retrn_PS3_1d(k_grid,k_grid,n_grid)/(1-beta);
crit = 1;
tol = 10^-2;

tic;
it = 0;
while crit > tol*(1-beta)
    it = it + 1;
    tag = 1;
    v1 = zeros(ngrid,1);
    p = zeros(ngrid,1);
    
    for i = 1:ngrid
        vstore = -999*ones(ngrid,1);
        for j = tag:ngrid
            n = fsolve(@(n) fsys_raw(n,k_grid(i),k_grid(j)),0.1,optimset('MaxIter',50,'MaxFunEvals',50,'Display','Off'));
            c = prodfunc_PS3_1d(k_grid(i),n)-k_grid(j);
            if c > 0
                vstore(j) = retrn_PS3_1d(k_grid(i),k_grid(j),n) + beta*v0(j);
            end
            if j > 1 && vstore(j) < vstore(j-1)                            % exploiting concavity of the value function
                break;
            end
        end
        [vjj,jj] = max(vstore);
        v1(i) = vjj;
        p(i) = jj;
        tag = jj;
    end
    
    crit = max(abs(v1-v0));
    disp(crit);
    v0 = v1;
end
toc;
et = toc;
vf = v0;
% finding the policy functions for n (hours worked) and c (consumption)
for i = 1:ngrid
    n(i) = fsolve(@(n) fsys_raw(n,k_grid(i),k_grid(p(i))),0.1,optimset('MaxIter',50,'MaxFunEvals',50,'Display','Off'));
    c(i) = prodfunc_PS3_1d(k_grid(i),n(i))-k_grid(p(i));
end
n = n';
c = c';
end