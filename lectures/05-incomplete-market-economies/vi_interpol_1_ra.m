function [et,it,vf,p] = vi_interpol_1_ra(k_grid,q)
%%
% VI_INTERPOL_1_RA.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% Created on 10.27.11
% Modified on 06.21.13
%
% PURPOSE   Performs an improved (exploiting the monotonicity of the policy
%           function) value function iteration method using a linear interpolation to
%           refine the policy function. The refinement in the optimization step uses
%           a golden search method
% USAGE     [et,it,vf,p] = vi_interpol_1_ra(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : policy function for tomorrow's capital
% USES      retrn.m
%           golden1_ra.m
%
global beta

ngrid = length(k_grid);
%--------------------------------------------------------------------------
% discretizing the stochastic process
%--------------------------------------------------------------------------
zngrid = 2;
z_grid = [1;0.1];
PI = [0.925 0.075;0.5 0.5];
%--------------------------------------------------------------------------
% starting value function interation
%--------------------------------------------------------------------------
v0 = zeros(ngrid,zngrid);
crit = 1;
tol = 10^-2;
tic;
it = 0;
while crit > tol*(1-beta);
    it = it + 1;
    v1 = zeros(ngrid,zngrid);
    p = zeros(ngrid,zngrid);
    for h = 1:zngrid;
        tag = 1;
        for i = 1:ngrid;
            vstore = -inf*ones(ngrid,zngrid);
            for j = tag:ngrid;
                c = k_grid(i)+z_grid(h)-k_grid(j)*q;
                if c > 0;
                    vstore(j,h) = retrn(c) + beta*( PI(h,:)*v0(j,:)' ) - 0.5*1000*( min([k_grid(j)-k_grid(1);0]).^2 );
                end;
            end;
            [~,jj] = max(vstore(:,h));
            tag = jj;
            if jj == 1;
                vss = [v0(jj,:);v0(jj+1,:)];
                jss = [jj;jj+1];
            elseif jj == ngrid;
                vss = [v0(jj-1,:);v0(jj,:)];
                jss = [jj-1;jj];
            else
                vss = [v0(jj-1,:);v0(jj,:);v0(jj+1,:)];
                jss = [jj-1;jj;jj+1];
            end;
            [kiter viter] = golden1_ra(k_grid(i),k_grid(jss),vss,z_grid,PI,h,q);
            v1(i,h) = viter;
            p(i,h) = kiter;
        end;
    end;
    crit = max(max(abs(v1-v0)));
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn(k_grid,p);
vf = v0;
end