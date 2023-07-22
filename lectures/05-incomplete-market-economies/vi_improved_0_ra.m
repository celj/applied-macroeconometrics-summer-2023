function [et,it,vf,p] = vi_improved_0_ra(k_grid,q)
%%
% VI_IMPROVED_0_RA.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 07.08.13
%
% PURPOSE   Performs an improved value function iteration method by
%           exploiting both the monotonicity of the policy function and the concavity
%           of the value function as in Judd(1998)
% USAGE     [et,it,vf,p] = vi_improved_0_ra(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : policy function index for tomorrow's capital
% USES      retrn.m
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
tol = 10^-6;
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
            for j = tag:ngrid;                                              % exploiting monotonicity of the policy function
                c = k_grid(i)+z_grid(h)-k_grid(j)*q;
                if c > 0;
                    vstore(j,h) = retrn(c) + beta*( PI(h,:)*v0(j,:)') - 0.5*1000*( min([k_grid(j)-k_grid(1);0]).^2 );
                end;
                if j > 1 && vstore(j,h) < vstore(j-1,h);                    % exploiting concavity of the value function
                    break;
                end;
            end;
            [vjj,jj] = max(vstore(:,h));
            v1(i,h) = vjj;
            p(i,h) = jj;
            tag = jj;                                                       % exploiting monotonicity of the policy function
        end
    end;
    crit = max(max(abs(v1-v0)));
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn(k_grid,k_grid(p));
vf = v0;
end