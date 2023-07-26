function [et,it,vf,p] = vi_improved_0_ra(k_grid,z_grid,PI)
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
%           prodfunc.m
%    
global beta

ngrid = length(k_grid);
zngrid = length(z_grid);
%--------------------------------------------------------------------------
% starting value function interation 
%--------------------------------------------------------------------------
v0 =  kron(retrn(k_grid,k_grid,0)/(1-beta),ones(1,zngrid));
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
            vstore = -999*ones(ngrid,zngrid);
            for j = tag:ngrid;                                              % exploiting monotonicity of the policy function
                c = prodfunc(k_grid(i),z_grid(h))-k_grid(j);
                if c > 0;
                    vstore(j,h) = retrn(k_grid(i),k_grid(j),z_grid(h)) + beta*( PI(h,:)*v0(j,:)');
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