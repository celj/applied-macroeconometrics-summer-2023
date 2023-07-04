function [et,it,vf,p] = pi_improved(k_grid)
%%
% PI_IMPROVED.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 07.10.13
%
% PURPOSE   Performs an improved policy function iteration method. It exploits the
%           monotonicity of the policy function and the concavity of the
%           value function. Instead of solving a large system of linear
%           equations it uses an iterative approach
% USAGE     [et,it,vf,p] = pi_improved(k_grid)
% INPUTS    k_grid  : capital grid (Nx1)
% OUTPUTS   et      : elapsed time
%           it      : number of iterations
%           vf      : value function
%           p       : policy function for tomorrow's capital
% USES      retrn.m
%           prodfunc.m
% NOTE      Remains to include the case delta = 1 together with ngrid-500 (07.10.13)
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
    Q = zeros(ngrid,ngrid);
    it = it + 1;
    tag = 1;
    p = zeros(ngrid,1);
    
    for i = 1:ngrid;
        vstore = -999*ones(ngrid,1);
        for j = tag:ngrid;                                                  % exploiting monotonicity of the policy function
            c = prodfunc(k_grid(i))-k_grid(j);
            if c > 0;
                vstore(j) = retrn(k_grid(i),k_grid(j)) + beta*v0(j);
            end;
            if j > 1 && vstore(j) < vstore(j-1);                            % exploiting concavity of the value function
                break;
            end;
        end;
        [~,jj] = max(vstore);
        p(i) = jj;
        tag = jj;                                                           % exploiting monotonicity of the policy function
        Q(i,jj) = 1;
    end
    %----------------------------------------------------------------------
    % iterative procedure
    %----------------------------------------------------------------------
    k = 0;
    v0k = v0;
    while k < 1
        v1k = retrn(k_grid,k_grid(p)) + beta*Q*v0k;
        v0k = v1k;
        k = k + 1;
    end;
    %----------------------------------------------------------------------
    % solving the system of linear equations
    % v1 = ( eye(ngrid)-beta*Q )\retrn(k_grid,k_grid(p));
    %----------------------------------------------------------------------
    v1 = v0k;
    crit = max(abs(v1-v0));
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn(k_grid,k_grid(p));
vf = v0;

end