function [et,it,vf,p] = vi_improved_1(k_grid)
%%
% VI_IMPROVED_1.M
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
%           of the value function as in Heer and Maussner (2004)
% USAGE     [et,it,vf,p] = vi_improved_1(k_grid)
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
    
    for i = 1:ngrid;                                                        % over states
        jmin = 1;
        jmax = ngrid;
        %------------------------------------------------------------------
        % binary search
        %------------------------------------------------------------------
        while (jmax-jmin) > 2;
            jl = floor( 0.5*(jmin+jmax) );
            ju = jl + 1;
            cl = prodfunc(k_grid(i))-k_grid(jl);
            cu = prodfunc(k_grid(i))-k_grid(ju);
            if cu > 0 && cl > 0;
                vu = retrn(k_grid(i),k_grid(ju)) + beta*v0(ju);
                vl = retrn(k_grid(i),k_grid(jl)) + beta*v0(jl);
            else
                vu = -999;
                vl = -999;
            end;
            
            if vu > vl;
                jmin = jl;
            else
                jmax = ju;
            end;
        end
        jmin1 = jmin + 1;
        vmin = retrn(k_grid(i),k_grid(jmin)) + beta*v0(jmin);
        vmin1 = retrn(k_grid(i),k_grid(jmin1)) + beta*v0(jmin1);
        vmax = retrn(k_grid(i),k_grid(jmax)) + beta*v0(jmax);
        jstore = [ jmin;jmin1;jmax ];
        vstore = [ vmin;vmin1;vmax ];
        [vjj,jj] = max(vstore);
        v1(i) = vjj;
        p(i) = jstore(jj);
        
    end
    crit = max(abs(v1-v0));
    v0 = v1;
end
toc;
et = toc;
%vfs = 1/(1-beta)*retrn(k_grid,k_grid(p));
vf = v0;

end