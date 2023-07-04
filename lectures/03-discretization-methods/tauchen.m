function [y_grid,p] = tauchen(rho,sig_u,m,ngrid)
%%
% TAUCHEN.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in Fall 2011
% Modified on 07.10.13
%
% PURPOSE   Discretize an AR(1) of the form z = rho*z(-1)+u
% USAGE     [y_grid,p] = tauchen(rho,sig_u,m,ngrid)
% INPUTS    rho     : first-order serial autocorrelation (1x1)
%           sig_u   : standard deviation of u (1x1)
%           m       : scale factor to set the amplitud of the discretization
%                     on z (arbitrary) (1x1)
%           ngrid   : number of states (1x1)
% OUTPUTS   y_grid  : discretized values
%           p       : Markov transition matrix
% REFERENCES
%           1.Tauchen, G. (1986): "Finite State Markov-Chain Approximations to
%           Univariate and Vector Autoregressions." Economics Letters 20,
%           177-181.
%
sig_y = sqrt(sig_u^2/(1-rho^2));
% contructing the grid
y_end = m*sig_y;
y_beg = -y_end;
step = (y_end-y_beg)/(ngrid-1);
y_grid = (y_beg:step:y_end)';

if ngrid == 1;
    y_grid = 0;
    p = 1;
else
    p = zeros(ngrid,ngrid);
    % contructing the Markov matrix
    for k = 1:ngrid;
        if k == 1;
            for j = 1:ngrid;
                p(j,k) = cdf('norm',(y_grid(1)-rho*y_grid(j)+step/2)/sig_u,0,1);
            end;
        elseif k == ngrid;
            for j = 1:ngrid;
                p(j,k) = 1-cdf('norm',(y_grid(ngrid)-rho*y_grid(j)-step/2)/sig_u,0,1);
            end;
        else
            for j = 1:ngrid;
                p(j,k) = cdf('norm',(y_grid(k)-rho*y_grid(j)+step/2)/sig_u,0,1)-cdf('norm',(y_grid(k)-rho*y_grid(j)-step/2)/sig_u,0,1);
            end;
        end;
    end;
end;

end