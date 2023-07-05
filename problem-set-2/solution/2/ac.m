function [y_grid,p] = ac(rho,sig_u,ngrid)
%%
% (incomplete)
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
% 
%
sig_y = sqrt(sig_u^2/(1-rho^2));
for i = 1:ngrid+1
    cutoffs(i) = sig_y*norminv((i-1)/ngrid,0,1);
end

% Defining mean values of the intervals
% By definition, the expected value of x_{t} given that x_{t} lies ...
% in the interval [x_{n} x_{n+1}] is
b = 1;
for i = 1:ngrid
    y_grid(i) = -sig_y*( ( pdf('norm',cutoffs(i+1)/sig_y,0,1)-pdf('norm',cutoffs(i)/sig_y,0,1) ) / ( cdf('norm',cutoffs(i+1)/sig_y,0,1)-cdf('norm',cutoffs(i)/sig_y,0,1) ));
end
p = zeros(ngrid,ngrid);
end