function [ms] = weibull(m,theta,k)
%%
% PS2_3.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
% 
% PURPOSE   Generates random numbers following the weibull distribution
% USAGE     ms = weibull(m,theta,k)
% INPUTS    m       : length of the simulated weibull series
%           theta   : rate parameter
%           k       : shape parameter
% OUTPUT    ms      : vector of random numbers following the weibull distribution
% COMMENTS  when k = 1, this script works like exponential.m
%
ms = zeros(m,1);
for i = 1:m;
    u = rand(1,1);
    ms(i,:) = (1/theta)*(-log(1-u))^(1/k); 
end;

end