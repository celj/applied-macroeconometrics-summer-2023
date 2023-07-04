function [ms] = exponential(m,theta)
%%
% EXPONENTIAL.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in January 2007
% Modified on 07.17.13
%
% PURPOSE   Generates random numbers following the exponential distribution
% USAGE     ms = exponential(m,theta)
% INPUTS    m       : length of the simulated exponential series
%           theta   : rate parameter
% OUTPUT    ms      : vector of random numbers following the exponential distribution
%
ms = zeros(m,1);
for i = 1:m;
    u = rand(1,1);
    ms(i,:) = -(1/theta)*log(1-u);
end;

end