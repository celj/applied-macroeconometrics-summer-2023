function [ms] = geometric(m,p)
%%
% GEOMETRIC.M
% Numerical methods class
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in January 2007
% Modified on 07.17.13
% 
% PURPOSE   Generates random numbers following the geometric distribution
% USAGE     ms = geometric(m,p)
% INPUTS    m   : length of the simulated geometric series
%           p   : probability of success
% OUTPUT    ms  : vector of random numbers following the geometric distribution
%
ms = zeros(m,1);
for i = 1:m;
    u = rand(1,1);
    ms(i,:) = ceil(log(1-u)/log(1-p));
end;

end