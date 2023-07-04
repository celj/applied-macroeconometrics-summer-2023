%%
% STATIONARY_BOOTSTRAP.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in January, 2007
% Modified on 07.17.13

% DESCRIPTION:
% This program performs the stationary bootstrap on the original sample.

% USAGE:
% store = stationary(y,p)

% where:
% y       =   original data
% p       =   success probability of the geometric distribution
% store   =   bootstrap sample

% REFERENCES:
% Politis D. and J. Romano (1994). "The Stationary Bootstrap."
% Journal of the American Statistical Association. Vol. 89, No.
% 428, pp. 1303-1313.
%
function [store] = stationary_bootstrap(y,p)

n = length(y);
blockl = [];

while length(blockl) < length(y);
    geo = geometric(1,p);
    li = fix(1+(n-geo+1)*rand(1,1));
    index = li:li+geo-1;
    block = y(index,:);
    blockl = [blockl;block];
end;
store = blockl(1:n,:);

end