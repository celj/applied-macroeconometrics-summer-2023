%%
% SIMPLE_BOOTSTRAP.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in December, 2006
% Modified on 07.17.13

% DESCRIPTION:
% This program performs a simple i.i.d resampling on the original
% sample.

% USAGE:
% store = simple_bootstrap(y)

% where:
% y       =   original sample
% store   =   the resampled version of the original data
%

function[store] = simple_bootstrap(y)

n = length(y);
index = fix(1+(n)*rand(n,1));
store = y(index,:);

end