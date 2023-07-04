%%
% BLOCK_BOOTSTRAP.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in December, 2006
% Modified on 07.17.13

% DESCRIPTION:
% This program performs a block bootstrap on the original sample.

% USAGE:
% store = block(y,l)

% where:
% y       =   original sample
% l       =   block's length
% store   =   block bootstrap sample

% REFERENCES:
% B�hlmann, P. (2002). "Bootstraps for Time Series." Statistical
% Science. Vol. 17, No. 1, pp. 52-72.
%
function [store] = block_bootstrap(y,l)

n = length(y);
blockl=[];

while length(blockl) < length(y);
    li = fix(1+(n-l+1)*rand(1,1));
    index = li:li+l-1;
    block = y(index,:);
    blockl = [blockl;block];
end;
store = blockl(1:n,:);

end
