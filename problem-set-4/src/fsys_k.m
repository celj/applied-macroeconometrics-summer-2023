function [F] = fsys_k(x,y)
%%
% PURPOSE   Solves for today's capital in the market resources function
% USAGE     F = fsys_k(x,y)
% INPUTS    x   : unknown today's capital (1x1)
%           y   : today's market resources (1x1)
% OUTPUT    F   : market resources equation
%
global alpha delta

F = y-x^alpha-(1-delta)*x;

end