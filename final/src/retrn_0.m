function [u] = retrn_0(x,y)
%%
% retrn_0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.18.13
%
% PURPOSE   Evaluates the return function
% USAGE     u = retrn_0(x,y)
% INPUTS    x : today's capital (1x1)
%           y : tomorrow's capital (1x1)
% OUTPUTS   u : evaluated return function
% USES      prodfunc.m
%

c = prodfunc(x) - y;
u = log(c);
end
