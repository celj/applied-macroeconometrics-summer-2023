function [u] = retrn_stoch(x,y,z)
%%
% RETRN.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 07.08.13
%
% PURPOSE   Evaluates the return function
% USAGE     u = retrn(x,y,z)
% INPUTS    x : today's capital (1x1)
%           y : tomorrow's capital (1x1)
%           z : today's shock (1x1)
% OUTPUTS   u : evaluated return function
% USES      prodfunc.m
%

c = prodfunc_stoch(x,z) - y;
u = log(c);

end