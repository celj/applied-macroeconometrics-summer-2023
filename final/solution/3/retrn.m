function [u] = retrn(x,y,z)
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
global sig

c = prodfunc(x,z)-y;
if sig == 1;
    u = log(c);
else
    u = (1/(1-sig))*(c).^(1-sig);
end;
    
end