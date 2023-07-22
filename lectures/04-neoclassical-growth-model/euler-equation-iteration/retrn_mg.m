function [u] = retrn_mg(x,y)
%%
% RETRN_MG.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.18.13
%
% PURPOSE   Evaluates the first derivative of the return function
% USAGE     u = retrn_mg(x,y)
% INPUTS    x : today's capital (1x1)
%           y : tomorrow's capital (1x1)
% OUTPUTS   u : evaluated return function
% USES      prodfunc.m
%
global sig

if sig == 1;
    u = (prodfunc(x)-y).^(-1);
else
    u = (prodfunc(x)-y).^(-sig);
end;

end