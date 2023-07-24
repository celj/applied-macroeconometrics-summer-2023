function [F] = prodfunc(x,z)
%%
% PRODFUNC.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.18.13
%
% PURPOSE   Evaluates the market resources function
% USAGE     F = prodfunc(x,z)
% INPUTS    x : today's capital (1x1)
%           z : today's shock (1x1)
% OUTPUTS   F : evaluated market resources function
%
global alpha delta fac

F = exp(z).*fac*x.^alpha+(1-delta).*x;

end