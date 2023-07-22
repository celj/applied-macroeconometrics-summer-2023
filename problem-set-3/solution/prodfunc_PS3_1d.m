function [F] = prodfunc_PS3_1d(k,n)
%%
% PRODFUNC_PS3_1D.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

global alpha delta fac

F = fac*(k.^alpha).*(n.^(1-alpha))+(1-delta).*k;

end