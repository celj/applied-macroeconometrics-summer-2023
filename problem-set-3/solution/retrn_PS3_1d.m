function [u] = retrn_PS3_1d(k,y,n)
%%
% RETRN_PS3_1D.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

global sig sh

c = prodfunc_PS3_1d(k,n)-y;
if sig == 1;
    u = sh.*log(c) + (1-sh).*log(1-n);
end;
end




