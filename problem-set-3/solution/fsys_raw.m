function [F] = fsys_raw(n,ki,kj)
%%
% FSYS_RAW.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

global alpha sh
y = (ki.^alpha)*(n.^(1-alpha));
phi = (1-sh)./(sh);
F = (1-alpha).*y./n - (1-alpha+phi).*y + phi*kj;
end








