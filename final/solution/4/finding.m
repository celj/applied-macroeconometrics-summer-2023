function [k_star] = finding(k_grid,gpol)
%%
% FINDING.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.14.13
% Modified on 07.14.13
% 
% PURPOSE   Finds a k_star on the k_grid such that k_star is the max over all x such
%           that g(x)=k_bar

if gpol(1) > min(k_grid);
    k_star = min(k_grid);
else
    i = 1;
    while gpol(i) <= min(k_grid)
        i = i+1;        
    end;
    k_star = k_grid(i-1);

end