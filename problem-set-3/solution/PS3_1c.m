%%
% PS3_1c.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

close all
clc
global beta sig
v = 0;
for i = 1:T+1
    if sig == 1
        v = v + beta^(i-1)*log(c(i));
    else
        v = v + beta^(i-1)*(1/(1-sig))*(c(i))^(1-sig);
    end
end
disp('compare')
disp('value function evaluated at k0, which is')
disp(vf(51))
disp('to')
disp(' ')
disp('discounted sum of utilities, when starting at k0, which is')
disp(v)