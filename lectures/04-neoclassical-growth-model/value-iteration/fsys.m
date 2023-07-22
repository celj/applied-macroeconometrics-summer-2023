function [F] = fsys(x,k,v,ki)
%%
% PURPOSE   Solves for tomorrow's capital in the Euler equation. It uses k
%           and v to compute the first derivative of v evaluated in kj (see
%           below)
% USAGE     F = fsys(x,k,ki,v)
% INPUTS    x   : consumption (1x1)
%           ki  : today's capital (1x1)
%           k   : capital grid (Nx1)
%           v   : value function (Nx1)
% OUTPUT    F   : Euler equation
%
global beta

kj = prodfunc(ki) - x;
s = slope(k,v,kj);
F = 1/x-(beta*s);
end








