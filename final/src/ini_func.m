function [s0,ind] = ini_func(u,p,s)
%%
% INI_FUNC.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.15.13
% Modified on 07.17.13
%
% PURPOSE   Computes the initial realization of the Markov process
% USAGE     [psi,ind] = ini_func(p,s,u)
% INPUTS    u   : rand(1,1) (1x1)
%           p   : ergodic distribution as initialization
%           s   : states values
% OUTPUTS   s0  : initial draw from the p
%           ind : index of psi in s
% REFERENCES
%           1.Haggstrom, O. (2003). "Finite Markov Chains and Algorithmic
%           Applications." London Mathematical Society. Student Texts 00.
%
ps = cumsum(p);
i = 1;
while u >= ps(i);
    i = i+1;
end;
ind = i;
s0 = s(i);

end