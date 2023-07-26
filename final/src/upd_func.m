function [sj,ind] = upd_func(u,p,s,si)
%%
% UPD_FUNC.M
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
% USAGE     [psi,ind] = ini_func(u,p,s,si)
% INPUTS    u   : rand(1,1) (1x1)
%           p   : stochastic kernel (Markov transition matrix)
%           si  : current state
% OUTPUTS   sj  : next state
%           ind : index of sj in s
% REFERENCES
%           1.H�ggstr�m, O. (2003). "Finite Markov Chains and Algorithmic
%           Applications." London Mathematical Society. Student Texts 00.
%
%% finding the corresponding index of si in s
i = 1;
while si >= s(i) && i <= length(s)-1;
    i = i+1;
end;
index = i-1;
%% computing sj
ps = cumsum(p(index,:))';
i = 1;
while u >= ps(i);
    i = i+1;
end;
ind = i;
sj = s(i);

end