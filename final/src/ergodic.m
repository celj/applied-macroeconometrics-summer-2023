function [p0] = ergodic(p)
%%
% ERGODIC.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created in Fall 2011
% Modified on 07.16.13
%
% PURPOSE   Computes the ergodic vector of probabilities associated with
%           Markov transition matrix p
% USAGE     p0 = ergodic(p)
% INPUTS    p   : Markov transition matrix (rows-summing-to-one setup) (kxk)
% OUTPUTS   p0  : ergodic probabilities (kx1)
%
[ve,va] = eig(p');
[~,i] = max(diag(va));
p0 = abs(ve(:,i));                                                          % the eigvector of ergodic probabilities corresponding to the largest (1) eigvalue
p0 = p0/sum(p0);

end