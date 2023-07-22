function [Q] = quadfunc(rho,w)
%%
% QUADFUNC.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.15.13
% Modified on 07.22.13
%
% PURPOSE   Constructs the quadratic form for the standard GMM
% USAGE     Q = quadfunc(rho,w)
% INPUTS    rho : initial value for the structural parameter
%           w   : weighting matrix
% OUTPUTS   Q   : value function
%
mms = mm(rho);
Q = mean(mms)*w*mean(mms)';

end