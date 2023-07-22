function [Q] = quadfunc_smith(rho,w)
%%
% QUADFUNC_SMITH.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.16.13
% Modified on 07.22.13
%
% PURPOSE   Constructs the quadratic form for the Smith(1993)'s
%           implementation of Indirect Inference.
% USAGE     Q = quadfunc_smith(rho,w)
% INPUTS    rho : initial value for the structural parameter
%           w   : weighting matrix
% OUTPUTS   Q   : value function
%
mms = mm_smith(rho);
Q = (mms)'*w*(mms);

end