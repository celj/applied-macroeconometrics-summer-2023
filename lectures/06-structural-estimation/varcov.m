function [vm] = varcov(rho)
%%
% VARCOV.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.15.13
% Modified on 07.22.13
%
% PURPOSE   Computes the covariance matrix of the GMM residuals
% USAGE     mm = mm(rho)
% INPUT     rho : structural parameter
% OUTPUT    vm  : covariance matrix (qxq)
%
T = length(mm(rho));
vm = mm(rho)'*mm(rho)/T;

end