function [mm] = mm(rho)
%%
% MM.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.15.13
% Modified on 07.22.13
%
% PURPOSE   Computes the moment condition vector (qx1) as in GMM
% USAGE     mm = mm(rho)
% INPUT     rho : structural parameter
% OUTPUT    mm  : vector of moments (nxq)
%
global data nq
T = length(data);
mm = zeros(T,nq);
mm(:,1) = (((rho(1)*(data(:,1).^-rho(2)).*data(:,2))-1).*data(:,3));
mm(:,2) = (((rho(1)*(data(:,1).^-rho(2)).*data(:,2))-1).*data(:,4));
mm(:,3) = (((rho(1)*(data(:,1).^-rho(2)).*data(:,2))-1).*data(:,5));

end