function [mm] = mm_gt(rho)
%%
% MM_GT.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.16.13
% Modified on 07.22.13
%
% PURPOSE   Computes the moment condition vector (qx1) as the scores of the
%           log-likelihood function evaluated with simulated data and 'true'
%           auxiliary model parameters
% USAGE     mm = mm_gt(rho)
% INPUT     rho : structural parameter
% OUTPUT    mm  : vector of moments (nxq)
%
global r nq N
[theta,~,sig2,~] = estimavar(r,1);
theta = [theta';sig2];

rng(25);
zngrid = 10;
lambda = 0.178;
sig_u = sqrt(0.0144);
[zs_grid,PI] = rouwenhorst(lambda,sig_u,zngrid);
zs_grid = zs_grid + (0.013)/(1-lambda);
mcstores = mcsimul(zs_grid,PI,N);
ys_index = mcstores(:,2);
% solving for q in the Euler equation
qs = zeros(N,1);
for i = 1:N;
    qs(i) = PI(ys_index(i),:)*(rho(1)*exp(zs_grid).^-rho(2));
end;
rs = 1./qs;                                                                 %rs_t+1 (simulated)
xs = [ones(N-1,1) rs(1:N-1,:)];
us = rs(2:N,:)-xs*theta(1:2,:);

mm = zeros(N-1,nq);
mm(:,1) = ((1/theta(3))*(xs(:,1).*us));
mm(:,2) = ((1/theta(3))*(xs(:,2).*us));
mm(:,3) = ((-0.5/theta(3))+((0.5/(theta(3)^2))*(us.^2)));
end