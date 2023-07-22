function [mm] = mm_smith(rho)
%%
% MM_SMITH.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.16.13
% Modified on 07.22.13
%
% PURPOSE   Computes the moment condition vector (qx1) as distance between vector
%           of parameters of the auxiliary model and its simulated counterpart
% USAGE     mm = mm_smith(rho)
% INPUT     rho : structural parameter
% OUTPUT    mm  : vector of moments (qx1)
%
global r N
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
% estimation of the axuliary model based on simulated series
[thetas,~,sig2s,~] = estimavar(rs,1);
thetas = [thetas';sig2s];

mm = (theta-thetas);
end