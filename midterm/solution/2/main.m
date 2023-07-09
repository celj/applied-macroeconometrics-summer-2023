%%
% M_2.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%
close all;
clear all;
clc;
%%
% AR(1) parameters
rho = 0.9;
sig_u = sqrt(0.06);
sig_z = sqrt(sig_u^2/(1-rho^2));
ngrid = 50;
n = 250;
[grid0,p0] = rouwenhorst(rho,sig_u,ngrid);
z_simul = mcsimul(grid0,p0,n);
z_simul = z_simul(:,1);
subplot(2,2,1)
plot(z_simul)
%%
mu = 0.01;
sig = sqrt(0.25);
l = -0.5*z_simul + (mu+sig*randn(n,1));
subplot(2,2,2)
plot(l)
%%
R = corrcoef([z_simul l]);
s = R(2,1);
%%
m = 2000;
p = 1/12;
for i = 1:m
    d = stationary_bootstrap([z_simul l],p);
    R = corrcoef(d);
    R = R(2,1);
    store(i) = R;
end
subplot(2,2,3)
plot(d(:,1))
subplot(2,2,4)
plot(d(:,2))
%
a = 0.10;   % level of significance
store = sort(store);
store = store';
store_l = quantile(store,a/2);
store_u = quantile(store,1-a/2);
disp('The Efron confidence interval is:');
disp([store_l s store_u]);

