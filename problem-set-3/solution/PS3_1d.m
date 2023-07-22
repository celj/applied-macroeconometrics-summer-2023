%%
% PS3_1d.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
% Notes: this script works for delta = 1.
%        For a different delta, you have to modify it
%

clear all
clc;
global alpha beta delta sig fac sh z_ss
%--------------------------------------------------------------------------
% Setting parameters
%--------------------------------------------------------------------------
alpha = 0.3;
beta = 0.95;
sig = 1;
delta = 1;
sh = 0.6;
ngrid = 250;
fac = 1;
rho = (1/beta)-1;
z_ss = fsolve(@(z) ss(z),[0.5;0.6],optimset('MaxIter',50,'MaxFunEvals',50,'Display','Off'));
k_0 = 0.01;
k_T = (alpha*beta)^(1/(1-alpha));
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
n_grid = k_grid*(1/(alpha*beta)^(1/(1-alpha)));
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et,it,vf,p,n,c] = vi_PS3_1d_raw(k_grid,n_grid);
disp('Raw VI done ...');

polk = k_grid(p);
%% figures for vi
d = 5;
figure(1);
subplot(2,2,1);
hold on;
plot(k_grid(1:d:end),polk(1:d:end),'ko');
plot(k_grid,polk,'k--');
l1 = legend('raw VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;
subplot(2,2,2);
hold on;
plot(k_grid(1:d:end),vf(1:d:end),'ko');
plot(k_grid,vf,'k--');
l1 = legend('raw VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;
subplot(2,2,3);
hold on;
plot(k_grid(1:d:end),c(1:d:end),'ko');
plot(k_grid,c,'k--');
l1 = legend('raw VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for consumption');
hold off;
subplot(2,2,4);
hold on;
plot(k_grid(1:d:end),n(1:d:end),'ko');
plot(k_grid,n,'k--');
ylim([0 1])
l1 = legend('raw VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for hours worked');
hold off;