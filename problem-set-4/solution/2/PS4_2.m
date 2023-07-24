%%
% PS4_2.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

clear all
clc;
global alpha beta delta sig fac
%--------------------------------------------------------------------------
% Setting parameters
%--------------------------------------------------------------------------
alpha = 0.3;
beta = 0.95;
sig = 1;
delta = 1;
ngrid = 250;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
k_0 = 0.01;
k_T = 5;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et,it,vf,p] = vi_improved_0_ra(k_grid);
polk = k_grid(p);
%%
zngrid = 3;
rho = 0.95;
sig_u = 0.01;
[z_grid,~] = rouwenhorst(rho,sig_u,zngrid);
polc(:,1) = prodfunc(k_grid,z_grid(1)) - k_grid(p(:,1));
polc(:,2) = prodfunc(k_grid,z_grid(2)) - k_grid(p(:,2));
polc(:,3) = prodfunc(k_grid,z_grid(3)) - k_grid(p(:,3));

kstore1 = (alpha*beta).*exp(z_grid(1)).*k_grid.^alpha;
kstore2 = (alpha*beta).*exp(z_grid(2)).*k_grid.^alpha;
kstore3 = (alpha*beta).*exp(z_grid(3)).*k_grid.^alpha;
kstore = [kstore1 kstore2 kstore3];

cstore1 = (1-alpha*beta).*exp(z_grid(1)).*k_grid.^alpha;
cstore2 = (1-alpha*beta).*exp(z_grid(2)).*k_grid.^alpha;
cstore3 = (1-alpha*beta).*exp(z_grid(3)).*k_grid.^alpha;
cstore = [cstore1 cstore2 cstore3];

a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);
%% figures for vi
d = 5;
figure(1);
subplot(2,2,1);
hold on;
plot(k_grid(1:d:end),kstore(1:d:end,1),'ko',k_grid(1:d:end),kstore(1:d:end,2),'ko',k_grid(1:d:end),kstore(1:d:end,3),'ko',k_grid(1:d:end),polk(1:d:end,1),'b+',k_grid(1:d:end),polk(1:d:end,2),'b+',k_grid(1:d:end),polk(1:d:end,3),'b+');
plot(k_grid,kstore(:,1),'k--',k_grid,kstore(:,2),'k--',k_grid,kstore(:,3),'k--',k_grid,polk(:,1),'b:',k_grid,polk(:,2),'b:',k_grid,polk(:,3),'b:');
l1 = legend('exact','exact','exact','improved stochastic VI','improved stochastic VI','improved stochastic VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;
subplot(2,2,2);
hold on;
plot(k_grid(1:d:end),vstore(1:d:end),'ko',k_grid(1:d:end),vf(1:d:end,1),'b+',k_grid(1:d:end),vf(1:d:end,2),'b+',k_grid(1:d:end),vf(1:d:end,3),'b+');
plot(k_grid,vstore,'k--',k_grid,vf(:,1),'b:',k_grid,vf(:,2),'b:',k_grid,vf(:,3),'b:');
l1 = legend('exact','improved stochastic VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;
subplot(2,2,3);
hold on;
plot(k_grid(1:d:end),cstore(1:d:end,1),'ko',k_grid(1:d:end),cstore(1:d:end,2),'ko',k_grid(1:d:end),cstore(1:d:end,3),'ko',k_grid(1:d:end),polc(1:d:end,1),'b+',k_grid(1:d:end),polc(1:d:end,2),'b+',k_grid(1:d:end),polc(1:d:end,3),'b+');
plot(k_grid,cstore(:,1),'k--',k_grid,cstore(:,2),'k--',k_grid,cstore(:,3),'k--',k_grid,polc(:,1),'b:',k_grid,polc(:,2),'b:',k_grid,polc(:,3),'b:');
l1 = legend('exact','exact','exact','improved stochastic VI','improved stochastic VI','improved stochastic VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for consumption');
hold off;