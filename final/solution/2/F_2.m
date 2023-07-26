%%
% F_2.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

clear all
clc;
global alpha beta delta sig fac

%--------------------------------------
% Setting parameters                 
%--------------------------------------
alpha = 0.3;
beta = 0.95;
sig = 1;
delta = 0.05;
ngrid = 1000;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
k_0 = 0.01;
k_T = 20;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
%--------------------------------------
% value function implementation
%--------------------------------------
[et,it,vf,p] = vi_improved_0(k_grid);
[et0,it0,vf0,p0] = vi_improved_0_co(k_grid);
polk = k_grid(p);
polk0 = k_grid(p0);
poli = polk-(1-delta)*k_grid;
poli0 = polk0-(1-delta)*k_grid;
polc = prodfunc(k_grid) - polk;
polc0 = prodfunc(k_grid) - polk0;

kstore = (alpha*beta)*k_grid.^alpha-(1-delta)*k_grid;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);
%% figures for vi
d = 5;
figure(1);
subplot(2,2,1);
hold on;
plot(k_grid(1:d:end),poli0(1:d:end),'ko',k_grid(1:d:end),poli(1:d:end),'b+');
plot(k_grid,poli0,'k--',k_grid,poli,'b:');
l1 = legend('improved VI constrained','improved VI unconstrained.','Location','NorthEast');
set(l1,'box','off');
xlabel('capital grid');
ylabel('investment');
hold off;
subplot(2,2,2);
hold on;
plot(k_grid(1:d:end),vf0(1:d:end),'ko',k_grid(1:d:end),vf(1:d:end),'b+');
plot(k_grid,vf0,'k--',k_grid,vf,'b:');
l1 = legend('improved VI constrained','improved VI unconstrained','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;
subplot(2,2,3);
hold on;
plot(k_grid(1:d:end),polk0(1:d:end),'ko',k_grid(1:d:end),polk(1:d:end),'b+');
plot(k_grid,polk0,'k--',k_grid,polk,'b:');
l1 = legend('improved VI constrained','improved VI unconstrained','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('tomorrows capital');
hold off;
subplot(2,2,4);
hold on;
plot(k_grid(1:d:end),polc0(1:d:end),'ko',k_grid(1:d:end),polc(1:d:end),'b+');
plot(k_grid,polc0,'k--',k_grid,polc,'b:');
l1 = legend('improved VI constrained','improved VI unconstrained','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('consumption');
hold off;