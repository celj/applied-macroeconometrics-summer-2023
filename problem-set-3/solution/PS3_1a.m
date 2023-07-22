%%
% PS3_1a.M
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
sig = 1; % log
delta = 1;
ngrid = 250;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
c_ss = prodfunc(k_ss)- k_ss;
k_0 = 0.01;
k_T = 5;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
y_0 = 0.1;
y_T = 2;
ystep = (y_T-y_0)/(ngrid-1);
y_grid = (y_0:ystep:y_T)';
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et,it,vf,p] = vi(k_grid);
disp('Raw VI done ...');
%%
kstore = (alpha*beta)*k_grid.^alpha;
cstore = (1-alpha*beta)*k_grid.^alpha;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);

polk = k_grid(p);
%% figures for vi
d = 5;
figure(1);
subplot(1,2,1);
hold on;
plot(k_grid(1:d:end),polk(1:d:end),'ko');
plot(k_grid,polk,'k--');
l1 = legend('raw VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;
subplot(1,2,2);
hold on;
plot(k_grid(1:d:end),vf(1:d:end),'ko');
plot(k_grid,vf,'k--');
l1 = legend('raw VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;