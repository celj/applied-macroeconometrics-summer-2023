%%
% PS4_1.M
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
y_0 = 0.1;
y_T = 2;
ystep = (y_T-y_0)/(ngrid-1);
y_grid = (y_0:ystep:y_T)';
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et1,it1,vf1,polk1] = vi_interpol_2(k_grid);
disp('Improved VI with quadratic interpolation done ...');
polc1 = prodfunc(k_grid)-polk1;
%% figures for vi-interpolation
d = 5;
figure(1);
subplot(2,2,1);
hold on;
plot(k_grid(1:d:end),polk1(1:d:end),'b+');
plot(k_grid,polk1,'b:');
l1 = legend('quad. intp.','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;
subplot(2,2,2);
hold on;
plot(k_grid(1:d:end),vf1(1:d:end),'b+');
plot(k_grid,vf1,'b:');
l1 = legend('quad. intp.','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;
subplot(2,2,3);
hold on;
plot(k_grid(1:d:end),polc1(1:d:end),'b+');
plot(k_grid,polc1,'b:');
l1 = legend('quad. intp.','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for consumption');
hold off;