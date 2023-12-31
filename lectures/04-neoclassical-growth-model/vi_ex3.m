%%
% VI_EX0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 06.21.13
%
% PURPOSE   Compares the performance of several methods to solve for the
%           deterministic neoclassical growth model
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
delta = 0.05;
ngrid = 100;
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
[et,it,vf,p] = vi(k_grid);
disp('Raw VI done ...');
[et0,it0,vf0,p0] = vi_improved_0(k_grid);
disp('Improved VI - Judd(1998) done ...');
[et1,it1,vf1,p1] = pi_improved(k_grid);
disp('Improved PI done ...');
%%
kstore = (alpha*beta)*k_grid.^alpha;
cstore = (1-alpha*beta)*k_grid.^alpha;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);

polk = k_grid(p);
polk0 = k_grid(p0);
polk1 = k_grid(p1);
%% figures for vi&pi
d = 5;
figure(1);
subplot(1,2,1);
hold on;
plot(k_grid(1:d:end),polk(1:d:end),'ko',k_grid(1:d:end),polk0(1:d:end),'b+',k_grid(1:d:end),polk1(1:d:end),'r*');
plot(k_grid,polk,'k--',k_grid,polk0,'b:',k_grid,polk1,'r--');
l1 = legend('raw VI','improved VI','improved PI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;
subplot(1,2,2);
hold on;
plot(k_grid(1:d:end),vf(1:d:end),'ko',k_grid(1:d:end),vf0(1:d:end),'b+',k_grid(1:d:end),vf1(1:d:end),'r*');
plot(k_grid,vf,'k--',k_grid,vf0,'b:',k_grid,vf1,'r--');
l1 = legend('raw VI','improved VI','improved PI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;