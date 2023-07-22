%%
% EI_EX0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.02.13
% Modified on 07.10.13
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
sig = 1.5;
delta = 0.05;
ngrid = 250;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
k_0 = 0.01;
k_T = 5;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
%--------------------------------------------------------------------------
% policy function implementation
%--------------------------------------------------------------------------
[et,it,vf,p] = vi_improved_0(k_grid);
disp('Improved VI - Judd(1998) done ...');
[et0,it0,polk0] = ei_baxter(k_grid);
disp('EI - Baxter(1990) done ...');
%%
kstore = (alpha*beta)*k_grid.^alpha;
cstore = (1-alpha*beta)*k_grid.^alpha;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);

polk = k_grid(p);
%% figures for vi&ei
d = 5;
figure(1);
hold on;
plot(k_grid(1:d:end),polk(1:d:end),'ko',k_grid(1:d:end),polk0(1:d:end),'r+');
plot(k_grid,polk,'k--',k_grid,polk0,'r:');
l1 = legend('Improved VI','EI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;