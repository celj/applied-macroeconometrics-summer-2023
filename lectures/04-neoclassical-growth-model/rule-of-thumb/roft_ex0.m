%%
% ROFT_EX0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.03.13
% Modified on 07.04.13
%
% PURPOSE   Solves the deterministic neoclassical growth model using rules
%           of thumb as in Smith(1990) and Smith(1991).
%
close all;
clear all;
clc;
global alpha beta delta sig fac
%--------------------------------------------------------------------------
% setting parameters
%--------------------------------------------------------------------------
alpha = 0.3;
beta = 0.95;
sig = 1;
delta = 1;
ngrid = 100;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
k_0 = 0.01;
k_T = 3;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';

% analytical solution
kstore = (alpha*beta).*k_grid.^alpha;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);
% using the value function iteration
[et,it,vf,k] = vi_improved_0(k_grid);
%--------------------------------------------------------------------------
% starting the parametric approximation
%--------------------------------------------------------------------------
k0 = k_ss;                                                                  % initial state value
T = 800;                                                                    % approximate infinite horizon
% using a linear policy rule
x0 = [0.01;0.1];
[et0,v0,k0,x] = smith_linear(x0,k0,T,k_grid);
%--------------------------------------------------------------------------
% figures
%--------------------------------------------------------------------------
d = 5;
figure(1);
subplot(1,2,1);
hold on;
plot(k_grid(1:d:end),v0(1:d:end),'bo',k_grid(1:d:end),vstore(1:d:end),'k+',k_grid(1:d:end),vf(1:d:end),'r+')
plot(k_grid,v0,'b--',k_grid,vstore,'k-',k_grid,vf,'r-');
xlabel('capital grid');
ylabel('value function');
l1 = legend('linear','exact','VI','Location','NorthWest');
set(l1,'box','off');
hold off;
kpol = k_grid(k);
subplot(1,2,2);
hold on;
plot(k_grid(1:d:end),k0(1:d:end),'bo',k_grid(1:d:end),kstore(1:d:end),'k+',k_grid(1:d:end),kpol(1:d:end),'r+')
plot(k_grid,k0,'b--',k_grid,kstore,'k-',k_grid,kpol,'r-');
xlabel('capital grid');
ylabel('policy function for tomorrow''s capital');
l2 = legend('linear','exact','VI','Location','NorthWest');
set(l2,'box','off');
hold off;