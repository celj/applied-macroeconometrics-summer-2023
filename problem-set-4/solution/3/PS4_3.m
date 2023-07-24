%%
% PS4_3.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
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
k_analy = (alpha*beta).*k_grid.^alpha;
c_analy = (1-alpha*beta).*k_grid.^alpha;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);
% using the value function iteration
[et2,it2,vf2,k2] = vi_improved_0(k_grid);
%--------------------------------------------------------------------------
% starting the parametric approximation
%--------------------------------------------------------------------------
k0 = k_ss;                                                                  % initial state value
T = 800;                                                                    % approximate infinite horizon
% using a linear policy rule
x0 = [0.01;0.1];
[et_1,v_1,k1_1,x_1] = smith_linear(x0,k0,T,k_grid);
c1_1 = prodfunc(k_grid) - k1_1;
% using a quadratic policy rule
x0 = [0.001;0.13;-0.02];
[et_2,v_2,k1_2,x_2] = smith_quad(x0,k0,T,k_grid);
c1_2 = prodfunc(k_grid) - k1_2;
%--------------------------------------------------------------------------
% figures
%--------------------------------------------------------------------------
d = 5;
figure(1);
subplot(2,2,1);
hold on;
plot(k_grid(1:d:end),v_1(1:d:end),'bo',k_grid(1:d:end),vstore(1:d:end),'k+',k_grid(1:d:end),v_2(1:d:end),'b*',k_grid(1:d:end),vf2(1:d:end),'r+')
plot(k_grid,v_1,'b--',k_grid,vstore,'k-',k_grid,v_2,'b:',k_grid,vf2,'r-');
xlabel('capital grid');
ylabel('value function');
l1 = legend('linear','exact','quad','VI','Location','NorthWest');
set(l1,'box','off');
hold off;
kpol = k_grid(k2);
cpol = prodfunc(k_grid) - kpol;
subplot(2,2,2);
hold on;
plot(k_grid(1:d:end),k1_1(1:d:end),'bo',k_grid(1:d:end),k_analy(1:d:end),'k+',k_grid(1:d:end),k1_2(1:d:end),'b*',k_grid(1:d:end),kpol(1:d:end),'r+')
plot(k_grid,k1_1,'b--',k_grid,k_analy,'k-',k_grid,k1_2,'b:',k_grid,kpol,'r-');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
legend('linear','analy','quad','VI');
l2 = legend('linear','exact','quad','VI','Location','NorthWest');
set(l2,'box','off');
hold off;
subplot(2,2,3);
hold on;
plot(k_grid(1:d:end),c1_1(1:d:end),'bo',k_grid(1:d:end),c_analy(1:d:end),'k+',k_grid(1:d:end),c1_2(1:d:end),'b*',k_grid(1:d:end),cpol(1:d:end),'r+')
plot(k_grid,c1_1,'b--',k_grid,c_analy,'k-',k_grid,c1_2,'b:',k_grid,cpol,'r-');
xlabel('capital grid');
ylabel('policy function for consumption');
legend('linear','analy','quad','VI');
l2 = legend('linear','exact','quad','VI','Location','NorthWest');
set(l2,'box','off');
hold off;