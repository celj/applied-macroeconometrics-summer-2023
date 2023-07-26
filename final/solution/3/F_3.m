%%
% F_3.M
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
delta = 0.05;
ngrid = 1000;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
c_ss = prodfunc(k_ss,0)-k_ss;
k_0 = 0.01;
k_T = 5;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
%--------------------------------------------------------------------------
% discretizing the stochastic process
%--------------------------------------------------------------------------
zngrid = 5;
rho = 0.95;
sig_u = 0.01;
[z_grid,PI] = rouwenhorst(rho,sig_u,zngrid);
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et,it,vf,p] = vi_improved_0_ra(k_grid,z_grid,PI);
polk = k_grid(p);
%% figures for vi
d = 5;
figure(1);
subplot(1,2,1);
hold on;
plot(k_grid(1:d:end),polk(1:d:end,:),'b--');
l1 = legend('improved stochastic VI','improved stochastic VI','improved stochastic VI','improved stochastic VI','improved stochastic VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for tomorrows capital');
hold off;
subplot(1,2,2);
hold on;
plot(k_grid(1:d:end),vf(1:d:end,:),'b--');
l1 = legend('improved stochastic VI','improved stochastic VI','improved stochastic VI','improved stochastic VI','improved stochastic VI','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;
%% simulations
clc
T = 20;
c = zeros(T+1,1);
k1 = zeros(T+1,1);
z_simul = mcsimul(z_grid,PI,T+1);
k0 = k_ss; % starts at the steady state value
for i = 1:T+1
    k1(i) = pw_linear(k_grid,k_grid(p(:,z_simul(i,2))),k0);
    c(i) = prodfunc(k0,z_simul(i,1)) - k1(i);
    k0 = k1(i);
end
k1 = [k0;k1(1:T,:)];
figure(2);
subplot(3,1,1)
plot((1:T+1)',z_simul(:,1),(1:T+1)',ones(T+1,1).*zeros(T+1,1),'k--')
l1 = legend('productivity shock time series','steady-state','Location','NorthWest');
xlim([0 19])
set(l1,'box','off');
xlabel('time');
ylabel('productivity shock');
hold off;
subplot(3,1,2)
hold on;
plot((1:T+1)',c,(1:T+1)',ones(T+1,1).*c_ss,'k--')
xlim([0 19])
l1 = legend('consumption time series','steady-state','Location','NorthWest');
set(l1,'box','off');
xlabel('time');
ylabel('consumption');
hold off;
subplot(3,1,3)
plot((1:T+1)',k1,(1:T+1)',ones(T+1,1).*k_ss,'k--')
xlim([0 19])
l1 = legend('capital stock time series','steady-state','Location','NorthWest');
set(l1,'box','off');
xlabel('time');
ylabel('capital');
hold off;
