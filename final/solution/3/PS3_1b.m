%%
% PS3_1b.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

close all
clc
T = 20;
c = zeros(T+1,1);
k1 = zeros(T+1,1);
z_simul = mcsimul(z_grid,PI,T+1);
k0 = k_ss;
for i = 1:T+1
    k1(i) = pw_linear(k_grid,k_grid(p(:,z_simul(i,2))),k0);
    c(i) = prodfunc(k0,z_simul(i,1)) - k1(i);
    k0 = k1(i);
end
k1 = [k0;k1(1:T,:)];
subplot(3,1,1)
plot((1:T+1)',z_simul(:,1),(1:T+1)',ones(T+1,1).*zeros(T+1,1),'k--')
l1 = legend('productivity shock series','steady-state','Location','NorthWest');
set(l1,'box','off');
xlabel('time');
ylabel('productivity shock');
hold off;
subplot(3,1,2)
hold on;
plot((1:T+1)',c,(1:T+1)',ones(T+1,1).*c_ss,'k--')
l1 = legend('consumption times series','steady-state','Location','NorthWest');
set(l1,'box','off');
xlabel('time');
ylabel('consumption');
hold off;
subplot(3,1,3)
plot((1:T+1)',k1,(1:T+1)',ones(T+1,1).*k_ss,'k--')
l1 = legend('capital stock time series','steady-state','Location','NorthWest');
set(l1,'box','off');
xlabel('time');
ylabel('capital');
hold off;