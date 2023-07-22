%%
% PS3_1b.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

close all
clc
T = 10000;
c = zeros(T+1,1);
k1 = zeros(T+1,1);
index = 51;
k0 = k_grid(index); % taking any initial value for k0 (could be k0 = 1.0120)
for i = 1:T+1
    k1(i) = pw_linear(k_grid,k_grid(p),k0);
    c(i) = prodfunc(k0) - k1(i);
    k0 = k1(i);
end
k1 = [k_grid(index);k1(1:T,:)];
%% c and k' in steady-state (see PS3_1a.m)
c_ss_v = ones(T+1,1).*c_ss;
k_ss_v = ones(T+1,1).*k_ss;
subplot(1,2,1)
plot((1:1:T+1)',[c c_ss_v])
l1 = legend('consumption','steady-state','Location','NorthWest');
set(l1,'box','off');
xlabel('period');
ylabel('value of sequences');
hold off;
subplot(1,2,2)
hold on;
plot((1:1:T+1)',[k1 k_ss_v])
l2 = legend('savings','steady-state','Location','NorthWest');
set(l2,'box','off');
xlabel('period');
ylabel('value of sequences');
%%
disp('consumption appears to converge to c_ss')
disp(c_ss)
disp('savings appears to converge to k_ss')
disp(k_ss)