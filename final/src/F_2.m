% Final Exam
% Applied Macroeconometrics
% ITAM
% Summer 2023
% Written by Carlos Lezama

% The investment constraint pulls the consumption policy function down
% due to the fact that a portion of this is reserved for investment.
% Note that c = y - i, where y = f(k) = k^alpha.

setup;

global alpha beta delta

alpha = 0.3;
beta = 0.95;
delta = 0.05;

k_grid = linspace(0.01, 20, 1e3)';
% v(k) and c are reduced,
% but the capital policy function increases.
% If we change 20 to 200, we can observe greater differences.
% We may also want to increase the number of grid points
% to smooth out the policy functions.

[et0,it0,v0,p0] = vi_improved_0(k_grid);
[et1,it1,v1,p1] = vi_improved_1(k_grid);

policy_k0 = k_grid(p0);
policy_k1 = k_grid(p1);

policy_c0 = prodfunc(k_grid) - policy_k0;
policy_c1 = prodfunc(k_grid) - policy_k1;

policy_i = (k_grid .^ alpha) - policy_c1;
% slightly goes below zero, barely noticeable

for i = 1:length(policy_i)
    if policy_i(i) < 0
        policy_i(i) = 0;
    end
end

figure;

subplot(2,2,1);
hold on;
plot(k_grid, policy_k0,'DisplayName','With no interest constraint');
plot(k_grid, policy_k1,'DisplayName','With interest constraint');
xlabel('$k$');
ylabel('$k^\prime$');
grid on;
legend;
hold off;

subplot(2,2,2);
hold on;
plot(k_grid, policy_c0,'DisplayName','With no interest constraint');
plot(k_grid, policy_c1,'DisplayName','With interest constraint');
xlabel('$k$');
ylabel('$c$');
grid on;
legend;
hold off;

subplot(2,2,3);
hold on;
plot(k_grid, v0,'DisplayName','With no interest constraint');
plot(k_grid, v1,'DisplayName','With interest constraint');
xlabel('$k$');
ylabel('$v(k)$');
grid on;
legend;
hold off;

subplot(2,2,4);
hold on;
plot(k_grid, policy_i);
xlabel('$k$');
ylabel('$i$');
grid on;
hold off;