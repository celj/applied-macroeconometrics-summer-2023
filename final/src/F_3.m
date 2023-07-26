% Final Exam
% Applied Macroeconometrics
% ITAM
% Summer 2023
% Written by Carlos Lezama

setup;

global alpha beta delta z_states

T = 20;
alpha = 0.3;
beta = 0.95;
delta = 0.05;
k_grid = linspace(0.01, 5, 1000)';
rho = 0.95;
sig_e = 0.01;
z_states = 5;

kss = (alpha / (((1 / beta) - 1) + delta)) ^ (1 / (1 - alpha));

[z_grid,p_z] = rouwenhorst(rho,sig_e,z_states);

[et,it,vf,p] = vi_improved_0_ra(k_grid,z_grid,p_z);

policy_k = k_grid(p);

for i=1:z_states
    policy_c(:,i) = prodfunc_stoch(k_grid,z_grid(i)) - policy_k(:,i);
end

figure
subplot(1,3,1)
hold on
for i=1:z_states
    plot(k_grid,vf(:,i),'DisplayName',sprintf('$z = %0.4f$',z_grid(i)))
end
legend
xlabel('$k$')
ylabel('$v(k)$')
hold off
subplot(1,3,2)
hold on
for i=1:z_states
    plot(k_grid,policy_k(:,i),'DisplayName',sprintf('$z = %0.4f$',z_grid(i)))
end
legend
xlabel('$k$')
ylabel('$k^\prime$')
hold off
subplot(1,3,3)
hold on
for i=1:z_states
    plot(k_grid,policy_c(:,i),'DisplayName',sprintf('$z = %0.4f$',z_grid(i)))
end
legend
xlabel('$k$')
ylabel('$c$')
hold off
sgtitle("Policy functions")

zt = mcsimul(z_grid,p_z,T);

k0 = kss;

kt = zeros(T + 1, 1);
ct = zeros(T, 1);

kt(1) = k0;
[~, idx] = min(abs(k_grid - k0));
kt(1) = k_grid(idx);

for t=1:T
    kt(t + 1) = policy_k(k_grid == kt(t), zt(t,2));
    ct(t) = policy_c(k_grid == kt(t), zt(t,2));
end

figure
subplot(1,3,1)
plot(0:T,kt)
xlabel('$t$')
ylabel('$k_t$')
subplot(1,3,2)
plot(1:T,ct)
xlabel('$t$')
ylabel('$c_t$')
subplot(1,3,3)
plot(1:T,zt(:,1))
xlabel('$t$')
ylabel('$z_t$')
sgtitle("Simulation")

% Better results for time series can be seen with greater T.
% Both capital and consumption randomly oscillate around their steady state.
% For T = 20, we cannot easily observe this behavior.