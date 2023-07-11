global A alpha beta delta eta grid_size k0 rho T

T = 1e4;

grid_size = 250;

A = 1; % total factor productivity
alpha = 0.3; % capital share
beta = 0.95; % discount factor
delta = 1; % depreciation rate
eta = 1; % relative risk aversion
k0 = 1.0120; % initial capital stock

rho = (1 / beta) - 1;

kss = ((alpha * A) / (rho + delta)) ^ (1 / (1 - alpha));
css = prodfunc(kss) - delta * kss;

kmin = 0.5 * kss;
kmax = 1.5 * kss;

k_grid = linspace(kmin, kmax, grid_size)';

[it,vf,p] = vi(k_grid);

policy = k_grid(p);

figure;
hold on;
subplot(1,2,1);
plot(k_grid,policy);
xlabel('$k$');
ylabel('Policy function for the capital of tomorrow');
subplot(1,2,2);
plot(k_grid,vf);
xlabel('$k$');
ylabel('Value function');
hold off;

k = zeros(T,1);
c = zeros(T,1);

k(1) = k0;
c(1) = prodfunc(k(1)) - delta * k(1);

for t = 2:T
    k(t) = pw_linear(k_grid,policy,k(t - 1));
    c(t) = prodfunc(k(t)) - delta * k(t);
end

figure;
hold on;

subplot(1,2,1);
semilogx(1:T,k);
yline(kss,'k--');
legend('','Steady state');
title('Capital');
xlabel('$t$');
ylabel('$k_t$');

subplot(1,2,2);
semilogx(1:T,c);
yline(css,'k--');
legend('','Steady state');
title('Consumption');
xlabel('$t$');
ylabel('$c_t$');

hold off;

c_positive = c(c > 0);

v0 = pw_linear(k_grid,vf,k0); % v(k0)
c_agg = sum(beta .^ (1:length(c_positive))' .* log(c_positive)); % ignores the first negative consumption

fprintf('The initial value of the value function is %f.\n',v0);
fprintf('The aggregate consumption under the optimal plan is %f.\n',c_agg);
fprintf('The difference between them is %f.\n',abs(v0 - c_agg));
% which is very close to zero, as expected