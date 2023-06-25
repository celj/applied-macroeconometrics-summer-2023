%%
% Bootstrap

close all;
clear all;
clc;

% available sample
beta = [0.5;2];
n = 10;                                                                     % sample size
x = [ones(n,1) randn(n,1)];
u = randn(n,1);
y = x*beta + u;

m = 100;                                                                   % number of Monte Carlo simulations
beta_store = zeros(2,m);

% Computing the OLS estimator
beta_hat = (x'*x)\x'*y;
u_hat = y-x*beta_hat;

% Bootstrap starts...
for i = 1:m
    u_res = simple_bootstrap(u_hat);
    y_res = x*beta_hat+u_res;
    beta_store(:,i) = (x'*x)\x'*y_res;
end;

% estimated density (kernel) for beta(i)
i = 2;
l = 1;
u = 3;
ngrid = 10000;
grid1 = (l:(u-l)/(ngrid-1):u)';

dens = kernel_est(grid1,beta_store(i,:)');
dens = dens/sum(dens);
figure(1);
plot(grid1,dens);

% Bias correction for betas
mean_beta_store = mean(beta_store');
mean_beta_store = mean_beta_store';
bias_bootstrap = mean_beta_store-beta_hat;
beta_corrected = beta_hat - (bias_bootstrap);

% Confidence intervals for beta(2,1)
beta2 = beta_store(2,:);
% Efron
a = 0.05;   % level of significance
% have to sort beta2
beta2 = sort(beta2);
beta2 = beta2';
beta2_l = quantile(beta2,a/2);
beta2_u = quantile(beta2,1-a/2);
disp('The Efron confidence interval is:');
disp([beta2_l beta_hat(2,1) beta2_u]);
% Hall
e = beta_store(2,:)'-beta_hat(2,1);
el = quantile(e,a/2);
eu = quantile(e,1-a/2);
disp('The Hall confidence interval is:');
disp([beta_hat(2,1)-eu beta_hat(2,1) beta_hat(2,1)-el]);



