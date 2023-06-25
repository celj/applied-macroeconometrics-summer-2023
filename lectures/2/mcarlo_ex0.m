%%
% Monte Carlo simulation

close all;
clear all;
clc;

% true model: y = beta*x=u
beta = [0.5;2];
n = 100;                                                                    % sample size
x = [ones(n,1) randn(n,1)];

m = 1000;                                                                   % number of Monte Carlo simulations
beta_store = zeros(2,m);

for i = 1:m
    u = sqrt(0.01)*randn(n,1);
    y = x*beta+u;
    beta_store(:,i) = (x'*x)\x'*y;
end;

% estimated density (kernel) for beta(i)
i = 2;
l = 1.9;
u = 2.1;
ngrid = 10000;
grid1 = (l:(u-l)/(ngrid-1):u)';

dens = kernel_est(grid1,beta_store(i,:)');
dens = dens/sum(dens);
figure(1);
plot(grid1,dens);