%%
% PS2_1.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

close all;
clear all;
clc;
global uu ll mm hmax p
data = xlsread('passthrough.xlsx','base','B28:C207');
p = 1;                  % lag of the local projection model
hmax = 20;              % months of horizon for the impulse response
y = data(:,1);
x = data(:,2);
m = 2000;               % bootstrap replications
for j = 1:m
    for h = 0:hmax
        y_b = simul_proj(data(:,1),data(:,2),h,p);
        x_b = simul_proj(data(:,2),data(:,2),h,p);
        logy = log(y);
        logx = log(x);
        T = length(logx);
        for i=13:T-h
            d12y(i)   = logy(i)  -logy(i-12);
            d12x(i)   = logx(i)  -logx(i-12);
        end
        t = 13;
        g = [d12y(t:T-h)' d12x(t:T-h)'];
        n = length(g);
        X1 = [ones(n-p,1) g(p+1:n,2) g(1:n-p,1) g(1:n-p,2)];
        X2 = [ones(n-p,1) g(p+1:n,2) g(1:n-p,2)];
        beta_h_1 = (X1'*X1)\(X1'*y_b);
        beta_h_2 = (X2'*X2)\(X2'*x_b);
        n_b(h+1) = beta_h_1(2,:);
        d_b(h+1) = beta_h_2(2,:);
    end
    pass_b(:,j) = cumsum(n_b')./cumsum(d_b');
end
uu = quantile(sort(pass_b'),0.95);
ll = quantile(sort(pass_b'),0.05);
mm = quantile(sort(pass_b'),0.5);
run localproj_ex.m
plot((0:hmax)',[ll' pass' mm' uu'])