%%
% PS1_3.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

clear all;
clc;
data = xlsread('LM_grossflows','Hoja1','B3:J49');
[n,k] = size(data);
%%
h = [0;0;100];
for i = 1:n
    %creating system of equations
    a = data(i,4) + data(i,6);
    b = data(i,2);
    c = data(i,8);
    d = data(i,2) + data(i,3);
    e = data(i,4);
    f = data(i,7);
    A = [-b a -c; d -e -f; 1 1 1];
    xj(i,:) = gauss(A,h);
end
%%
urate = xj(:,2)./(xj(:,1)+xj(:,2))*100;
%%
data_avg = sum(data)/n;
for i = 1:n
    %creating system of equations
    a = data(i,4) + data_avg(6);
    b = data_avg(2);
    c = data_avg(8);
    d = data_avg(2) + data_avg(3);
    e = data(i,4);
    f = data_avg(7);
    A = [-b a -c; d -e -f; 1 1 1];
    xj(i,:) = gaussj(A,h);
end
%%
urate_cf = xj(:,2)./(xj(:,1)+xj(:,2))*100;
%%
X = [ones(n,1) urate (1:1:n)'];
Y = urate_cf;
b_ols = (X'*X)\(X'*Y)*100;
plot([urate urate_cf])
h = legend('flows-based urate','counterfactual urate');
ylabel('% of the labor force');
set(h,'fontsize',12,'box','off');

