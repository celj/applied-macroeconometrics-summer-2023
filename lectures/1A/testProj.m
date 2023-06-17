%%
% Written by Yousef Saad

close all;
clear all;
clc;

A = lap2D(30,30);
n = size(A,1);
b = A*ones(n,1);
x = randn(n,1);
nits = 300;
tol = 1.e-06;
[x1,r1] = steep(A,x,b,nits,tol);
[x2,r2] = mr(A,x,b,nits,tol);

semilogy([1:length(r1)], r1, 'linewidth',2)
hold on
semilogy([1:length(r2)], r2,'r--', 'linewidth',2)
h=legend('Steep. desc.','min. res.');
%%
set(h,'fontsize',16)
