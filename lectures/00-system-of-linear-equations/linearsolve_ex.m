close all;
clear all;
clc;
%% Direct methods
global xsol
% Setting up the linear system
A = [10 0 1;0.5 7 1;1 0 6];
b = [21;9;8];
% Gauss-Jordan elimination
xj = gaussj(A,b);
% Gauss elimination
xg = gauss(A,b);
xsol = xj;
%% Iterative methods
x0 = zeros(3,1);
[xs,rs,er] = sor(A,x0,b,5,1);
% a glance of convergence
figure(1);
semilogy(1:length(rs), rs, 'linewidth',2);
ylabel('log of (the norm of) residuals');
xlabel('number of iterations');
h = legend('gauss-seidel');
set(h,'fontsize',12,'box','off');
% rate (and order) of convergence
figure(2);
ratio = er(2:length(er)) ./ er(1:length(er)-1);
plot(1:length(er)-1,ratio,'linewidth',2);
ylabel('ratio of (the norm of) errors (around the exact solution)');
xlabel('number of iterations');
h = legend('gauss-seidel');
set(h,'fontsize',12,'box','off');
figure(3);
% Constructing the G matrix (see slides)
G = (eye(length(b))-tril(A)\A);
[V,D] = eig(G);
d = diag(D);
plot(real(d),imag(d),'xb','linewidth',2);
legend('eigs of G');
