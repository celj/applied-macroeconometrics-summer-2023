%%
% Steepest Descent and condition number
% On how the condition number of matrix A in Ax=b affects the convergence
% of the Steepest Descent method

close all;
clear all;
clc;

% system
A = hilb(2);
b = [0.1;0.2];
% grid for x
x = (-20:0.1:20)';
n = length(x);
f = zeros(n,n);
% evaluating f(x) = x'Ax-b'x
for i = 1:n;
    for j = 1:n;
        f(i,j) = 0.5*[ x(i) x(j) ]*A*[ x(i);x(j) ]-b'*[ x(i);x(j) ];
    end;
end;

figure(1);
contour(x,x,f);
axis equal
