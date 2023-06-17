%%
% BISECTION_EX.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% Created in July, 2013
%

close all;
clear all;
clc;

xl = 0;
xr = 1;
f = '27*x^4+162*x^3-180*x^2+62*x-7';
[iter_store,sol_store,error_store] = bisection(xl,xr,f);
disp([iter_store sol_store error_store [nan;error_store(2:length(error_store))./error_store(1:length(error_store)-1)]]);

% figure
x = 0:0.01:1;
y = 27*x.^4+162*x.^3-180*x.^2+62*x-7;
figure(1);
plot(x,y);
