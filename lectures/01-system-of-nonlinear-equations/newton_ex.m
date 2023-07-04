%%
% NEWTON_EX.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% Created in July, 2013
%

close all;
clear all;
clc;

tol = 10^-5;
x0 = 0;
[iter_store,sol_store,error_store] = newton(x0);
disp([iter_store sol_store error_store [nan;error_store(2:length(error_store))./error_store(1:length(error_store)-1)]]);

% when using func1.m
%x=(-5:0.1:5)';
%plot(x,x.^1/3.*exp(-x.^2))
