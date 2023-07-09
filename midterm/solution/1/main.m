%%
% M_1.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%
clear all
close all
clc

x0 = [1.2;1.2];
%x0 = [-1.2;1];
[iter_store,sol_store,funct_store] = newton(x0);
[iter_store_,sol_store_,funct_store_] = chord(x0);

figure(1); 
semilogy(iter_store, funct_store, 'linewidth',2);
hold on
semilogy(iter_store_, funct_store_, 'linewidth',2);
ylabel('log of (the norm of) residuals');
xlabel('number of iterations');
h=legend('newton','chord');
set(h,'fontsize',12,'box','off');