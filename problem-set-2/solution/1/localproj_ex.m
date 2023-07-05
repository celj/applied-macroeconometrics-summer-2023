%%
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

close all;
%clear all;
clc;
% this script runs after executing simul_proj_ex

global hmax p

data = xlsread('passthrough.xlsx','base','B28:C207');
%p = 1;                  % lag of the local projection model
%hmax = 20;              % months of horizon for the impulse response
for h = 0:hmax
    [proj_coef0,proj_res0] = localproj(data(:,1),data(:,2),h,p);
    n(h+1) = proj_coef0(2);
    [proj_coef1,proj_res1] = localproj(data(:,2),data(:,2),h,p);
    d(h+1) = proj_coef1(2);
end
pass = cumsum(n)./cumsum(d);
plot((0:hmax)',pass)
