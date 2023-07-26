%%
% F_4.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%
clear all
clc;
global beta sig
data = xlsread('F_4','Hoja1','B13:J13');
%--------------------------------------------------------------------------
% setting parameters                 
%--------------------------------------------------------------------------
beta = 0.96;
sig = 2;
ngrid = 500;
k_0 = -2;
k_T = 6;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
q = 1.0338;
r = 1/q;
z_grid = [1.0;0.5;0.1];
PI = reshape(data,3,3);
PI = PI';
zngrid = length(z_grid);
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et1,it1,vf1,p1] = vi_interpol_1_ra(k_grid,q,z_grid,PI);
c1 = k_grid+z_grid(1)-p1(:,1)*q;
c2 = k_grid+z_grid(2)-p1(:,2)*q;
c3 = k_grid+z_grid(3)-p1(:,3)*q;
polk1 = p1;
polc1 = [c1 c2 c3];
figure(1);
subplot(2,2,1);
hold on;
plot(k_grid,k_grid,'k-',k_grid,polk1(:,1),'r.-',k_grid,polk1(:,2),'b.-',k_grid,polk1(:,3),'g.-');
l1 = legend('45°','employed','unemployed','inactive','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy functio for tomorrows capital');
hold off;
subplot(2,2,2);
hold on;
plot(k_grid,vf1(:,1),'r.-',k_grid,vf1(:,2),'b.-',k_grid,vf1(:,3),'g.-');
l1 = legend('employed','unemployed','inactive','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('value function');
hold off;
subplot(2,2,3);
hold on;
plot(k_grid,polc1(:,1),'r.-',k_grid,polc1(:,2),'b.-',k_grid,polc1(:,3),'g.-');
l1 = legend('employed','unemployed','inactive','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('policy function for consumption');
%--------------------------------------------------------------------------
% iterating over the dstribution F(a,e)
%--------------------------------------------------------------------------
% first find the ergodic distribution of the shocks
ergo = ergodic(PI);
%--------------------------------------------------------------------------
% initializing the distribution
%--------------------------------------------------------------------------
ngrid = 1000;
step = (k_T-k_0)/(ngrid-1);
kf_grid = (k_0:step:k_T)';
F0 = zeros(ngrid,zngrid);

for h = 1:zngrid;
    for i = 1:ngrid;
        F0(i,h) = ( (kf_grid(i)-kf_grid(1))/(kf_grid(ngrid)-kf_grid(1)) )*ergo(h);
    end;
end;
% finding in case the policy function being horizontal at initial values in
% k_grid
s = zeros(zngrid,1);
for h = 1:zngrid;
    if abs(polk1(1,h)-k_0) < 10^-6;
        s(h) = finding(k_grid,polk1(:,h));    
    else
        s(h) = -999999;
    end;
end
%--------------------------------------------------------------------------
% computing the inverse of a'(a,e)
%--------------------------------------------------------------------------
af = zeros(ngrid,zngrid);
for h = 1:zngrid;
    for j = 1:ngrid;
        if abs(kf_grid(j)- k_0) < 10^-6 && s(h) ~= -999999;
            af(j,h) = s(h);
        else
            af(j,h) = pw_linear(polk1(:,h),k_grid,kf_grid(j));        
        end;
    end;
end;
%--------------------------------------------------------------------------
% finding the stationary distribution
%--------------------------------------------------------------------------
crit = 1;
tol = 10^-3;
it = 0;
while it < 150;
    it = it+1;
    %disp(it);
    % inteporlating F0 in af
    F0i = zeros(ngrid,zngrid);
    for h = 1:zngrid;
        for j = 1:ngrid;
            if af(j,h) < k_0;
                F0i(j,h) = 0;
            elseif af(j,h) > k_T; 
                F0i(j,h) = ergo(h);
            else
                F0i(j,h) = pw_linear(kf_grid,F0(:,h),af(j,h));
            end;
        end;
    end;
    % computation of F1
    F1i = zeros(ngrid,zngrid);
    for i = 1:ngrid;
        F1i(i,:) = F0i(i,:)*PI;
    end;
    crit = max(abs(F1i(:,1)-F0i(:,1)));
    F0 = F1i;
    %disp(crit);
end;
hold off;
subplot(2,2,4);
hold on;
plot(kf_grid,F0(:,1),'r.-',kf_grid,F0(:,2),'b.-',kf_grid,F0(:,3),'g.-');
l1 = legend('employed','unemployed','inactive','Location','NorthWest');
set(l1,'box','off');
xlabel('capital grid');
ylabel('distribution function');
axis([k_0 k_T 0 1]);
%--------------------------------------------------------------------------
% computing the MC condition for given q
%--------------------------------------------------------------------------
% fix h=1;
count = ngrid-2+1;
store_a = zeros(count,1);
store_e = zeros(zngrid,1);
for h = 1:zngrid;
    for i = 2:ngrid;
        store_a(i,h) = 0.5*( F0(i,h)-F0(i-1,h) )*( kf_grid(i)+kf_grid(i-1) );
    end;
   store_e(h) = sum(store_a(:,h)) + F0(1,h)*kf_grid(1);
end;
sum_mc = sum(store_e);