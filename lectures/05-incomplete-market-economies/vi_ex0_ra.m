%%
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 10.27.11
% Modified on 07.14.13
%
% Notes: with input from Heer & Maussner (2005)
%

clear all
clc;
global beta sig
%--------------------------------------------------------------------------
% setting parameters
%--------------------------------------------------------------------------
beta = 0.96;
sig = 1;
ngrid = 100;
k_0 = -2;
k_T = 4;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';
q = 1.0129;
%--------------------------------------------------------------------------
% discretizing the stochastic process
%--------------------------------------------------------------------------
zngrid = 2;
z_grid = [1;0.1];
PI = [0.925 0.075;0.5 0.5];
%--------------------------------------------------------------------------
% value function implementation
%--------------------------------------------------------------------------
[et0,it0,vf0,p0] = vi_improved_0_ra(k_grid,q);
[et1,it1,vf1,p1] = vi_interpol_1_ra(k_grid,q);
polk0 = k_grid(p0);
p1f = p1;                                                                   % check this!!!
%p1f(1:3,2) = -2;
polk = p1f;
c1 = k_grid+z_grid(1)-polk0(:,1)*q;
c2 = k_grid+z_grid(2)-polk0(:,2)*q;
polc0 = [c1 c2];
figure(1);
plot(k_grid,k_grid,'k-',k_grid,polk0(:,1),'r.-',k_grid,polk0(:,2),'b.-');
axis([k_0 5 k_0 5]);
figure(2);
plot(k_grid,vf0);
figure(3);
plot(k_grid,polc0);
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
    if abs(polk(1,h)-k_0) < 10^-6;
        s(h) = finding(k_grid,polk(:,h));
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
            af(j,h) = pw_linear(polk(:,h),k_grid,kf_grid(j));
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
plot(kf_grid,F0(:,1),kf_grid,F0(:,2))
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