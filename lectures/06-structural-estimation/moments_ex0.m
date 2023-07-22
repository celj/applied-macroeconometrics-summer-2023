%%
% MOMENTS_EX0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.15.13
% Modified on 07.22.13
%
% PURPOSE   Computes three moments estimators: GMM, EMM, and II
%
close all;
clear all;
clc;

global data nq r N
%--------------------------------------------------------------------------
% Generating artificial time series for endowment and prices - 'the data'
%--------------------------------------------------------------------------
N = 5000;                                                                   %size of the simulation series
rng(250944);                                                                % seed
rng(1);
beta = 0.97;
gama = 1.3;
% Discretizing the (growth of) the endowment process
zngrid = 10;                                                                % number of states
lambda = 0.178;
sig_u = sqrt(0.0144);
[z_grid,PI] = rouwenhorst(lambda,sig_u,zngrid);
z_grid = z_grid + (0.013)/(1-lambda);
% simulation of the (growth of) the endowment process
T = 100;                                                                    % size of the sample
mcstore = mcsimul(z_grid,PI,T);
y = mcstore(:,1);                                                           %y_t+1
y_index = mcstore(:,2);
% solving for the price q in the Euler equation
q = zeros(T,1);
for i = 1:T;
    q(i) = PI(y_index(i),:)*(beta*exp(z_grid).^-gama);
end;
r = 1./q;                                                                   %r_t+1
%% Generating the final data (after lags) including instruments
nq = 3;                                                                     % number of instruments
data = [(exp(y(2:T))) 1./q(2:T) ones(T-1,1) exp(y(1:T-1)) 1./q(1:T-1)];
%%
%--------------------------------------------------------------------------
% Estimations
%--------------------------------------------------------------------------
%% GMM - 2step
rho0 = [beta;gama];

w0 = eye(3);
rho0 = fminunc(@(rho) quadfunc(rho,w0),rho0,optimset('MaxIter',1000,'MaxFunEvals',1000,'LargeScale','off','HessUpdate','bfgs','Display','Off'));
w0 = inv(varcov(rho0));
rho0 = fminunc(@(rho) quadfunc(rho,w0),rho0,optimset('MaxIter',1000,'MaxFunEvals',1000,'LargeScale','off','HessUpdate','bfgs','Display','Off'));
disp('GMM')
disp('')
disp(rho0);
%% Indirect Inference (score as moments conditions - Gallant and Tauchen)
% Constructing the weighting matrix using the auxiliary model based on
% actual data
[theta,~,sig2,u] = estimavar(r,1);
theta = [theta';sig2];
x = [ones(T-1,1) r(1:T-1,:)];
mm = zeros(T-1,nq);
mm(:,1) = ((1/theta(3))*(x(:,1).*u));
mm(:,2) = ((1/theta(3))*(x(:,2).*u));
mm(:,3) = ((-0.5/theta(3))+((0.5/(theta(3)^2))*(u.^2)));
vm = mm'*mm/(T-1);

rho0_gt = [beta;gama];
w0_gt = inv(vm);
rho0_gt = fminunc(@(rho) quadfunc_gt(rho,w0_gt),rho0_gt,optimset('MaxIter',1000,'MaxFunEvals',1000,'LargeScale','off','HessUpdate','bfgs','Display','Off'));
disp('II GT')
disp('')
disp(rho0_gt);
%% Indirect Inference (param distance as moments conditions - Smith)
% Constructing the weighting matrix using the auxiliary model based on
% actual data
rho0_s = [beta;gama];
w0_s = inv(vm);
rho0_s = fminunc(@(rho) quadfunc_smith(rho,w0_s),rho0_s,optimset('MaxIter',1000,'MaxFunEvals',1000,'LargeScale','off','HessUpdate','bfgs','Display','Off'));
disp('II Smith')
disp('')
disp(rho0_s);











