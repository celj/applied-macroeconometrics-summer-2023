%%
% ROFT_EX0.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.03.13
% Modified on 07.04.13
%
% PURPOSE   Solves the deterministic neoclassical growth model using rules
%           of thumb as in Smith(1990) and Smith(1991).
%
global alpha beta delta sig fac
%--------------------------------------------------------------------------
% setting parameters
%--------------------------------------------------------------------------
alpha = 0.3;
beta = 0.95;
sig = 1;
delta = 1;
ngrid = 100;
fac = 1;
rho = (1/beta)-1;
k_ss = (fac*alpha/(rho+delta))^(1/(1-alpha));
k_0 = 0.01;
k_T = 3;
step = (k_T-k_0)/(ngrid-1);
k_grid = (k_0:step:k_T)';

% analytical solution
kstore = (alpha*beta).*k_grid.^alpha;
a = 1/(1-beta)*( log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta) );
b = alpha/(1-alpha*beta);
vstore = a+b*log(k_grid);
% using the value function iteration
[et,it,vf,k] = vi_improved_0(k_grid);
%--------------------------------------------------------------------------
% starting the parametric approximation
%--------------------------------------------------------------------------
k0 = k_ss;                                                                  % initial state value
T = 800;                                                                    % approximate infinite horizon
% using a linear policy rule
x0 = [0.01;0.1;0.1];
[~,v0_q,k0_q,~] = smith_quadratic(x0,k0,T,k_grid);
x0 = [0.01;0.1];
[et0,v0_l,k0_l,x] = smith_linear(x0,k0,T,k_grid);

c0_l = prodfunc(k0_l(1:end-1)) - k0_l(2:end);
c0_q = prodfunc(k0_q(1:end-1)) - k0_q(2:end);
cstore = prodfunc(kstore(1:end-1)) - kstore(2:end);

%--------------------------------------------------------------------------
% figures
%--------------------------------------------------------------------------
d = 5;
figure;
subplot(1,3,1);
hold on;
plot(k_grid(1:d:end),v0_l(1:d:end),'go',k_grid(1:d:end),v0_q(1:d:end),'bo',k_grid(1:d:end),vstore(1:d:end),'k+',k_grid(1:d:end),vf(1:d:end),'r+')
plot(k_grid,v0_l,'g--',k_grid,v0_q,'b--',k_grid,vstore,'k-',k_grid,vf,'r-');
xlabel('capital grid');
ylabel('value function');
l1 = legend('linear','quadratic','exact','VI');
set(l1,'box','off');
hold off;
kpol = k_grid(k);
subplot(1,3,2);
hold on;
plot(k_grid(1:d:end),k0_l(1:d:end),'go',k_grid(1:d:end),k0_q(1:d:end),'bo',k_grid(1:d:end),kstore(1:d:end),'k+',k_grid(1:d:end),kpol(1:d:end),'r+')
plot(k_grid,k0_l,'g--',k_grid,k0_q,'b--',k_grid,kstore,'k-',k_grid,kpol,'r-');
xlabel('capital grid');
ylabel('policy function for tomorrow''s capital');
l2 = legend('linear','quadratic','exact','VI');
set(l2,'box','off');
hold off;
subplot(1,3,3);
hold on;
plot(k_grid(1:d:end-1),c0_l(1:d:end),'go',k_grid(1:d:end-1),c0_q(1:d:end),'bo',k_grid(1:d:end-1),cstore(1:d:end),'k+')
plot(k_grid(1:end-1),c0_l,'g--',k_grid(1:end-1),c0_q,'b--',k_grid(1:end-1),cstore,'k-');
xlabel('capital grid');
ylabel('policy function for consumption');
l3 = legend('linear','quadratic','exact');
set(l3,'box','off');
hold off;