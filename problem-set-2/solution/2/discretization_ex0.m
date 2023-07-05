%%
% PS2_2.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
% 
% Purpose       Evaluates relative performance of three discretization
%               methods
% References
%               1. Tauchen, G. (1986): "Finite State Markov-Chain Approximations to 
%               Univariate and Vector Autoregressions." Economics Letters 20, 177-181
%               2. Tauchen, G. and R. Hussey (1991): "Quadrature-Based Methods for 
%               Obtaining Approximate Solutions to Nonlinear Asset Pricing Models."
%               Econometrica, Vol. 59, No. 2, pp. 371-396
%               3. Rouwenhorst, K.G (1995): "Asset Pricing Implications of Equilibrium 
%               Business Cycle Models." In Frontiers of Business Cycle Research, Chapter 10
%               4. Floden, M. (2008): "A note on the accuracy of Markov-chain approximations 
%               to highly persistent AR(1) processes." Economics Letters 99,
%               516-520

close all;
clear all;
clc;
%%
% AR(1) parameters
rho = 0.99;
sig_u = sqrt(0.020);
sig_z = sqrt(sig_u^2/(1-rho^2));
ngrid = 5;
m = sqrt(ngrid-1);                                                          % to make Tauchen's comparable to Rouwenhorst's method
%% http://www.efunda.com/math/num_integration/findgausshermite.cfm
%% n_grid = 5
nodes = [-2.02018287046;-0.958572464614;0;0.958572464614;2.02018287046];
weights = [0.019953242059;0.393619323152;0.945308720483;0.393619323152;0.019953242059];
%% n_grid = 9
%nodes = [-3.19099320178;-2.26658058453;-1.46855328922;-0.723551018753;0;0.723551018753;1.46855328922;2.26658058453;3.19099320178];
%weights = [3.96069772633E-005;0.00494362427554;0.0884745273944;0.432651559003;0.720235215606;0.432651559003;0.0884745273944;0.00494362427554;3.96069772633E-005];
%% n_grid = 15
%nodes = [-4.49999070731;-3.6699503734;-2.96716692791;-2.32573248617;-1.71999257519;-1.13611558521;-0.565069583256;0;0.565069583256;1.13611558521;1.71999257519;2.32573248617;2.96716692791;3.6699503734;4.49999070731];
%weights = [1.52247580425E-009;1.05911554771E-006;0.000100004441232;0.00277806884291;0.0307800338725;0.158488915796;0.412028687499;0.564100308726;0.412028687499;0.158488915796;0.0307800338725;0.00277806884291;0.000100004441233;1.05911554771E-006;1.52247580425E-009];
%% 
[grid0,p0] = rouwenhorst(rho,sig_u,ngrid);
[grid1,p1] = tauchen(rho,sig_u,m,ngrid);
[grid2,p2] = tauchen_hussey(rho,sig_u,nodes,weights,1);
[grid3,p3] = ac(rho,sig_u,ngrid);
%%
[rhos0,s_grid0] = mkv_moments(grid0,p0);
[rhos1,s_grid1] = mkv_moments(grid1,p1);
[rhos2,s_grid2] = mkv_moments(grid2,p2);
%%
disp('moments of the continuous model');
disp('     rho       sig_z');
disp([rho sig_z]);
disp('Rouwenhorst(1995) discrete approximation moments');
disp('     rho       sig_z');
disp([rhos0 s_grid0]);
disp('Tauchen(1986) discrete approximation moments');
disp('     rho       sig_z');
disp([rhos1 s_grid1]);
disp('Tauchen&Hussey(1991) discrete approximation moments');
disp('     rho       sig_z');
disp([rhos2 s_grid2]);