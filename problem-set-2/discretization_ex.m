%%
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.14.13
% Modified on 07.17.13
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

% Updated by Carlos Lezama, 03.07.2023

% AR(1) parameters
rho = 0.98;
sig_u = sqrt(0.020);
sig_z = sqrt(sig_u^2/(1-rho^2));
ngrid = 5;
m = sqrt(ngrid-1); % to make Tauchen's comparable to Rouwenhorst's method

nodes = [-2.02018287046;-0.958572464614;0;0.958572464614;2.02018287046];
weights = [0.019953242059;0.393619323152;0.945308720483;0.393619323152;0.019953242059];

[grid0,p0] = rouwenhorst(rho,sig_u,ngrid);
[grid1,p1] = tauchen(rho,sig_u,m,ngrid);
[grid2,p2] = tauchen_hussey(rho,sig_u,nodes,weights,1);
[grid3,p3] = adda_cooper(rho,sig_u,ngrid);

[rhos0,s_grid0] = mkv_moments(grid0,p0);
[rhos1,s_grid1] = mkv_moments(grid1,p1);
[rhos2,s_grid2] = mkv_moments(grid2,p2);
[rhos3,s_grid3] = mkv_moments(grid3,p3);

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
disp('Adda&Cooper(2003) discrete approximation moments');
disp('     rho       sig_z');
disp([rhos3 s_grid3]);