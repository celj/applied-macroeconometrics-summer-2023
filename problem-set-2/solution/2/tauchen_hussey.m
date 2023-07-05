function [y_grid,p] = tauchen_hussey(rho,sig_u,abcs,wegs,a)
%%
% TAUCHEN_HUSSEY.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.14.13
% Modified on 07.16.13
%
% PURPOSE   Discretize an AR(1) of the form z = rho*z(-1)+u
% USAGE     [y_grid,p] = tauchen_hussey(rho,sig_u,m,ngrid,a)
% INPUTS    rho     : first-order serial autocorrelation (1x1)
%           sig_u   : standard deviation of u (1x1)
%           abcs    : Gauss-Hermite abscissas (kx1)     
%           wegs    : Gauss-Hermite weights (kx1)     
%           a       : weight in linear combination between conditional and 
%                     unconditional dev. est. a*sig_u+(1-a)*sig_y. When
%                     a fifth argument is not provided, the code accomodates 
%                     the recommendation of Floden(2008)
% OUTPUTS   y_grid  : discretized values
%           p       : Markov transition matrix
% REFERENCES
%           1.Tauchen, G. and R. Hussey (1991): "Quadrature-Based Methods for 
%           Obtaining Approximate Solutions to Nonlinear Asset Pricing Models."
%           Econometrica, Vol. 59, No. 2, pp. 371-396.
%           2.Floden, M. (2008): "A note on the accuracy of Markov-chain approximations 
%           to highly persistent AR(1) processes." Economics Letters 99,
%           516-520.
%         
ngrid = length(abcs);
if nargin < 5;
    a = 0.5+0.25*rho;
end;
sig_y = sqrt(sig_u^2/(1-rho^2));
sig = (a*sig_u+(1-a)*sig_y);
y_grid = sqrt(2)*sig*abcs;
w = wegs/sqrt(pi);
p = zeros(ngrid,ngrid);

for i = 1:ngrid;
    for j = 1:ngrid;
        p(i,j) = ( pdf('norm',(y_grid(j)-rho*y_grid(i))/sig_u,0,1)/pdf('norm',(y_grid(j)-0)/sig,0,1) )*w(j);        
    end
    s = sum(p(i,:));
    p(i,:) = p(i,:)/s;
end

end