function [rho,s_grid] = mkv_moments(grid,p)
%%
% MKV_MOMENTS.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.14.13
% Modified on 06.18.13
%
% PURPOSE   Computes the moments of a Markov chain
% USAGE     [rho,s_grid] = mkv_moments(grid,p)
% INPUTS    grid    : grid points
%           p       : Markov matrix
% OUTPUTS   rho     : AR(1) coefficient
%           s_grid  : standard deviation of the discretized approximation of the
%           process
% USES: ergodic.m
%
n = length(grid);
p0 = ergodic(p);
m_grid = p0'*grid;                                                          % mean
s_grid = sqrt(p0'*(grid.^2)-(m_grid)^2);                                    % standard deviation
cv_grid = (p.*(p0*(ones(n,1)'))).*(grid*grid');
cv_grid = sum(sum(cv_grid))-m_grid^2;                                       % covariance
rho = cv_grid/s_grid^2;

end