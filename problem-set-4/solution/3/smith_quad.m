function[et,v,k1,x] = smith_quad(x0,k0,T,k_grid)
%%
% SMITH_QUAD.M
% Numerical methods class
% Summer 2013
% Written by Gustavo Leyva
% Created on 07.04.13
% Modified on 07.16.13
%
% PURPOSE   Solves the neoclassical growth model using a quadratic policy rule
% USAGE     [et,v,k1,x] = smith_quadratic(x0,k0,T,k_grid)
% INPUTS    x0      : initial values for the shock
%           k0      : initial values for capital
%           t       : length of the horizon
%           k_grid  : capital grid
% OUTPUTS   et      : elapsed time
%           v       : value function evaluated in k_grid
%           k1      : quadratic policy rule
%           x       : optimal parameters of the quadratic policy rule
%
tic;
x = fminunc(@(x) myfun_quad(x,k0,T),x0,optimset('MaxIter',1000,'MaxFunEvals',1000,'LargeScale','off','HessUpdate','bfgs'));
% constructing the value function
ngrid = length(k_grid);
v = zeros(ngrid,1);
for i = 1:ngrid;
    v(i) = -myfun_quad(x,k_grid(i),T);
end;
toc;
et = toc;
k1 = x(1) + x(2)*k_grid + x(3)*k_grid.^2;

end