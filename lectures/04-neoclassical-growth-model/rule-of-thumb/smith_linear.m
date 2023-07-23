function [et,v,k1,x] = smith_linear(x0,k0,t,k_grid)
%%
% SMITH_LINEAR.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.04.13
% Modified on 07.16.13
%
% PURPOSE   Solves the neoclassical growth model using a linear policy rule
% USAGE     [et,v,k1,x] = smith_linear(x0,k0,T,k_grid)
% INPUTS    x0      : initial guess for the parameters
%           k0      : initial values for capital
%           t       : length of the horizon
%           k_grid  : capital grid
% OUTPUTS   et      : elapsed time
%           v       : value function evaluated in k_grid
%           k1      : linear policy rule
%           x       : optimal parameters of the linear policy rule
%
tic;
x = fminunc(@(x) myfun_linear(x,k0,t),x0,optimset('MaxIter',1000,'MaxFunEvals',1000,'LargeScale','off','HessUpdate','bfgs'));
% constructing the value function
ngrid = length(k_grid);
v = zeros(ngrid,1);
for i = 1:ngrid;
    v(i) = -myfun_linear(x,k_grid(i),t);
end;
k1 = x(1) + x(2)*k_grid;
toc;
et = toc;

end