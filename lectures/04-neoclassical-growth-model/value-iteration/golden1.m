function [t1,f1] = golden1(xi,x,fx)
%%
% GOLDEN1.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 06.21.13
% Modified on 07.18.13
%
% PURPOSE   Performs a golden search optimization method for VI with linear
%           interpolation
% USAGE     [t1,f1] = golden1(xi,x,fx)
% INPUTS    xi  : today's capital state value
%           x   : indices on optimal capital
%           fx  : values on rhs value function
% OUTPUTS   t1  : solution
%           f1  : value function
%
global beta

L = x(1);
U = x(length(x));
p = 0.5*(sqrt(5)-1);
t1 = p*L + (1-p)*U;
t2 = (1-p)*L + p*U;
tol= 10^-8;
f1 = retrn(xi,t1) + beta*pw_linear(x,fx,t1);
f2 = retrn(xi,t2) + beta*pw_linear(x,fx,t2);

it = 0;
while abs(U-L) > tol*max([1;(abs(t1)+abs(t2))]);
    
    if f1 > f2;
        U = t2;
        t2 = t1;
        f2 = f1;
        t1 = p*L + (1-p)*U;
        f1 = retrn(xi,t1) + beta*pw_linear(x,fx,t1);
    else
        L = t1;
        t1 = t2;
        f1 = f2;
        t2 = (1-p)*L + p*U;
        f2 = retrn(xi,t2) + beta*pw_linear(x,fx,t2);
    end;
    it = it+1;
    
end;

end