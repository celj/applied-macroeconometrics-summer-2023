function [t1,f1] = golden1_ra(xi,x,fx,z,P,hz,q)
%%
% GOLDEN1_RA.M
% Numerical methods class
% Summer 2013
% Written by Gustavo Leyva
% Created on 06.21.13
% Modified on 06.21.13
% 
% PURPOSE   Performs a golden search optimization method for VI with linear
%           interpolation
% USAGE     [t1,f1] = golden1(xi,x,fx,z,P,hz,q) 
% INPUTS    xi  : today's capital state value
%           x   : indices on optimal capital
%           xi  : values on rhs value function
%           z   : shocks
%           P   : Markov transition matrix
%           hz  : index of today's shock
%           q   : price
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

zn = length(z);
pw1 = zeros(1,zn);
pw2 = zeros(1,zn);
for hi = 1:zn;
    pw1(:,hi) = pw_linear(x,fx(:,hi),t1);
    pw2(:,hi) = pw_linear(x,fx(:,hi),t2);
end;    
c1 = xi+z(hz)-t1*q;
c2 = xi+z(hz)-t2*q;
f1 = retrn(c1) + beta*(P(hz,:)*pw1' );
f2 = retrn(c2) + beta*(P(hz,:)*pw2' );

it = 0;
while abs(U-L) > tol*max([1;(abs(t1)+abs(t2))]);

    if f1 > f2;
        U = t2;
        t2 = t1;
        f2 = f1;
        t1 = p*L + (1-p)*U;        
        pw1 = zeros(1,zn);        
        for hi = 1:zn;
            pw1(:,hi) = pw_linear(x,fx(:,hi),t1);            
        end;    
        c1 = xi+z(hz)-t1*q;
        f1 = retrn(c1) + beta*(P(hz,:)*pw1' );   
    else
        L = t1;
        t1 = t2;
        f1 = f2;
        t2 = (1-p)*L + p*U;
        pw2 = zeros(1,zn);
        for hi = 1:zn;
            pw2(:,hi) = pw_linear(x,fx(:,hi),t2);            
        end;    
        c2 = xi+z(hz)-t2*q;
        f2 = retrn(c2) + beta*(P(hz,:)*pw2' );
    end;
it = it+1;

end;

end