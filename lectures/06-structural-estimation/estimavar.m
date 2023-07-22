function [B,GAM,SUUB,U] = estimavar(data,p)
%%
% ESTIMAVAR.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 09.12.11
% Modified on 06.25.12
%
% PURPOSE   Computes a VAR by OLS
% USAGE     [B,GAM,SUUB,U] = estimavar(data,p)
%
% INPUTS    data : data set '[y1 y2 y3]' (Txk)
%           p    : VAR lag order (maximum = 12)
%           B    : matrix of VAR estimated coefficients (by OLS) (kx(3*p+1)).
%                  The format B for VAR(data,1) is as follows:
%                  y1 = cons1 a1*y1(-1) a2*y2(-1) a3*y3(-1) ...
%                  y2 = cons2 b1*y1(-1) b2*y2(-1) b3*y3(-1) ...
%                  y3 = cons3 c1*y1(-1) c2*y2(-1) c3*y3(-1) ...
%           GAM  : matrix X'X [(k+1)x(k+1)]
%           SUUB : the variance-covariance matrix of the residuals. This matrix is
%                  degrees-of-freedom adjuted. If want to work with the unadjusted
%                  version of that matrix choose SUB (kxk)
%           U    : the vector of estimated residuals [(T-p)xk]
% REFERENCES
%           1.L�tkepohl, H. (2005). "New Introduction to Multiple Time
%           Series Analysis." Springer.
%
[T,k] = size(data);
T = T-p;
Y = data';
cY = size(Y,2);
YY = Y(:,p+1:cY);
Z = zeros(k*p+1,T);

for i = 1:T;
    if     p == 1;
        Z(:,i) = [1;Y(:,i+p-1)];
    elseif p == 2;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2)];
    elseif p == 3;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3)];
    elseif p == 4;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4)];
    elseif p == 5;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5)];
    elseif p == 6;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6)];
    elseif p == 7;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7)];
    elseif p == 8;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8)];
    elseif p == 9;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9)];
    elseif p == 10;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10)];
    elseif p == 11;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);Y(:,i+p-11)];
    elseif p == 12;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);Y(:,i+p-11);Y(:,i+p-12)];
    end;
end;

B = YY*Z'/(Z*Z');
U = YY-(YY*Z'/(Z*Z'))*Z; U = U';
SUB = (1/T)*YY*(eye(T)-Z'/(Z*Z')*Z)*YY';
SUUB = (T/(T-k*p-1))*SUB;                                                   %bias correction
GAM = Z*Z';