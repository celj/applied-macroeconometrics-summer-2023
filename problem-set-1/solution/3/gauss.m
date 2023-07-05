function x = gauss(A,b)
%%
% Written by Yousef Saad
%---------------------------------------------------
% function x = gauss (A, b)
% solves A x = b by Gauss elimination
%---------------------------------------------------
n = size(A,1);
A = [A,b];

for k = 1:n-1
    for i = k+1:n
        piv = A(i,k) / A(k,k);
        A(i,k+1:n+1) = A(i,k+1:n+1) - piv*A(k,k+1:n+1);
    end
end

x = backsolv(A,A(:,n+1));