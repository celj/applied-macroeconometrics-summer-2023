function [x,nrmHist,error] = sor(A,x,b,nsteps,om) 
%%
% Written by Yousef Saad
%---------------------------------------------------
%%function [x,nrmHist] = sor (A, x, b, nsteps, om) 
%---------------------------------------------------
global xsol
D = diag(diag(A));
E = -tril(A,-1);
F = -triu(A,1);
L = D-om*E;
U = om*F + (1-om)*D;
%%-------------------- record res. norm
res0 = norm(b - A*x);
nrmHist(1) = res0; 
error(1) = norm(x-xsol); 
%%-------------------- iterate 
for i=1: nsteps
    x = L \ (U*x + om*b);
    nrmHist(i+1) = norm(b - A*x);
    error(i+1) = norm(x-xsol);
end

