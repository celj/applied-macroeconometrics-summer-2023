function [x,nrmHist] = steep(A,x,b,nsteps,tol)
%%
% Written by Yousef Saad
%%    function [x,nrmHist] = steep(A, x, b, nsteps, tol)
% steepest descent
% function [x,r] = steep(A, x, b, nsteps)
% A = matrix =
% x = initial guess on entry -- app sol on return
% b = rhs,
% nsteps = number of steps --
r = b - A * x;
nrmHist(1) = norm(r);
if (nargin<5);
    tol = eps;
end
tol1 = tol*nrmHist(1);
%%-------------------- loop
for i = 1:nsteps
    ar = A * r;
    alp = r'*r / ( r' * ar);
    x = x + alp * r;
    r = r - alp* ar;
    ro = norm(r);
    %%--------------------  save norms for plotting
    nrmHist(i+1) = ro;
    if (ro < tol1), break, end
end
