function [x,nrmHist] = mr(A,x,b,nsteps,tol)
%%
% Written by Yousef Saad
%% function [x,nrmHist] = mr(A,x,b, nsteps,tol)
% minimal residual method
% function [x,r] = mr(A,x,b, nsteps)
% A = matrix =
% x = initial guess on entry -- app sol on return
% b = rhs, m = number of steps --
%
r = b - A * x;
nrmHist(1) = norm(r);
if (nargin<5);
    tol = eps;
end
tol1 = tol*nrmHist(1);
%%
for i = 1:nsteps
    ar  = A * r;
    alp = (ar'*r) / ( ar' * ar);
    x   = x + alp * r;
    r   = r - alp* ar;
    ro = norm(r);
    nrmHist(i+1) = ro;
    if (ro < tol1), break, end
end
