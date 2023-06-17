function A = lap2D(nx,ny)
%%
% Written by Yousef Saad
% function A = lap2D(nx,ny)
%--- generates a 2-D laplacean -----------------
%%
tx = sptridiag(-1, 2, -1, nx) ;
ty = sptridiag(-1, 2, -1,  ny);
A  = kron(speye(ny,ny),tx) + kron(ty,speye(nx,nx));
