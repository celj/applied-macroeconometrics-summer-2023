function [iter_store,sol_store,funct_store] = chord(x0)

crit = 100;
tol = 10^-5;
i = 0;
sol_store(1,:) = x0';
iter_store(1) = i;
funct_store(1) = norm(func0(x0),'inf');
[f0,df0] = func0(x0);
while crit > tol
    i = i+1;
    iter_store(i+1) = i;
    [fx,dfx] = func0(x0);
    A = df0;
    b = fx;
    dirtn = -gaussj(A,b);
    x1 = x0+dirtn;
    crit = norm(x1-x0);
    sol_store(i+1,:) = x1';
    funct_store(i+1) = norm(func0(x1),'inf')./funct_store(1);
    x0 = x1;
end
end