function [fx,dfx] = func5(x)
syms z
f = z.^2 + 1;
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end