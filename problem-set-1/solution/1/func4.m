function [fx,dfx] = func4(x)
syms z
f = z.^2;
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end