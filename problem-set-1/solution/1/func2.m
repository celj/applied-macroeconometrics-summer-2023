function [fx,dfx] = func2(x)
syms z
f = atan(z);
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end