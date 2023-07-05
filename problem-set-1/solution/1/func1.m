function [fx,dfx] = func1(x)
syms z
f = cos(z)-z;
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end