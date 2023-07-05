function [fx,dfx] = func3(x)
syms z
f = sin(z);
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end