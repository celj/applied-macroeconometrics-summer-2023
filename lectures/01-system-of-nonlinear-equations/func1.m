function [fx,dfx] = func1(x)
syms z
f = (z^1/3)*exp(-z^2);
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end
