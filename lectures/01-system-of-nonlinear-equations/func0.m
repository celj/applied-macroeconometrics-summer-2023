function [fx,dfx] = func0(x)
syms z
f = 27*z^4+162*z^3-180*z^2+62*z-7;
df = diff(f,z);
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x);
dfx = dfm(x);

end
