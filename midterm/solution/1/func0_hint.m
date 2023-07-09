function [fx,dfx] = func0(x)
syms z1 z2

h = ... + (1-z1)^2;
dh1 = diff(h,z1);
dh2 = diff(h,...);
f = [...;dh2];
df11 =diff(f(1),z1);
df12 =diff(f(1),...);
df21 =diff(...,...);
df22 =diff(f(2),z2);
df = [df11 df12;... df22];
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x(1),...);
dfx = dfm(...,x(2));

end