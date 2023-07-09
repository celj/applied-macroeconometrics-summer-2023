function [fx,dfx] = func0(x)
syms z1 z2

h = 100*(z2-z1^2)^2 + (1-z1)^2;
dh1 = diff(h,z1);
dh2 = diff(h,z2);
f = [dh1;dh2];
df11 =diff(f(1),z1);
df12 =diff(f(1),z2);
df21 =diff(f(2),z1);
df22 =diff(f(2),z2);
df = [df11 df12;df21 df22];
fm = matlabFunction(f);
dfm = matlabFunction(df);
fx = fm(x(1),x(2));
dfx = dfm(x(1),x(2));

end