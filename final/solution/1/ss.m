function F = ss(z)

global alpha sh beta delta
k = z(1);
n = z(2);
phi = (1-sh)/sh;
y = (k.^alpha)*(n.^(1-alpha));
F(1) = 1/beta - alpha*(k.^(alpha-1))*(n.^(1-alpha)) - (1-delta);
F(2) = (1-alpha).*y./n - (1-alpha+phi).*y + phi*k - phi*(1-delta)*k;

end