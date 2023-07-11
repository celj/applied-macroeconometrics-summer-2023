function val=valfun2(x)

global vlast beta delta alpha k0 kt phi theta

k = x(1);
h = x(2);
g = interp1(k0,vlast,k,'linear');
kk = 1.75 * kt ^ alpha * h ^ (1 - alpha) - k + (1 - delta) * kt; % add technology A_t = 1.75 to improve efficiency and convergence

if kk <= .001
    val = log(.001) + phi * log(1 - h) + beta * g + (kk - .001);
else
    val = log(kk) + phi * log(1 - h) + beta * g;
end

val = -theta * val;