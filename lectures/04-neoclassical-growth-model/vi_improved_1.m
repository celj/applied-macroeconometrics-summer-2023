global vlast beta delta theta k0 kt alpha grid_size phi rho

grid_size = 250;

alpha = 0.3; % capital share
beta = 0.95; % discount factor
delta = 1; % depreciation rate
theta = 0.6; % utility weight on leisure

rho = (1 / beta) - 1;

phi = (1 - theta) / theta;

h=.5*ones(1,grid_size);
k0=linspace(0.2, 10, grid_size);
kt1=k0;
vlast=ones(1,grid_size);

xmax =[9.99 .99];
xmin =[.21 .01];

options = optimset('Display','off','LargeScale','off');

crit = 1;
tol = 1e2;

while crit > tol * (1 - beta)
    for j = 1:grid_size
        kt = k0(j);
        z0 = [kt,h(j)];
        z = fmincon(@valfun2,z0,[],[],[],[],xmin,xmax,[],options);
        v(j) = -valfun2(z);
        kt1(j) = z(1);
        h(j) = z(2);
    end
    crit = max(abs(vlast - v));
    vlast = v;
end

figure;
plot(k0,vlast);
title('Model with labor choice');
xlabel('$k_t$');
ylabel('Value function');