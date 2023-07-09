tic;

setup;

tol = 10^-5;
crit = 100;

x0 = [1.2, 1.2];
x1 = [-1.2, 1];

[f, dfx] = func0(x0);

[it0,sol0,err0] = newton(x0,tol,crit);
[it1,sol1,err1] = newton(x1,tol,crit);

[it2,sol2,err2] = chord(x0,tol,crit);
[it3,sol3,err3] = chord(x1,tol,crit); % does not converge

figure;
semilogy(it0, err0, it1, err1);
grid on;
legend('Newton $x_0$', 'Newton $x_1$');
title('$x_0 = (1.2, 1.2)$, $x_1 = (-1.2, 1)$');
xlabel('Iterations');
ylabel('Error');

figure;
semilogy(it2, err2);
grid on;
legend('Chord $x_0$');
title('$x_0 = (1.2, 1.2)$');
xlabel('Iterations');
ylabel('Error');

figure;
hold on;
plot(sol0(:,1), sol0(:,2));
plot(sol1(:,1), sol1(:,2));
plot(sol2(:,1), sol2(:,2));
plot(1,1,'k.','MarkerSize',30);
grid on;
legend('Newton $x_0$', 'Newton $x_1$', 'Chord $x_0$','True solution');
title('$x_0 = (1.2, 1.2)$, $x_1 = (-1.2, 1)$');
xlabel('$x_1$');
ylabel('$x_2$');
hold off;

ngrid = 50;
rho = 0.9;
sig_u = sqrt(0.06);
sig_z = sqrt(sig_u^2 / (1 - rho^2));

[grid0,p0] = rouwenhorst(rho,sig_u,ngrid);

% sims = 250;

% initial_state = 1;

toc;