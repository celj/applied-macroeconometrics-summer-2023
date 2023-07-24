% T = 1e4;
alpha = 0.3;
beta = 0.95;
delta = 1;
k_max = 5;
k_min = 0.01;
p = 250;
q = 3;
rho = 0.95;
sig_e = 0.01;

[grid_z,p_z] = rouwenhorst(rho,sig_e,q);

grid_k = linspace(k_min,k_max,p)';

% ergodic_z = ergodic(p_z);

% F = zeros(p,q);

% for i = 1:q
%     for j = 1:p
%         F(j,i) = ((grid_k(j) - grid_k(1)) / (grid_k(p) - grid_k(1))) * ergodic_z(i);
%     end
% end

% F = cumsum(F) ./ sum(F);

u = zeros(p,q,p);

for i=1:p % today's capital
    for j=1:q % today's productivity
        for l=1:p % tomorrow's capital
            c = max((exp(grid_z(j)) * (grid_k(i) ^ alpha)) + ((1 - delta) * grid_k(i)) - grid_k(l),1e-200);
            u(i,j,l) = log(c);
        end
    end
end

v0 = zeros(p,q);
v = zeros(p,q);

crit = 1e-5;
termin = 0;

policy = zeros(p,q); % policy function

while (termin == 0)
    for i=1:p % today's capital
        for j=1:q % today's productivity
            aux = zeros(p,1);
            for l=1:p % tomorrow's capital
                aux(l) = u(i,j,l);
                for m=1:q % tomorrow's productivity
                    aux(l) = aux(l) + beta * p_z(j,m) * u(l,m,j);
                end
            end
            [v(i,j),policy(i,j)] = max(aux);
        end
    end
    if norm(v0 - v) < crit
        termin = 1;
    end
    v0 = v;
end

k = zeros(p,q);
c = zeros(p,q);

for i=1:p % today's capital
    for j=1:q % today's productivity
        k(i,j) = grid_k(policy(i,j));
        c(i,j) = (exp(grid_z(j)) * (grid_k(i) ^ alpha)) + ((1 - delta) * grid_k(i)) - k(i,j);
    end
end

figure
subplot(1,3,1)
hold on
for i=1:q
    plot(grid_k,v(:,i),'DisplayName',sprintf('$z = %0.4f$',grid_z(i)))
end
legend
xlabel('$k$')
ylabel('$v(k)$')
hold off
subplot(1,3,2)
hold on
for i=1:q
    plot(grid_k,k(:,i),'DisplayName',sprintf('$z = %0.4f$',grid_z(i)))
end
legend
xlabel('$k$')
ylabel('$k^\prime$')
hold off
subplot(1,3,3)
hold on
for i=1:q
    plot(grid_k,c(:,i),'DisplayName',sprintf('$z = %0.4f$',grid_z(i)))
end
legend
xlabel('$k$')
ylabel('$c$')
hold off
sgtitle("Policy functions")
