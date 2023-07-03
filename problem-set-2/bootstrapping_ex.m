tic

data = readmatrix('data/passthrough.xlsx','Sheet','base','Range','B28:C207');

cpi = data(:,1);
exr = data(:,2);

h_max = 20;
p = 1;

simulations = 2000;

projections = 0:h_max;

[exchange_rate, price_level] = deal(zeros(h_max+1,1)');

for h = projections
    exchange_rate(h+1) = localproj(exr,exr,h,p);
    price_level(h+1) = localproj(cpi,exr,h,p);
end

pass_through = cumsum(price_level)./cumsum(exchange_rate);

exchange_rate_store = bootstrap(exr,exr,h_max,p,exchange_rate,simulations);
price_level_store = bootstrap(cpi,exr,h_max,p,price_level,simulations);

pass_through_store = cumsum(price_level_store,2)./cumsum(exchange_rate_store,2);

[lower, upper] = efron_ci(pass_through_store,0.1);

figure
hold on
plot(projections,pass_through,"Color","black")
plot(projections,lower,'--',"Color","black")
plot(projections,upper,'--',"Color","black")
title('Cumulative pass-through elasticity')
xlabel('Months')
hold off

disp('The Efron confidence interval for every h is:');
disp([fix(0:h_max)' lower' pass_through' upper']);

toc