tol = 1e-6;
kmax = 1e4;

data = readtable('data/LM_grossflows.xls');
data(any(ismissing(data), 2), :) = [];

data.a = data.flowsEU;
data.b = data.flowsUE + data.flowsUO;
data.c = data.flowsOU;
data.d = data.flowsEU + data.flowsEO;
data.e = data.flowsUE;
data.f = data.flowsOE;

k = 100;
b = [0; 0; k];
quarters = size(data, 1);

results = array2table(zeros(quarters, size(b, 1) + 1));
results.Properties.VariableNames = {'quarter','employed','unemployed','out'};
results.quarter = string(results.quarter);

for q = 1:quarters
    row = data(q, :);
    A = [-row.a row.b -row.c; row.d -row.e -row.f; 1 1 1];
    x = linsolve(A, b);
    results(q, 1).quarter = row.date;
    results(q, 2).employed = x(1);
    results(q, 3).unemployed = x(2);
    results(q, 4).out = x(3);
end

results.rateU = results.unemployed ./ (results.employed + results.unemployed) * 100;

results.dt = datetime(results.quarter, 'InputFormat', 'yyyyQQQ');

disp(results(:, 1:4));

figure;
plot(results.dt, results.rateU);
title('Unemployment rate in Mexico');
xlabel('Quarter');
ylabel('$\%$');

aux = ones(quarters, 1);

data.flowsEO = aux * mean(data.flowsEO);
data.flowsEU = aux * mean(data.flowsEU);
data.flowsOE = aux * mean(data.flowsOE);
data.flowsOU = aux * mean(data.flowsOU);
data.flowsUO = aux * mean(data.flowsUO);

data.a = data.flowsEU;
data.b = data.flowsUE + data.flowsUO;
data.c = data.flowsOU;
data.d = data.flowsEU + data.flowsEO;
data.e = data.flowsUE;
data.f = data.flowsOE;

for q = 1:quarters
    row = data(q, :);
    A = [-row.a row.b -row.c; row.d -row.e -row.f; 1 1 1];
    x = linsolve(A, b);
    results(q, 1).quarter = row.date;
    results(q, 2).employed = x(1);
    results(q, 3).unemployed = x(2);
    results(q, 4).out = x(3);
end

results.rateU = results.unemployed ./ (results.employed + results.unemployed) * 100;

results.dt = datetime(results.quarter, 'InputFormat', 'yyyyQQQ');

disp(results(:, 1:4));

figure;
plot(results.dt, results.rateU);
title('Counterfactual unemployment rate in Mexico');
xlabel('Quarter');
ylabel('$\%$');