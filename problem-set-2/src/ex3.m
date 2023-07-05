k = 3;
lambda = 1;
samples = 1e3;

U = rand(samples,1); % U ~ Uniform(0,1)
W = inverse_weibull(U,lambda,k); % W ~ Weibull(lambda,k)

% true distribution
x = linspace(min(W),max(W),samples);
X = wblcdf(x,lambda,k);

% fit test
weibullpd = fitdist(W, 'Weibull');

figure
plot(weibullpd,'PlotType',"probability")
title("Weibull Distribution with $\lambda=$" + lambda + " and $k=$" + k)

figure
hold on
cdfplot(W)
plot(x,X)
title("Weibull CDF with $\lambda=$" + lambda + " and $k=$" + k)
legend('Empirical CDF','Theoretical CDF')
hold off

% theoretical values
meanW = lambda * gamma(1 + 1/k);
medianW = lambda * ((log(2))^(1/k));
varW = (lambda^2) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2);

fprintf('Absolute error between theoretical and empirical values:\n')
fprintf('Mean:     |%f - %f| = %f\n', meanW, mean(W), abs(meanW - mean(W)))
fprintf('Median:   |%f - %f| = %f\n', medianW, median(W), abs(medianW - median(W)))
fprintf('Variance: |%f - %f| = %f\n', varW, var(W), abs(varW - var(W)))