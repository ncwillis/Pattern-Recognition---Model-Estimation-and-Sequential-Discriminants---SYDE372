clear all;
close all;
load('lab2_1.mat');

% Part 2: Model Estimation 1-D Case

x_Values = linspace(0,10,100);

a_muTrue = 5;
a_sigmaTrue = 1;
b_lambdaTrue = 1;
gaussian_trueA = normpdf(x_Values, a_muTrue, a_sigmaTrue);
exp_trueB =exppdf(x_Values,1/b_lambdaTrue);

% Parametric Estimation: Gaussian
[a_mu, a_sigma] = gaussianMLE(a);
[b_mu, b_sigma] = gaussianMLE(b);

gaussian_estA = normpdf(x_Values, a_mu, a_sigma);
gaussian_estB = normpdf(x_Values, b_mu, b_sigma);

figure(1);
plot(x_Values,gaussian_estA);
hold on;
plot(x_Values,gaussian_trueA);
hold off;
title('Set A');
xlabel('Samples');
ylabel('Probability');
legend('Estimated', 'True');

figure(2);
plot(x_Values, gaussian_estB);
hold on;
plot(x_Values, exp_trueB);
hold off;
title('Set B');
xlabel('Samples');
ylabel('Probability');
legend('Estimated', 'True');


% Parametric Estimation: Exponential
a_lambdaEst = lambdaGen(a);
b_lambdaEst = lambdaGen(b);

exp_estA = exppdf(x_Values, 1/a_lambdaEst);
exp_estB = exppdf(x_Values, 1/b_lambdaEst);

figure(3);
plot(x_Values, exp_estA);
hold on;
plot(x_Values, gaussian_trueA);
hold off;
title('Set A');
xlabel('Samples');
ylabel('Probability');
legend('Estimated', 'True');

figure(4);
plot(x_Values, exp_estB);
hold on;
plot(x_Values, exp_trueB);
hold off;
title('Set B');
xlabel('Samples');
ylabel('Probability');
legend('Estimated', 'True');


% Parametric Estimation: Uniform
[a_lowerBound, a_upperBound] = uniformMLE(a);
[b_lowerBound, b_upperBound] = uniformMLE(b);
x_Values_Unif = linspace(0,10,500);

unif_estA = unifpdf(x_Values_Unif, a_lowerBound, a_upperBound);
unif_estB = unifpdf(x_Values_Unif, b_lowerBound, b_upperBound);

figure(5);
plot(x_Values_Unif, unif_estA);
hold on;
plot(x_Values, gaussian_trueA);
hold off;
title('Set A');
xlabel('Samples');
ylabel('Probability');
legend('Estimated', 'True');

figure(6);
plot(x_Values_Unif, unif_estB);
hold on;
plot(x_Values, exp_trueB);
hold off;
title('Set B');
xlabel('Samples');
ylabel('Probability');
legend('Estimated', 'True');


% Non-Parametric Estimation: Parzen Window
a_parzen_1 = parzen(x_Values, a, 0.1);
a_parzen_4 = parzen(x_Values, a, 0.4);
a_norm = normpdf(x_Values, 5, 1);

b_parzen_1 = parzen(x_Values, b, 0.1);
b_parzen_4 = parzen(x_Values, b, 0.4);
b_norm = exppdf(x_Values, 1);

% Set A, StdDev = 0.1
figure(7);
plot(x_Values, a_parzen_1)
hold on;
plot(x_Values, a_norm)
xlabel('Samples');
ylabel('Probability');
legend('Estimated - StdDev: 0.1', 'True');
title('Set A');

% Set A, StdDev = 0.4
figure(8);
plot(x_Values, a_parzen_4)
hold on;
plot(x_Values, a_norm)
xlabel('Samples');
ylabel('Probability');
legend('Estimated - StdDev: 0.4', 'True');
title('Set A');

% Set B, StdDev = 0.1
figure(9);
plot(x_Values, b_parzen_1)
hold on;
plot(x_Values, b_norm)
xlabel('Samples');
ylabel('Probability');
legend('Estimated - StdDev: 0.1', 'True');
title('Set B');

% Set B, StdDev = 0.4
figure(10);
plot(x_Values, b_parzen_4)
hold on;
plot(x_Values, b_norm)
xlabel('Samples');
ylabel('Probability');
legend('Estimated - StdDev: 0.4', 'True');
title('Set B');

function [mu, sigma] = gaussianMLE(dataSet)
    sampleSize = length(dataSet);
    mu = sum(dataSet)/sampleSize;
    sigma = sqrt(sum((dataSet - mu).^2)/sampleSize);
end

function lambda = lambdaGen(dataSet)
    sampleSize = length(dataSet);
    lambda = sampleSize/sum(dataSet);
end

function [lowerBound, upperBound] = uniformMLE(dataSet)
    lowerBound = min(dataSet);
    upperBound = max(dataSet);
end

function distribution = parzen(x_Values, dataSet, stDev)
    sampleSize = length(dataSet);
    distribution = zeros(size(x_Values));
    
    for index = 1:sampleSize
        x_dist = normpdf(x_Values, dataSet(index), stDev);
        distribution = distribution + x_dist;
    end
    
    distribution = (1/sampleSize)*distribution;
end