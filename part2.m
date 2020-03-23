clear all;
load('lab2_1.mat');

% Part 2: Model Estimation 1-D Case

% Parametric Estimation: Gaussian
[a_mu, a_sigma] = gaussianMLE(a);
[b_mu, b_sigma] = gaussianMLE(b);

% Parametric Estimation: Exponential
a_lambda = lambdaGen(a);
b_lambda = lambdaGen(b);

% Parametric Estimation: Uniform
[a_lowerBound, a_upperBound] = uniformMLE(a);
[b_lowerBound, b_upperBound] = uniformMLE(b);

% Non-Parametric Estimation: Parzen Window
x_Values = linspace(0,10,100);

a_parzen_1 = parzen(x_Values, a, 0.1);
a_parzen_4 = parzen(x_Values, a, 0.4);
a_norm = normpdf(x_Values, 5, 1);

figure(1)
plot(x_Values, a_parzen_1)
hold on;
plot(x_Values, a_parzen_4)
hold on;
plot(x_Values, a_norm)

b_parzen_1 = parzen(x_Values, b, 0.1);
b_parzen_4 = parzen(x_Values, b, 0.4);
b_norm = exppdf(x_Values, 1);

figure(2)
plot(x_Values, b_parzen_1)
hold on;
plot(x_Values, b_parzen_4)
hold on;
plot(x_Values, b_norm)

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