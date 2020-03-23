clear all;
load('lab2_2.mat');
minX = min([min(al(:,1)), min(bl(:,1)), min(cl(:,1))]);
minY = min([min(al(:,2)), min(bl(:,2)), min(cl(:,2))]);
maxX = max([max(al(:,1)), max(bl(:,1)), max(cl(:,1))]);
maxY = max([max(al(:,2)), max(bl(:,2)), max(cl(:,2))]);

res = [1 minX minY maxX maxY];

[a_mu, a_coV] = gaussianMLE(al);
[b_mu, b_coV] = gaussianMLE(bl);
[c_mu, c_coV] = gaussianMLE(cl);

% Create multivariate normal PDF, sigma = 400
[X1,Y1] = meshgrid(1:1:400);
win = mvnpdf([X1(:) Y1(:)], [200, 200], [400 0; 0 400]);
win = reshape(win,length(Y1),length(X1));

[parzen_A, X_A_Vals, Y_A_Vals] = parzen2D(al, res, win);
[parzen_B, X_B_Vals, Y_B_Vals] = parzen2D(bl, res, win);
[parzen_C, X_C_Vals, Y_C_Vals] = parzen2D(cl, res, win);

grid = gridMaker(minX, maxX, minY, maxY);
for i = 1:size(grid.grid, 1)
    for j = 1:size(grid.grid, 2)
        [~, classification_Parzen] = max([parzen_A(i,j), parzen_B(i,j) ,parzen_C(i,j)]);
        [~, classification_Gaussian] = max([ML_classifier(i,j,a_mu,a_coV),ML_classifier(i,j,b_mu,b_coV),ML_classifier(i,j,c_mu,c_coV)]);
        grid.grid(i,j,4) = classification_Parzen;
        grid.grid(i,j,3) = classification_Gaussian;
    end
end

[X2,Y2] = meshgrid(X_A_Vals, Y_A_Vals);
ML = zeros(size(X2));
for i = 1:size(X2,1)
   for j = 1:size(Y2,2)
       [~, class] = max([parzen_A(i,j), parzen_B(i,j), parzen_C(i,j)]);
       ML(i,j) = class;
   end
end

%Plot
figure(1)
hold on;
%contour(X_A_Vals, Y_A_Vals, classGauss);
contour(grid.grid(:,:,1), grid.grid(:,:,2), grid.grid(:,:,3));
%contour(grid.grid(:,:,1), grid.grid(:,:,2), grid.grid(:,:,4));
contour(X_A_Vals, Y_A_Vals, ML);
scatter(at(:,1), at(:,2), 20, 'r','filled');
scatter(bt(:,1), bt(:,2), 20, 'g','filled');
scatter(ct(:,1), ct(:,2), 20, 'b','filled');
hold off;

function [mu, coV] = gaussianMLE(dataSet)
    sampleSize = size(dataSet,1);
    mu = [sum(dataSet(:,1))/sampleSize sum(dataSet(:,2))/sampleSize];
    coV = cov(dataSet);
end

function probML = ML_classifier(x, y, u, s)
    u = u';
    X = [x; y];
    h = X - u;
    r = h';
    probML = (exp(-0.5*((r)*inv(s)*(h))))/(2*pi*sqrt(det(s)));
end