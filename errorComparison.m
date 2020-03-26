clear all;
close all;
load('lab2_3.mat');

% Make grid
minX = min([min(a(:,1)), min(b(:,1))]);
minY = min([min(a(:,2)), min(b(:,2))]);
maxX = max([max(a(:,1)), max(b(:,1))]);
maxY = max([max(a(:,2)), max(b(:,2))]);
grid = gridMaker(minX, maxX, minY, maxY);

aData = a;
bData = b;
maxJ = 5;
errorMatrix = [];

% Want to make 20 classifiers
for u = 1:20
    
    %Initialize G (discriminant holder), nbAArray, naBArray
    G = zeros(400);
    nbAArray = [];
    naBArray = [];
    naB = 1;
    nbA = 1;
    a = aData;
    b = bData;
    
    % Construct a sequential classifier with 5 discriminants
    j = 1;
    while((size(a,1) ~= 0) && (size(b,1) ~= 0) && (j<maxJ+1))    
        [disc, naB, nbA, randPA, randPB] = getDisc(a, b, grid);
        grid.grid(:,:,3) = disc;
        G(:,:,j) = disc;
        naBArray = [naBArray naB];
        nbAArray = [nbAArray nbA];
        
        if naB == 0
            %remove all points from b that G classifies as B
            newA = a;
            newB = [];
            for i = 1:size(b,1)
                point = [b(i,1) b(i,2)];
                closestIndex = errorFunctions.findClosestIndex(point(1),point(2), grid.grid(:,:,1), grid.grid(:,:,2));
                xInd = closestIndex(1);
             	yInd = closestIndex(2);
             	if grid.grid(xInd, yInd, 3) == 1
                 	newB = [newB; point];
             	end
            end
        end
        if nbA == 0
            %remove all points from a that G classifies as A
            newA = [];
            newB = b;
            for i = 1:size(a,1)
                point = [a(i,1) a(i,2)];
                closestIndex = errorFunctions.findClosestIndex(point(1),point(2), grid.grid(:,:,1), grid.grid(:,:,2));
                xInd = closestIndex(1);
             	yInd = closestIndex(2);
             	if grid.grid(xInd, yInd, 3) == 2
                 	newA = [newA; point];
             	end
            end
        end
        a = newA;
        b = newB;
        j = j+1;
    end
    
    % Classify points for each of J= 1, 2, ...
    errorArray = [];
    for v = 1:maxJ
        %classify points
        errCnt = 0;
        for i = 1:size(aData,1)
            k = 1;
            classified = 0;
            xA = aData(i,1);
            yA = aData(i,2);
            closestIndex = errorFunctions.findClosestIndex(xA, yA, grid.grid(:,:,1), grid.grid(:,:,2));
            minErr = [];
            while(classified == 0)
                classifier = G(:,:,k);
                nABi = naBArray(k) + nbAArray(k);
                minErr = [minErr nABi];
                if (k==v)
                    [~, minK] = min(minErr);
                    classification = G(closestIndex(1,1),closestIndex(1,2), minK);
                    classified = 1;
                elseif ((classifier(closestIndex(1,1),closestIndex(1,2)) == 2) && (naBArray(k) == 0))
                    classification = 2;
                    classified = 1;
                elseif ((classifier(closestIndex(1,1),closestIndex(1,2)) == 1) && (nbAArray(k) == 0))
                    classification = 1;
                    classified = 1;
                else
                    k = k+1;
                end
            end
            if classification ~= 1
                errCnt = errCnt+1;
            end
        end
        for i = 1:size(bData,1)
            k = 1;
            classified = 0;
            xB = bData(i,1);
            yB = bData(i,2);
            closestIndex = errorFunctions.findClosestIndex(xB, yB, grid.grid(:,:,1), grid.grid(:,:,2));
            minErr = [];
            while(classified == 0)
                classifier = G(:,:,k);
                nABi = naBArray(k) + nbAArray(k);
                minErr = [minErr nABi];
                if (k==v)
                    [~, minK] = min(minErr);
                    classification = G(closestIndex(1,1),closestIndex(1,2), minK);
                    classified = 1;
                elseif ((classifier(closestIndex(1,1),closestIndex(1,2)) == 2) && (naBArray(k) == 0))
                    classification = 2;
                    classified = 1;
                elseif ((classifier(closestIndex(1,1),closestIndex(1,2)) == 1) && (nbAArray(k) == 0))
                    classification = 1;
                    classified = 1;
                else
                    k = k+1;
                end
            end
            if classification ~= 2
                errCnt = errCnt+1;
            end
        end
        %calculate error
        N = 400;
        errorRate = errCnt/N;
        errorArray = [errorArray errorRate]
    end
    errorMatrix = [errorMatrix; errorArray]
end

% Calculate max and min errors
avgError = [];
maxError = [];
minError = [];
stdDev = [];

for i = 1:5
    maxi = [max(errorMatrix(:,i));i];
    maxError = [maxError maxi];
    mini = [min(errorMatrix(:,i));i];
    minError = [minError mini];
    avgi = [mean(errorMatrix(:,i));i];
    avgError = [avgError avgi];
    stdDevi = [std(errorMatrix(:,i));i];
    stdDev = [stdDev stdDevi];
end

%make error rate plots
figure(1)
hold on;
plot(maxError(2,:), maxError(1,:), 'o-b','linewidth',2,'markersize',5,'markerfacecolor','b');
plot(minError(2,:), minError(1,:), 'o-r','linewidth',2,'markersize',5,'markerfacecolor','r');
plot(avgError(2,:), avgError(1,:), 'o-g','linewidth',2,'markersize',5,'markerfacecolor','g');
plot(stdDev(2,:), stdDev(1,:), 'o-m','linewidth',2,'markersize',5,'markerfacecolor','m');
hold off;
xlabel('Max J');
ylabel('Error');
legend('Max Error Rate', 'Min Error Rate', 'Average Error Rate', 'Standard Deviation');


function [discriminant, naB, nbA, randPA, randPB]  = getDisc(a_data, b_data, grid)
    confMatrixMED1 = [1 1; 1 1];
    while((confMatrixMED1(1,2) ~= 0) && (confMatrixMED1(2,1) ~= 0))
        indA = randi([1 size(a_data, 1)],1,1);
        indB = randi([1 size(b_data, 1)],1,1);
        
        pA = [a_data(indA,1) a_data(indA,2)];
        pB = [b_data(indB,1) b_data(indB,2)];
        
        for i = 1:size(grid.grid, 1)
            for j = 1:size(grid.grid, 2)
                classification = eDist(grid.grid(i,j,1),grid.grid(i,j,2),pA,pB);
                grid.grid(i,j,3) = classification;
            end
        end
        confMatrixMED1 = [errorFunctions.getErrorCount(a_data, 1, grid, 3) errorFunctions.getErrorCount(a_data, 2, grid, 3); errorFunctions.getErrorCount(b_data, 1, grid, 3) errorFunctions.getErrorCount(b_data, 2, grid, 3)];
    end
    discriminant = grid.grid(:,:,3);
    randPA = pA;
    randPB = pB;
    naB = confMatrixMED1(1,2);
    nbA = confMatrixMED1(2,1);
end

function euclideanDistance = eDist(x,y,muA,muB)
            % calculates euclidean distance
            distA = (x-muA(1))^2 + (y-muA(2))^2;
            distB = (x-muB(1))^2 + (y-muB(2))^2;
            if distA > distB
                %choose distB
                euclideanDistance = 2;
            else
                %choose distA
                euclideanDistance = 1;
            end
        end