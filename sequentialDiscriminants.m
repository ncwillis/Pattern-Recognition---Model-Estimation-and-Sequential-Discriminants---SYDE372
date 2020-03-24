clear all;
close all;
load('lab2_3.mat');

j = 1;

minX = min([min(a(:,1)), min(b(:,1))]);
minY = min([min(a(:,2)), min(b(:,2))]);
maxX = max([max(a(:,1)), max(b(:,1))]);
maxY = max([max(a(:,2)), max(b(:,2))]);
grid = gridMaker(minX, maxX, minY, maxY);

aData = a;
bData = b;

G = zeros(400);
nbAArray = [];
naBArray = [];
naB = 1;
nbA = 1;

while((size(a,1) ~= 0) && (size(b,1) ~= 0))    
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
    figure(j)
    hold on;
    contour(grid.grid(:,:,1), grid.grid(:,:,2), disc);
    scatter(a(:,1), a(:,2), 20, 'r','filled');
    scatter(b(:,1), b(:,2), 20, 'g','filled');
    scatter(randPA(1), randPA(2), 'black', 'filled');
    scatter(randPB(1), randPB(2), 'black', 'filled');
    hold off;
    a = newA;
    b = newB;
    j = j+1;
    figNum = j;
end

k = 1;
for i = 1:size(grid.grid, 1)
    for j = 1:size(grid.grid, 2)
        classified = false;
        while(classified == false)
            classifier = G(:,:,k);
            if (classifier(i,j) == 2) && (naBArray(k) == 0)
                classification = 2;
                classified = true;
            elseif (classifier(i,j) == 1) && (nbAArray(k) == 0)
                classification = 1;
                classified = true;
            else
                k = k + 1;
            end
        end
        grid.grid(i,j,4) = classification;
    end
end

figure(figNum+1)
hold on;
contour(grid.grid(:,:,1), grid.grid(:,:,2), grid.grid(:,:,4));
scatter(aData(:,1), aData(:,2), 20, 'r','filled');
scatter(bData(:,1), bData(:,2), 20, 'g','filled');
hold off;



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