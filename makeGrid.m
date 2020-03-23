  
function [xValues, yValues, grid] = makeGrid (NumberDataPoints, varargin)

    numberClasses = length(varargin);
    minXValues = zeros(numberClasses);
    minYValues = zeros(numberClasses);
    maxXValues = zeros(numberClasses);
    maxYValues = zeros(numberClasses);
    
    for i = 1:numberClasses
        xValues = varargin{i}(:,1);
        yValues = varargin{i}(:,2);
        minXValues(i) = min(xValues);
        minYValues(i) = min(yValues);
        maxXValues(i) = max(xValues);
        maxYValues(i) = max(yValues);
    end
    
    xValues = min(minXValues)-1:NumberDataPoints:max(maxXValues)+1;
    yValues = min(minYValues)-1:NumberDataPoints:max(maxYValues)+1;
    
    grid = zeros(length(yValues), length(xValues));   
end