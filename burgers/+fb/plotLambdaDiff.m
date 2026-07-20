function plotLambdaDiff(varargin)
    [~, markerSize] = fb.getMarkersAndSize(0);


    valueList = varargin{1};
    residualList = varargin{2};
    errorList = varargin{3};
    residualQuadraticList = varargin{4};
    

    for i = 1:length(valueList)
        lambda = valueList(i);
        bestLambda = 0;
        for j = 1:10
            testLambda = 1/(2*j);
            if abs(testLambda-lambda)<abs(bestLambda-lambda)
                bestLambda = testLambda;
                minDiff = abs(testLambda-lambda);
            end
        end

    end

    title = strcat("lambda - " + bestLambda + " (min diff = " + minDiff +")");
    figure("Name",title);

    legendLambda = "";
    legendResidual = "";
    legendError = "";
    legendResidualQuadratic = "";
    diff = zeros(length(valueList));
    for i = 1:length(valueList)
        lambda = valueList(i);
        residual = residualList(i);
        error = errorList(i);
        residualQuadratic = residualQuadraticList(i);
        diff(i) = abs(lambda-bestLambda);
        if(i~=1)
            legendLambda = legendLambda + newline;
            legendResidual = legendResidual + newline;
            legendError = legendError + newline;
            legendResidualQuadratic = legendResidualQuadratic + newline;
        end
        legendLambda = legendLambda + lambda;
        legendResidual = legendResidual + residual;
        legendError = legendError + error;
        legendResidualQuadratic = legendResidualQuadratic + residualQuadratic;
        %semilogy(i, abs(lambda-bestLambda), 'x', 'MarkerSize', 5);
    end
    semilogy(1:length(valueList), diff(1:length(valueList)), 'x', 1:length(residualList), residualList(1:length(residualList)), 'o', 1:length(errorList), errorList(1:length(errorList)), '*',  1:length(residualQuadraticList), residualQuadraticList(1:length(residualQuadraticList)), '*', 'MarkerSize', markerSize);


    legend("lambdas" + newline + legendLambda + newline + newline, "residuals" + newline + legendResidual, "errors" + newline + legendError, "residual quadratic" + newline + legendResidualQuadratic, 'Location', 'southwest');
end