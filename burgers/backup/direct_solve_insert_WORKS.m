clear all;
dom = [0 2];

y = chebfun('y', dom);
N = chebop(dom);
epsilon=1e-12;  % precision

N.op = @(y, u) -(sum((y+u).*diff(u))/sum(u-y.*diff(u))).*u + ((1+sum((y+u).*diff(u))/sum(u-y.*diff(u))).*y+u).*diff(u);
N.lbc = @(u) u;
N.rbc = @(u) u + 1;


uInit = chebfun(@(y) -sin(pi*y/4), dom);   % WORKS!!
%uInit = chebfun(@(y) -sin(5*pi*y/4), dom);
%uInit = chebfun(@(y) -1/2*y, dom);
%uInit = chebfun(@(y) -sqrt(y/2), dom);
uInitBest = chebfun(@(y) -sin(pi*y/4), dom); % seems close to converge but newtons method fails


N.init = uInit;
u = N\0;



lambda = sum((y+u).*diff(u))/sum(u-y.*diff(u));


fprintf('\nsolved finished \n\t lambda = %f\n\n\n', lambda);

residual = y + u + u.^(1+1/lambda);
residual2 = -lambda.*u + ((1+lambda).*y+u).*diff(u);


%lambdaFormula=2^(-5);
lambdaFormula=lambda;

% use positive values, by point symmetry to avoid complex things
y_of_minusU = chebfun(@(u) u + u.^(1+1/lambdaFormula), [0 1]);
minusU_of_y = inv(y_of_minusU)
if abs(minusU_of_y.domain(1)-dom(1))<epsilon
    minusU_of_y.domain(1)=dom(1);
else
    plot(-minusU_of_y);
    legend('y(U)');
    error("error\n\tminusU_of_y.domain(1) = %f is not close to dom(1) = %f\n\t |minusU_of_y.domain(1) - dom(1)| ~ 10^(%f)",minusU_of_y.domain(1), dom(1), log10(abs(minusU_of_y.domain(1)-dom(1))));
end
if abs(minusU_of_y.domain(end)-dom(end))<epsilon
    minusU_of_y.domain(end)=dom(end);
else
    plot(-minusU_of_y);
    legend('y(U)');
    error("error\n\tminusU_of_y.domain(end) = %f is not close to dom(end) = %f\n\t |minusU_of_y.domain(end) - dom(end)| ~ 10^(%f)",minusU_of_y.domain(end), dom(end), log10(abs(minusU_of_y.domain(end)-dom(end))));
end
u_formula = - minusU_of_y
residualFormula = max(y + u_formula + u_formula.^(1+1/lambdaFormula))

myPlot = plot(y, u, y, uInit, y, residual, y, residual2);
myPlot(1).LineWidth = 1.5;

if true
    myPlot = plot(y, u, y, u_formula, y, uInit, y, residual, y, residual2);
    markers = {'x','o','+','*','s','d','^','v','>','<','.','p','h'};
    for i = 1:numel(myPlot)
        myPlot(i).MarkerIndices = 1:7:numel(myPlot(i).XData);
        myPlot(i).Marker = markers(i);
    end
    legend(myPlot,strcat('solve (lambda=',string(lambda),')'), strcat('formula (lambda=',string(lambdaFormula),')'), 'uInit', 'residual solve formula', 'residual solve burgers');
end