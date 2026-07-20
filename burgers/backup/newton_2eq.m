clear all;
close all;
dom = [0 2];
y = chebfun('y', dom);
N = chebop(dom);
epsilon=1e-12;  % precision

splitting on;

initLin = chebfun(@(y) -1/2*y, dom);
initLinPer = initLin*0.00001*cos(2*pi*y);
init2Sin = chebfun(@(y) -sin(5*pi*y/4), dom);
initSin = chebfun(@(y) -sin(pi*y/4), dom);
initSqrt = chebfun(@(y) -sqrt(y/2), dom);
uInit = initSin;
lambdaOld = 0.5;
uOld = uInit;


%N.op = @(y, u , lambda) -lambda.*u + ((1+lambda).*y+u).*diff(u);
N.lbc = @(u,lambda) u;
N.rbc = @(u,lambda) u + 1;

tau=1;


for i = 1:1
    N.op = @(y, u , lambda) [((2-tau)*lambdaOld-lambda).*uOld-lambdaOld.*u+((lambda-(2-tau)*lambdaOld).*y+u-(2-tau)*uOld).*diff(uOld)+((1+lambdaOld).*y+uOld).*diff(u)-(1-tau)*y.*diff(uOld);
        sum(((2-tau)*lambdaOld-lambda).*uOld-lambdaOld.*u+((lambda-(2-tau)*lambdaOld).*y+u-(2-tau)*uOld).*diff(uOld)+((1+lambdaOld).*y+uOld).*diff(u)-(1-tau)*y.*diff(uOld))];
    %N.op = @(y, u , lambda) ((2-tau)*lambdaOld-lambda).*uOld-lambdaOld.*u+((lambda-(2-tau)*lambdaOld).*y+u-(2-tau)*uOld).*diff(uOld)+((1+lambdaOld).*y+uOld).*diff(u)-(1-tau)*y.*diff(uOld);
    N.init = [uOld; lambdaOld];
    uLambda = N\0;
    u = uLambda{1};
    lambda = uLambda{2};
    uOld=u;
    lambdaOld=lambda;
    %disp(strcat("new lambda = ",string(lambda)));
    fprintf('\niteration %d finished\n\tnew lambda = %f\n\n\n', i, lambda');
end


fprintf('\nnewton iteration finished \n\t new lambda = %f\n\n\n', lambda);

%lambdaFormula=2^(-5);
lambdaFormula=1;

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

if false

legend({'u(y)', 'uInit(y)', 'formula residual(y)', 'operator residual (y)'});
    formulaPlot = plot(y, u_formula, y, u_formula + u_formula.^(1+1/lambdaFormula), y, y + u_formula + u_formula.^(1+1/lambdaFormula));
    legend(formulaPlot, 'u_formula', 'u_formula + u_formula.^(1+1/lambdaFormula)', 'residual = y + u_formula + u_formula.^(1+1/lambdaFormula)');
    title(gca, 'ONLY REAL PARTS!!!')
end


if true
    myPlot = plot(y, u, y, uInit, y, u_formula);
    myPlot(1).LineWidth=2;
    legend(myPlot,strcat('solve (lambda = ',string(lambda),')'), 'uInit', strcat('formula (lambda = ',string(lambdaFormula),')'));
end

