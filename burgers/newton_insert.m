clear all;
dom = [0 2];

y = chebfun('y', dom);
%epsilon=1e-12;  % precision
epsilon=1e-6;  % precision

[f, fPrime, lambda_of_U, lambdaPrime] = fb.def_f_fPrime_lambda_lambdaPrime(dom);

uInitGood = chebfun(@(y) -sin(pi*y/4), dom);
uInit2 = chebfun(@(y) -sin(pi*y/4)*sin((1-cos(pi*y/4))*pi/2), dom);
%uInit = chebfun(@(y) -sin((1-cos(pi*y/4))*pi/2), dom);
%uInit = chebfun(@(y) -sin(pi*y/4)-0.1*sin((1-cos(pi*y/4))*pi), dom);
%uInit = chebfun(@(y) -sin(pi*y/4)-0.01*sin((1-cos(pi*y/4))*pi), dom);
%uInit = chebfun(@(y) -sin(pi*y/4)-0.09*sin(pi*y/2), dom);
%uInit = chebfun(@(y) -0.9*sin(pi*y/4)-0.09*sin(pi*y/2), dom);
%uInit = chebfun(@(y) -sin(pi*y/4)-0.095*sin(pi*y/2), dom);
%uInit = chebfun(@(y) -sin(pi*sqrt(y/4)), dom);
%uInit = chebfun(@(y) -1/2*y, dom);
%uInit = chebfun(@(y) -sqrt(y/2), dom);
%uInit = chebfun(@(y) -sin(5*pi*y/4), dom);
uRandom = randnfun(0.1, dom);
uArbitrary1= chebfun(@(x) sin(x) + x.^2, dom);
uArbitrary2= chebfun(@(x) exp(-x.^2).*cos(5*x) + log(1+x.^2), dom);
uArbitrary3= chebfun(@(x) exp(-x.^2).*cos(5*x) + log(1+x.^2) + x.^2 + sqrt(x), dom);
uTest = chebfun(@(y) -sin(5*pi*y/4), dom);

load('input/u_NewtonOptimized.mat'); %u_NewtonOptimized
load('input/u_lambda05.mat'); %u_lambda05
load('input/u_lambda00908.mat');

pert0 = chebfun(@(y) sin(20*pi*y), dom)*chebfun(@(y) sin(pi*y), dom);
%pert0 = chebfun(@(y) sin(pi*y), dom);
epsPert = 0.0000000001;
pert = 1*epsPert*pert0;
uInit = u_lambda05 + pert;


u = uInit;
lambdaInit = lambda_of_U(uInit);

tau=0.0001;

converged = false;

uBase=uInit;

N_test = chebop(dom);
N_test.op = @(dir) fPrime(uBase,dir);
N_test.init = 0;
N_test.lbc = @(dir) dir;
N_test.rbc = @(dir) dir;
test = N_test\pert;
test=0.5*test
zero = chebfun(@(x) 0, dom);

fb.plot(test,pert);
fb.plotFunctionsAndDiff(test,pert);
%fb.plot(test,pert,uBase);
%fb.plotFunctionsAndDiff(test,pert,uBase);

i=0;
%while (~converged)
while (i<0)
    i = i + 1;
    uOld=u;

    N_solve = chebop(dom);
    N_solve.op = @(u) fPrime(uOld,u) - fPrime(uOld,uOld) + f(uOld);
    N_solve.init = uOld;
    N_solve.lbc = @(u) u;
    N_solve.rbc = @(u) u + 1;
    u_step = N_solve\0;


    N_solve2 = chebop(dom);
    N_solve2.op = @(u) fPrime(uOld,u-uOld) + f(uOld);
    N_solve2.init = uOld;
    N_solve2.lbc = @(u) u;
    N_solve2.rbc = @(u) u + 1;
    u_step2 = N_solve2\0;


    N_solve3 = chebop(dom);
    N_solve3.op = @(u) fPrime(uOld,u-uOld);
    N_solve3.init = uOld;
    N_solve3.lbc = @(u) u;
    N_solve3.rbc = @(u) u + 1;
    u_step3 = -N_solve3\f(uOld);

    dir_uStep = (u_step - uOld);
    dir_uStep2 = (u_step2 - uOld);
    dir_uStep3 = (u_step3 - uOld);

    N_direction = chebop(dom);
    N_direction.op = @(dir) fPrime(uOld,dir);
    N_direction.init = 0;
    N_direction.lbc = @(dir) dir;
    N_direction.rbc = @(dir) dir;
    direction = -N_direction\f(uOld);


    %u = uOld - tau*direction;
    %tau = tau;

    pert_=-pert;
    direction_=direction;

    %figure();
    %plot(y,u);
    fPrime_Uold_Direction = fPrime(uOld, direction);
    fPrime_Uold_diruStep = fPrime(uOld, dir_uStep);
    fPrime_Uold_diruStep2 = fPrime(uOld, dir_uStep2);
    fPrime_Uold_diruStep3 = fPrime(uOld, dir_uStep3);
    
    fb.plot(direction_, dir_uStep, dir_uStep2, dir_uStep3);
    
    fb.plot(direction_, dir_uStep, dir_uStep2, dir_uStep3, pert_);
    fb.plotFunctionsAndDiff(direction_, dir_uStep, dir_uStep2, dir_uStep3, pert_);


    fb.plot(fPrime_Uold_Direction, fPrime_Uold_diruStep, fPrime_Uold_diruStep2, fPrime_Uold_diruStep3);
    fb.plotFunctionsAndDiff(fPrime_Uold_Direction, fPrime_Uold_diruStep, fPrime_Uold_diruStep2, fPrime_Uold_diruStep3);
    
    fb.plot(fPrime_Uold_Direction);

    %if(max(abs(u-uOld))<10e-10)
    %    disp("convergence criterion fulfilled");
    %    converged = true;
    %end
    %tau = tau*0.8;
    
    

    fb.kappaTest(uArbitrary2, uRandom, f, fPrime);
    fb.kappaTest(uArbitrary2, uRandom, lambda_of_U, lambdaPrime);
    fb.kappaTest(uOld, uRandom, f, fPrime);
    fb.kappaTest(uOld, uRandom, lambda_of_U, lambdaPrime);
    fb.kappaTest(uOld, direction, f, fPrime);
    fb.kappaTest(uOld, direction, lambda_of_U, lambdaPrime);
    lambda = lambda_of_U(u)
    fprintf('\niteration %d finished\n\tnew lambda = %e\n\n\n', i, lambda');
end



if max(abs(imag(u))) < epsilon
    fprintf('u is real');
else
    fprintf('u is complex');    
end
fprintf('\t(max(abs(imag(u))) = %f)\n',max(abs(imag(u))));



lambda = lambda_of_U(u);
fprintf('lambda(uInit) = %.20f\n', lambda_of_U(uInit));
fprintf('lambda(u) = %.20f\n', lambda);



fprintf('\nsolved finished \n\t lambda = %f\n\n\n', lambda);



%lambdaFormula=2^(-5);
lambdaFormula=lambda;

plotFormula = false;
if(plotFormula)
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
    %residualFormula = max(y + u_formula + u_formula.^(1+1/lambdaFormula))
else
    u_formula = chebfun(@(y) 0, dom);
end

if false
    residual_u_to_formula = y + u + u.^(1+1/lambda);
    residual_u_to_equation = -lambda.*u + ((1+lambda).*y+u).*diff(u);
    residual_u_to_newtonstep = fPrime(uOld,u) - fPrime(uOld,uOld) + tau*f(uOld);
    residual_u_minus_formula= u-u_formula;
    residual_uFormual_formula = y + u_formula + u_formula.^(1+1/lambdaFormula);
    residual_uFormual_equation = -lambdaFormula.*u_formula + ((1+lambdaFormula).*y+u_formula).*diff(u_formula);
end



if false
    myPlot = plot(y, u, y, u_formula, y, uInit, y, residual_u_to_formula, y, residual_u_to_equation);
    markers = {'x','o','+','*','s','d','^','v','>','<','.','p','h'};
    for i = 1:numel(myPlot)
        myPlot(i).MarkerIndices = 1:7:numel(myPlot(i).XData);
        myPlot(i).Marker = markers(i);
    end
    legend(myPlot,strcat('solve (lambda=',string(lambda),')'), strcat('formula (lambda=',string(lambdaFormula),')'), 'uInit', 'residual solve formula', 'residual solve burgers');
end


if false
    figure();
    myLogPlot = semilogy(y, abs(residual_u_to_equation), y, abs(residual_u_to_newtonstep), y, abs(residual_u_minus_formula), y, abs(residual_uFormual_formula), y, abs(residual_uFormual_equation));
    legend({"residual u to equation (lambda=" + lambda + ")","residual u to newton step (lambda=" + lambda + ")","u("+ lambda + ") - u\_formula(" + lambdaFormula + ")","residual u\_formula to formula (" + lambdaFormula + ")","residual u\_formula to equation (" + lambdaFormula + ")"}, 'Location', 'southeast');
end

%u_NewtonOptimized = u;
%save("u_NewtonOptimized.mat","u_NewtonOptimized");