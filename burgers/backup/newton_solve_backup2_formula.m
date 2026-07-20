clear all;
dom = [0 2];
y = chebfun('y', dom);
N = chebop(dom);

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

tau=1e-04;

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

residualNewton = ((2-tau)*lambdaOld-lambda).*uOld-lambdaOld.*u+((lambda-(2-tau)*lambdaOld).*y+u-(2-tau)*uOld).*diff(uOld)+((1+lambdaOld).*y+uOld).*diff(u)-(1-tau)*y.*diff(uOld);
residualBurgers = -lambda.*u + ((1+lambda).*y+u).*diff(u);
residualFromFormula = - u - u.^(1+1/lambda) - y;

lambdaFormula=0.5;
N_formula = chebop(@(y,u) - u - u.^(1+1/lambdaFormula) - y, dom);
N_formula.lbc = 0;
N_formula.rbc = -1;

N_formula2 = chebop(@(y,u) 1+(1+(1+1/lambdaFormula)*u.^(1/lambdaFormula)).*diff(u), dom);

N_formula2.lbc = 0;
N_formula2.rbc = -1;
N_formula2.init = uInit;
%u_formula2 = N_formula 2 \ 0;

bounds=1000;
dom2=[dom -4 4]
f = chebfun2(@(y,u) - u - u.^(1+1/lambdaFormula) - y, dom2);
c = roots(f) %The zero contours of a function are computed by Chebfun2 to plotting accuracy and they are typically not accurate to machine precision.


y_of_U = chebfun(@(u) -u-u.^(1+1/lambdaFormula), [-1 0]);
U_of_y = inv(y_of_U);
plot(U_of_y);
legend(strcat('formule solve ( lambda = ',string(lambdaFormula),')'));

figure;

%u_formula3 = roots(@(u) - u - y2, dom);
%N_formula.init = uInit;
%u_formula = N_formula \ 0;
%myPlot = plot(y, u, y, uInit, y, c, y, residualNewton);

myPlot(1).LineWidth=2;

R=2.2;
X=40;
%T = fzero(@(t) exp(t)-exp(1))

legend(strcat('solve ( lambda = ',string(lambda),' )'), 'uInit', strcat('formula ( lambda = ',string(lambdaFormula),' )'), 'residualNewton');