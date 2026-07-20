clear all;
dom = [0 2];
y = chebfun('y', dom);
N = chebop(dom);

N.op = @(y, u , lambda) [-lambda.*u + ((1+lambda).*y+u).*diff(u); sum(lambda.*(y.*diff(u)-y))+sum((y+u).*diff(u))];
N.lbc = @(u,lambda) u;
N.rbc = @(u,lambda) u + 1;

uInit = chebfun(@(y) -sin(pi*y/4), dom);
%uInit = chebfun(@(y) -sin(5*pi*y/4), dom);
%uInit = chebfun(@(y) -1/2*y, dom);
%uInit = chebfun(@(y) -sqrt(y/2), dom);
uInitBest = chebfun(@(y) -sin(pi*y/4), dom); % seems close to converge but newtons method fails
lambdaInit = 0.5;


N.init = [uInit; lambdaInit];
uLambda = N\0;
u = uLambda{1};
lambda = uLambda{2}
lambdaResidual = sum(lambda.*(y.*diff(u)-y))+sum((y+u).*diff(u))

residual = y + u + u.^(1+1/lambda);
residual2 = -lambda.*u + ((1+lambda).*y+u).*diff(u);
myPlot = plot(y, u, y, uInit, y, residual, y, residual2);
myPlot(1).LineWidth = 1.5;
legend({'u(y)', 'uInit(y)', 'formula residual(y)', 'operator residual (y)'});