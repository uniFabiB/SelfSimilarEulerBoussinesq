clear all;
dom = [-2 2];
cheboppref.setDefaults('display','iter');
cheboppref.setDefaults('maxIter',50); % help cheboppref

x = chebfun('x', dom);
N = chebop(dom);
epsilon=1e-12;  % precision

[fOdd, fOddPrime, lambdaOdd, lambdaOddPrime] = fb.def_fOdd_fOddPrime_lambdaOdd_lambdaOddPrime(dom);
[fFull, fFullPrime, lambdaFull, lambdaFullPrime] = fb.def_fFull_fFullPrime_lambdaFull_lambdaFullPrime(dom);

uInit = chebfun(@(x) -sin(pi*x/4), dom);
%load("u_lambda05Odd");
%uInit = u_lambda05Odd;

N = chebop(dom);

N.op = @(x,u) [
    fOdd(u);
    fFull(u)
];

N.bc = @(x,u) [u(-2)-1; u(0)];

N.init = uInit;

u = N \ [0; 0];

load("input/u_lambda05","u_lambda05");
fb.plotFunctionsAndDiff(uInit, u, u_lambda05);
drawnow;


fprintf('max(abs(uInit-u)) = %.5e\n', max(abs(uInit-u)));
fprintf('abs(lambda(uInit)-lambda(u))) = %.5e\n', abs(lambdaOdd(uInit)-lambdaOdd(u)));
fprintf('lambda(u)-0.5 = %.5e\n', (lambdaOdd(u)-0.5));

%u_lambda05Odd=u;
%save('u_lambda05Odd.mat', 'u_lambda05Odd');