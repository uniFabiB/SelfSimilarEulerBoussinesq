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
%load("input/u_lambda05Odd");
%uInit = u_lambda05Odd;

N = chebop(dom);


deflatedLambdas = [0.5]; % [0.5, 0.25, ...]         expected: lambda_i = 1/(2i+2), i in N_0
w = @(u) 1;
wPrime = @(u,v) 0;
p = 1;
sigma = 1;
lambdaUsed = @(u) lambdaFull(u);
for i = 1:length(deflatedLambdas)
    w = @(u) w(u).*(1./(((lambdaUsed(u)-deflatedLambdas(i)).^2.0).^(p/2.0)) + sigma);
    wWithoutI = @(u) 1;
    for j = 1:length(deflatedLambdas)
        if i ~= j
            wWithoutI = @(u) wWithoutI(u) * (1./(((lambdaUsed(u)-deflatedLambdas(j)).^2.0).^(p/2.0)) + sigma);
        end
    end
    wPrime = @(u,v) wPrime(u,v) - p.*(lambdaUsed(u)-deflatedLambdas(i)).*lambdaUsed(v) ./ (((lambdaUsed(u)-deflatedLambdas(i)).^2.0).^((p-2.0)/(2.0))) .* wWithoutI(u);
end

gOdd = @(u) w(u).*fOdd(u);
gFull = @(u) w(u).*fFull(u);
gOddPrime = @(u,v) w(u).*fOddPrime(u,v) + fOdd(u).*wPrime(u,v);
gFullPrime = @(u,v) w(u).*fFullPrime(u,v) + fFull(u).*wPrime(u,v);

%N.op = @(x,u) [fOdd(u);fFull(u)];
N.op = @(x,u) [gOdd(u);gFull(u)];
%N.op = @(x,u) [fOdd(u)/(lambdaOdd(u)-0.5);fFull(u)/(lambdaFull(u)-0.5)];

N.bc = @(x,u) [u(-2)-1; u(0)];

N.init = uInit;

u = N \ [0; 0];

load("input/u_lambda05","u_lambda05");
residual = fFull(u);
fb.plot(uInit, u, u_lambda05);
fb.plotFunctionsAndDiff(uInit, u, u_lambda05);
drawnow;

fprintf('lambda = %.5e\n', lambdaFull(u));
fprintf('||f(u)||_infty = %.5e\n', max(abs(residual)));
fprintf('\n');

fprintf('max(abs(uInit-u)) = %.5e\n', max(abs(uInit-u)));
fprintf('abs(lambda(uInit)-lambda(u))) = %.5e\n', abs(lambdaOdd(uInit)-lambdaOdd(u)));
fprintf('lambda(u)-0.5 = %.5e\n', (lambdaOdd(u)-0.5));

u_temp=u;
save('u_temp.mat', 'u_temp');