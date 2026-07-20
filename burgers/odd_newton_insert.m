clear all;
dom = [-2 2];

%convergence
myTol = 1e-12;                   %~1e-13 seems to be machine precision
                                % after ~1e-4 often slow
convCount = 1;

%deflation
deflatedP = 2.0;
deflatedSigma = 0.0;

%maxNumCompThreads(6);


[fOdd, fOddPrime, lambdaOdd, lambdaOddPrime] = fb.def_fOdd_fOddPrime_lambdaOdd_lambdaOddPrime(dom);
[fFull, fFullPrime, lambdaFull, lambdaFullPrime] = fb.def_fFull_fFullPrime_lambdaFull_lambdaFullPrime(dom);

load("u05.mat");
load("u025.mat");
load("u01666.mat");
load("um15.mat");
load("utemp.mat");
load("utemp2.mat");


initialParamsFor05 = struct('deflatedUs',[], ...
    'deflatedLambdas',[], ...
    'uInit',chebfun(@(y) -sin(pi*y/4), dom));
    % 6 min runtime for <1e-12


initialParamsFor025 = struct('deflatedUs',u05, ...
    'deflatedLambdas',[], ...
    'uInit',u05 + 0.01* chebfun(@(y) -sin(pi*y/2), dom));
    % 1 min runtime for <1e-12


initialParamsFor01666 = struct('deflatedUs',u025, ...
    'deflatedLambdas',[], ...
    'uInit',u025 + 0.01* chebfun(@(y) -sin(pi*y/2), dom));
    % 1 min runtime for <1e-12


initialParams = initialParamsFor05;

deflatedLambdas = initialParams.deflatedLambdas;
deflatedUs = initialParams.deflatedUs;
uInit = initialParams.uInit;

%deflatedLambdas = []; % [0.5, 0.25, ...]          expected: lambda_i = 1/(2i), i in N => 0.5, 0.25, 0.1666, 0.125, 0.1, 0.08333, ...

%deflatedUs = [];   % [u05, u025, u01666];

%uInit = chebfun(@(y) -sin(pi*y/4), dom);
%uInit = u05;
%uInit = u025;
%uInit = u01666;
%uInit = utemp;
%uInit = utemp2;
%uInit = uInit + 0.01* chebfun(@(y) -sin(pi*y/2), dom);


%chebfunpref.setDefaults('splitting',true)
%cheboppref.setDefaults('display','iter');
xFunc = chebfun(@(x) x , dom);



iter=0;
converged=0;
tau=1;
dir=chebfun(@(x) 0, dom);
dBar=0;
u=uInit;


lambdaList = [];
residualList = [];
errorList = [];
residualQuadraticList = [];


preTime = datetime('now');
startTime = preTime;


lambdaUsed = @(u) lambdaFull(u);
lambdaPrimeUsed = @(u,v) lambdaFullPrime(u,v);

q_normQ = deflatedP;
normQ = @(f) norm(f,q_normQ);
normQPrime = @(f,g) normQ(f).^(1.0-q_normQ) .* sum( (f.^(2.0)).^((q_normQ-2.0)/2.0) .* f.* g );


if ~isempty(deflatedUs)
    [wUsed, wPrimeUsed] = fb.def_w_wPrime(deflatedUs, @(u,uTilde) normQ(u-uTilde), @(u, uTilde, v) normQPrime(u-uTilde,v), deflatedP, deflatedSigma);
elseif ~isempty(deflatedLambdas)
    [wUsed, wPrimeUsed] = fb.def_w_wPrime(deflatedLambdas, @(u,uTilde) lambdaUsed(u), @(u, v) lambdaPrimeUsed(u,v), deflatedP, deflatedSigma);
else
    [wUsed, wPrimeUsed] = fb.def_w_wPrime(deflatedUs, @(u,uTilde) normQ(u-uTilde), @(u, uTilde, v) normQPrime(u-uTilde,v), deflatedP, deflatedSigma);
end


gOdd = @(u) wUsed(u).*fOdd(u);
gFull = @(u) wUsed(u).*fFull(u);
gOddPrime = @(u,v) wUsed(u).*fOddPrime(u,v) + fOdd(u).*wPrimeUsed(u,v);
gFullPrime = @(u,v) wUsed(u).*fFullPrime(u,v) + fFull(u).*wPrimeUsed(u,v);


N_dir = chebop(dom);
N_dir.bc = @(x,dir) [dir(-2); dir(0); dir(2)];


fprintf('\nSTART VALUES\n');
fprintf('\tclock time          %12s\n', string(datetime('now','Format','HH:mm:ss')));
fprintf('\tlambda              %12.10f\n', lambdaUsed(u));
fprintf('\tresidual  f norm    %e\n', norm(fFull(u)));
fprintf('\tresidual wf norm    %e\n', norm(gFull(u)));
%fprintf('\terror       norm    %e\n', errorNorm);
fprintf('\n\n');

while (converged<convCount)
    iter = iter + 1;
    uOld=u;
    tauOld=tau;
    dirOld=dir;
    dBarOld=dBar;
    
    N_dir.init = dirOld;
    N_dir.op = @(x,dir) [gOddPrime(uOld,dir); gFullPrime(uOld,dir)];
    direction = N_dir \ [-gOdd(uOld); - gFull(uOld)];
    
    
    

    f_step = @(u) [gOdd(u); gFull(u)];
    fPrime_step = @(u,v) [gOddPrime(u,v); gFullPrime(u,v)];
    %f_step = @(u) [fOdd(u); fFull(u)];
    %fPrime_step = @(u,v) [fOddPrime(u,v); fFullPrime(u,v)];
    [tau, dir, dBar, cFactor] = fb.calcDampedStepSizeAndDirectionV1(uOld, direction, f_step, fPrime_step, dom, iter, tauOld, dirOld, dBarOld, false);

    u = uOld + tau*dir;

    tempLambda = lambdaUsed(u);
    lambdaGuess = 0;
    for i = 1:10
        testLambda = 1/(2*i);
        if abs(testLambda-tempLambda)<abs(lambdaGuess-tempLambda)
            lambdaGuess = testLambda;
            lambdaGuessReciprocal = 2*i;
        end
    end

    % solution formula y = -U-U^{1+\frac{1}{\lambda}}
    errorFunction = u + u.^(1+lambdaGuessReciprocal) + xFunc;
    errorNorm = norm(errorFunction);

    elapsedTime = seconds(datetime('now') - preTime);
    preTime = datetime('now');
    fprintf('\niteration %d finished\n', iter);
    fprintf('\tclock time          %12s\n', string(datetime('now','Format','HH:mm:ss')));
    fprintf('\tdelta clock time    %11.1fs\n', elapsedTime);
    fprintf('\tlambda              %12.10f\n', lambdaUsed(u));
    fprintf('\tresidual  f norm    %e\n', norm(fFull(u)));
    fprintf('\tresidual wf norm    %e\n', norm(gFull(u)));
    fprintf('\tu-uOld      norm    %e\n', norm(u-uOld));
    fprintf('\terror       norm    %e\n', errorNorm);
    fprintf('\n');

    lambdaList = [lambdaList, lambdaUsed(u)];
    residualList = [residualList, norm(fFull(u))];
    errorList = [errorList, errorNorm];
    residualQuadratic = norm(fFull(u))/((norm(fFull(uOld)))^2);
    residualQuadraticList = [residualQuadraticList, residualQuadratic];
    
    if norm(fFull(u))< myTol
        converged = converged + 1;
        if converged == 1
            fprintf('\nconvergence criterium fulfilled %d time', converged);
        else
            fprintf('\nconvergence criterium fulfilled %d times', converged);
        end
    end
end



fprintf('\n\nEND VALUES\n');
fprintf('\titerations          %12d\n', iter);
fprintf('\tclock time          %12s\n', string(datetime('now','Format','HH:mm:ss')));
fprintf('\ttotal clock time    %11.1fs\n',  seconds(datetime('now') - startTime));
fprintf('\tlambda              %12.10f\n', lambdaUsed(u));
fprintf('\tresidual  f norm    %e\n', norm(fFull(u)));
fprintf('\tresidual wf norm    %e\n', norm(gFull(u)));
fprintf('\terror       norm    %e\n', errorNorm);
fprintf('\n');

%fb.plotFunctionsAndDiff(lambdaList-0.25);
fb.plotLambdaDiff(lambdaList,residualList, errorList, residualQuadraticList);

%fb.plot(uInit, u, u_lambda05Odd);
%fb.plotFunctionsAndDiff(uInit, u, u_lambda05Odd, residual);

%u05=u;
%save("u05", u5);