clear all;
dom = [-2 2];
domPos= [0 2];
myTol = 10e-4;

x = chebfun('x', dom);


[fOdd, fOddPrime, lambdaOdd, lambdaOddPrime] = fb.def_fOdd_fOddPrime_lambdaOdd_lambdaOddPrime(dom);

[fFull, fFullPrime, lambdaFull, lambdaFullPrime] = fb.def_fFull_fFullPrime_lambdaFull_lambdaFullPrime(dom);

uInit = chebfun(@(y) -sin(pi*y/4), dom);

load("input/u_lambda05Odd");

i=0;
converged=0;
tau=1;
dir=chebfun(@(x) 0, dom);
dBar=0;
u=uInit;


N_solve = chebop(dom);
N_solve.op = @(x,u) [fOdd(u); fFull(u)];
N_solve.init = uInit;
%N_solve.bc = @(u) [u(-2)-1; u(0); u(1)+u(-1); u(0.5)+u(-0.5); u(0.1)+u(-0.1); u(0.2)+u(-0.2); u(0.3)+u(-0.3); u(0.4)+u(-0.4); u(2)+1; sum(u)];
%N_solve.bc = @(x,u) sum(u);
N_solve.bc = @(x,u) [u(-2)-1; u(0)];
%cheboppref.setDefaults('display','iter')
%uSolve = N_solve\[0;0];
%fb.plot(uSolve);

while (~converged)
    i = i + 1;
    uOld=u;
    tauOld=tau;
    dirOld=dir;
    dBarOld=dBar;
    
    N_dir = chebop(dom);

    N_dir.op = @(x,dir) [
        fOddPrime(uOld,dir) + fOdd(uOld);
        fFullPrime(uOld,dir) + fFull(uOld)
    ];
    
    N_dir.bc = @(x,dir) [dir(-2); dir(0)];   % if appropriate
    
    N_dir.init = dirOld;
    
    direction = N_dir \ [0; 0];

    [tau, newDir, dBar] = fb.calcDampedStepSizeAndDirectionV0(uOld, direction, N_dir, fFull, fFullPrime, dom, i, tauOld, dirOld, dBarOld);

    %newDir = direction;

    u = uOld + tau*newDir;
    lambda = fb.calcLambda(u);

    if norm(u-uOld)< myTol
        converged=1;
        fprintf('\nconvergence criterium fulfilled, exiting loop\n');

    end
    fprintf('\niteration %d finished\n\tnew lambda = %e\n\tnew norm(u-uOld) = %e\n\n\n', i, lambda, norm(u-uOld));
    %fb.plotFunctionsAndDiff(uInit, u, u_lambda05);
    %drawnow; % Updates the figure window immediately
end


lambda = fb.calcLambda(u);
fprintf('\nfinished after %d iterations\n\tfinal lambda = %e\n\n\n', i, lambda);
fb.plot(uInit, u, u_lambda05Odd);
fb.plotFunctionsAndDiff(uInit, u, u_lambda05Odd);