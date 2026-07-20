clear all;
dom = [0 2];

y = chebfun('y', dom);
N = chebop(dom);
epsilon=1e-12;  % precision



uInit = chebfun(@(y) -sin(pi*y/4), dom);
%uInit = chebfun(@(y) -sin(pi*y/4)*sin((1-cos(pi*y/4))*pi/2), dom);
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
uRandom = randnfun(15, dom);
uArbitrary1= chebfun(@(x) sin(x) + x.^2, dom);
uArbitrary2= chebfun(@(x) exp(-x.^2).*cos(5*x) + log(1+x.^2), dom);
uTest = chebfun(@(y) -sin(5*pi*y/4), dom);


%plot(uInit);

u = uInit;

tau=1;

for i = 1:2
    uOld=u;
    lambda_of_U = @(U) sum((y+U).*diff(U))./sum(U-y.*diff(U));
    lambdaPrime = @(U,U_tilde) ( ...
        (sum(U_tilde.*diff(U))+sum((y+U).*diff(U_tilde))).*sum(U-y.*diff(U)) ...
        -sum((y+U).*diff(U)).*sum(U_tilde-y.*diff(U_tilde)) ...
        ) / ( ...
        ((sum(U-y.*diff(U)))^2)...
        );
    f = @(U) - lambda_of_U(U).*U + ((1+lambda_of_U(U)).*y+U).*diff(U);
    fPrime = @(U,U_tilde) ( ...
        - lambdaPrime(U,U_tilde).*U ...
        - lambda_of_U(U) .* U_tilde ...
        + (lambdaPrime(U,U_tilde).*y+U_tilde).*diff(U) + ...
        ((1+lambda_of_U(U)).*y+U).*diff(U_tilde) ...
        );
    N.op = @(u) fPrime(uOld,u) - fPrime(uOld,uOld) + tau*f(uOld);
    N.init = uOld;
    N.lbc = @(u) u;
    N.rbc = @(u) u + 1;
    %N.lbc = @(u) diff(u);   % wrong bc but shows that it goes trough a couple of iterations before error
    %N.rbc = @(u) diff(u);   % wrong bc but shows that it goes trough a couple of iterations before error
    u = N\0;
    %kappaTest(uArbitrary2, uRandom, f, fPrime);
    lambda = lambda_of_U(u);
    disp(strcat("new lambda = ",string(lambda)));
    fprintf('\niteration %d finished\n\tnew lambda = %f\n\n\n', i, lambda');
end



function kappaTest(baseField, perturbation, func, directionalDerivativeExpression)
    tau = logspace(-20,0,21);
    fU = func(baseField);
    gradExpression = directionalDerivativeExpression(baseField,perturbation);
    denom = sum(abs(gradExpression));
    for i = 1:numel(tau)
        fUphPert = func(baseField+tau(i)*perturbation);
        enum = sum(abs(fU - fUphPert)/tau(i));
        frac = enum./denom;
        sum(frac);
        fprintf('tau %e,\tkappa = %e, \tlog(abs(1-kappa)) = %f\n', tau(i), frac, log10(abs(1-frac)));
    end
end

if max(abs(imag(u))) < epsilon
    fprintf('u is real');
else
    fprintf('u is complex');    
end
fprintf('\t(max(abs(imag(u))) = %f)\n',max(abs(imag(u))));



lambda = sum((y+u).*diff(u))/sum(u-y.*diff(u));


fprintf('\nsolved finished \n\t lambda = %f\n\n\n', lambda);

residual_u_to_formula = y + u + u.^(1+1/lambda);
residual_u_to_equation = -lambda.*u + ((1+lambda).*y+u).*diff(u);
residual_u_to_newtonstep = fPrime(uOld,u) - fPrime(uOld,uOld) + tau*f(uOld);


%lambdaFormula=2^(-5);
lambdaFormula=lambda;

plotFormula = true;
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


if true
    myPlot = plot(y, u, y, u_formula, y, uInit, y, residual_u_to_formula, y, residual_u_to_equation);
    markers = {'x','o','+','*','s','d','^','v','>','<','.','p','h'};
    for i = 1:numel(myPlot)
        myPlot(i).MarkerIndices = 1:7:numel(myPlot(i).XData);
        myPlot(i).Marker = markers(i);
    end
    legend(myPlot,strcat('solve (lambda=',string(lambda),')'), strcat('formula (lambda=',string(lambdaFormula),')'), 'uInit', 'residual solve formula', 'residual solve burgers');
end


if true
    figure();
    myLogPlot = semilogy(y, abs(residual_u_to_equation), y, abs(residual_u_to_newtonstep));
    legend({"residual u to equation (lambda=" + lambda + ")","residual u to newton step (lambda=" + lambda + ")"}, 'Location', 'southeast');
end