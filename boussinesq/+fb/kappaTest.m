function kappaTest(baseField, perturbation, func, directionalDerivativeExpression)
    tau = logspace(-20,0,21);
    kappa = zeros(size(tau));
    fU = func(baseField);
    gradExpression = directionalDerivativeExpression(baseField,perturbation);
    denom = sum(abs(gradExpression));

    % more pointwise test!!!!
    % see notes

   
    
    for i = 1:numel(tau)
        fUphPert = func(baseField+tau(i)*perturbation);
        enum = sum(abs(fU - fUphPert)/tau(i));  % L^1 type int|f(u+tau per)-f(u)|dx/tau
        frac = enum./denom;
        sum(frac);
        kappa(i)=frac;
        fprintf('tau %e,\tkappa = %e, \tlog(abs(1-kappa)) = %f\n', tau(i), frac, log10(abs(1-frac)));
        if i==1
            bestTau = tau;
            bestFrac = frac;
            bestLog = log10(abs(1-frac));
        else
            if abs(frac-1) < abs(bestFrac-1)
                bestTau = tau(i);
                bestFrac = frac;
                bestLog = log10(abs(1-frac));
            end
        end
    end
    fprintf('__________________________________________________________________\n');
    fprintf('best tau %e,\tkappa = %e, \tlog(abs(1-kappa)) = %f\n', bestTau, bestFrac, bestLog);
    %title = strcat(title , ", ", inputname(i)); % does not work if argument is not a single variable like a+b
    figure("name", strcat("kappa test ", inputname(3)));
    loglog(tau,abs(kappa-1),'x',MarkerSize=16);
end