clear all;
dom = [0 2];

y = chebfun('y', dom);
N_solve = chebop(dom);
N_direction = chebop(dom);
%epsilon=1e-12;  % precision
epsilon=1e-6;  % precision
convergence=1e-4; % convergence tolerance


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
uRandom = randnfun(0.1, dom);
uArbitrary1= chebfun(@(x) sin(x) + x.^2, dom);
uArbitrary2= chebfun(@(x) exp(-x.^2).*cos(5*x) + log(1+x.^2), dom);
uTest = chebfun(@(y) -sin(5*pi*y/4), dom);



%plot(uInit);

u = uInit;
direction = 0;
dBar = 0;
tau=1.0;

for i = 1:3
    uOld=u;
    dOld=direction;
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

    N_solve.op = @(u) fPrime(uOld,u) - fPrime(uOld,uOld) + tau*f(uOld);
    N_solve.init = uOld;
    N_solve.lbc = @(u) u;
    N_solve.rbc = @(u) u + 1;
    %u = N_solve\0;

    %N_direction.op = @(direction) fPrime(uOld,direction) + f(uOld);
    N_direction.op = @(direction) fPrime(uOld,direction);
    N_direction.init = uOld;
    N_direction.lbc = @(direction) direction;
    N_direction.rbc = @(direction) direction;
    direction = N_direction\f(uOld);
    direction = -direction;
    
    myU=uOld;
    d=direction;
    fuOld=f(uOld);
    save(strcat("myNewtonVars",string(i)),"f", "myU", "N_direction", "d", "fPrime", "tau", "fuOld")


    %%% Algorithm 4 GENIB: G on p. 70 in thesis %%%
    %%% Birkisson - Numerical Solution of Bo... %%%

    %% 
    [tau, direction, dBar] = calcDampedStepSizeAndDirectionMyFromCode(u, direction, N_direction, f, fPrime, dom, i, tau, dOld, dBar);




    tau;
    norm(direction);
    u = uOld + tau*direction;

    
    %N.lbc = @(u) diff(u);   % wrong bc but shows that it goes trough a couple of iterations before error
    %N.rbc = @(u) diff(u);   % wrong bc but shows that it goes trough a couple of iterations before error
    
    %kappaTest(uArbitrary2, uRandom, f, fPrime);
    lambda = lambda_of_U(u);
    
    %fprintf('\niteration %d finished\n\tnew lambda = %f\n\n\n', i, lambda');
end

function [tau, dOut, dBar] = calcDampedStepSizeAndDirectionMyFromThesis(u, d, N_d_lin, f, fPrime, dom, newtonIter, tauOld, dOld, dBarOld)
    tau = 1.0;
    dOut = d;
    dBar = d;
end


function [tau, dOut, dBar] = calcDampedStepSizeAndDirectionMyFromCode(myU, d, N_d_lin, f, fPrime, dom, newtonIter, tauOld, dOld, dBarOld)
    tau = tauOld;
    epsilon = 1e-6;
    tauMin = epsilon;
    endSearch = false;
    cFactor=1.0;
    dOut = d;
    initPrediction = 1;
    while (~ endSearch)


        if(newtonIter>1 && initPrediction)
            mu = norm(dOld)*norm(dBarOld)/(norm(dBarOld-d, 'fro')*norm(d))*tau;
            tau=min(1,mu);
            initPrediction = 0;
        end

        if(tau<tauMin)
            % tau too small, try full newton to try a new base point
            tau=1.0;
            dOut = d;
            dBar = d;
            break;
        end

        uTrial = myU + tau*d;
        resTrial = feval(f, uTrial);



        N_temp = chebop(dom);
        N_temp.op = @(d) fPrime(myU,d);
        N_temp.init = 0;
        N_temp.lbc = @(d) d;
        N_temp.rbc = @(d) d;
        dBar = N_temp\resTrial;

        %dBar = N_d_lin\0;
        dBar = -dBar;


        % contraction factor
        cFactor = norm(dBar)/norm(d);

        % correction factor
        muPrime = 1.0/2.0*norm(d)*tau^2/(norm(dBar - (1-tau)*d,'fro'));



        if(cFactor>1)
            tau = min(muPrime,1.0/2.0*tau);
            continue
        end

        tauPrime = min(muPrime,1.0);
        if(abs(tauPrime-1.0)<epsilon)
            tau = 1.0;
            %dOut = dBar;
            break;
        else
            if(tauPrime>4*tau)
                tau = tauPrime;
            else
                break
            end
        end
        
    end

    %norm(d)
    %norm(dBar)
    %norm(myU+tau*dOut)
    save(strcat("myDampVars",string(newtonIter)), "myU", "d", "dBar", "fPrime", "resTrial")
    fprintf('\t%8s   %8s   %8s\n','||d||','cFactor','tau');
    fprintf('\t%8.4f   %8.4f   %8.4f\n', norm(dOut), cFactor, tau);
end


function [tau, dOut] = calcDampedStepSizeAndDirectionCopied(u,d,N_d_lin)

    accept = false;
    tau = 1.0;
    epsilon = 1e-06;

    initPredictor = true;

    % Iterate until we find a step-size lambda that we accept:
    while ( ~accept )
        
        if ( tau < epsilon)
            % If LAMBDA falls below LAMBDAMIN, we try to take a full Newton step in
            % the hope that will put us in a different point in solution space that
            % we have a chance of converging from. If we don't observe convergence
            % in the following step, we then give up.
            dOut = d;
            tau = 1.0;
            accept = 1;
            continue
        end
        
        % Take a trial step
        uTrial = u + tau*d;


        resTrial = feval(N_d_lin, uTrial);
        
        % Compute a simplified Newton step using the current derivative of the
        % operator, but with a new right-hand side.
        dBar = N_d_lin\resTrial;
        
        
        % We had two output arguments above, need to negate deltaBar:
        dBar = -dBar;    
        
        % TODO: Do we need to update the values of the RHS for the BCs here?
          
        
        % Contraction factor:
        cFactor = norm(dBar)/norm(d);
        
        % Correction factor for the step-size:
        %%% implementation and thesis mismatch!!! %%%
        muPrime = (.5*norm(d)*tau^2) / norm(dBar*tau-(1-tau)*d, 'fro');
        
        % If we don't observe contraction, decrease LAMBDA
        if ( cFactor >= 1 )
            tau = min(muPrime, .5*tau);
            % Go back to the start of the loop.
            continue
        end
        
        % New potential candidate for LAMBDA
        tauPrime = min(1, muPrime);
        
        if ( tauPrime == 1 && norm(dBar) < convergence)
            % We have converged within the damped phase! 
            % solvebvpNonlinear() will find out about our success.
            dOut = dBar;
            tau = 1.0;
            break
        end
       
        
        % TODO: Document
        if ( tauPrime >= 4*tau )
            tau = tauPrime;
            continue
        end
        
        % If we get all the way here, accept iterate, and tell the Newton iteration
        % to keep up the good work!


        dOut = d;
        break;


    end
end

function kappaTest(baseField, perturbation, func, directionalDerivativeExpression)
    tau = logspace(-20,0,21);
    fU = func(baseField);
    gradExpression = directionalDerivativeExpression(baseField,perturbation);
    denom = sum(abs(gradExpression));

    % more pointwise test!!!!
    % see notes
    
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


residual_u_to_formula = y + u + u.^(1+1/lambda);
residual_u_to_equation = -lambda.*u + ((1+lambda).*y+u).*diff(u);
residual_u_to_newtonstep = fPrime(uOld,u) - fPrime(uOld,uOld) + tau*f(uOld);
residual_u_minus_formula= u-u_formula;
residual_uFormual_formula = y + u_formula + u_formula.^(1+1/lambdaFormula);
residual_uFormual_equation = -lambdaFormula.*u_formula + ((1+lambdaFormula).*y+u_formula).*diff(u_formula);


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
    myLogPlot = semilogy(y, abs(residual_u_to_equation), y, abs(residual_u_to_newtonstep), y, abs(residual_u_minus_formula), y, abs(residual_uFormual_formula), y, abs(residual_uFormual_equation));
    legend({"residual u to equation (lambda=" + lambda + ")","residual u to newton step (lambda=" + lambda + ")","u("+ lambda + ") - u\_formula(" + lambdaFormula + ")","residual u\_formula to formula (" + lambdaFormula + ")","residual u\_formula to equation (" + lambdaFormula + ")"}, 'Location', 'southeast');
end