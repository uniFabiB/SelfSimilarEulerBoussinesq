function [w, wPrime] = def_w_wPrime(deflatedVarsList, metricFunction, metricFunctionPrime, p, sigma, c)
    fprintf(['\nFB ERROR\n\def_w_wPrime\n\tNOT IMPLEMENTED YET']);
    % w(u) = c prod_i ( 1/(|f(u, tilde u_i)|^p) + sigma )
    % wPrime(u,v) = c sum_i ( p fPrime(u, tilde u_i; v) / |f(u - tilde
    % u_i)|^(p+1) prod_{i /= j} ( 1/(|f(u, tilde u_j)|^p) + sigma )
    % c = normalizing constant: c^{-1} = prod_i ( 1/(|avgDist|^p) + sigma )  ->
    averageDistance = 0.05; % assumes average distance f(u-random) ~ averageDistance

    % most likely f(u, tilde u_i) = ||u-tilde u_i||
    % but in case of lambda f(u, tilde u_i)= |lambda(u)-lambda(tilde u_i)|
    % so need to specify f(u,tilde u) and fPrime(u, tilde u; v)



    %%% USAGE %%%
    % u based %


    % lambda based %
    %[wUsed, wPrimeUsed] = fb.def_w_wPrime(deflatedLambdas, @(u,uTilde) lambdaUsed(u), @(u, v) lambdaPrimeUsed(u,v), p, sigma, deflatedC);


    if isempty(deflatedVarsList)
        % no previous solutions
        disp("nothing to deflate");
        w = @(u) 1.0;
        wPrime = @(u,v) 0.0;
        return;
    end


    allNummeric = true;
    for var = deflatedVarsList
        if isnumeric(var)
            %pass;
        elseif isa(var, 'chebfun')
            allNummeric = false;
        else
            fprintf('FB WARNING: "%d" is neither chebfun nor numeric', var);
            Error('var neither chebfun nor numeric');
        end
    end



   if allNummeric
        % deflated based on lambda
        % w(u) = c prod_i ( 1/(|lambda(u) - lambda(tilde u_i)|^p) + sigma )
        % usage 
        disp("deflate based on lambda");

        c = 1.0/( (averageDistance.^(2.0)).^(-p/2.0) + sigma ).^length(deflatedVarsList);
        c = 1.0;
        [w, wPrime] = def_wLambda_wLambdaPrime(deflatedVarsList, metricFunction, metricFunctionPrime, p, sigma, c);
   else
        % deflated based on u
        % w(u) = c prod_i ( 1/(|f(u - tilde u_i)|^p) + sigma )
        disp("deflate based on u");

        prodFactor = @(u, tildeu) (((metricFunction(u, tildeu)).^(2.0)).^(-p/2.0) + sigma);
        c = 1/( (averageDistance.^(2.0)).^(-p/2.0) + sigma ).^size(deflatedVarsList,2);
        c = 1.0;
        w = @(u) c;
        wPrime = @(u,v) 0.0;
    
        for i = 1:size(deflatedVarsList,2)
            deflatedVari = deflatedVarsList(:,i);
            w = @(u) w(u).*prodFactor(u,deflatedVari);
            wWithoutI = @(u) 1.0;
            for j = 1:size(deflatedVarsList,2)
                if i ~= j
                    deflatedVarj = deflatedVarsList(:,j);
                    wWithoutI = @(u) wWithoutI(u) .* prodFactor(u,deflatedVarj);
                end
            end
            wPrime = @(u,v) wPrime(u,v) - c.* p.* ( metricFunctionPrime(u, deflatedVari, v) ).* ( ((metricFunction(u, deflatedVari)).^(2.0)).^(-(p+1.0)/2.0) ) .* wWithoutI(u);
        end
   end


end


%%% old deflated lambda %%%
function [wLambda, wLambdaPrime] = def_wLambda_wLambdaPrime(deflatedLambdas, lambdaUsed, lambdaPrimeUsed, p, sigma, c)
    wLambda = @(u) c;
    wLambdaPrime = @(u,v) 0.0;
    for i = 1:length(deflatedLambdas)
        wLambda = @(u) wLambda(u).*( ((lambdaUsed(u)-deflatedLambdas(i)).^2.0).^(-p/2.0) + sigma);
        wLambdaWithoutI = @(u) 1.0;
        for j = 1:length(deflatedLambdas)
            if i ~= j
                wLambdaWithoutI = @(u) wLambdaWithoutI(u) .* ( ((lambdaUsed(u)-deflatedLambdas(j)).^2.0).^(-p/2.0) + sigma);
            end
        end
        wLambdaPrime = @(u,v) wLambdaPrime(u,v) - p.*(lambdaUsed(u)-deflatedLambdas(i)).*lambdaPrimeUsed(u,v) .* ( ((lambdaUsed(u)-deflatedLambdas(i)).^2.0).^(-(p+2.0)/(2.0)) ) .* wLambdaWithoutI(u);
    end
end