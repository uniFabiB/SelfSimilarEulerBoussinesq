function [fOdd, fOddPrime, lambdaOdd_used, lambdaOddPrime] = def_fOdd_fOddPrime_lambdaOdd_lambdaOddPrime(dom, varargin)
    x = chebfun('x', dom);
    domLambdaEnd = dom(2);
    domEndLambdaPrime = dom(2);

    %%% diff(U(-x)) = -U'(-x) = d/dx U(-x) %%%
    

    %lambdaOdd = @(U) sum( (x+0.5*(U(x)-U(-x))).*(diff(U(x)) - diff(U(-x))))./sum( U(x)-U(-x)-x.*(diff(U(x))-diff(U(-x))));
    lambdaOdd_halfDomain = @(U) sum( (x+0.5*(U(x)-U(-x))).*(diff(U(x)) - diff(U(-x))), 0, domLambdaEnd)./sum( U(x)-U(-x)-x.*(diff(U(x))-diff(U(-x))),  0, domLambdaEnd);
    lambdaOdd_used = lambdaOdd_halfDomain;
    lambdaOddPrime = @(U, U_tilde) ( sum( 0.5*(U_tilde(x)-U_tilde(-x)) .* (diff(U(x)) - diff(U(-x))) + (x + 0.5*(U(x)-U(-x))) .* (diff(U_tilde(x)) - diff(U_tilde(-x))), 0, domEndLambdaPrime) .* sum( U(x) - U(-x) - x.* (diff(U(x)) - diff(U(-x))), 0, domEndLambdaPrime))./( (sum(U(x)-U(-x) - x.* (diff(U(x)) - diff(U(-x))), 0, domEndLambdaPrime))^2 ) - ( sum( (x+ 0.5*(U(x)-U(-x))).* (diff(U(x))- diff(U(-x))) , 0, domEndLambdaPrime) .* sum( U_tilde(x)-U_tilde(-x) - x.* (diff(U_tilde(x)) - diff(U_tilde(-x)) ), 0, domEndLambdaPrime)) ./ ( (sum(U(x)-U(-x) - x.* (diff(U(x)) - diff(U(-x))), 0, domEndLambdaPrime))^2 );
    
    if nargin > 1
        lambdaOdd_used = @(u) varargin{1};
        lambdaOddPrime = @(u, v) 0;
        fprintf('\nFB WARNING\n\tin def_fOdd_fOddPrime_lambdaOdd_lambdaOddPrime\n\tusing fixed lambda = %e\n\n\n', varargin{1});
    end


    fOdd = @(U) - lambdaOdd_used(U).*(U(x)-U(-x)) + ((1+lambdaOdd_used(U)).*x+0.5*(U(x)-U(-x))).*(diff(U(x))-diff(U(-x)));
    fOddPrime = @(U, U_tilde) - lambdaOddPrime(U,U_tilde) .* (U(x)-U(-x)) - lambdaOdd_used(U) .* (U_tilde(x) - U_tilde(-x)) + (lambdaOddPrime(U,U_tilde) .* x + 0.5*(U_tilde(x) - U_tilde(-x)) ) .* (diff(U(x)) - diff(U(-x))) + ((1+lambdaOdd_used(U)) .* x + 0.5* (U(x)- U(-x))) .* (diff(U_tilde(x)) - diff(U_tilde(-x)));

end