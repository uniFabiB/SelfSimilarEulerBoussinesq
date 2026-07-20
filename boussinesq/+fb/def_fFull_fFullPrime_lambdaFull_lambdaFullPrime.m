function [fFull, fFullPrime, lambdaFull, lambdaFullPrime] = def_fFull_fFullPrime_lambdaFull_lambdaFullPrime(dom, varargin)
    fprintf(['\nFB ERROR\n\tdef_fFull_fFullPrime_lambdaFull_lambdaFullPrime\n\tNOT IMPLEMENTED YET']);

    r = chebfun2(@(r,phi) r);
    phi = chebfun2(@(r,phi) phi);


    % (210) in self_similar_euler_2026_06_01.pdf
    F1 =@() 1./r .* diffx(r.* )
    
    
    
    y = chebfun('y', dom);
    domLambdaEnd = dom(2);
    domEndLambdaPrime = dom(2);
    lambdaFull = @(U) sum((y+U).*diff(U), 0, domLambdaEnd)./sum(U-y.*diff(U), 0, domLambdaEnd);
    lambdaFullPrime = @(U,U_tilde) ( ...
        (sum(U_tilde.*diff(U), 0, domEndLambdaPrime)+sum((y+U).*diff(U_tilde), 0, domEndLambdaPrime)).*sum(U-y.*diff(U), 0, domEndLambdaPrime) ...
        -sum((y+U).*diff(U), 0, domEndLambdaPrime).*sum(U_tilde-y.*diff(U_tilde), 0, domEndLambdaPrime) ...
        ) / ( ...
        ((sum(U-y.*diff(U), 0, domEndLambdaPrime))^2)...
        );


    if nargin > 1
        lambdaFull = @(u) varargin{1};
        lambdaFullPrime = @(u, v) 0;
        fprintf('\nFB WARNING\n\tin def_fFull_fFullPrime_lambdaFull_lambdaFullPrime\n\tusing fixed lambda = %e\n\n\n', varargin{1});
    end


    fFull = @(U) - lambdaFull(U).*U + ((1+lambdaFull(U)).*y+U).*diff(U);
    fFullPrime = @(U,U_tilde) ( ...
        - lambdaFullPrime(U,U_tilde).*U ...
        - lambdaFull(U) .* U_tilde ...
        + (lambdaFullPrime(U,U_tilde).*y+U_tilde).*diff(U) + ...
        ((1+lambdaFull(U)).*y+U).*diff(U_tilde) ...
        );
end
