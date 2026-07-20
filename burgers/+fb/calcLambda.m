function lambda = calcLambda(U)
    dom = domain(U);
    y = chebfun('y', dom);
    domLambdaEnd = dom(2);
    lambda = sum((y+U).*diff(U), 0, domLambdaEnd)./sum(U-y.*diff(U), 0, domLambdaEnd);
end
