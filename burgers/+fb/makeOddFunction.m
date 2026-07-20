function oddFunction = makeOddFunction(expression, dom)
    posPart= chebfun(expression, [0 dom(2)]);
    oddFunction = chebfun(@(x) (x>=0).*posPart(x) - (x<0).*posPart(-x), dom, 'splitting', 'on');
    % "splitting on" splits the domain in multiple subdomains so that
    % singularities can be better represented by the series
end
