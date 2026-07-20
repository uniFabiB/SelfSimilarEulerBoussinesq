
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
