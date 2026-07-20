clear all;
dom = [0 2];
y = chebfun('y', dom);
N = chebop(dom);
N.op = @(y, u , lambda) -lambda.*u + ((1+lambda).*y+u).*diff(u);
N.lbc = @(u,lambda) u;
N.rbc = @(u,lambda) u + 1;

%N.init = [chebfun(@(y) -sin(pi*y/4), dom); 0.52];
N.init = [chebfun(@(y) -sin(5*pi*y/4), dom); 0.52];
uLambda = N\0;
u = uLambda{1};
lambda = uLambda{2}


residual = y + u + u.^(1+1/lambda);
residual2 = -lambda.*u + ((1+lambda).*y+u).*diff(u);
plot(u);