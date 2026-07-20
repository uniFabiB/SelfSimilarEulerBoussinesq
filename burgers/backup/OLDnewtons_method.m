clear all;


% domain
dom = [0 2];
y = chebfun('y', dom);


initSin = -sin(pi*y/4);

u_init = initSin;
lambda_init = 0.52;

u_old = u_init;
lambda_old = lambda_init;
% nonlinear operator
N = chebop(dom);
N.op = @(y,u,lambda) [ -lambda_old * u + (u-u_old)*diff(u_old,1) + ((1+lambda_old)*y+u_old)*diff(u,1);  -lambda*u_old + ((1+lambda)*y+u_old)*diff(u_old,1)];


% boundary conditions
N.lbc = 0;
N.rbc = 0;



%N.init = [u_init; lambda_init];
u0 = u_init;
lambda0 = lambda_init;


% differential operator
%N.op = @(y,u,lambda) (lambda_old-lambda)*u_old - lambda_old*u + ((lambda-lambda_old)*y+u-u_old)*diff(u_old,1) + ((1+lambda_old)*y+u_old)*diff(u,1);

[u_new lambda_new] = N \ 0;

