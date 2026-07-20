clear all;


% domain
dom = [0 2];
y = chebfun('y', dom);

% nonlinear operator
N_base = chebop(dom);


% boundary conditions
N_base.lbc = 0;
N_base.rbc = -1;


% different problems: formula solve and ode solve
A_burgers = N_base;
B_burgers = N_base;



A_burgers.op = @(y,u) (y + u) * diff(u,1);
B_burgers.op = @(y,u) u - y * diff(u,1);

% solve
[V, lambda] = eigs(A_burgers, B_burgers, 5, 'LR');

eigenvalues = diag(lambda);
eigenvectors = V;