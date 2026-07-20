clear all;

N = chebop2(@(u) diffx(u,2) + diffy(u,2));
N.domain = [-1 1 -1 1];
f = chebfun2(@(x,y) 0);
N.rbc = @(u) u + 1;

N.lbc = @(y) sin(pi*y);   % left boundary
N.rbc = @(y) exp(-y);   % right boundary
N.ubc = @(x) x*exp(-1);   % upper boundary
N.dbc = @(x) x*exp(1);   % lower boundary

u = N \ f;

plot(u)