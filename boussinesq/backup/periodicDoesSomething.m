clear;
domain = [0 10 -1 1];

N = chebop2(@(u) diffx(u,2) + diffy(u,2), domain);
N.lbc = 0;
N.rbc = @(y) 1-y^2;
N.ubc = 'periodic';

fprintf("\nFB ERROR\n\tclearly this periodic bc for 'up' (N.ubc = 'periodic';) changes something, but it is not at all periodic");

f = chebfun2(@(x,y) 0, domain);
u = N \ f;

plot(u)