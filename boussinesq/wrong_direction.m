clear; clf;
top = 1;
bot = -1;
left = 0;
right = 5;
domain = [left right bot top];
xIndex = 2;     % unintuitively see 13.6 diff, diffx, diffy in https://www.chebfun.org/docs/guide/guide13.html
yIndex = 1;     % unintuitively or https://github.com/chebfun/chebfun/issues/1710

N = chebop2(@(u) diff(u,2,xIndex) + diff(u,2,yIndex), domain);

N.lbc = 0;
N.rbc = 5;
N.ubc = @(x) x;
N.dbc = @(x) x;
u = N \ 0;

plot(diff(u,1,xIndex))
plot(diff(u,1,yIndex))