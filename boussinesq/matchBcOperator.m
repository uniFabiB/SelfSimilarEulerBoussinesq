clear; clf; close all;
eps = 1.01;

top = 1;
bot = -1;
left = 0;
right = 2;
domain = [left right bot top];
xIndex = 2;     % unintuitively see 13.6 diff, diffx, diffy in https://www.chebfun.org/docs/guide/guide13.html
yIndex = 1;     % unintuitively or https://github.com/chebfun/chebfun/issues/1710

x=chebfun2(@(x,y) x, domain);
y=chebfun2(@(x,y) y, domain);

N = chebop2(@(x,y,u,v) [diff(u,2,xIndex)+diff(u,2,yIndex); v], domain);

test = chebfun(1);
test2 = chebfun(2);

N.lbc = @(y) sin(pi*y);
%N.rbc = @(y,u) diff(u,1);
N.rbc = @(y) 0;
N.ubc = 0;
N.dbc = 0;

u = N \ 0;


uu = u(:,top);
ud = u(:,bot);

if norm(u(:,top)-u(:,bot)) > 1e-13
    warning("u(:,top) neq u(:,bot), so no periodic boundary conditions on top/bot\n\tnorm(u(:,top)-u(:,bot)) = %e",norm(u(:,top)-u(:,bot)));
end
ul = u(left,:);
ur = u(right,:);

%fb.plotFunctionsAndDiff(uu,ud);

figure("Name","u");
plot(u);

u_diffx = diff(u,1,xIndex);
u_diffy = diff(u,1,yIndex);

u_diffxu = u_diffx(:,top);
u_diffxd = u_diffx(:,bot);
u_diffyu = u_diffy(:,top);
u_diffyd = u_diffy(:,bot);

u_diffxl = u_diffx(left,:);
u_diffxr = u_diffx(right,:);
u_diffyl = u_diffy(left,:);
u_diffyr = u_diffy(right,:);

%fb.plot(ul, ur, u_diffxl, u_diffxr, u_diffyl, u_diffyr);
%fb.plot(uu, ud, u_diffxd, u_diffxu, u_diffyd, u_diffyu);

%fb.plot(ul, ur, u_diffxl, u_diffxr);
%fb.plot(uu, ud, u_diffyd, u_diffyu);


%fb.plot(u_diffyu);

%figure("Name","diff u");
%plot(diff(u));
%figure("Name","diff_x u");
%plot(diff(u,1,xIndex));
%figure("Name","diff_y u");
%plot(diff(u,1,yIndex));
