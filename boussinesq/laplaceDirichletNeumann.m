clear; clf; close all;
eps = 1.01;

top = 1;
bot = -1;
left = -1;
right = 1;
domain = [left right bot top];
xIndex = 2;     % unintuitively see 13.6 diff, diffx, diffy in https://www.chebfun.org/docs/guide/guide13.html
yIndex = 1;     % unintuitively or https://github.com/chebfun/chebfun/issues/1710

N = chebop2(@(u) lap(u), domain);

%works
N.lbc = @(y) sin(pi*y);
N.rbc = @(y,u) diff(u,1);
N.ubc = 0;
N.dbc = 0;


%N.bc = @(x,u) [feval(diff(u),0)-1];
%N.lbc = 'periodic';
%N.rbc = 'periodic';                        % u(1,y) = 0
%N.ubc = 'periodic';                        % u(x,1) = 0
%N.dbc = 'periodic';                        % u(x,-1) = 0
%N.bc = 'periodic'

%N.bc =  @(x, y, u) [u(-1), u(1)];
%N.bc = { @(y,u) u, @(y,u) u, @(x,u) u, @(x,u) u };

%N.lbc = @(x) sin(pi*x);
%N.rbc = @(x) sin(pi*x);
%N.ubc = @(x) sin(pi*x);
%N.dbc = @(x) sin(pi*x);


%N.bc = @(x) sin(pi*x); %works
%N.bc = @(x,u) u;



%N.bc = @(x,y,u) [feval(u, x, 'bot'), feval(u, x, 'right'), feval(u, x, 'top'), feval(u, x, 'left')];

%N.lbc = @(y,u) u - sin(pi*y);
%N.rbc = @(y,u) diff(u,1);
%N.ubc = @(x,u) diff(u,1);
%N.dbc = @(x,u) diff(u,1);

u = N \ 0;

f=chebfun2(@(x,y) sin(pi*x).*exp(y),domain,'trigx');
%plot(f)

uu = u(:,top);
ud = u(:,bot);
ul = u(left,:);
ur = u(right,:);

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

fb.plot(ul, ur, u_diffxl, u_diffxr);
fb.plot(uu, ud, u_diffyd, u_diffyu);


%fb.plot(u_diffyu);

%figure("Name","diff u");
%plot(diff(u));
%figure("Name","diff_x u");
%plot(diff(u,1,xIndex));
%figure("Name","diff_y u");
%plot(diff(u,1,yIndex));
