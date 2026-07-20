clear all;

lambda = 0.5;

% domain
dom = [0 2];
y = chebfun('y', dom);

% nonlinear operator
N = chebop(dom);

N.op = @(y,u) - lambda * u + ((1+lambda)*y+u) * diff(u,1) ;

% boundary conditions
N.lbc = 0;
N.rbc = -1;


initSin = -sin(pi*y/4);
initLin = -1/2*y;

N.init = initSin;

% Solve using Newton's method
u = N \ 0;

figure;

% Plot solution
plot(y, u), grid on, hold on
xlabel('y'), ylabel('u(y)')

%yLimBot = u(u.domain(end));
yLimBot = min(u);
yLimTop = max(u);
ylim([yLimBot yLimTop]);

residual = y + u + u.^(1+1/lambda);

title("0 =  - \lambda u + ((1+\lambda)*y+u) * U_y , \lambda="+num2str(lambda));
subtitle({"L^2-error="+num2str(norm(residual)), "   "});