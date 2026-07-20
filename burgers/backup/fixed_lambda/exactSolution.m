clear all;

lambda = 0.125;

% domain
dom = [0 2];
y = chebfun('y', dom);

% nonlinear operator
N_base = chebop(dom);


% boundary conditions
N_base.lbc = 0;
N_base.rbc = -1;


% different problems: formula solve and ode solve
N_burgers = N_base;
N_formula = N_base;

N_burgers.rbc = -1;



N_burgers.op = @(y,u) - lambda * u + ((1+lambda)*y+u) * diff(u,1);
N_formula.op = @(y,u) - u - u.^(1+1/lambda) - y;

initSin = -sin(pi*y/4);
initLin = -1/2*y;

N_burgers.init = initSin;
N_formula.init = initLin;


% solve
u_burgers = N_burgers \ 0;
u_formula = N_formula \ 0;

figure;

% Plot solution
plot(y, u_formula, u_burgers), grid on
xlabel('y'), ylabel('u')

%yLimBot = u(u.domain(end));
yLimBot = min(u_burgers);
yLimTop = max(u_burgers);
ylim([yLimBot yLimTop]);

residual_burgers = y + u_burgers + u_burgers.^(1+1/lambda);
residual_formula = y + u_formula + u_formula.^(1+1/lambda);

title("solving  -\lambda u + ((1+\lambda)y+u) u_y = 0 for \lambda="+num2str(lambda)+" using the solution formula and like it is");
subtitle({"L^2-error burgers="+num2str(norm(residual_burgers)), "L^2-error formula="+num2str(norm(residual_formula))});