clear all;
dom = [0 2];

y = chebfun('y', dom);
N = chebop(dom);
epsilon=1e-12;  % precision

[f, fPrime, lambda_of_U, lambdaPrime] = fb.def_half_f_fPrime_lambda_lambdaPrime(dom);
N.op = @(u) f(u);
%N.op = @(y, u) -(sum((y+u).*diff(u))/sum(u-y.*diff(u))).*u + ((1+sum((y+u).*diff(u))/sum(u-y.*diff(u))).*y+u).*diff(u);
N.lbc = @(u) u;
N.rbc = @(u) u + 1;


%load('input/u_lambda05.mat');
%uInit=u_lambda05;

%load('input/u_lambda03355.mat');
%uInit=u_lambda03355;

%load('input/u_lambda00908.mat');
%uInit=u_lambda00908;

uInit = chebfun(@(y) -sin(pi*y/4), dom);
    % WORKS!!
    % lambda = 0.5
    % says newton failed
%uInit = chebfun(@(y) -sin(pi*y/4)*sin((1-cos(pi*y/4))*pi/2), dom);
    % WORKS!! 
    % lambda = -1.5
%uInit = chebfun(@(y) -sin((1-cos(pi*y/4))*pi/2), dom);
    % WORKS
    % lambda = -2.0
%uInit = chebfun(@(y) -sin(pi*y/4)-0.1*sin((1-cos(pi*y/4))*pi), dom);
    % WORKS
    % lambda = -2.0
%uInit = chebfun(@(y) -sin(pi*y/4)-0.01*sin((1-cos(pi*y/4))*pi), dom);
    % WORKS
    % lambda = -0.498213
    % says newton failed
%uInit = chebfun(@(y) -sin(pi*y/4)-0.09*sin(pi*y/2), dom);
    % WORKS
    % lambda = 0.09078
    % says newton failed
%uInit = chebfun(@(y) -0.9*sin(pi*y/4)-0.09*sin(pi*y/2), dom);
    % WORKS after two runs
    % lambda = 0.500002
    % newton number of iter exceeded
%uInit = chebfun(@(y) -sin(pi*y/4)-0.095*sin(pi*y/2), dom);
    % WORKS
    % lambda = 0.335551
    % says newton failed
%uInit = chebfun(@(y) -sin(pi*sqrt(y/4)), dom);
%uInit = chebfun(@(y) -1/2*y, dom);
%uInit = chebfun(@(y) -sqrt(y/2), dom);
%uInit = chebfun(@(y) -sin(5*pi*y/4), dom);


%plot(uInit);
cheboppref.setDefaults('display','iter')
N.init = uInit;
u = N\0;

if max(abs(imag(u))) < epsilon
    %fprintf('u is real');
else
    %fprintf('u is complex');    
end
fprintf('\t(max(abs(imag(u))) = %f)\n',max(abs(imag(u))));



lambda = sum((y+u).*diff(u))/sum(u-y.*diff(u));


fprintf('\nsolved finished \n\t lambda = %f\n\n\n', lambda);

residual = y + u + u.^(1+1/lambda);
residual2 = -lambda.*u + ((1+lambda).*y+u).*diff(u);


%lambdaFormula=2^(-5);
lambdaFormula=lambda;

% use positive values, by point symmetry to avoid complex things
y_of_minusU = chebfun(@(u) u + u.^(1+1/lambdaFormula), [0 1]);
minusU_of_y = inv(y_of_minusU);
if abs(minusU_of_y.domain(1)-dom(1))<epsilon
    minusU_of_y.domain(1)=dom(1);
else
    plot(-minusU_of_y);
    legend('y(U)');
    error("error\n\tminusU_of_y.domain(1) = %f is not close to dom(1) = %f\n\t |minusU_of_y.domain(1) - dom(1)| ~ 10^(%f)",minusU_of_y.domain(1), dom(1), log10(abs(minusU_of_y.domain(1)-dom(1))));
end
if abs(minusU_of_y.domain(end)-dom(end))<epsilon
    minusU_of_y.domain(end)=dom(end);
else
    plot(-minusU_of_y);
    legend('y(U)');
    error("error\n\tminusU_of_y.domain(end) = %f is not close to dom(end) = %f\n\t |minusU_of_y.domain(end) - dom(end)| ~ 10^(%f)",minusU_of_y.domain(end), dom(end), log10(abs(minusU_of_y.domain(end)-dom(end))));
end
u_formula = - minusU_of_y
%residualFormula = max(y + u_formula + u_formula.^(1+1/lambdaFormula))

if true
    figure("Name","u, residuals");
    myPlot = plot(y, u, y, u_formula, y, uInit, y, residual, y, residual2);
    markers = {'x','o','+','*','s','d','^','v','>','<','.','p','h'};
    for i = 1:numel(myPlot)
        myPlot(i).MarkerIndices = 1:7:numel(myPlot(i).XData);
        myPlot(i).Marker = markers(i);
    end
    legend(myPlot,strcat('solve (lambda=',string(lambda),')'), strcat('formula (lambda=',string(lambdaFormula),')'), 'uInit', 'residual solve formula', 'residual solve burgers');
end


if true
    figure("Name","residual solve burgers");
    semilogy(y,abs(residual2));
end

%u_lambda05 = u;
%u_lambda03355 = u;
u_lambda00908 = u;

fprintf('max(abs(uInit-u)) = %.5e\n', max(abs(uInit-u)));
fprintf('abs(lambda(uInit)-lambda(u))) = %.5e\n', abs(lambda_of_U(uInit)-lambda_of_U(u)));
fprintf('lambda(u)-0.5 = %.5e\n', (lambda_of_U(u)-0.5));
save("u_lambda00908.mat","u_lambda00908");