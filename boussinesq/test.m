clear all;
% 1. Create 1D periodic chebop operators for x and y
Lx = chebop(-pi, pi); Lx.bc = 'periodic';
Ly = chebop(-pi, pi); Ly.bc = 'periodic';

% 2. Get the discrete differentiation matrices
Nx = 32; Ny = 32; % Grid resolutions
Dx2 = diffmat(Lx, 2, Nx); % Second derivative matrix in x
Dy2 = diffmat(Ly, 2, Ny); % Second derivative matrix in y
Ix = eyemat(Lx, Nx);
Iy = eyemat(Ly, Ny);

% 3. Construct the 2D Laplacian using Kronecker products
L_2D = kron(Iy, Dx2) + kron(Dy2, Ix);

% 4. Set up your periodic RHS on the grid and solve via backslash
[X, Y] = meshgrid(chebpts(Nx, [-pi, pi], 'trig'), chebpts(Ny, [-pi, pi], 'trig'));
F = sin(X).*cos(Y);

% Flatten, solve the linear system, and reshape back to 2D
U_vec = L_2D \ F(:);
U = reshape(U_vec, Ny, Nx);

% Turn back into a smooth Chebfun2 representation
u = chebfun2(U, [-pi, pi, -pi, pi], 'trig');
surf(u)