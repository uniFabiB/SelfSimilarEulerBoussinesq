# Description

This a repo for Matlab/[Chebfun](https://www.chebfun.org/) to solve the question of self-similar blow-up for the 3d axisymmetric Euler equations or (almost, see Majda, Bertozzi - Vorticity and incompressible flow 5.4.1) equivalently the 2d Boussinesq Equations.
As a proof of concept, the 1d Burgers' equation is solved with an implementation of a damped Newton's method that can use deflation techniques to get different self-similarity exponents ($\lambda$ in the following).

For more mathematical background see Wang, Lai, G\'omez-Serrano, Buckmaster 2023 ([Journal](https://doi.org/10.1103/PhysRevLett.130.244002), [arXiv](https://arxiv.org/abs/2201.06780))

The idea behind self-similar blow-up is to look for solutions of the form
$$u(x, t) = (1 - t)^\lambda  U\left(\frac{x}{(1 - t)^{\lambda+1}} \right)=(1-t)^\lambda U(y).$$

For the Boussinesq equations one gets velocity formulation
$$
    -\lambda U + \big((\lambda+1) y + U\big)\cdot \nabla U + \nabla P = \rho e_2
    \\
    (1-\lambda) \rho + \left[(\lambda+1) y + U\right]\cdot \nabla \rho = 0
    \\
    \nabla \cdot U = 0
$$
or the vorticity formulation
$$
    \omega + \big((\lambda+1) y + U\big)\cdot \nabla \omega = \Phi
    \\
    (2 + \partial_1 U_1) \Phi + \left[(\lambda+1) y + U\right]\cdot \nabla \Phi = -\partial_1 U_2 \Psi
    \\
    (2+ \partial_2 U_2 ) \Psi + \left[(\lambda+1) y + U\right]\cdot \nabla \Psi = - \partial_2 U_1 \Phi
    \\
    \nabla \cdot U = 0
$$
where $\Phi(y)=\partial_{y_1}\rho(y)$, $\Psi(y)=\partial_{y_2} \rho(y)$, and $\omega =\nabla^\perp_y \cdot U$
and for Burgers' equation one gets
$$
\begin{aligned}
    -\lambda U + ((1+\lambda) y + U) U_y = 0
\end{aligned}
$$

## Usage

TODO    
