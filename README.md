# Heat-Diffusion-LED
Model heat diffusion on thin rod in real time on addressable RGBs.

## Demonstration
Below is a couple of simulations using colormaps from matplotlib displayed on addressable RGBs.
Gnuplot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

## The model
The heat conduction is modeled using the one dimensional heat equation to solve for the temperature $T$, on the domain $\Omega = [0, L]$:

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T - \beta \frac{\partial T}{\partial x} - \gamma T + q(x, t), \quad x \in \Omega.$$

The first term on the right-hand side (RHS) represents heat spread through diffusion, with $\alpha$ being the diffusion coefficient. The second term on the RHS is an advection term, with $\beta$ being the advection velocity. The third term models the heat exchange with the surrounding air, with $\gamma$ being the thermal conductivity between the air and the rod. Lastly $q$ is a forcing term representing an external source heating the rod. This could for for example be a lighter centered in the middle of the rod, or even radioactive decay.

For boundary conditions insulated ends are used:
$$\frac{\partial T}{\partial x} = 0, \quad x=0 or x=1$$,
meaning no heat escapes through the boundary.

The initial temperature is set to T = 0.

### Spatial discretization
The problem is solved numerically using the finite element method (FEM). To discretizate we first define the space of possible solutions (also known as a trial space) $$\mathcal{U}$$, so that $T \in \mathcal{U}$. Next we multiply by a test function $v \in \mathcal{V}$, where is $\mathcal{V}$ is the space of test functions (the so-called test space), and integratinging over the domain. This leads to the time-continuous weak formulation of the problem: Find $$T \in \mathcal{U}$$ such that:

$$ \left(v, \frac{\partial T}{\partial t} \right)_{\Omega} = -\alpha \left(\frac{\partial v}{\partial x}, \frac{\partial u}{\partial x} \right)_{\Omega} - \beta \left( v, \frac{\partial u}{\partial x} \right)_{\Omega} - \gamma \left(v, u \right) + \left(v, q \right)_{\Omega}.$$

Here the notation $$ \left( \cdot, \cdot \right) = \int_{\Omega} \cdot d \Omega $$ is used to denote an inner product.

With the weak form set up, we perform a discretization by restricting test and trial spaces to functions that are piece wise polynomial over each element of the domain.


### Temporal discretization
To discretize in time we use the Crank-Nicolson discretization which for the problem is unconditionally stable and second order accurate w.r.t. to the time-step size $$\Delta t$$.
