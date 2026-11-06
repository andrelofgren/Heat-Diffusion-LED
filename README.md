# Heat-Diffusion-LED
Model heat diffusion on thin rod in real time on addressable RGBs.

## Demonstration
Gnuplot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

## The model
The heat conduction is modeled using the one dimensional heat equation over the domain $\Omega = [0, L]$, where $L$ is the length of the rod:
$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T - \beta \frac{\partial T}{\partial x} - \gamma T, \quad x \in [0, L]$$,
where $T$ is the sought temperature field.

The first term in the right-hand side (RHS) is diffusive term, with $\alpha$ being the diffusion coefficient. The second term in the RHS is advection term, with $\beta$ being the advection velocity. The last term models the heat transfer with the air surrounding the rod, with $\gamma$ being the thermal conductivity between the air and the rod.

For boundary conditions insulated ends are used:
$$\frac{\partial T}{\partial x} = 0, \quad x=0~or~x=1$$,
meaning no heat escapes through the boundary.

The initial temperature is set to T = 0.

### The finite element discretization
The problem is solved numerically using the finite element method (FEM). To discretizate we first define the space of possible solutions $$\mathcal{U}$$, so that $T \in \mathcal{U}$. Next we multiply by a test function $v \in \matcal{V}$, where is $\mathcal{V}$ is the space of test functions, and integratinging over the domain. This leads to the time continuous weak formulation of the problem: Find $$T \in \mathcal{U} s.t:$$
$$ \left (v, \frac{\partial T}{\partial t} = - \alpha (\frac{\partial v}{}) \right)
