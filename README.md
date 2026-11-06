# Heat-Diffusion-LED
Model heat diffusion on thin rod in real time on addressable RGBs.

## Demonstration
Below is a couple of simulations using colormaps from matplotlib displayed on addressable RGBs.
Gnuplot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

## The model
The heat conduction is modeled using the one dimensional heat equation to solve for the temperature $$T$$, on the domain $$\Omega = [0, L]$$:

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T - \beta \frac{\partial T}{\partial x} - \gamma T + q(x, t), \quad x \in \Omega.$$

The first term on the right-hand side (RHS) represents heat spread through diffusion, with $$\alpha$$ being the diffusion coefficient. The second term on the RHS is an advection term, with $$\beta$$ being the advection velocity. The third term models the heat exchange with the surrounding air, with $$\gamma$$ being the thermal conductivity between the air and the rod. Lastly $$q$$ is a forcing term representing an external source heating the rod. This could for for example be a lighter centered in the middle of the rod, or really anything that can generate heat such as radioactive decay.

For boundary conditions insulated ends are used:
$$\frac{\partial T}{\partial x} = 0, \quad x=0$$ or $$x=1$$,
meaning no heat escapes through the boundary.

The initial temperature is set to $$T=0$$.

### Spatial discretization
The problem is solved numerically using the finite element method (FEM). To discretizate we first define the space of possible solutions (the so-called trial space) $$\mathcal{U}$$, so that $$T \in \mathcal{U}$$. Next we multiply by a test function $$v \in \mathcal{V}$$, where $$\mathcal{V}$$ is the space of test functions (the so-called test space), and integratinging over the domain. This leads to the time-continuous weak formulation of the problem: Find $$T \in \mathcal{U}$$ such that:

$$ \left(v, \frac{\partial T}{\partial x} \right) =$$

Here the notation $$ \left( \cdot, \cdot \right) = \int_{\Omega} \cdot d \Omega $$ is used to denote an inner product.

With the weak form in place, the finite element discretization can now be performed by restricting trial- and test spaces to functions that are piecewise polynomial over each element of the domain, which we denote by $$\mathcal{U}_h$$ and $$\mathcal{V}_h$$, respectively. To discretize we make use of the fact that any function $$u \in \mathcal{U}_h$$ can be written as a linear combination $$T(x, t) = \sum_{j=1}^{n} u_j(t) \varphi_j(x)$$, and consequently its spatial derivative is $$\frac{\partial T}{\partial x} = u_j(t) \frac{d \varphi}{dx}$$. Here $$\varphi_j$$ is a polynomial of arbitrary order centered on the node with index $$j$$. Inserting the expression for $$T$$ into the above weak form and letting $$v = \varphi_i(x)$$, we for each $$i \in \{1,...,n\}$$ obtain the following ordinary differential equations:

$$\sum_{j=1}^{n} \frac{d u_j}{dt} (\varphi_i, \varphi_j) = -\alpha \sum_{j=1}^{n} \left(\frac{\varphi_i}{dx}, \frac{\varphi_j}{dx} \right)$$

This system can be written more compactly in terms vectors and matrices:

$$M\frac{d\mathbf{u}}{dt} = -\alpha A \mathbf{u} - \beta S \mathbf{u} - \gamma M \mathbf{u} + \mathbf{f}$$,

where $$M$$ is the mass matrix, $$A$$ the stiffness matrix, $$S$$ is a mixed matrix and $$\mathbf{f}$$ the load vector.

### Temporal discretization
To discretize in time we use the Crank-Nicolson discretization which for the problem is unconditionally stable and second order accurate w.r.t. to the time-step size $$\Delta t$$.
