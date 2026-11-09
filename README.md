# Heat-Diffusion-LED
Model heat diffusion on thin rod in real time on addressable RGBs.

## Demonstration
Below is a couple of simulations using colormaps from matplotlib displayed on addressable RGBs.

Gnuplot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

Hot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

Jet:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

## The model
![Some text](assets/rod_drawing.png)

The heat conduction is modeled using the one dimensional heat equation to solve for the temperature $$T$$, on the domain $$\Omega = [x_0, x_1]$$:

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T - \beta \frac{\partial T}{\partial x} - \gamma T + q(x, t), \quad x \in \Omega.$$

The first term on the right-hand side (RHS) represents heat spread through diffusion, with $$\alpha$$ being the diffusion coefficient. The second term on the RHS is an advection term, with $$\beta$$ being the advection velocity. The third term models the heat exchange with the surrounding air, with $$\gamma$$ being the thermal conductivity between the air and the rod. Lastly $$q$$ is a forcing term representing an external source heating the rod. This could for for example be a lighter centered in the middle of the rod, or really anything that can generate heat such as radioactive decay.

For boundary conditions insulated ends are used:

$$\frac{\partial T}{\partial x} = 0, \quad x=x_0$$ or $$x=x_1$$,

meaning no heat is allowed to escape through the boundary.

The initial temperature is set to $$T=0$$.

### Spatial discretization
The problem is solved numerically using the finite element method (FEM). To discretizate one first has to define the space of possible solutions (the so-called trial space) $$\mathcal{U}$$, so that $$T \in \mathcal{U}$$. Next we multiply by a test function $$v \in \mathcal{V}$$, where $$\mathcal{V}$$ is the space of test functions (the so-called test space), and integratining over the domain. This leads to the time-continuous weak formulation of the problem: Find $$T \in \mathcal{U}$$ such that:

$$ \left(v, \frac{\partial T}{\partial x} \right)_{\Omega} =$$

Here the notation $$ \left( \cdot, \cdot \right) = \int_{\Omega} \cdot d \Omega $$ is used to denote an inner product.

With the weak form in place, the finite element discretization can now be performed by restricting trial- and test spaces to functions that are piecewise polynomial over each element of the discretized domain, which we denote by $$\mathcal{U}_h$$ and $$\mathcal{V}_h$$, respectively. To discretize the PDE we make use of the fact that any function $$u \in \mathcal{U}_h$$ can be written as a linear combination

$$T(x, t) = \sum_{j=1}^{n} u_j(t) \varphi_j(x)$$,

and consequently its spatial derivative is

$$\frac{\partial T}{\partial x} = \sum_{j=1}^{n} u_j(t) \frac{d \varphi}{dx}$$.

Here $$\varphi_j$$ is a polynomial of arbitrary order centered on the node with index $$j$$, shown in Fig. (2) below. Inserting the expression for $$T$$ into the above weak form and letting $$v = \varphi_i(x)$$, we for each $$i \in \{1,...,n\}$$ obtain the following ordinary differential equations:

$$\sum_{j=1}^{n} \frac{d u_j}{dt} (\varphi_i, \varphi_j) = -\alpha \sum_{j=1}^{n} \left(\frac{\varphi_i}{dx}, \frac{\varphi_j}{dx} \right)$$

By defining the following the following matrices and vector with elements given by:

$$M_{ij} = (\varphi_i, \varphi_j)_{\Omega}$$,

$$A_{ij} = (\frac{d \varphi_i}{dx}, \frac{d\varphi_j}{dx})_{\Omega}$$,

$$S_{ij} = (\varphi_i, \frac{d\varphi_j}{dx})_{\Omega}$$,

$$f_i = (\varphi_i, q(x, t))_{\Omega}$$,

this system can be written more compactly:

$$M\frac{d\mathbf{u}}{dt} = -\alpha A \mathbf{u} - \beta S \mathbf{u} - \gamma M \mathbf{u} + \mathbf{f}$$.

### Temporal discretization
To discretize in time we use the Crank-Nicolson (CN) discretization. This method is very efficient as it is second order accurate w.r.t. to the time-step size $$\Delta t$$ and unconditionally stable.

Applying the CN discretization we obtain

$$M\frac{\mathbf{u}^{k+1} - \mathbf{u}^{k}}{\Delta t} = -\alpha A \frac{\mathbf{u}^{k+1} + \mathbf{u}^{k}}{2} - \beta S \frac{\mathbf{u}^{k+1} + \mathbf{u}^{k}}{2} - \gamma M \frac{\mathbf{u}^{k+1} + \mathbf{u}^{k}}{2} + \frac{\mathbf{f}^{k+1} + \mathbf{f}^k}{2}$$,

where $k$ denots the time step. Rearranging term we then arrive at the following linear system for the unknown temperature coefficients, $$\mathbf{u}^{k+1}$$, at time step $k+1$

$$\left ( \left ( 1 + \frac{\gamma \Delta t}{2} \right) M + \frac{\alpha \Delta t}{2} A + \frac{\beta \Delta t}{2} S \right ) \mathbf{u}^{k+1} = \left ( \left ( 1 - \frac{\gamma \Delta t}{2} \right) M - \frac{\alpha \Delta t}{2} A - \frac{\beta \Delta t}{2} S \right ) \mathbf{u}^{k+1} + \Delta t \frac{\mathbf{f}^{k+1} + \mathbf{f}^k}{2}$$.

Defining $$L = \left ( 1 + \frac{\gamma \Delta t}{2} \right) M + \frac{\alpha \Delta t}{2} A + \frac{\beta \Delta t}{2} S$$ and $$R=\left ( 1 - \frac{\gamma \Delta t}{2} \right) M - \frac{\alpha \Delta t}{2} A - \frac{\beta \Delta t}{2} S$$, this system becomes

$$L \mathbf{u}^{k+1} = R \mathbf{u}^k + \Delta t \frac{\mathbf{f}^{k+1} + \mathbf{f}^{k}}{2}.$$

Because matrix $$L$$ is time independent this system can be intgrated in time efficiently since is only needs to be inverted once for the whole simulation.

## Software and libraries used
* Gnu Scientific Library for dense linear algebra
* Suitesparse for sparse linear algebra
* cJSON for parsing the configuration file

## License
All code-related files are realsed under MIT license while all media (found in the assets directory) under Creative Commons license CC-BY-4.0. Please see the respective LICENSE files for terms of use.
