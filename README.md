# Heat-Diffusion-LED (WIP)
![Some text](assets/rod_drawing.png)
Model heat diffusion on a thin rod in real time and display the result using addressable RGBs.

## Demonstration
Below is a couple of simulations using colormaps from matplotlib displayed on addressable RGBs.

Gnuplot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

Hot:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

Jet:
https://github.com/user-attachments/assets/866e1d23-5a75-4f3f-b477-7837770ec363

## Wiring
![Some text](assets/circuit.png)

## The model
The heat conduction is modeled using the one dimensional heat equation to solve for the temperature $$T$$, on the spatial domain $$\Omega = [x_0, x_1]$$ over a time interval $$[t_0, t_1]$$:

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T - \beta \frac{\partial T}{\partial x} - \gamma T + q(x, t), \quad (x, t) \in \Omega \times [t_0, t_1]$$.

The first term on the right-hand side (RHS) represents heat spread through diffusion, with $$\alpha$$ being the diffusion coefficient. The second term on the RHS is an advection term, with $$\beta$$ being the advection velocity. The third term models the heat exchange with the surrounding air, with $$\gamma$$ being the thermal conductivity between the air and the rod. Lastly $$q$$ is a forcing term representing an external source heating the rod. This could for for example be a lighter centered in the middle of the rod, or really anything that can generate heat such as radioactive decay.

For boundary conditions insulated ends are used:

$$\frac{\partial T}{\partial x} = 0, \quad x=x_0$$ or $$x=x_1$$,

meaning no heat is allowed to escape through the boundary.

The initial temperature is set to $$T=0$$.

### Spatial discretization
The problem is solved numerically using the finite element method (FEM). To discretizate one first has to define the space of possible solutions (the so-called trial space) $$\mathcal{U}$$, so that $$T \in \mathcal{U}$$. Next we multiply by a test function $$v \in \mathcal{V}$$, where $$\mathcal{V}$$ is the space of test functions (the so-called test space), and integratining over the domain. This leads to the time-continuous weak formulation of the problem: Find $$T \in \mathcal{U}$$ such that:

$$ \left(v, \frac{\partial T}{\partial t} \right)_{\Omega} = \left$$,

where

$$ \left(v, u \right) = \int_{\Omega} v u d \Omega $$

denotes an inner product between functions $$v$$ and $$u$$.

With the weak form in place, the finite element discretization can now be performed by restricting trial- and test spaces to functions that are piecewise polynomial over each element of the discretized domain. These spaces we denote by $$\mathcal{U}_h$$ and $$\mathcal{V}_h$$, respectively. To discretize the PDE we then make use of the fact that any function $$u \in \mathcal{U}_h$$ can be written as a linear combination:

$$T(x, t) = \sum_{j=1}^{n} u_j(t) \varphi_j(x)$$:
Here $$\varphi_j$$ is a polynomial of arbitrary order centered on the node with index $$j$$. Order 1 and order 2 polynomials are shown in the figure below.

Taking the derivative w.r.t. x then gives:

$$\frac{\partial T}{\partial x} = \sum_{j=1}^{n} u_j(t) \frac{d \varphi}{dx}$$.

![Lagrange basis of (a) order 1 and (b) order 2](assets/lagrange_basis.png)

Inserting the expression for $$T$$ into the above weak form and letting $$v = \varphi_i(x)$$, we for each $$i \in \{1,...,n\}$$ obtain the following ordinary differential equations:

$$\sum_{j=1}^{n} \frac{d u_j}{dt} \left (\varphi_i, \varphi_j \right) = -\alpha \sum_{j=1}^{n} \left(\frac{d \varphi_i}{dx}, \frac{d \varphi_j}{dx} \right) u_j - \beta \sum_{j=1}^{n} \left(\varphi_i, \frac{d\varphi_j}{dx} \right ) u_j - \gamma \sum_{j=1}^{n} \left(\varphi_i, \varphi_j \right) u_j + \left(\varphi_i, q(x, t) \right)_{\Omega}$$

By defining the following the following matrices and vector with elements given by:

$$M_{ij} = (\varphi_i, \varphi_j)_{\Omega}$$,

$$A_{ij} = (\frac{d \varphi_i}{dx}, \frac{d\varphi_j}{dx})_{\Omega}$$,

$$S_{ij} = (\varphi_i, \frac{d\varphi_j}{dx})_{\Omega}$$,

$$f_i = (\varphi_i, q(x, t))_{\Omega}$$,

the system can be written more compactly as:

$$M\frac{d\mathbf{u}}{dt} = -\alpha A \mathbf{u} - \beta S \mathbf{u} - \gamma M \mathbf{u} + \mathbf{f}$$.

### Temporal discretization
To discretize in time we use the Crank-Nicolson (CN) discretization. This method is highly efficient as it is second order accurate w.r.t. to the time-step size $$\Delta t$$ and unconditionally stable.

Applying the CN discretization we obtain

$$M\frac{\mathbf{u}^{k+1} - \mathbf{u}^{k}}{\Delta t} = -\alpha A \frac{\mathbf{u}^{k+1} + \mathbf{u}^{k}}{2} - \beta S \frac{\mathbf{u}^{k+1} + \mathbf{u}^{k}}{2} - \gamma M \frac{\mathbf{u}^{k+1} + \mathbf{u}^{k}}{2} + \frac{\mathbf{f}^{k+1} + \mathbf{f}^k}{2}$$,

where $k$ denots the time step. Rearranging terms we then arrive at the following linear system for the unknown temperature coefficients, $$\mathbf{u}^{k+1}$$, at time step $k+1$

$$\left ( \left ( 1 + \frac{\gamma \Delta t}{2} \right) M + \frac{\alpha \Delta t}{2} A + \frac{\beta \Delta t}{2} S \right ) \mathbf{u}^{k+1} = \left ( \left ( 1 - \frac{\gamma \Delta t}{2} \right) M - \frac{\alpha \Delta t}{2} A - \frac{\beta \Delta t}{2} S \right ) \mathbf{u}^k + \Delta t \frac{\mathbf{f}^{k+1} + \mathbf{f}^k}{2}$$.

Defining $$L = \left ( 1 + \frac{\gamma \Delta t}{2} \right) M + \frac{\alpha \Delta t}{2} A + \frac{\beta \Delta t}{2} S$$ and $$R=\left ( 1 - \frac{\gamma \Delta t}{2} \right) M - \frac{\alpha \Delta t}{2} A - \frac{\beta \Delta t}{2} S$$, the system can be written as

$$L \mathbf{u}^{k+1} = R \mathbf{u}^k + \frac{\Delta t}{2} \left ( \mathbf{f}^{k+1} + \mathbf{f}^{k} \right )$$.

Because the matrix $$L$$ is time independent, it only needs to be inverted once for the whole simulation. The solution procedure is summarized by the following (pseudo) C algorithm:

```C
Vector run_simulation(const double t0, const double t1, const size_t nt)
{
    const double dt = (t1 - t0)/nt;

    const Matrix L = assemble_lhs_mat(dt);
    const Matrix R = assemble_rhs_mat(dt);
    Vector u0 = 0;

    Matrix_invert(L_inv);
    for (size_t step = 0; step < num_time_steps; step++) {
        Vector f = assemble_rhs_vec(dt);
        Vector u = L_inv*R*u0 + L_inv*f;

        u0 = u;
        t += dt;
    }

    return u;
}
```

Note that in the code the inverse is never formed explicitly, rather it relies on its (sparse) LU factorization to calculate its action on vectors and matrices.

## Software and libraries used
* GNU Scientific Library for dense linear algebra
* Suitesparse for sparse linear algebra
* cJSON for parsing the configuration file
* Inkscape for creating figures
* Fritzing for creating the circuit diagram
* Matplotlib for color schemes

## License
All code-related files are realsed under MIT license while all media (found in the assets directory) under Creative Commons license CC-BY-4.0. Please see the respective LICENSE files for terms of use.
