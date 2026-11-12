#include <suitesparse/cs.h>
#include <suitesparse/klu.h>

#include "../include/heat_problem_1d.h"
#include "../include/fem_1d.h"

static void assemble_mass_matrix(HeatProblem1D *hp);
static void assemble_mixed_matrix(HeatProblem1D *hp);
static void assemble_stiffness_matrix(HeatProblem1D *hp);

void HeatProblem1D_init(HeatProblem1D *hp, Mesh1D *mesh)
{
    hp->mesh = mesh;
    assemble_mass_matrix(hp);
    assemble_mixed_matrix(hp);
    assemble_stiffness_matrix(hp);
    hp->lhs_mat = NULL;
    hp->rhs_mat = NULL;
    hp->klu_symbolic = NULL;
    hp->klu_numeric = NULL;

    const size_t num_dofs = Mesh1D_get_num_vertices(hp->mesh);
    hp->sol_vec = gsl_vector_calloc(num_dofs);
    hp->source_vec = gsl_vector_calloc(num_dofs);
    hp->rhs_vec = gsl_vector_calloc(num_dofs);
}

void HeatProblem1D_finalize(HeatProblem1D *hp)
{
    cs_spfree(hp->stiff_mat);
    cs_spfree(hp->mixed_mat);
    cs_spfree(hp->mass_mat);

    HeatProblem1D_reset_lhs_mat(hp);
    HeatProblem1D_reset_rhs_mat(hp);
    HeatProblem1D_reset_klu(hp);

    gsl_vector_free(hp->source_vec);
    gsl_vector_free(hp->rhs_vec);
    gsl_vector_free(hp->sol_vec);
}

void HeatProblem1D_setup(HeatProblem1D *hp)
{
    HeatProblem1D_assemble_lhs_mat(hp);
    HeatProblem1D_assemble_rhs_mat(hp);
    HeatProblem1D_factorize_lhs_mat(hp);
}

static inline void get_element_coords(
    const gsl_vector_uint *element,
    const gsl_vector *pvec,
    gsl_vector *element_node_coords_ret
) {
    for (size_t i = 0; i < element->size; i++) {
        const size_t dof = gsl_vector_uint_get(element, i);
        const double x = gsl_vector_get(pvec, dof);
        gsl_vector_set(element_node_coords_ret, i, x);
    }
}

static void assemble_mass_matrix(HeatProblem1D *hp)
{
    // The polynomial degree of basis functions
    const size_t deg = Mesh1D_get_degree(hp->mesh);

    // Gauss precision
    const size_t gp_mass_mat = deg+1;

    // Allocate memory for quadrature
    gsl_vector *qpoints_r = gsl_vector_calloc(gp_mass_mat);
    gsl_vector *qpoints_x = gsl_vector_calloc(qpoints_r->size);
    gsl_vector *qweights = gsl_vector_calloc(gp_mass_mat);

    // Calculate quadrature points
    FEM1D_gauss_legendre_quadrature(gp_mass_mat, qpoints_r, qweights);

    // Allocate memory for basis functions
    const size_t num_dofs_in_cell = hp->mesh->cmat->size2;
    gsl_matrix *phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_matrix *grad_phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_vector *detJ_r = gsl_vector_calloc(qpoints_r->size);
    gsl_matrix *grad_phi_x = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);

    // Calculate reference cell basis functions
    FEM1D_lagrange_basis(deg, qpoints_r, phi_r, grad_phi_r);

    // Allocate memory for storing cell vertex information:
    gsl_vector_uint *element = gsl_vector_uint_calloc(num_dofs_in_cell);  // Connectivity
    gsl_vector *element_coords = gsl_vector_calloc(num_dofs_in_cell);  // Vertex coordinates

    // Matrix for storing local contributions
    gsl_matrix *M_block = gsl_matrix_calloc(num_dofs_in_cell, num_dofs_in_cell);

    const size_t num_dofs = Mesh1D_get_num_vertices(hp->mesh);
    const size_t max_nnz = (2*deg + 1)*num_dofs;

    // Sparse matrix coordinate format
    cs_di *lhs_triplet = cs_di_spalloc(num_dofs, num_dofs, max_nnz, true, true);

    // Assembly:
    const size_t num_cells = hp->mesh->cmat->size1;
    // Loop over mesh cells
    for (size_t ci = 0; ci < num_cells; ci++) {
        gsl_matrix_uint_get_row(element, hp->mesh->cmat, ci);
        get_element_coords(element, hp->mesh->pvec, element_coords);
        FEM1D_map_from_reference_cell(
            element_coords, phi_r, grad_phi_r, detJ_r, qpoints_x, grad_phi_x
        );
        // Loop over quadrature points
        for (size_t q = 0; q < qpoints_r->size; q++) {
            const double detJ = gsl_vector_get(detJ_r, q);
            const double Wq = gsl_vector_get(qweights, q);
            const double detJxW = detJ*Wq;
            // Calculate cell local matrix contributions
            for (size_t li = 0; li < element->size; li++) {
                const double phi_li = gsl_matrix_get(phi_r, li, q);
                for (size_t lj = 0; lj < element->size; lj++) {
                    const double phi_lj = gsl_matrix_get(phi_r, lj, q);
                    const double val = gsl_matrix_get(M_block, li, lj) + phi_li*phi_lj*detJxW;
                    gsl_matrix_set(M_block, li, lj, val);
                }
            }
        }

        // Insert values into the global matrix
        for (size_t li = 0; li < element->size; li++) {
            size_t gi = gsl_vector_uint_get(element, li);
            for (size_t lj = 0; lj < element->size; lj++) {
                size_t gj = gsl_vector_uint_get(element, lj);
                double val = gsl_matrix_get(M_block, li, lj);
                cs_di_entry(lhs_triplet, gi, gj, val);
            }
        }
        gsl_matrix_set_zero(M_block);
    }

    // Compress triplet matrix into CSC format
    hp->mass_mat = cs_compress(lhs_triplet);
    // Add up duplicate entries
    cs_dupl(hp->mass_mat);
    // We may now release the triplet matrix
    cs_di_spfree(lhs_triplet);

    gsl_matrix_free(M_block);
    gsl_vector_free(element_coords);
    gsl_vector_uint_free(element);
    gsl_matrix_free(grad_phi_x);
    gsl_vector_free(qpoints_x);
    gsl_vector_free(detJ_r);
    gsl_matrix_free(grad_phi_r);
    gsl_matrix_free(phi_r);
    gsl_vector_free(qpoints_r);
    gsl_vector_free(qweights);
}

static void assemble_mixed_matrix(HeatProblem1D *hp)
{
    // The polynomial degree of basis functions
    const size_t deg = Mesh1D_get_degree(hp->mesh);

    // Gauss precision
    const size_t gp_mixed_mat = deg+1;

    // Allocate memory for quadrature
    gsl_vector *qpoints_r = gsl_vector_calloc(gp_mixed_mat);
    gsl_vector *qpoints_x = gsl_vector_calloc(qpoints_r->size);
    gsl_vector *qweights = gsl_vector_calloc(gp_mixed_mat);

    // Calculate quadrature points
    FEM1D_gauss_legendre_quadrature(gp_mixed_mat, qpoints_r, qweights);

    // Allocate memory for basis functions
    const size_t num_dofs_in_cell = hp->mesh->cmat->size2;
    gsl_matrix *phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_matrix *grad_phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_vector *detJ_r = gsl_vector_calloc(qpoints_r->size);
    gsl_matrix *grad_phi_x = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);

    // Calculate reference cell basis functions
    FEM1D_lagrange_basis(deg, qpoints_r, phi_r, grad_phi_r);

    // Allocate memory for storing cell vertex information:
    gsl_vector_uint *element = gsl_vector_uint_calloc(num_dofs_in_cell);  // Connectivity
    gsl_vector *element_coords = gsl_vector_calloc(num_dofs_in_cell);  // Vertex coordinates

    // Matrix for storing local contributions
    gsl_matrix *S_block = gsl_matrix_calloc(num_dofs_in_cell, num_dofs_in_cell);

    const size_t num_dofs = Mesh1D_get_num_vertices(hp->mesh);
    const size_t max_nnz = (2*deg + 1)*num_dofs;

    // Sparse matrix coordinate format
    cs_di *lhs_triplet = cs_di_spalloc(num_dofs, num_dofs, max_nnz, true, true);

    // Assembly:
    const size_t num_cells = hp->mesh->cmat->size1;
    // Loop over mesh cells
    for (size_t ci = 0; ci < num_cells; ci++) {
        gsl_matrix_uint_get_row(element, hp->mesh->cmat, ci);
        get_element_coords(element, hp->mesh->pvec, element_coords);
        FEM1D_map_from_reference_cell(
            element_coords, phi_r, grad_phi_r, detJ_r, qpoints_x, grad_phi_x
        );
        // Loop over quadrature points
        for (size_t q = 0; q < qpoints_r->size; q++) {
            const double detJ = gsl_vector_get(detJ_r, q);
            const double Wq = gsl_vector_get(qweights, q);
            const double detJxW = detJ*Wq;
            // Calculate cell local matrix contributions
            for (size_t li = 0; li < element->size; li++) {
                const double phi_li = gsl_matrix_get(phi_r, li, q);
                for (size_t lj = 0; lj < element->size; lj++) {
                    const double ddx_phi_lj = gsl_matrix_get(grad_phi_x, lj, q);
                    const double val = gsl_matrix_get(S_block, li, lj) + phi_li*ddx_phi_lj*detJxW;
                    gsl_matrix_set(S_block, li, lj, val);
                }
            }
        }

        // Insert values into the global matrix
        for (size_t li = 0; li < element->size; li++) {
            const size_t gi = gsl_vector_uint_get(element, li);
            for (size_t lj = 0; lj < element->size; lj++) {
                const size_t gj = gsl_vector_uint_get(element, lj);
                const double val = gsl_matrix_get(S_block, li, lj);
                cs_di_entry(lhs_triplet, gi, gj, val);
            }
        }
        gsl_matrix_set_zero(S_block);
    }

    // Compress triplet matrix into CSC format
    hp->mixed_mat = cs_compress(lhs_triplet);
    // Add up duplicate entries
    cs_dupl(hp->mixed_mat);
    // We may now release the triplet matrix
    cs_di_spfree(lhs_triplet);

    gsl_matrix_free(S_block);
    gsl_vector_free(element_coords);
    gsl_vector_uint_free(element);
    gsl_matrix_free(grad_phi_x);
    gsl_vector_free(qpoints_x);
    gsl_vector_free(detJ_r);
    gsl_matrix_free(grad_phi_r);
    gsl_matrix_free(phi_r);
    gsl_vector_free(qpoints_r);
    gsl_vector_free(qweights);
}

static void assemble_stiffness_matrix(HeatProblem1D *hp)
{
    // The polynomial degree of basis functions
    const size_t deg = Mesh1D_get_degree(hp->mesh);

    // Gauss precision
    const size_t gp_stiff_mat = deg;

    // Allocate memory for quadrature
    gsl_vector *qpoints_r = gsl_vector_calloc(gp_stiff_mat);
    gsl_vector *qpoints_x = gsl_vector_calloc(qpoints_r->size);
    gsl_vector *qweights = gsl_vector_calloc(gp_stiff_mat);

    // Calculate quadrature points
    FEM1D_gauss_legendre_quadrature(gp_stiff_mat, qpoints_r, qweights);

    // Allocate memory for basis functions
    const size_t num_dofs_in_cell = hp->mesh->cmat->size2;
    gsl_matrix *phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_matrix *grad_phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_vector *detJ_r = gsl_vector_calloc(qpoints_r->size);
    gsl_matrix *grad_phi_x = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);

    // Calculate reference cell basis functions
    FEM1D_lagrange_basis(deg, qpoints_r, phi_r, grad_phi_r);

    // Allocate memory for storing cell vertex information:
    gsl_vector_uint *element = gsl_vector_uint_calloc(num_dofs_in_cell);  // Connectivity
    gsl_vector *element_coords = gsl_vector_calloc(num_dofs_in_cell);  // Vertex coordinates

    // Matrix for storing local contributions
    gsl_matrix *A_block = gsl_matrix_calloc(num_dofs_in_cell, num_dofs_in_cell);

    const size_t num_dofs = Mesh1D_get_num_vertices(hp->mesh);
    const size_t max_nnz = (2*deg + 1)*num_dofs;

    // Sparse matrix coordinate format
    cs_di *lhs_triplet = cs_di_spalloc(num_dofs, num_dofs, max_nnz, true, true);

    // Assembly:
    const size_t num_cells = hp->mesh->cmat->size1;
    // Loop over mesh cells
    for (size_t ci = 0; ci < num_cells; ci++) {
        gsl_matrix_uint_get_row(element, hp->mesh->cmat, ci);
        get_element_coords(element, hp->mesh->pvec, element_coords);
        FEM1D_map_from_reference_cell(
            element_coords, phi_r, grad_phi_r, detJ_r, qpoints_x, grad_phi_x
        );
        // Loop over quadrature points
        for (size_t q = 0; q < qpoints_r->size; q++) {
            const double detJ = gsl_vector_get(detJ_r, q);
            const double Wq = gsl_vector_get(qweights, q);
            const double detJxW = detJ*Wq;
            // Calculate cell local matrix contributions
            for (size_t li = 0; li < element->size; li++) {
                const double ddx_phi_li = gsl_matrix_get(grad_phi_x, li, q);
                for (size_t lj = 0; lj < element->size; lj++) {
                    const double ddx_phi_lj = gsl_matrix_get(grad_phi_x, lj, q);
                    const double val =
                        gsl_matrix_get(A_block, li, lj) + ddx_phi_li*ddx_phi_lj*detJxW;
                    gsl_matrix_set(A_block, li, lj, val);
                }
            }
        }

        // Insert values into the global matrix
        for (size_t li = 0; li < element->size; li++) {
            const size_t gi = gsl_vector_uint_get(element, li);
            for (size_t lj = 0; lj < element->size; lj++) {
                const size_t gj = gsl_vector_uint_get(element, lj);
                const double val = gsl_matrix_get(A_block, li, lj);
                cs_di_entry(lhs_triplet, gi, gj, val);
            }
        }
        gsl_matrix_set_zero(A_block);
    }

    // Compress triplet matrix into CSC format
    hp->stiff_mat = cs_compress(lhs_triplet);
    // Add up duplicate entries
    cs_dupl(hp->stiff_mat);
    // We may now release the triplet matrix
    cs_di_spfree(lhs_triplet);

    gsl_matrix_free(A_block);
    gsl_vector_free(element_coords);
    gsl_vector_uint_free(element);
    gsl_matrix_free(grad_phi_x);
    gsl_vector_free(qpoints_x);
    gsl_vector_free(detJ_r);
    gsl_matrix_free(grad_phi_r);
    gsl_matrix_free(phi_r);
    gsl_vector_free(qpoints_r);
    gsl_vector_free(qweights);
}

// Assemble the LHS matrix as obtained by the Crank-Nicolson discretization
void HeatProblem1D_assemble_lhs_mat(HeatProblem1D *hp)
{
    const double alpha = hp->alpha;
    const double beta = hp->beta;
    const double gamma = hp->gamma;
    const double dt = hp->dt;
    cs_di *tmp_mat = cs_di_add(hp->stiff_mat, hp->mixed_mat, 0.5*alpha*dt, 0.5*beta*dt);
    hp->lhs_mat = cs_di_add(
        hp->mass_mat, tmp_mat, 1.0 + 0.5*gamma*dt, 1.0
    );
    cs_di_spfree(tmp_mat);
}

// Assemble the RHS matrix as obtained by the Crank-Nicolson discretization
void HeatProblem1D_assemble_rhs_mat(HeatProblem1D *hp)
{
    const double alpha = hp->alpha;
    const double beta = hp->beta;
    const double gamma = hp->gamma;
    const double dt = hp->dt;
    cs_di *tmp_mat = cs_di_add(hp->stiff_mat, hp->mixed_mat, -0.5*alpha*dt, -0.5*beta*dt);
    hp->rhs_mat = cs_di_add(
        hp->mass_mat, tmp_mat, 1.0 - 0.5*gamma*dt, 1.0
    );
    cs_di_spfree(tmp_mat);
}

// Assemble the RHS vector as obtained by the Crank-Nicolson discretization
void HeatProblem1D_assemble_rhs_vec(HeatProblem1D *hp)
{
    const double t = hp->t;
    const double dt = hp->dt;

    // The polynomial degree of basis functions
    const size_t deg = Mesh1D_get_degree(hp->mesh);

    // Gauss precision
    const size_t gp_mass_mat = deg+1;

    // Allocate memory for quadrature
    gsl_vector *qpoints_r = gsl_vector_calloc(gp_mass_mat);
    gsl_vector *qpoints_x = gsl_vector_calloc(qpoints_r->size);
    gsl_vector *qweights = gsl_vector_calloc(gp_mass_mat);
    FEM1D_gauss_legendre_quadrature(gp_mass_mat, qpoints_r, qweights);

    // Allocate memory for basis functions
    const size_t num_dofs_in_cell = hp->mesh->cmat->size2;
    gsl_matrix *phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_matrix *grad_phi_r = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);
    gsl_vector *detJ_r = gsl_vector_calloc(qpoints_r->size);
    gsl_matrix *grad_phi_x = gsl_matrix_calloc(num_dofs_in_cell, qpoints_r->size);

    // Calculate reference cell basis functions
    FEM1D_lagrange_basis(deg, qpoints_r, phi_r, grad_phi_r);

    // Allocate memory for storing cell vertex information:
    // Connectivity
    gsl_vector_uint *element = gsl_vector_uint_calloc(num_dofs_in_cell);
    // Vertex coordinates
    gsl_vector *element_coords = gsl_vector_calloc(num_dofs_in_cell);

    // Matrix for storing local contributions
    gsl_vector *f_block = gsl_vector_calloc(num_dofs_in_cell);

    // Assembly loop
    // Loop over mesh cells
    const size_t num_cells = hp->mesh->cmat->size1;
    for (size_t ci = 0; ci < num_cells; ci++) {
        gsl_matrix_uint_get_row(element, hp->mesh->cmat, ci);
        get_element_coords(element, hp->mesh->pvec, element_coords);
        FEM1D_map_from_reference_cell(
            element_coords, phi_r, grad_phi_r, detJ_r, qpoints_x, grad_phi_x
        );
        // Loop over quadrature points
        for (size_t q = 0; q < qpoints_r->size; q++) {
            const double detJ = gsl_vector_get(detJ_r, q);
            const double Wq = gsl_vector_get(qweights, q);
            const double xq = gsl_vector_get(qpoints_x, q);
            const double f0_val = hp->source(xq, t);
            const double f1_val = hp->source(xq, t+dt);
            const double dtxfxdetJxW = dt*0.5*(f0_val+f1_val)*detJ*Wq;
            // Calculate cell local vector contributions
            for (size_t li = 0; li < element->size; li++) {
                const double phi_li = gsl_matrix_get(phi_r, li, q);
                const double val =
                    gsl_vector_get(f_block, li) + dtxfxdetJxW*phi_li;
                gsl_vector_set(f_block, li, val);
            }
        }

        // Insert values into the global vector
        for (size_t li = 0; li < element->size; li++) {
            const size_t gi = gsl_vector_uint_get(element, li);
            const double val = gsl_vector_get(hp->source_vec, gi) + gsl_vector_get(f_block, li);
            gsl_vector_set(hp->source_vec, gi, val);
        }
        gsl_vector_set_zero(f_block);
    }

    // Calculate the rhs vector: rhs_vec = rhs_mat*u0 + f_vec
    cs_di_gaxpy(hp->rhs_mat, hp->sol_vec->data, hp->rhs_vec->data);
    gsl_vector_add(hp->rhs_vec, hp->source_vec);

    // Finalize
    gsl_vector_free(qweights);
    gsl_vector_free(qpoints_r);
    gsl_matrix_free(phi_r);
    gsl_matrix_free(grad_phi_r);
    gsl_vector_free(detJ_r);
    gsl_vector_free(qpoints_x);
    gsl_matrix_free(grad_phi_x);
    gsl_vector_uint_free(element);
    gsl_vector_free(element_coords);
    gsl_vector_free(f_block);
}

// Apply symmetric homogenous Dirichlet BC to LHS matrix
void HeatProblem1D_apply_zero_dirichlet_bc_lhs_mat(HeatProblem1D *hp)
{
    const size_t num_nnz = hp->lhs_mat->nzmax;
    hp->lhs_mat->x[0] = 1.0;
    hp->lhs_mat->x[1] = 0.0;
    hp->lhs_mat->x[2] = 0.0;
    hp->lhs_mat->x[num_nnz-3] = 0.0;
    hp->lhs_mat->x[num_nnz-2] = 0.0;
    hp->lhs_mat->x[num_nnz-1] = 1.0;
}

// Apply symmetric homogenous Dirichlet BC to RHS vector
void HeatProblem1D_apply_zero_dirichlet_bc_rhs_vec(HeatProblem1D *hp)
{
    const size_t num_dofs = hp->rhs_vec->size;
    hp->rhs_vec->data[0] = 0.0;
    hp->rhs_vec->data[num_dofs-1] = 0.0;
}

// Perform LU factorization on LHS matrix
void HeatProblem1D_factorize_lhs_mat(HeatProblem1D *hp)
{
    // Analyze the sparsity pattern of the matrix
    hp->klu_symbolic = klu_analyze(
        hp->lhs_mat->n, hp->lhs_mat->p, hp->lhs_mat->i, &hp->klu_common
    );
    // Perform the actual LU factorization
    hp->klu_numeric = klu_factor(
        hp->lhs_mat->p, hp->lhs_mat->i, hp->lhs_mat->x, hp->klu_symbolic, &hp->klu_common
    );
}

// Solve linear system using sparse LU factorization
void HeatProblem1D_solve_linsys(HeatProblem1D *hp)
{
    const size_t num_rhs = 1;
    klu_solve(
        hp->klu_symbolic,
        hp->klu_numeric,
        hp->rhs_vec->size,
        num_rhs,
        hp->rhs_vec->data,
        &hp->klu_common
    );
    gsl_vector_memcpy(hp->sol_vec, hp->rhs_vec);
}

// Advance the system in time
void HeatProblem1D_step(HeatProblem1D *hp)
{
    HeatProblem1D_assemble_rhs_vec(hp);
    HeatProblem1D_solve_linsys(hp);
}

void HeatProblem1D_reset_lhs_mat(HeatProblem1D *hp)
{
    if (hp->lhs_mat)
        cs_di_spfree(hp->lhs_mat);
    hp->lhs_mat = NULL;
}

void HeatProblem1D_reset_rhs_mat(HeatProblem1D *hp)
{
    if (hp->rhs_mat)
        cs_di_spfree(hp->rhs_mat);
    hp->rhs_mat = NULL;
}

void HeatProblem1D_reset_rhs(HeatProblem1D *hp)
{
    gsl_vector_set_zero(hp->source_vec);
    gsl_vector_set_zero(hp->rhs_vec);
}

void HeatProblem1D_reset_klu(HeatProblem1D *hp)
{
    if (hp->klu_numeric) {
        klu_free_numeric(&hp->klu_numeric, &hp->klu_common);
        hp->klu_numeric = NULL;
    }
    if (hp->klu_symbolic) {
        klu_free_symbolic(&hp->klu_symbolic, &hp->klu_common);
        hp->klu_symbolic = NULL;
    }
}

void HeatProblem1D_reset_solution(HeatProblem1D *hp)
{
    gsl_vector_set_zero(hp->sol_vec);
}

void HeatProblem1D_reset_system(HeatProblem1D *hp)
{
    HeatProblem1D_reset_lhs_mat(hp);
    HeatProblem1D_reset_rhs_mat(hp);
    HeatProblem1D_reset_rhs(hp);
    HeatProblem1D_reset_klu(hp);
    HeatProblem1D_reset_solution(hp);
}
