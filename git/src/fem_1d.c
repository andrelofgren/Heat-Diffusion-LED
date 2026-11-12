#include "../include/fem_1d.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>

void FEM1D_map_from_reference_cell(
    const gsl_vector *element_coords,
    const gsl_matrix *phi_r,
    const gsl_matrix *grad_phi_r,
    gsl_vector *detJ_r_ret,
    gsl_vector *qpoints_x_ret,
    gsl_matrix *grad_phi_x_ret
) {
    double dF;  // Jacobian
    for (size_t q = 0; q < qpoints_x_ret->size; q++) {
        // Quadrature points
        gsl_vector_const_view phi_r_col_q = gsl_matrix_const_column(phi_r, q);
        gsl_blas_ddot(&phi_r_col_q.vector, element_coords, &qpoints_x_ret->data[q]);

        // Jacobian determinant
        gsl_vector_const_view grad_phi_r_col_q = gsl_matrix_const_column(grad_phi_r, q);
        gsl_blas_ddot(&grad_phi_r_col_q.vector, element_coords, &dF);
        gsl_vector_set(detJ_r_ret, q, dF);

        // Gradient of basis functions
        gsl_vector_view grad_phi_x_col_q = gsl_matrix_column(grad_phi_x_ret, q);
        gsl_vector_axpby(1.0/dF, &grad_phi_r_col_q.vector, 0.0, &grad_phi_x_col_q.vector);
    }
}

static void get_reference_element_vertex_coordinates(gsl_vector *nodes_ret)
{
    // Calculate vertex coordinates of the reference element
    double r = 0.0;
    double dr = 1.0/(nodes_ret->size-1);
    for (size_t i = 0; i < nodes_ret->size; i++) {
        gsl_vector_set(nodes_ret, i, r);
        r += dr;
    }
}

void FEM1D_lagrange_basis(
    const size_t deg,
    const gsl_vector *qpoints_r,
    gsl_matrix *phi_r_ret,
    gsl_matrix *grad_phi_r_ret
) {
    // Assemble the basis function coefficient matrix, R
    gsl_matrix *R = gsl_matrix_alloc(deg+1, deg+1);
    gsl_vector *verts_ref_element = gsl_vector_alloc(deg+1);
    get_reference_element_vertex_coordinates(verts_ref_element);
    for (size_t i = 0; i < R->size1; i++) {
        double r = gsl_vector_get(verts_ref_element, i);
        double R_ji = 1.0;
        for (int j = 0; j < R->size2; j++) {
            gsl_matrix_set(R, j, i, R_ji);
            R_ji *= r;
        }
    }

    // Initialize LU decomp. to invert for coefficients
    int signum;
    gsl_permutation *P = gsl_permutation_alloc(R->size1);
    gsl_linalg_LU_decomp(R, P, &signum);

    // Calculate phi and grad_phi at quadrature points, where:
    // R * phi = [1, r, r^2, ..., r^d]^T,
    // R * grad_phi = [0, 1, 2r, ..., d r^d-1]^T
    gsl_vector *r_vec = gsl_vector_alloc(deg+1);  // Vector [1, r, r^2, ..., r^d]^T
    gsl_vector *dr_vec = gsl_vector_alloc(deg+1);  // Vector [0, 1, 2r, ..., d r^d-1]^T 
    for (size_t q = 0; q < qpoints_r->size; q++) {
        double r = gsl_vector_get(qpoints_r, q);
        for (size_t i = 0; i < deg+1; i++) {
            double r_val, dr_val;
            if (i == 0) {
                r_val = 1.0;
                dr_val = 0.0;
            } else {
                // Obtain r_vec and dr_vec through recurrence relation
                r_val = r*gsl_vector_get(r_vec, i-1);
                dr_val = i*gsl_vector_get(r_vec, i-1);
            }
            gsl_vector_set(r_vec, i, r_val);
            gsl_vector_set(dr_vec, i, dr_val);
        }
        gsl_vector_view phi_r_col = gsl_matrix_column(phi_r_ret, q);
        gsl_vector_view grad_phi_r_col = gsl_matrix_column(grad_phi_r_ret, q);

        // Invert for phi and grad_phi at point r
        gsl_linalg_LU_solve(R, P, r_vec, &phi_r_col.vector);
        gsl_linalg_LU_solve(R, P, dr_vec, &grad_phi_r_col.vector);
    }

    // Finalize
    gsl_vector_free(dr_vec);
    gsl_vector_free(r_vec);
    gsl_permutation_free(P);
    gsl_vector_free(verts_ref_element);
    gsl_matrix_free(R);
}

void FEM1D_gauss_legendre_quadrature(
    const size_t precision, gsl_vector *points_ret, gsl_vector *weights_ret
) {
    // Gauss-Legendre quadrature on [0, 1]
    gsl_integration_glfixed_table *gauss_table =
        gsl_integration_glfixed_table_alloc(precision);
    for (size_t i = 0; i < points_ret->size; i++)
        gsl_integration_glfixed_point(
            0.0, 1.0, i, &points_ret->data[i], &weights_ret->data[i], gauss_table
        );
    gsl_integration_glfixed_table_free(gauss_table);
}
