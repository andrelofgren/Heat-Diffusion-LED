#ifndef __FEM_1D_H__
#define __FEM_1D_H__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void FEM1D_map_from_reference_cell(
    const gsl_vector *vertex_coords,
    const gsl_matrix *phi_r,
    const gsl_matrix *grad_phi_r,
    gsl_vector *detJ_r_ret,
    gsl_vector *qpoints_x_ret,
    gsl_matrix *grad_phi_x_ret
);

void FEM1D_map_points_to_reference_cell(
    const gsl_vector *points_x,
    const gsl_vector *node_coords,
    gsl_vector *points_r_ret
);

void FEM1D_lagrange_basis(
    const size_t deg,
    const gsl_vector *qpoints_r,
    gsl_matrix *phi_r_ret,
    gsl_matrix *grad_phi_r_ret
);

void FEM1D_gauss_legendre_quadrature(
    const size_t precision,
    gsl_vector *points_ret,
    gsl_vector *weights_ret
);

#endif  // __FEM_1D_H__
