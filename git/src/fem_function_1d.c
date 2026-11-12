#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_blas.h>

#include "../include/fem_function_1d.h"
#include "../include/fem_1d.h"

void FEMFunction1D_init(FEMFunction1D *fem_func, const Mesh1D *mesh)
{
    fem_func->mesh = mesh;
    const size_t num_dofs = Mesh1D_get_num_vertices(fem_func->mesh);
    fem_func->vals = gsl_vector_calloc(num_dofs);
}

void FEMFunction1D_finalize(FEMFunction1D *fem_func)
{
    gsl_vector_free(fem_func->vals);
}

void FEMFunction1D_set_vals_from_vec(FEMFunction1D *fem_func, const gsl_vector *vec)
{
    gsl_vector_memcpy(fem_func->vals, vec);
}

void FEMFunction1D_set_vals_from_func(FEMFunction1D *fem_func, double (*func)(double))
{
    const size_t num_dofs = FEMFunction1D_get_num_dofs(fem_func);
    for (size_t dof = 0; dof < num_dofs; dof++) {
        const double x = FEMFunction1D_get_dof_coordinate(fem_func, dof);
        const double fval = func(x);
        gsl_vector_set(fem_func->vals, dof, fval);
    }
}

void FEMFunction1D_set_vals_from_array(
    FEMFunction1D *fem_func, const double *array, const size_t size
) {
    memcpy(fem_func->vals->data, array, size * sizeof(double));
}

// Evaluate the fem function on points inside specified cell
static void eval_cell(
    const FEMFunction1D *fem_func,
    const size_t cell_id,
    const gsl_vector *cell_points_x,
    gsl_vector *cell_points_r,
    gsl_vector_uint *cell_vertex_ids,
    gsl_vector *cell_vertex_coords,
    gsl_matrix *phi_r,
    gsl_matrix *grad_phi_r,
    gsl_vector *f_cell,
    gsl_vector *cell_values_ret
) {
    // Get cell vertices
    Mesh1D_get_cell_vertex_coordinates(fem_func->mesh, cell_id, cell_vertex_coords);

    // Get cell indices
    Mesh1D_get_cell_vertex_ids(fem_func->mesh, cell_id, cell_vertex_ids);

    // Calculate cell extent
    const double x_min = gsl_vector_min(cell_vertex_coords);
    const double x_max = gsl_vector_max(cell_vertex_coords);

    // Map points_x in cell to point_r in reference cell, according to
    // r = (x - x_min)/(x_max - x_min)
    gsl_vector_memcpy(cell_points_r, cell_points_x);
    gsl_vector_add_constant(cell_points_r, -x_min);
    gsl_vector_scale(cell_points_r, 1.0/(x_max - x_min));

    // Evaluate lagrange basis on points
    const size_t deg = Mesh1D_get_degree(fem_func->mesh);
    FEM1D_lagrange_basis(deg, cell_points_r, phi_r, grad_phi_r);

    // Get the fem function values on cell
    for (size_t i = 0; i < f_cell->size; i++) {
        size_t dof = gsl_vector_uint_get(cell_vertex_ids, i);
        double f_val = gsl_vector_get(fem_func->vals, dof);
        gsl_vector_set(f_cell, i, f_val);
    }
    // Evaluate fem function on points: f(x) = phi(r)^T*f
    gsl_blas_dgemv(CblasTrans, 1.0, phi_r, f_cell, 0.0, cell_values_ret);
}

// Evaluate the fem function on a sorted set of points
void FEMFunction1D_eval(
    const FEMFunction1D *fem_func, const gsl_vector *points_x, gsl_vector *values_ret
) {
    const size_t num_dofs_in_cell = Mesh1D_get_num_vertices_in_cell(fem_func->mesh);
    gsl_vector_uint *cell_vertex_ids = gsl_vector_uint_alloc(num_dofs_in_cell);
    gsl_vector *cell_vertex_coords = gsl_vector_alloc(num_dofs_in_cell);

    gsl_vector_uint *point_cell_ids = gsl_vector_uint_alloc(points_x->size);
    Mesh1D_find_cell_ids(fem_func->mesh, points_x, point_cell_ids);

    // All points may be in a single cell
    const size_t max_num_points_in_cell = points_x->size;

    gsl_vector *cell_points_x = gsl_vector_alloc(max_num_points_in_cell);
    gsl_vector *cell_points_r = gsl_vector_alloc(max_num_points_in_cell);

    gsl_matrix *phi_r = gsl_matrix_alloc(num_dofs_in_cell, max_num_points_in_cell);
    gsl_matrix *grad_phi_r = gsl_matrix_alloc(num_dofs_in_cell, max_num_points_in_cell);

    gsl_vector *f_cell = gsl_vector_alloc(num_dofs_in_cell);

    // Loop over all points
    size_t k = 0;
    while (k < points_x->size) {
        // Get the cell index of the point
        size_t ci = gsl_vector_uint_get(point_cell_ids, k);

        // Find all points inside the cell, so we for efficiency
        // reasons may evaluate them all at once
        size_t num_points_in_cell = 0;
        do {
            double x = gsl_vector_get(points_x, k);
            gsl_vector_set(cell_points_x, num_points_in_cell, x);
            num_points_in_cell += 1;
            k++;
        } while (k < points_x->size && gsl_vector_uint_get(point_cell_ids, k) == ci);

        cell_points_x->size = num_points_in_cell;
        cell_points_r->size = num_points_in_cell;
        phi_r->size2 = num_points_in_cell;
        grad_phi_r->size2 = num_points_in_cell;

        // Create a vector view for cell local function values
        gsl_vector_view cell_vals = gsl_vector_subvector(
            values_ret,
            k-num_points_in_cell,
            num_points_in_cell
        );

        // Evaluate the fem function at all points inside the cell
        eval_cell(
            fem_func, ci, cell_points_x, cell_points_r,
            cell_vertex_ids, cell_vertex_coords,
            phi_r, grad_phi_r, f_cell, &cell_vals.vector
        );

        cell_points_x->size = max_num_points_in_cell;
        cell_points_r->size = max_num_points_in_cell;
        phi_r->size2 = max_num_points_in_cell;
        grad_phi_r->size2 = max_num_points_in_cell;
    }

    // Finalize
    gsl_vector_free(f_cell);
    gsl_matrix_free(grad_phi_r);
    gsl_matrix_free(phi_r);
    gsl_vector_free(cell_points_r);
    gsl_vector_free(cell_points_x);
    gsl_vector_uint_free(point_cell_ids);
    gsl_vector_free(cell_vertex_coords);
    gsl_vector_uint_free(cell_vertex_ids);
}

size_t FEMFunction1D_get_num_dofs(const FEMFunction1D *fem_func)
{
    return Mesh1D_get_num_vertices(fem_func->mesh);
}

double FEMFunction1D_get_dof_coordinate(const FEMFunction1D *fem_func, size_t const dof_id)
{
    return Mesh1D_get_vertex_coordinate(fem_func->mesh, dof_id);
}
