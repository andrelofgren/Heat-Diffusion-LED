#include "../include/mesh_1d.h"

static void Mesh1D_assemble_pvec(Mesh1D *mesh, const double x0, const double x1);
static void Mesh1D_assemble_cmat(Mesh1D *mesh);

void Mesh1D_init(
    Mesh1D *mesh,
    const double x0,
    const double x1,
    const size_t n,
    const size_t deg
) {
    const size_t num_cells = n;
    const size_t num_verts = deg*n + 1;
    const double dx = (x1 - x0)/(num_verts - 1);

    mesh->pvec = gsl_vector_alloc(num_verts);
    mesh->cmat = gsl_matrix_uint_alloc(num_cells, deg+1);
    Mesh1D_assemble_pvec(mesh, x0, dx);
    Mesh1D_assemble_cmat(mesh);
}

static void Mesh1D_assemble_pvec(Mesh1D *mesh, const double x0, const double dx)
{
    size_t N = Mesh1D_get_num_vertices(mesh);

    double x = x0;
    for (size_t i = 0; i < N; i++) {
        gsl_vector_set(mesh->pvec, i, x);
        x += dx;
    }
}

static void Mesh1D_assemble_cmat(Mesh1D *mesh)
{
    const size_t n = Mesh1D_get_num_cells(mesh);
    const size_t deg = Mesh1D_get_degree(mesh);
    for (size_t ci = 0; ci < n; ci++)
        for (size_t cj = 0; cj < deg+1; cj++)
            gsl_matrix_uint_set(mesh->cmat, ci, cj, ci*deg+cj);
}


size_t Mesh1D_get_num_vertices(const Mesh1D *mesh)
{
    return mesh->pvec->size;
}

size_t Mesh1D_get_num_cells(const Mesh1D *mesh)
{
    return mesh->cmat->size1;
}

size_t Mesh1D_get_degree(const Mesh1D *mesh)
{
    return mesh->cmat->size2-1;
}

size_t Mesh1D_get_num_vertices_in_cell(const Mesh1D *mesh) {
    return mesh->cmat->size2;
}

void Mesh1D_get_cell_vertex_ids(
    const Mesh1D *mesh, const size_t cell_id, gsl_vector_uint *cell_vertex_ids_ret
) {
    gsl_matrix_uint_get_row(cell_vertex_ids_ret, mesh->cmat, cell_id);
}

double Mesh1D_get_vertex_coordinate(const Mesh1D *mesh, const size_t vertex_id)
{
    return gsl_vector_get(mesh->pvec, vertex_id);
}

void Mesh1D_get_cell_vertex_coordinates(
    const Mesh1D *mesh,
    const size_t cell_id,
    gsl_vector *cell_vertex_coords_ret
) {

    const size_t num_verts_in_cell = mesh->cmat->size2;
    for (size_t i = 0; i < num_verts_in_cell; i++) {
        size_t vertex_id = gsl_matrix_uint_get(mesh->cmat, cell_id, i);
        double x = Mesh1D_get_vertex_coordinate(mesh, vertex_id);
        gsl_vector_set(cell_vertex_coords_ret, i, x);
    }
}

static int bin_search(
    const Mesh1D *mesh,
    gsl_vector *cell_vertex_coords,
    int ci,
    size_t ci_min,
    size_t ci_max,
    const double x
) {
    Mesh1D_get_cell_vertex_coordinates(mesh, ci, cell_vertex_coords);
    const double cell_x_min = gsl_vector_min(cell_vertex_coords);
    const double cell_x_max = gsl_vector_max(cell_vertex_coords);

    if (x > cell_x_max) {
        ci_min = ci;
        ci = ((ci_min+ci_max) >> 1);
        return bin_search(mesh, cell_vertex_coords, ci, ci_min, ci_max, x);
    } else if (x < cell_x_min) {
        ci_max = ci;
        ci = ((ci_min+ci_max) >> 1);
        return bin_search(mesh, cell_vertex_coords, ci, ci_min, ci_max, x);
    } else {
        return ci;
    }
}

void Mesh1D_find_cell_ids(
    const Mesh1D *mesh, const gsl_vector *points, gsl_vector_uint *cell_ids_ret
) {
    const double mesh_x_min = gsl_vector_min(mesh->pvec);
    const double mesh_x_max = gsl_vector_max(mesh->pvec);
    const size_t num_cells = Mesh1D_get_num_cells(mesh);
    const size_t num_verts_in_cell = Mesh1D_get_num_vertices_in_cell(mesh);
    gsl_vector *cell_vertex_coords = gsl_vector_alloc(num_verts_in_cell);

    for (size_t i = 0; i < points->size; i++) {
        double x = gsl_vector_get(points, i);
        size_t ci;
        if (x < mesh_x_min)
            ci = 0;
        else if (x > mesh_x_max)
            ci = num_cells-1;
        else {
            ci = bin_search(
                mesh, cell_vertex_coords, num_cells >> 1, 0, num_cells, x
            );
        }
        gsl_vector_uint_set(cell_ids_ret, i, ci);
    }

    gsl_vector_free(cell_vertex_coords);
}

void Mesh1D_finalize(Mesh1D *mesh)
{
    gsl_vector_free(mesh->pvec);
    gsl_matrix_uint_free(mesh->cmat);
}
