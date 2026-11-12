#ifndef __MESH_1D_H__
#define __MESH_1D_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct Mesh1D {
    gsl_vector *pvec;
    gsl_matrix_uint *cmat;
} Mesh1D;

void Mesh1D_init(
    Mesh1D *mesh,
    const double x0,
    const double x1,
    const size_t n,
    const size_t deg
);
void Mesh1D_finalize(Mesh1D *mesh);

size_t Mesh1D_get_degree(const Mesh1D *mesh);
size_t Mesh1D_get_num_cells(const Mesh1D *mesh);
size_t Mesh1D_get_num_vertices(const Mesh1D *mesh);
size_t Mesh1D_get_num_vertices_in_cell(const Mesh1D *mesh);
double Mesh1D_get_vertex_coordinate(const Mesh1D *mesh, const size_t vertex_id);
void Mesh1D_find_cell_ids(
    const Mesh1D *mesh, const gsl_vector *points, gsl_vector_uint *cell_ids_ret
);
void Mesh1D_get_cell_vertex_ids(
    const Mesh1D *mesh, const size_t cell_id, gsl_vector_uint *cell_vertex_ids_ret
);
void Mesh1D_get_cell_vertex_coordinates(
    const Mesh1D *mesh,
    const size_t cell_id,
    gsl_vector *cell_vertex_coords_ret
);

#endif  // __MESH_1D_H__
