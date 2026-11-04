#ifndef __FEM_FUNCTION_1D_H__
#define __FEM_FUNCTION_1D_H__

#include "mesh_1d.h"

typedef struct FEMFuction1D {
    const Mesh1D *mesh;
    gsl_vector *vals;
} FEMFunction1D;

void FEMFunction1D_init(FEMFunction1D *fem_func, const Mesh1D *mesh);
void FEMFunction1D_finalize(FEMFunction1D *fem_func);

void FEMFunction1D_set_vals_from_array(
    FEMFunction1D *fem_func, const double *array, const size_t size
);
void FEMFunction1D_set_vals_from_func(FEMFunction1D *fem_func, double (*func)(double));
void FEMFunction1D_set_vals_from_vec(FEMFunction1D *fem_func, const gsl_vector *vec);

void FEMFunction1D_eval(
    const FEMFunction1D *fem_func, const gsl_vector *points_x, gsl_vector *values_ret
);

size_t FEMFunction1D_get_num_dofs(const FEMFunction1D *fem_func);
double FEMFunction1D_get_dof_coordinate(
    const FEMFunction1D *fem_func, const size_t dof_id
);

#endif  // __FEM_FUNCTION_1D_H__
