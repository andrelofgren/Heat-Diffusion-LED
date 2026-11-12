#ifndef __HEAT_PROBLEM_H__
#define __HEAT_PROBLEM_H__

#include "../include/mesh_1d.h"
#include <suitesparse/cs.h>
#include <suitesparse/klu.h>

typedef struct HeatProblem1D {
    Mesh1D *mesh;
    double alpha;
    double beta;
    double gamma;
    double (*source)(double, double);
    double t;
    double dt;
    cs_di *mass_mat;
    cs_di *mixed_mat;
    cs_di *stiff_mat;
    cs_di *lhs_mat;
    cs_di *rhs_mat;
    gsl_vector *source_vec;
    gsl_vector *rhs_vec;
    gsl_vector *sol_vec;
    klu_common klu_common;
    klu_symbolic *klu_symbolic;
    klu_numeric *klu_numeric;
} HeatProblem1D;

void HeatProblem1D_init(HeatProblem1D *hp, Mesh1D *mesh);
void HeatProblem1D_finalize(HeatProblem1D *hp);

void HeatProblem1D_setup(HeatProblem1D *hp);

void HeatProblem1D_assemble_lhs_mat(HeatProblem1D *hp);
void HeatProblem1D_assemble_rhs_mat(HeatProblem1D *hp);
void HeatProblem1D_assemble_rhs_vec(HeatProblem1D *hp);

void HeatProblem1D_apply_zero_dirichlet_bc_lhs_mat(HeatProblem1D *hp);
void HeatProblem1D_apply_zero_dirichlet_bc_rhs_vec(HeatProblem1D *hp);

void HeatProblem1D_factorize_lhs_mat(HeatProblem1D *hp);
void HeatProblem1D_solve_linsys(HeatProblem1D *hp);

void HeatProblem1D_step(HeatProblem1D *hp);

void HeatProblem1D_reset_lhs_mat(HeatProblem1D *hp);
void HeatProblem1D_reset_rhs_mat(HeatProblem1D *hp);
void HeatProblem1D_reset_rhs(HeatProblem1D *hp);
void HeatProblem1D_reset_klu(HeatProblem1D *hp);
void HeatProblem1D_reset_solution(HeatProblem1D *hp);
void HeatProblem1D_reset_system(HeatProblem1D *hp);

#endif  // __HEAT_PROBLEM_H__
