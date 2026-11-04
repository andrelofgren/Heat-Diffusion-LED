#include <stddef.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <suitesparse/cs.h>

#include "../include/heat_problem_1d.h"

#define NUM_CASES 10
#define CASES {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560}

double heat_source(const double x, const double t)
{
    const double A0 = 50.0;
    const double sigma = 0.05;
    const double mu = 0.4;
    const double X = sqrt(0.5)*((x-mu)/sigma);
    const double tau = 2.0;
    return A0*exp(-t/tau)*exp(-X*X)/(sigma*sqrt(2*M_PI));
}

void solve_heatproblem(HeatProblem1D *hp, const size_t nt)
{
    HeatProblem1D_setup(hp);
    for (int k = 0; k < nt; k++) {
        HeatProblem1D_step(hp);
        HeatProblem1D_reset_rhs(hp);
        hp->t += hp->dt;
    }
}

int main(int argc, char *argv[])
{
    const double T = 10.0;
    const double x0 = 0.0;
    const double x1 = 0.8;
    const size_t nx = 100;
    const size_t nt_ref = 5000;
    const size_t deg = 1;

    Mesh1D mesh;
    Mesh1D_init(&mesh, x0, x1, nx, deg);
    HeatProblem1D hp_ref;
    hp_ref.alpha = 0.0025;
    hp_ref.beta = 0.1;
    hp_ref.gamma = 0.2;
    hp_ref.source = heat_source;
    hp_ref.t = 0.0;
    hp_ref.dt = T/nt_ref;
    klu_defaults(&hp_ref.klu_common);

    HeatProblem1D_init(&hp_ref, &mesh);
    solve_heatproblem(&hp_ref, nt_ref);

    const size_t num_dofs = Mesh1D_get_num_vertices(&mesh);
    gsl_vector *u_vec_ref = gsl_vector_alloc(num_dofs);
    gsl_vector_memcpy(u_vec_ref, hp_ref.sol_vec);

    const size_t nt_vec_data[NUM_CASES] = CASES;
    gsl_vector *u_vec = gsl_vector_alloc(num_dofs);
    double error_vec[NUM_CASES+1];
    error_vec[0] = 1.0;
    double error_norm;
    HeatProblem1D hp;
    hp.alpha = 0.0025;
    hp.beta = 0.1;
    hp.gamma = 0.2;
    hp.source = heat_source;
    HeatProblem1D_init(&hp, &mesh);
    for (size_t k = 0; k < NUM_CASES; k++){
        hp.t = 0.0;
        hp.dt = T/nt_vec_data[k];
        klu_defaults(&hp.klu_common);
        solve_heatproblem(&hp, nt_vec_data[k]);
        gsl_vector_memcpy(u_vec, hp.sol_vec);
        gsl_vector_sub(u_vec, u_vec_ref);
        gsl_blas_ddot(u_vec, u_vec, &error_norm);
        error_vec[k+1] = 1.0/num_dofs * sqrt(error_norm);
        printf(
            "case %d: \t dt = %g, \t e_k = %g, \t e_{k-1}/e_k = %g\n",
            (int) k+1, hp.dt, error_vec[k+1], error_vec[k]/error_vec[k+1]
        );
        HeatProblem1D_reset_system(&hp);
    }

    // Finalize
    HeatProblem1D_finalize(&hp);
    HeatProblem1D_finalize(&hp_ref);
    gsl_vector_free(u_vec_ref);
    gsl_vector_free(u_vec);
    Mesh1D_finalize(&mesh);

    return 0;
}
