#include <t8_cmesh.h>

double              t8dg_scalar1d_hat_function (const double x[3], const double t);

double              t8dg_scalar2d_hat_function (const double x[3], const double t);

double              t8dg_scalar3d_norm_function (const double x[3], const double t);

double              t8dg_scalar2d_step_function (const double x[3], const double t);

double              t8dg_scalar2d_triangle_step_function (const double x[3], const double t);

t8_cmesh_t          t8dg_cmesh_new_periodic_diagonal_line_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_half_moebius_more_trees (sc_MPI_Comm comm);
