#include <t8_cmesh.h>
#include <sc_mpi.h>

#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

double              t8dg_scalar1d_hat_function (const double x[3], const double t);

double              t8dg_scalar2d_hat_function (const double x[3], const double t);

double              t8dg_scalar3d_norm_function (const double x[3], const double t);

double              t8dg_scalar2d_step_function (const double x[3], const double t);

double              t8dg_scalar2d_triangle_step_function (const double x[3], const double t);

t8_cmesh_t          t8dg_cmesh_new_periodic_diagonal_line_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_half_moebius_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_more_trees_different_size (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_moebius (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_tilted (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm, int num_trees);

T8DG_EXTERN_C_END ();
