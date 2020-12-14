#ifndef SRC_T8DG_COMMON_H_
#define SRC_T8DG_COMMON_H_

#include <t8_cmesh.h>
#include <sc_mpi.h>

#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

t8dg_scalar_function_3d_time_fn t8dg_common_initial_cond_fn (int initial_cond_arg);

double              t8dg_scalar1d_hat_function (const double x[3], const double t);

double              t8dg_scalar2d_hat_function (const double x[3], const double t);

double              t8dg_scalar3d_norm_function (const double x[3], const double t);

double              t8dg_scalar2d_step_function (const double x[3], const double t);

double              t8dg_scalar2d_triangle_step_function (const double x[3], const double t);

double              t8dg_scalar3d_step_function (const double x[3], const double t);

T8DG_EXTERN_C_END ();

#endif
