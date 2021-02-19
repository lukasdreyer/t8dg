#ifndef SRC_T8DG_COMMON_H_
#define SRC_T8DG_COMMON_H_

#include <t8_cmesh.h>
#include <sc_mpi.h>

#include "t8dg.h"

/**A timedependent scalar function f:R^3 x R^+ -> R*/
typedef double      (*t8dg_scalar_function_3d_time_fn) (const double x[DIM3], const double t, void *fn_data);
typedef double      (*t8dg_scalar_function_3d_fn) (const double x[DIM3], void *scalar_fn_data);

typedef struct t8dg_scalar3d_cos_product_data
{
  int                 dim;
  double              diffusion_coefficient;
} t8dg_scalar3d_cos_product_data_t;

T8DG_EXTERN_C_BEGIN ();

t8dg_scalar_function_3d_time_fn t8dg_common_initial_cond_fn (int initial_cond_arg);

t8dg_scalar_function_3d_time_fn t8dg_common_analytic_solution_fn (int initial_cond_arg, double diffusion_coefficient);

double              t8dg_scalar1d_hat_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar1d_step_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar2d_hat_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar3d_norm_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar2d_step_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar2d_triangle_step_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar3d_step_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar3d_constant_one (const double x[3], const double t, void *fn_data);

double              t8dg_scalar3d_cos_product (const double x[3], const double t, void *fn_data);

double              t8dg_scalar3d_constant_zero (const double x[3], const double t, void *fn_data);

double              t8dg_circle_ring_step_function (const double x[3], const double t, void *fn_data);

double              t8dg_scalar2d_angle (const double x[3], const double t, void *fn_data);

double              t8dg_cylinder_ring_sin_product_fn (const double x[3], const double t, void *fn_data);

double              t8dg_cylinder_ring_step_function (const double x[3], const double t, void *fn_data);

double              t8dg_cylinder_ring_source_fn (const double x[3], const double t, void *fn_data);

double              t8dg_smooth_indicator1Dfn (const double x[3], const double t, void *fn_data);

double              t8dg_smooth_indicator2Dfn (const double x[3], const double t, void *fn_data);

double              t8dg_smooth_indicator3Dfn (const double x[3], const double t, void *fn_data);

double              t8dg_circle_ring_sin_product_fn (const double x[3], const double t, void *fn_data);

T8DG_EXTERN_C_END ();

#endif
