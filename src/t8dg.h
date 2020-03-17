/*
 * t8dg.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_H_
#define SRC_T8DG_H_

#include <t8.h>
#define T8DG_ASSERT T8_ASSERT
#define T8DG_ALLOC T8_ALLOC
#define T8DG_FREE T8_FREE

#define DIM3 3
#define MAX_FACES 6
#if 0
#define MAX_SUBFACES (1<<(MAX_DIM - 1))
#endif

typedef t8_locidx_t t8dg_locidx_t;

typedef double      (*t8dg_scalar_function_3d_fn) (const double x[DIM3]);

typedef void (*t8dg_matrix_application) (sc_array_t *dest, const sc_array_t *src, const void *application_data);
typedef void (*t8dg_time_matrix_application) (sc_array_t *dest, const sc_array_t *src, double t, const void *application_data);

typedef double t8dg_square_3D_matrix_t[DIM3][DIM3];


typedef struct t8dg_LGL_vertex_set t8dg_LGL_vertex_set_t;
typedef struct t8dg_LGL_quadrature t8dg_LGL_quadrature_t;
typedef struct t8dg_LGL_functionbasis t8dg_LGL_functionbasis_t;


typedef struct t8dg_coarse_geometry_3D t8dg_coarse_geometry_3D_t;
typedef struct t8dg_element_fine_to_coarse_geometry_data t8dg_element_fine_to_coarse_geometry_data_t;
typedef void (*t8dg_jacobian_fn_3D)(t8dg_square_3D_matrix_t jacobian,const double vertex[DIM3], void *geometry_data);
typedef void (*t8dg_geometry_fn_3D)(double image_vertex[DIM3], const double vertex[DIM3], void *geometry_data);


typedef struct t8dg_advect_flux_mortar t8dg_advect_flux_mortar_t;
typedef double (*t8dg_numerical_flux_1D_fn)(const double u_minus,const double u_plus, const void *application_data);



//typedef struct t8dg_advect_problem_linear_1D t8dg_advect_problem_linear_1D_t;




#endif /* SRC_T8DG_H_ */
