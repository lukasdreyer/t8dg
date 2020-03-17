/*
 * t8dg.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_H_
#define SRC_T8DG_H_

#define MAX_DIM 3
#define MAX_FACES (2*MAX_DIM)
#if 0
#define MAX_SUBFACES (1<<(MAX_DIM - 1))
#endif

typedef double      (*t8dg_scalar_function_MAX_DIMd_fn) (const double x[MAX_DIM]);

typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;


typedef void (*t8dg_matrix_application) (sc_array_t *dest, const sc_array_t *src, const void *application_data);
typedef void (*t8dg_time_matrix_application) (sc_array_t *dest, const sc_array_t *src, double t, const void *application_data);


typedef struct t8dg_LGL_vertex_set t8dg_LGL_vertex_set_t
typedef struct t8dg_LGL_quadrature t8dg_LGL_quadrature_t;
typedef struct t8dg_LGL_functionbasis t8dg_LGL_functionbasis_t;
typedef double (*t8dg_numerical_flux)(const double u_minus,const double u_plus, const void *application_data);

typedef struct t8dg_1D_advect_problem t8dg_1D_advect_problem_t;
typedef struct t8dg_1D_advect_element_precomputed_values t8dg_1D_advect_element_precomputed_values_t;
typedef struct t8dg_1D_advect_advance_element_data t8dg_1D_advect_advance_element_data_t;


typedef double      (*t8dg_scalar_function_MAX_DIMd_fn) (const double x[MAX_DIM]);
typedef double t8dg_square_MAX_DIM_matrix[MAX_DIM][MAX_DIM];
typedef t8dg_square_MAX_DIM_matrix t8dg_jacobian_matrix_t;

typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;
/*TODO: change sc_dmatrix_t* jacobian to double jacobian[3][3], maybe also change 3 to MAX_DIM to optimize when only 1 or 2 dims are needed. */
typedef void (*t8dg_jacobian_fn)(t8dg_jacobian_matrix_t jacobian,const double vertex[MAX_DIM], void *geometry_data);
typedef void (*t8dg_geometry_fn)(double image_vertex[MAX_DIM], const double vertex[MAX_DIM], void *geometry_data);



#endif /* SRC_T8DG_H_ */
