/*
 * t8dg_geometry.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_GEOMETRY_H_
#define SRC_T8DG_GEOMETRY_H_

#include <sc_containers.h>
#include <sc_dmatrix.h>
#include "global.h"


typedef double      (*t8dg_scalar_function_MAX_DIMd_fn) (const double x[MAX_DIM]);
typedef double square_MAX_DIM_matrix[MAX_DIM][MAX_DIM];
typedef square_MAX_DIM_matrix jacobian_matrix_t;

typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;
/*TODO: change sc_dmatrix_t* jacobian to double jacobian[3][3], maybe also change 3 to MAX_DIM to optimize when only 1 or 2 dims are needed. */
typedef void (*jacobian_fn)(jacobian_matrix_t jacobian,const double vertex[MAX_DIM], void *geometry_data);
typedef void (*geometry_fn)(double image_vertex[MAX_DIM], const double vertex[MAX_DIM], void *geometry_data);

struct t8dg_coarse_geometry
{
  jacobian_fn		jacobian;
  geometry_fn		geometry;
};

t8dg_coarse_geometry_t *t8dg_1D_linear_geometry();

void t8dg_coarse_geometry_destroy(t8dg_coarse_geometry_t **pgeometry);

void refined_to_coarse_geometry(double coarse_element_vertex[MAX_DIM], double refined_element_vertex[MAX_DIM],
				t8dg_1D_advect_element_precomputed_values_t *element_values);

void invert_jacobian_matrix(jacobian_matrix_t jacobian_invers, jacobian_matrix_t jacobian_matrix, int dim);

void determinant_jacobian_matrix(double *det, jacobian_matrix_t jacobian_matrix, int dim);

#endif /* SRC_T8DG_GEOMETRY_H_ */
