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
#define MAX_DIM 1

typedef struct jacobian_matrix{
  double matrix[MAX_DIM][MAX_DIM];
  int dim;
}jacobian_matrix_t;

typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;
/*TODO: change sc_dmatrix_t* jacobian to double jacobian[3][3], maybe also change 3 to MAX_DIM to optimize when only 1 or 2 dims are needed. */
typedef void (*jacobian_fn)(jacobian_matrix_t jacobian,const double vertex[3], void *geometry_data);
typedef void (*geometry_fn)(double image_vertex[3], const double vertex[3], void *geometry_data);

struct t8dg_coarse_geometry
{
  jacobian_fn		jacobian;
  geometry_fn		geometry;
};

t8dg_coarse_geometry_t *t8dg_1D_linear_geometry();

void t8dg_coarse_geometry_destroy(t8dg_coarse_geometry_t **pgeometry);

void refined_to_coarse_geometry(double coarse_element_vertex[MAX_DIM], const double refined_element_vertex[MAX_DIM],
				double scaling_factor,int idx_rotation_reflection, const double translation_vertex[MAX_DIM]);


#endif /* SRC_T8DG_GEOMETRY_H_ */
