/*
 * t8dg_geometry.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_GEOMETRY_H_
#define SRC_T8DG_GEOMETRY_H_

#include <sc_containers.h>
#include "t8dg.h"


struct t8dg_coarse_geometry_3D
{
  t8dg_jacobian_fn_3D		jacobian;
  t8dg_geometry_fn_3D		geometry;
};


struct t8dg_element_fine_to_coarse_geometry_data
{
/* those are needed for the transformation from reference element to the location in the coarse tree*/
  double scaling_factor;
  double translation_vector[DIM3];
  int idx_rotation_reflection;
};




t8dg_coarse_geometry_3D_t *t8dg_1D_linear_geometry();

void t8dg_coarse_geometry_destroy(t8dg_coarse_geometry_3D_t **pgeometry);

void t8dg_fine_to_coarse_geometry(double coarse_element_vertex[DIM3], double refined_element_vertex[DIM3],
				t8dg_element_fine_to_coarse_geometry_data_t *element_data);

void t8dg_invert_sub_square_matrix(t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim);

void t8dg_determinant_sub_square_matrix(double *det, t8dg_square_3D_matrix_t matrix, int dim);

#endif /* SRC_T8DG_GEOMETRY_H_ */
