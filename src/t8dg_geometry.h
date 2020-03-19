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
#include "t8.h"
#include "t8_element.h"


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

void
t8dg_vec_print (const double vec_x[3]);

void t8dg_element_set_geometry_data(t8dg_element_fine_to_coarse_geometry_data_t *geometry_data ,t8_element_t *element,t8dg_locidx_t idata,t8_eclass_scheme_c *scheme);


t8dg_coarse_geometry_3D_t *t8dg_coarse_geometry_new_1D_linear();

void t8dg_coarse_geometry_destroy(t8dg_coarse_geometry_3D_t **pgeometry);

void t8dg_fine_to_coarse_geometry(double coarse_element_vertex[DIM3], double refined_element_vertex[DIM3],
				t8dg_element_fine_to_coarse_geometry_data_t *element_data);

void t8dg_square3D_matrix_invert_sub_matrix(t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim);

void t8dg_square3D_matrix_determinant_sub_matrix(double *det, t8dg_square_3D_matrix_t matrix, int dim);

void t8dg_square3D_matrix_scale(t8dg_square_3D_matrix_t matrix, double scaling_factor, int dim);

void t8dg_square3D_matrix_apply_rotation_reflection_matrix_to_matrix(t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int idx_rotation_reflection, int dim);

#endif /* SRC_T8DG_GEOMETRY_H_ */
