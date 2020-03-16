/*
 * t8dg_geometry.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include <t8.h>
#include "t8dg_geometry.h"

//#include <t8_vec.h>


/*geometry_data are tree_vertices*/
static void constant_1D_jacobian_fn(jacobian_matrix_t jacobian, const double vertex[MAX_DIM], void *geometry_data){
  T8_ASSERT(jacobian.dim == 1);
  T8_ASSERT(geometry_data!=NULL);
  T8_ASSERT(vertex!=NULL);
  double *tree_vertices = (double*) geometry_data;

  /*tree-vertices is 3*number_of_vertices*/
  double x_0 = tree_vertices[0];
  double x_1 = tree_vertices[3];

  jacobian.matrix[0][0] = x_1 - x_0; /*this is the lenght of the coarse line*/
}

static void linear_1D_geometry_fn(double image_vertex[MAX_DIM], const double vertex[MAX_DIM], void *geometry_data){
  T8_ASSERT(geometry_data!=NULL);
  T8_ASSERT(vertex!=NULL&&image_vertex!=NULL);
  double *tree_vertices = (double*) geometry_data;
  double x_0 = tree_vertices[0];
  double x_1 = tree_vertices[3];
  double h = x_1 - x_0;

  image_vertex[0] = x_0 + vertex[0] * h;
}


t8dg_coarse_geometry_t *t8dg_1D_linear_geometry(){
  t8dg_coarse_geometry_t *geometry = T8_ALLOC(t8dg_coarse_geometry_t,1);
  geometry->geometry = linear_1D_geometry_fn;
  geometry->jacobian = constant_1D_jacobian_fn;
  return geometry;
}

void t8dg_coarse_geometry_destroy(t8dg_coarse_geometry_t **pgeometry){
  T8_FREE(*pgeometry);
  *pgeometry = NULL;
}

/*TODO: implement*/
void refined_to_coarse_geometry(double coarse_element_vertex[MAX_DIM], const double refined_element_vertex[MAX_DIM],
				double scaling_factor,int idx_rotation_reflection, const double translation_vertex[MAX_DIM]){



}

