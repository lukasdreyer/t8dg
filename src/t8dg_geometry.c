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
  T8_ASSERT(geometry_data!=NULL);
  T8_ASSERT(vertex!=NULL);
  double *tree_vertices = (double*) geometry_data;

  /*tree-vertices is 3*number_of_vertices*/
  double x_0 = tree_vertices[0];
  double x_1 = tree_vertices[3];

  jacobian[0][0] = x_1 - x_0; /*this is the lenght of the coarse line*/
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
void refined_to_coarse_geometry(double coarse_element_vertex[MAX_DIM], double refined_element_vertex[MAX_DIM],
				t8dg_1D_advect_element_precomputed_values_t *element_values){
  T8_ASSERT(element_values->idx_rotation_reflection >=0 && element_values->idx_rotation_reflection <=1);
  T8_ASSERT(coarse_element_vertex!=NULL && refined_element_vertex != NULL);
  int idim;

  /*x=Rx*/
  // apply_rotation_reflection(coarse_element_vertex, refined_element_vertex, element_values->idx_rotation_reflection);

  /*hx+x_0*/
  /*TODO: implement vector functions!*/
  for(idim = 0; idim < element_values->dim ; idim++){
      coarse_element_vertex[idim] = element_values->scaling_factor * coarse_element_vertex[idim]
		+ element_values->translation_vector[idim];
  }
  /*not needed if rot/reflection already sets unused values to 0 */
  for(idim = element_values->dim ; idim < MAX_DIM; idim++){
      coarse_element_vertex[idim] = 0;
  }
}

void invert_jacobian_matrix(jacobian_matrix_t jacobian_invers, jacobian_matrix_t jacobian_matrix, int dim){
  T8_ASSERT(dim > 0 && dim <= MAX_DIM);
  double det;
  determinant_jacobian_matrix(&det,jacobian_matrix,dim);
  if(dim ==1){
      jacobian_invers[0][0] = 1. / det;
  }
  else if(dim ==2){
      jacobian_invers [0][0] = 1./ det * jacobian_matrix[1][1];
      jacobian_invers [1][1] = 1./ det * jacobian_matrix[0][0];
      jacobian_invers [0][1] = - 1./ det * jacobian_matrix[0][1];
      jacobian_invers [1][0] = - 1./ det * jacobian_matrix[1][0];
  }
  else if(dim ==3){
      SC_ABORT("not yet implemented");
  }
}

void determinant_jacobian_matrix(double *det, jacobian_matrix_t jacobian_matrix, int dim){
  T8_ASSERT(dim > 0 && dim <= MAX_DIM);
  if(dim ==1){
      *det = jacobian_matrix[0][0];
  }
  else if(dim ==2){
      *det = jacobian_matrix[0][0]*jacobian_matrix[1][1] - jacobian_matrix[0][1] * jacobian_matrix [1][0];
  }
  else if(dim ==3){
      SC_ABORT("not yet implemented");
  }
}
