/*
 * t8dg_geometry.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include <t8.h>
#include <t8_vec.h>
#include "t8dg.h"
#include "t8dg_geometry.h"



/*geometry_data are tree_vertices*/
static void t8dg_constant_1D_jacobian_fn(t8dg_square_3D_matrix_t jacobian, const double vertex[DIM3], void *geometry_data){
  T8_ASSERT(geometry_data!=NULL);
  T8_ASSERT(vertex!=NULL);
  double *tree_vertices = (double*) geometry_data;

  /*tree-vertices is 3*number_of_vertices*/
  double *x_0 = tree_vertices;
  double *x_1 = tree_vertices + 3;

  double h = t8_vec_dist(x_0,x_1);

  jacobian[0][0] = h; /*this is the lenght of the coarse line*/
}

static void t8dg_linear_1D_geometry_fn(double image_vertex[DIM3], const double vertex[DIM3], void *geometry_data){
  T8_ASSERT(geometry_data!=NULL);
  T8_ASSERT(vertex!=NULL&&image_vertex!=NULL);
  double *tree_vertices = (double*) geometry_data;

  /*TODO: use index function for tree vertices*/
  double *x_0 = tree_vertices;
  double *x_1 = tree_vertices + 3;

  double h = t8_vec_dist(x_0,x_1);
  t8_vec_axpyz(vertex, x_0, image_vertex, h);
}


t8dg_coarse_geometry_3D_t *t8dg_coarse_geometry_new_1D_linear(){
  t8dg_coarse_geometry_3D_t *geometry = T8_ALLOC(t8dg_coarse_geometry_3D_t,1);
  geometry->geometry = t8dg_linear_1D_geometry_fn;
  geometry->jacobian = t8dg_constant_1D_jacobian_fn;
  return geometry;
}

void t8dg_coarse_geometry_destroy(t8dg_coarse_geometry_3D_t **pgeometry){
  T8_FREE(*pgeometry);
  *pgeometry = NULL;
}

/*TODO: implement*/
void t8dg_fine_to_coarse_geometry(double coarse_element_vertex[DIM3], double fine_element_vertex[DIM3],
				t8dg_element_fine_to_coarse_geometry_data_t *element_data){
  T8_ASSERT(element_data->idx_rotation_reflection ==0);
  T8_ASSERT(coarse_element_vertex!=NULL && fine_element_vertex != NULL);
  int idim;

  /*x=Rx*/
  // apply_rotation_reflection(coarse_element_vertex, refined_element_vertex, element_values->idx_rotation_reflection);

  /*hx+x_0*/
  /*TODO: implement vector functions! -> can now be done by t8_vec*/

}

void t8dg_invert_sub_square_matrix(t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim){
  T8_ASSERT(dim > 0 && dim <= DIM3);
  double det;
  t8dg_determinant_sub_square_matrix(&det,matrix,dim);
  if(dim ==1){
      matrix_invers[0][0] = 1. / det;
  }
  else if(dim ==2){
      matrix_invers [0][0] = 1./ det * matrix[1][1];
      matrix_invers [1][1] = 1./ det * matrix[0][0];
      matrix_invers [0][1] = - 1./ det * matrix[0][1];
      matrix_invers [1][0] = - 1./ det * matrix[1][0];
  }
  else if(dim ==3){
      SC_ABORT("not yet implemented");
  }
}

void t8dg_determinant_sub_square_matrix(double *det, t8dg_square_3D_matrix_t matrix, int dim){
  T8_ASSERT(dim > 0 && dim <= DIM3);
  if(dim ==1){
      *det = matrix[0][0];
  }
  else if(dim ==2){
      *det = matrix[0][0]*matrix[1][1] - matrix[0][1] * matrix [1][0];
  }
  else if(dim ==3){
      SC_ABORT("not yet implemented");
  }
}
