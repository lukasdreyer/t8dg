/*
 * t8dg_geometry.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include <t8.h>
#include <t8_vec.h>
#include <t8_element_cxx.hxx>

#include "t8dg.h"
#include "t8dg_geometry.h"

void
t8dg_vec_print (const double vec_x[3])
{
  int                 i;

  for (i = 0; i < 3; i++) {
    printf ("%f ", vec_x[i]);
  }
  printf ("\n");
}

void
t8dg_element_set_geometry_data (t8dg_element_fine_to_coarse_geometry_data_t * geometry_data, t8_element_t * element,
                                t8_locidx_t idata, t8_eclass_scheme_c * scheme)
{
  int                 level;
  double              length;
  int                 tree_int_coords[3];
  int                 idim;
  level = scheme->t8_element_level (element);
  geometry_data->idx_rotation_reflection = 0;
  geometry_data->scaling_factor = pow (2, -level);

  length = 1. / scheme->t8_element_root_len (element);
  scheme->t8_element_vertex_coords (element, 0, tree_int_coords);
  tree_int_coords[1] = 0;
  tree_int_coords[2] = 0;

  for (idim = 0; idim < DIM3; idim++) {
    geometry_data->translation_vector[idim] = length * tree_int_coords[idim];
  }
}

/*geometry_data are tree_vertices*/
static void
t8dg_constant_1D_jacobian_fn (t8dg_square_3D_matrix_t jacobian, const double vertex[DIM3], void *geometry_data)
{
  T8_ASSERT (geometry_data != NULL);
  T8_ASSERT (vertex != NULL);
  double             *tree_vertices = (double *) geometry_data;

  /*tree-vertices is 3*number_of_vertices */
  double             *x_0 = tree_vertices;
  double             *x_1 = tree_vertices + 3;

  double              h = t8_vec_dist (x_0, x_1);

  jacobian[0][0] = h;           /*this is the lenght of the coarse line */
}

static void
t8dg_linear_1D_geometry_fn (double image_vertex[DIM3], const double vertex[DIM3], void *geometry_data)
{
  T8_ASSERT (geometry_data != NULL);
  T8_ASSERT (vertex != NULL && image_vertex != NULL);
  double             *tree_vertices = (double *) geometry_data;

  double             *x_0 = tree_vertices;
  double             *x_1 = tree_vertices + DIM3;

  double              h = t8_vec_dist (x_0, x_1);
  t8_vec_axpyz (vertex, x_0, image_vertex, h);
}

t8dg_coarse_geometry_3D_t *
t8dg_coarse_geometry_new_1D_linear ()
{
  t8dg_coarse_geometry_3D_t *geometry = T8_ALLOC (t8dg_coarse_geometry_3D_t, 1);
  geometry->geometry = t8dg_linear_1D_geometry_fn;
  geometry->jacobian = t8dg_constant_1D_jacobian_fn;
  return geometry;
}

void
t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_3D_t ** pgeometry)
{
  T8_FREE (*pgeometry);
  *pgeometry = NULL;
}

void
t8dg_fine_to_coarse_geometry (double coarse_element_vertex[DIM3], double fine_element_vertex[DIM3],
                              t8dg_element_fine_to_coarse_geometry_data_t * element_data)
{
  T8_ASSERT (element_data->idx_rotation_reflection == 0);
  T8_ASSERT (coarse_element_vertex != NULL && fine_element_vertex != NULL);
  T8_ASSERT (element_data->scaling_factor != 0);

  /*TODO: implement */
  /*x=Rx */

  t8_vec_axpyz (fine_element_vertex, element_data->translation_vector, coarse_element_vertex, element_data->scaling_factor);
  /*hx+x_0 */

}
