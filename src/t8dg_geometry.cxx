/*
 * t8dg_geometry.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include <t8.h>
#include <t8_vec.h>
#include <t8_element_cxx.hxx>
#include <t8_element.h>
#include <t8_forest.h>

#include "t8dg.h"
#include "t8dg_geometry.h"

/*geometry_data are tree_vertices*/
static void
t8dg_constant_1D_jacobian_fn (const double vertex[DIM3], t8dg_square_3D_matrix_t jacobian, t8_forest_t forest, t8_locidx_t itree)
{
  T8DG_ASSERT (vertex != NULL);

  double             *tree_vertices = t8_forest_get_tree_vertices (forest, itree);

  /*tree-vertices is 3*number_of_vertices */
  double             *x_0 = tree_vertices;
  double             *x_1 = tree_vertices + 3;

  double              h = t8_vec_dist (x_0, x_1);

  jacobian[0][0] = h;           /*this is the length of the coarse line */
}

static void
t8dg_linear_1D_geometry_fn (const double vertex[DIM3], double image_vertex[DIM3], t8_forest_t forest, t8_locidx_t itree)
{
  T8DG_ASSERT (vertex != NULL && image_vertex != NULL);

  double             *tree_vertices = t8_forest_get_tree_vertices (forest, itree);

  double             *x_0 = tree_vertices;
  double             *x_1 = tree_vertices + DIM3;

  t8_debugf ("x_0\n");
  t8dg_vec_print (x_0);
  t8_debugf ("x_1\n");
  t8dg_vec_print (x_1);

  double              h = t8_vec_dist (x_0, x_1);
  t8_vec_axpyz (vertex, x_0, image_vertex, h);
}

t8dg_coarse_geometry_3D_t *
t8dg_coarse_geometry_new_1D_linear ()
{
  t8dg_coarse_geometry_3D_t *geometry = T8_ALLOC (t8dg_coarse_geometry_3D_t, 1);
  geometry->geometry = t8dg_linear_1D_geometry_fn;
  geometry->jacobian = t8dg_constant_1D_jacobian_fn;
  geometry->data_type = T8DG_TREE_VERTICES;
  return geometry;
}

void
t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_3D_t ** pgeometry)
{
  T8DG_FREE (*pgeometry);
  *pgeometry = NULL;
}

void
t8dg_fine_to_coarse_geometry (const double refined_element_vertex[DIM3],
                              double coarse_element_vertex[DIM3], t8_eclass_scheme_c * scheme, const t8_element_t * element)
{
  T8DG_ASSERT (coarse_element_vertex != NULL && refined_element_vertex != NULL);
  int                 idim;
  int                 vertex_int_coords[3];
  int                 level = scheme->t8_element_level (element);
  double              translation_vector[3] = { 0, 0, 0 };
  double              scaling_factor = pow (2, -level);
  double              length_inv = 1. / scheme->t8_element_root_len (element);
  scheme->t8_element_vertex_coords (element, 0, vertex_int_coords);
  for (idim = 0; idim < DIM3; idim++) {
    translation_vector[idim] = length_inv * vertex_int_coords[idim];
  }

  t8_vec_axpyz (refined_element_vertex, translation_vector, coarse_element_vertex, scaling_factor);
  /*hx+x_0 */

}
