/*
 * t8dg_geometry.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include "t8dg_coarse_geometry.h"

#include <t8.h>
#include <t8_vec.h>
#include <t8_forest.h>

#include "t8dg.h"

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

  double              h = t8_vec_dist (x_0, x_1);
  t8_vec_axpyz (vertex, x_0, image_vertex, h);
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_1D_linear ()
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_linear_1D_geometry_fn;
  geometry->jacobian = t8dg_constant_1D_jacobian_fn;
  geometry->data_type = T8DG_TREE_VERTICES;
  return geometry;
}

void
t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_t ** pgeometry)
{
  T8DG_FREE (*pgeometry);
  *pgeometry = NULL;
}
