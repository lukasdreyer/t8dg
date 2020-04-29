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

typedef double      (*t8dg_coarse_geometry_gram_determinant_fn) (const t8_forest_t forest, const t8_locidx_t itree, void *data,
                                                                 const double coarse_vertex[3]);

typedef void        (*t8dg_coarse_geometry_differential_fn) (const t8_forest_t forest, const t8_locidx_t itree, void *data,
                                                             const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                             double transformed_gradient_tangential_vector[3]);

/** returns the image vertex of the geometry function*/
typedef void        (*t8dg_coarse_geometry_fn) (const t8_forest_t forest, const t8_locidx_t itree, void *data,
                                                const double vertex[DIM3], double image_vertex[DIM3]);

typedef enum t8dg_coarse_geometry_data
{
  T8DG_TREE_VERTICES
} t8dg_coarse_geometry_data_t;

/** F_CE and DF_CE for the coarse geometry
 * TODO: How to chose unused dimensions so that DF_CE can be used to calculate the differential on the submanifold*/
struct t8dg_coarse_geometry
{
  t8dg_coarse_geometry_gram_determinant_fn gram_det; /**< Jacobian function*/
  t8dg_coarse_geometry_differential_fn differential_invers_transpose;
  t8dg_coarse_geometry_fn geometry; /**< Geometry function*/
  t8dg_coarse_geometry_data_t attribute_data_type;
  void               *data;
};

/*TODO: Document In and Output*/

double
t8dg_coarse_geometry_calculate_gram_determinant (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                                                 const t8_locidx_t itree, const double coarse_vertex[3])
{
  return coarse_geometry->gram_det (forest, itree, coarse_geometry->data, coarse_vertex);
}

void
t8dg_coarse_geometry_apply (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                            const t8_locidx_t itree, const double coarse_vertex[3], double image_vertex[3])
{
  coarse_geometry->geometry (forest, itree, coarse_geometry->data, coarse_vertex, image_vertex);
}

void
t8dg_coarse_geometry_calculate_gradient_tangential_vector (const t8dg_coarse_geometry_t * coarse_geometry,
                                                           const t8_forest_t forest, const t8_locidx_t itree,
                                                           const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                           double transformed_gradient_tangential_vector[3])
{
  coarse_geometry->differential_invers_transpose (forest, itree, coarse_geometry->data, coarse_vertex, coarse_tangential_vector,
                                                  transformed_gradient_tangential_vector);
}

void
t8dg_coarse_geometry_transform_normal_vector (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                                              const t8_locidx_t itree, const double coarse_vertex[3], const double coarse_normal_vector[3],
                                              double image_normal_vector[3])
{
  /*TODO: Check for nonlinear geometry! */
  coarse_geometry->differential_invers_transpose (forest, itree, coarse_geometry->data, coarse_vertex, coarse_normal_vector,
                                                  image_normal_vector);
}

static void
t8dg_linear_1D_gram_differential_invers_transpose_fn (const t8_forest_t forest, const t8_locidx_t itree, void *data,
                                                      const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                      double transformed_gradient_tangential_vector[3])
{
  double             *tree_vertices;
  double              image_coarse_element_length_vector[3];
  double              image_coarse_element_length_squared;

  tree_vertices = t8_forest_get_tree_vertices (forest, itree);

  /*tangential vector of the coarse element */
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector, -1);

  image_coarse_element_length_squared = t8_vec_norm (image_coarse_element_length_vector);

  t8_vec_axb (coarse_tangential_vector, transformed_gradient_tangential_vector, 1. / image_coarse_element_length_squared, 0);

}

static double
t8dg_linear_1D_gram_determinant_fn (const t8_forest_t forest, const t8_locidx_t itree, void *data, const double coarse_vertex[3])
{
  double             *tree_vertices;
  double              image_coarse_element_length_vector[3];
  double              image_coarse_element_length;

  tree_vertices = t8_forest_get_tree_vertices (forest, itree);

  /*tangential vector of the coarse element */
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector, -1);

  image_coarse_element_length = t8_vec_norm (image_coarse_element_length_vector);
  return image_coarse_element_length;
}

static void
t8dg_linear_1D_geometry_fn (const t8_forest_t forest, const t8_locidx_t itree, void *data, const double vertex[DIM3],
                            double image_vertex[DIM3])
{
  T8DG_ASSERT (vertex != NULL && image_vertex != NULL);

  double             *tree_vertices = t8_forest_get_tree_vertices (forest, itree);

  double              image_coarse_element_length_vector[3];

  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector, -1);
  t8_vec_axpyz (image_coarse_element_length_vector, tree_vertices, image_vertex, vertex[0]);
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_1D_linear ()
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_linear_1D_geometry_fn;
  geometry->gram_det = t8dg_linear_1D_gram_determinant_fn;
  geometry->differential_invers_transpose = t8dg_linear_1D_gram_differential_invers_transpose_fn;
  geometry->attribute_data_type = T8DG_TREE_VERTICES;
  geometry->data = NULL;
  return geometry;
}

void
t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_t ** pgeometry)
{
  T8DG_FREE (*pgeometry);
  *pgeometry = NULL;
}
