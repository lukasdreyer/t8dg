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

typedef double      (*t8dg_coarse_geometry_sqrt_gram_determinant_fn) (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                                      const double coarse_vertex[3]);

typedef double      (*t8dg_coarse_geometry_sqrt_face_gram_determinant_fn) (const t8_forest_t forest, const t8_gloidx_t iglobaltree,
                                                                           void *data, const int iface, const double coarse_vertex[3]);

typedef void        (*t8dg_coarse_geometry_differential_fn) (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                             const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                             double transformed_gradient_tangential_vector[3]);

/** returns the image vertex of the geometry function*/
typedef void        (*t8dg_coarse_geometry_fn) (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                const double vertex[DIM3], double image_vertex[DIM3]);

typedef enum t8dg_coarse_geometry_data
{
  T8DG_TREE_VERTICES
} t8dg_coarse_geometry_data_t;

struct t8dg_coarse_geometry
{
  t8dg_coarse_geometry_sqrt_gram_determinant_fn sqrt_gram_det; /**< Weighting factor for integrals*/
  t8dg_coarse_geometry_sqrt_face_gram_determinant_fn sqrt_face_gram_det; /**< Weighting factor for integrals*/
  t8dg_coarse_geometry_differential_fn differential_invers_transpose; /**< Transformation for gradient and normal vector*/
  t8dg_coarse_geometry_fn geometry; /**< Geometry function*/
  t8dg_coarse_geometry_data_t attribute_data_type;  /**< Determines which attribute from the cmesh should be taken*/
  void               *data; /**< Additional data provided to the coarse geometry*/
};

/*TODO: Document In and Output*/
double
t8dg_coarse_geometry_calculate_sqrt_face_gram_determinant (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                                                           const t8_gloidx_t iglobaltree, const int iface, const double coarse_vertex[3])
{
  return coarse_geometry->sqrt_face_gram_det (forest, iglobaltree, coarse_geometry->data, iface, coarse_vertex);
}

double
t8dg_coarse_geometry_calculate_sqrt_gram_determinant (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                                                      const t8_gloidx_t iglobaltree, const double coarse_vertex[3])
{
  return coarse_geometry->sqrt_gram_det (forest, iglobaltree, coarse_geometry->data, coarse_vertex);
}

void
t8dg_coarse_geometry_apply (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                            const t8_gloidx_t iglobaltree, const double coarse_vertex[3], double image_vertex[3])
{
  coarse_geometry->geometry (forest, iglobaltree, coarse_geometry->data, coarse_vertex, image_vertex);
}

void
t8dg_coarse_geometry_calculate_gradient_tangential_vector (const t8dg_coarse_geometry_t * coarse_geometry,
                                                           const t8_forest_t forest, const t8_gloidx_t iglobaltree,
                                                           const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                           double transformed_gradient_tangential_vector[3])
{
  coarse_geometry->differential_invers_transpose (forest, iglobaltree, coarse_geometry->data, coarse_vertex, coarse_tangential_vector,
                                                  transformed_gradient_tangential_vector);
}

void
t8dg_coarse_geometry_transform_normal_vector (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
                                              const t8_gloidx_t iglobaltree, const double coarse_vertex[3],
                                              const double coarse_normal_vector[3], double image_normal_vector[3])
{
  /*TODO: Check for nonlinear geometry! */
  coarse_geometry->differential_invers_transpose (forest, iglobaltree, coarse_geometry->data, coarse_vertex, coarse_normal_vector,
                                                  image_normal_vector);
}

static void
t8dg_linear_1D_differential_invers_transpose_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                 const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                 double transformed_gradient_tangential_vector[3])
{
  double             *tree_vertices;
  double              image_coarse_element_length_vector[3];
  double              image_coarse_element_length_squared;

  tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);
  T8DG_ASSERT (tree_vertices != NULL);
  /*tangential vector of the coarse element */
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector, -1);

  image_coarse_element_length_squared = t8_vec_dot (image_coarse_element_length_vector, image_coarse_element_length_vector);

  t8_vec_axb (image_coarse_element_length_vector, transformed_gradient_tangential_vector,
              coarse_tangential_vector[0] / image_coarse_element_length_squared, 0);
}

static double
t8dg_linear_1D_sqrt_face_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const int iface,
                                              const double coarse_vertex[3])
{
  return 1;
}

static double
t8dg_linear_1D_sqrt_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double coarse_vertex[3])
{
  double             *tree_vertices;
  double              image_coarse_element_length_vector[3];
  double              image_coarse_element_length;

  tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  /*tangential vector of the coarse element */
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector, -1);

  image_coarse_element_length = t8_vec_norm (image_coarse_element_length_vector);
  return image_coarse_element_length;
}

static void
t8dg_linear_1D_geometry_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double vertex[DIM3],
                            double image_vertex[DIM3])
{
  T8DG_ASSERT (vertex != NULL && image_vertex != NULL);

  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector[3];

  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector, -1);
  t8_vec_axpyz (image_coarse_element_length_vector, tree_vertices, image_vertex, vertex[0]);
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_1D_linear ()
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_linear_1D_geometry_fn;
  geometry->sqrt_gram_det = t8dg_linear_1D_sqrt_gram_determinant_fn;
  geometry->sqrt_face_gram_det = t8dg_linear_1D_sqrt_face_gram_determinant_fn;
  geometry->differential_invers_transpose = t8dg_linear_1D_differential_invers_transpose_fn;
  geometry->attribute_data_type = T8DG_TREE_VERTICES;
  geometry->data = NULL;
  return geometry;
}

static void
t8dg_linear_2D_geometry_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double vertex[DIM3],
                            double image_vertex[DIM3])
{

  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);

  t8_vec_axpyz (image_coarse_element_length_vector_x, tree_vertices, image_vertex, vertex[0]);
  t8_vec_axpy (image_coarse_element_length_vector_y, image_vertex, vertex[1]);
}

static void
t8dg_linear_2D_differential_invers_transpose_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                 const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                 double transformed_gradient_tangential_vector[3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);

  double              g11, g12, g22, det;

  g11 = t8_vec_dot (image_coarse_element_length_vector_x, image_coarse_element_length_vector_x);
  g12 = t8_vec_dot (image_coarse_element_length_vector_x, image_coarse_element_length_vector_y);
  g22 = t8_vec_dot (image_coarse_element_length_vector_y, image_coarse_element_length_vector_y);
  det = g11 * g22 - g12 * g12;

  t8_vec_axb (image_coarse_element_length_vector_x, transformed_gradient_tangential_vector,
              (g22 * coarse_tangential_vector[0] - g12 * coarse_tangential_vector[1]) / det, 0);
  t8_vec_axpy (image_coarse_element_length_vector_y, transformed_gradient_tangential_vector,
               (-g12 * coarse_tangential_vector[0] + g11 * coarse_tangential_vector[1]) / det);
}

static double
t8dg_linear_2D_sqrt_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double coarse_vertex[3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);

  double              g11, g12, g22;

  g11 = t8_vec_dot (image_coarse_element_length_vector_x, image_coarse_element_length_vector_x);
  g12 = t8_vec_dot (image_coarse_element_length_vector_x, image_coarse_element_length_vector_y);
  g22 = t8_vec_dot (image_coarse_element_length_vector_y, image_coarse_element_length_vector_y);

  return sqrt (g11 * g22 - g12 * g12);
}

static double
t8dg_linear_2D_sqrt_face_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const int iface,
                                              const double coarse_vertex[3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);

  switch (iface / 2) {
  case 0:
    return sqrt (t8_vec_dot (image_coarse_element_length_vector_y, image_coarse_element_length_vector_y));
    break;

  case 1:
    return sqrt (t8_vec_dot (image_coarse_element_length_vector_x, image_coarse_element_length_vector_x));
    break;

  default:
    T8DG_ABORT ("Facenumber too big");
    break;
  }
}

static void
t8dg_linear_3D_geometry_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double vertex[DIM3],
                            double image_vertex[DIM3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  double              image_coarse_element_length_vector_z[3];

  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 4 * DIM3, image_coarse_element_length_vector_z, -1);

  t8_vec_axpyz (image_coarse_element_length_vector_x, tree_vertices, image_vertex, vertex[0]);
  t8_vec_axpy (image_coarse_element_length_vector_y, image_vertex, vertex[1]);
  t8_vec_axpy (image_coarse_element_length_vector_z, image_vertex, vertex[2]);
}

static void
t8dg_linear_3D_differential_invers_transpose_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                 const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                 double transformed_gradient_tangential_vector[3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  double              image_coarse_element_length_vector_z[3];

  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 4 * DIM3, image_coarse_element_length_vector_z, -1);

  int                 idx;
/*TODO: improve geometry!!*/
  for (idx = 0; idx < DIM3; idx++) {
    transformed_gradient_tangential_vector[idx] = coarse_tangential_vector[idx];
  }
  return;
  T8DG_ABORT ("Not implemented \n ");
}

static double
t8dg_linear_3D_sqrt_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double coarse_vertex[3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  double              image_coarse_element_length_vector_z[3];
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 4 * DIM3, image_coarse_element_length_vector_z, -1);
  return 1;
  T8DG_ABORT ("Not implemented \n ");
//  return det;
}

static double
t8dg_linear_3D_sqrt_face_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const int iface,
                                              const double coarse_vertex[3])
{
  double             *tree_vertices = t8dg_forest_get_tree_vertices_gloidx (forest, iglobaltree);

  double              image_coarse_element_length_vector_x[3];
  double              image_coarse_element_length_vector_y[3];
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, image_coarse_element_length_vector_x, -1);
  t8_vec_axpyz (tree_vertices, tree_vertices + 2 * DIM3, image_coarse_element_length_vector_y, -1);
  return 1;
  T8DG_ABORT ("Not implemented \n ");

  switch (iface / 2) {
  case 0:
    return sqrt (t8_vec_dot (image_coarse_element_length_vector_y, image_coarse_element_length_vector_y));
    break;

  case 1:
    return sqrt (t8_vec_dot (image_coarse_element_length_vector_x, image_coarse_element_length_vector_x));
    break;

  default:
    T8DG_ABORT ("Facenumber too big");
    break;
  }
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_linear (int dim)
{
  switch (dim) {
  case 1:
    return t8dg_coarse_geometry_new_1D_linear ();
    break;
  case 2:
    return t8dg_coarse_geometry_new_2D_linear ();
    break;
  case 3:
    return t8dg_coarse_geometry_new_3D_linear ();

  default:
    T8DG_ABORT ("Not implemented \n ");
    break;
  }
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_3D_linear ()
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_linear_3D_geometry_fn;
  geometry->sqrt_gram_det = t8dg_linear_3D_sqrt_gram_determinant_fn;
  geometry->sqrt_face_gram_det = t8dg_linear_3D_sqrt_face_gram_determinant_fn;
  geometry->differential_invers_transpose = t8dg_linear_3D_differential_invers_transpose_fn;
  geometry->attribute_data_type = T8DG_TREE_VERTICES;
  geometry->data = NULL;
  return geometry;
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_2D_linear ()
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_linear_2D_geometry_fn;
  geometry->sqrt_gram_det = t8dg_linear_2D_sqrt_gram_determinant_fn;
  geometry->sqrt_face_gram_det = t8dg_linear_2D_sqrt_face_gram_determinant_fn;
  geometry->differential_invers_transpose = t8dg_linear_2D_differential_invers_transpose_fn;
  geometry->attribute_data_type = T8DG_TREE_VERTICES;
  geometry->data = NULL;
  return geometry;
}

void
t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_t ** pgeometry)
{
  T8DG_FREE ((*pgeometry)->data);
  T8DG_FREE (*pgeometry);
  *pgeometry = NULL;
}

typedef struct t8dg_circle_ring_data
{
  double              inner_radius;
  double              outer_radius;
} t8dg_circle_ring_data_t;

static double
t8dg_circle_ring_radius (double x, t8dg_circle_ring_data_t * data)
{
  return data->inner_radius + x * (data->outer_radius - data->inner_radius);
}

static double
t8dg_circle_ring_radius_derivative (double x, t8dg_circle_ring_data_t * data)
{
  return (data->outer_radius - data->inner_radius);
}

static double
t8dg_circle_ring_angle (double y, t8_gloidx_t iglobaltree, t8dg_circle_ring_data_t * data)
{
  return (y + iglobaltree) * M_PI_2;
}

static double
t8dg_circle_ring_angle_derivative (double y, t8_gloidx_t iglobaltree, t8dg_circle_ring_data_t * data)
{
  return M_PI_2;
}

static void
t8dg_circle_ring_2D_geometry_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double vertex[DIM3],
                                 double image_vertex[DIM3])
{
  image_vertex[0] = t8dg_circle_ring_radius (vertex[0], data) * cos (t8dg_circle_ring_angle (vertex[1], iglobaltree, data));
  image_vertex[1] = t8dg_circle_ring_radius (vertex[0], data) * sin (t8dg_circle_ring_angle (vertex[1], iglobaltree, data));
}

static double
t8dg_circle_ring_2D_sqrt_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                              const double coarse_vertex[3])
{
  return t8dg_circle_ring_radius (coarse_vertex[0], data) * t8dg_circle_ring_radius_derivative (coarse_vertex[0],
                                                                                                data) *
    t8dg_circle_ring_angle_derivative (coarse_vertex[1], iglobaltree, data);
}

static void
t8dg_circle_ring_2D_differential_invers_transpose_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                      const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                      double transformed_gradient_tangential_vector[3])
{
  double              radius, radius_der, angle, angle_der;
  double              a, b, c, d, x, y;
  double              det;
  x = coarse_vertex[0];
  y = coarse_vertex[1];
  radius_der = t8dg_circle_ring_radius_derivative (x, data);
  radius = t8dg_circle_ring_radius (x, data);
  angle = t8dg_circle_ring_angle (y, iglobaltree, data);
  angle_der = t8dg_circle_ring_angle_derivative (y, iglobaltree, data);
  a = radius_der * cos (angle);
  b = -radius * sin (angle) * angle_der;
  c = radius_der * sin (angle);
  d = radius * cos (angle) * angle_der;
  det = a * d - b * c;
  T8DG_ASSERT (fabs (det - t8dg_circle_ring_2D_sqrt_gram_determinant_fn (forest, iglobaltree, data, coarse_vertex)) < 1e-12);

  transformed_gradient_tangential_vector[0] = (d * coarse_tangential_vector[0] - c * coarse_tangential_vector[1]) / det;
  transformed_gradient_tangential_vector[1] = (-b * coarse_tangential_vector[0] + a * coarse_tangential_vector[1]) / det;
}

static double
t8dg_circle_ring_2D_sqrt_face_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const int iface,
                                                   const double coarse_vertex[3])
{
  switch (iface / 2) {
  case 0:
    return t8dg_circle_ring_angle_derivative (coarse_vertex[1], iglobaltree, data) * t8dg_circle_ring_radius (coarse_vertex[0], data);
    break;

  case 1:
    return t8dg_circle_ring_radius_derivative (coarse_vertex[0], data);
    break;

  default:
    T8DG_ABORT ("Facenumber too big");
    break;
  }
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_2D_circle_ring (double inner_radius, double outer_radius)
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_circle_ring_2D_geometry_fn;
  geometry->sqrt_gram_det = t8dg_circle_ring_2D_sqrt_gram_determinant_fn;
  geometry->sqrt_face_gram_det = t8dg_circle_ring_2D_sqrt_face_gram_determinant_fn;
  geometry->differential_invers_transpose = t8dg_circle_ring_2D_differential_invers_transpose_fn;
  geometry->attribute_data_type = T8DG_TREE_VERTICES;
  t8dg_circle_ring_data_t *circle_ring_data;
  circle_ring_data = T8DG_ALLOC_ZERO (t8dg_circle_ring_data_t, 1);
  circle_ring_data->inner_radius = inner_radius;
  circle_ring_data->outer_radius = outer_radius;

  geometry->data = circle_ring_data;
  return geometry;
}

static void
t8dg_cylinder_ring_3D_geometry_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const double vertex[DIM3],
                                   double image_vertex[DIM3])
{
  image_vertex[0] = t8dg_circle_ring_radius (vertex[0], data) * cos (t8dg_circle_ring_angle (vertex[1], iglobaltree, data));
  image_vertex[1] = t8dg_circle_ring_radius (vertex[0], data) * sin (t8dg_circle_ring_angle (vertex[1], iglobaltree, data));
  image_vertex[2] = vertex[2];
}

static double
t8dg_cylinder_ring_3D_sqrt_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                const double coarse_vertex[3])
{
  return t8dg_circle_ring_radius (coarse_vertex[0], data) * t8dg_circle_ring_radius_derivative (coarse_vertex[0],
                                                                                                data) *
    t8dg_circle_ring_angle_derivative (coarse_vertex[1], iglobaltree, data);
}

static void
t8dg_cylinder_ring_3D_differential_invers_transpose_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data,
                                                        const double coarse_vertex[3], const double coarse_tangential_vector[3],
                                                        double transformed_gradient_tangential_vector[3])
{
  double              radius, radius_der, angle, angle_der;
  double              a, b, c, d, x, y;
  double              det;
  x = coarse_vertex[0];
  y = coarse_vertex[1];
  radius_der = t8dg_circle_ring_radius_derivative (x, data);
  radius = t8dg_circle_ring_radius (x, data);
  angle = t8dg_circle_ring_angle (y, iglobaltree, data);
  angle_der = t8dg_circle_ring_angle_derivative (y, iglobaltree, data);
  a = radius_der * cos (angle);
  b = -radius * sin (angle) * angle_der;
  c = radius_der * sin (angle);
  d = radius * cos (angle) * angle_der;
  det = a * d - b * c;
  T8DG_ASSERT (fabs (det - t8dg_circle_ring_2D_sqrt_gram_determinant_fn (forest, iglobaltree, data, coarse_vertex)) < 1e-12);

  transformed_gradient_tangential_vector[0] = (d * coarse_tangential_vector[0] - c * coarse_tangential_vector[1]) / det;
  transformed_gradient_tangential_vector[1] = (-b * coarse_tangential_vector[0] + a * coarse_tangential_vector[1]) / det;
  transformed_gradient_tangential_vector[2] = (coarse_tangential_vector[2]);
}

static double
t8dg_cylinder_ring_3D_sqrt_face_gram_determinant_fn (const t8_forest_t forest, const t8_gloidx_t iglobaltree, void *data, const int iface,
                                                     const double coarse_vertex[3])
{
  switch (iface / 2) {
  case 0:
    return t8dg_circle_ring_angle_derivative (coarse_vertex[1], iglobaltree, data) * t8dg_circle_ring_radius (coarse_vertex[0], data);
    break;

  case 1:
    return t8dg_circle_ring_radius_derivative (coarse_vertex[0], data);
    break;
  case 2:
    return t8dg_circle_ring_2D_sqrt_gram_determinant_fn (forest, iglobaltree, data, coarse_vertex);
    break;

  default:
    T8DG_ABORT ("Facenumber too big");
    break;
  }
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_3D_cylinder_ring (double inner_radius, double outer_radius)
{
  t8dg_coarse_geometry_t *geometry = T8_ALLOC (t8dg_coarse_geometry_t, 1);
  geometry->geometry = t8dg_cylinder_ring_3D_geometry_fn;
  geometry->sqrt_gram_det = t8dg_cylinder_ring_3D_sqrt_gram_determinant_fn;
  geometry->sqrt_face_gram_det = t8dg_cylinder_ring_3D_sqrt_face_gram_determinant_fn;
  geometry->differential_invers_transpose = t8dg_cylinder_ring_3D_differential_invers_transpose_fn;
  geometry->attribute_data_type = T8DG_TREE_VERTICES;
  t8dg_circle_ring_data_t *circle_ring_data;
  circle_ring_data = T8DG_ALLOC_ZERO (t8dg_circle_ring_data_t, 1);
  circle_ring_data->inner_radius = inner_radius;
  circle_ring_data->outer_radius = outer_radius;

  geometry->data = circle_ring_data;
  return geometry;
}

t8dg_coarse_geometry_t *
t8dg_coarse_geometry_new_arg (int geometry_arg)
{
  switch (geometry_arg) {
  case 0:
    return t8dg_coarse_geometry_new_1D_linear ();
  case 1:
    return t8dg_coarse_geometry_new_2D_linear ();
  case 2:
    return t8dg_coarse_geometry_new_3D_linear ();
  case 3:
    return t8dg_coarse_geometry_new_2D_circle_ring (1, 2);
  case 4:
    return t8dg_coarse_geometry_new_3D_cylinder_ring (1, 2);

  default:
    break;
  }
  return NULL;
}
