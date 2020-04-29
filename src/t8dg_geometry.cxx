#include "t8dg.h"
#include "t8dg_geometry.h"
#include <t8_forest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>

static void
t8dg_geometry_transform_reference_vertex_to_coarse_vertex (const t8dg_geometry_transformation_data_t * data,
                                                           const double reference_vertex[DIM3], double coarse_vertex[DIM3])
{
  int                 idim;
  int                 coarse_vertex_int_translation_coords[3] = { 0, 0, 0 };
  int                 level;
  double              translation_vector[3] = { 0, 0, 0 };
  double              scaling_factor;
  double              length_inv;

  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;
  t8_element_t       *element;

  eclass = t8_forest_get_eclass (data->forest, data->itree);
  scheme = t8_forest_get_eclass_scheme (data->forest, eclass);
  element = t8_forest_get_element_in_tree (data->forest, data->itree, data->ielement);

  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -level);
  length_inv = 1. / scheme->t8_element_root_len (element);

  scheme->t8_element_vertex_coords (element, 0, coarse_vertex_int_translation_coords);
  for (idim = 0; idim < DIM3; idim++) {
    translation_vector[idim] = length_inv * coarse_vertex_int_translation_coords[idim];
  }

  /* For triangle reflection about x=y also needed */

  t8_vec_axpyz (reference_vertex, translation_vector, coarse_vertex, scaling_factor);
  /*hx+x_0 */
}

void
t8dg_geometry_transform_reference_vertex_to_image_vertex (const t8dg_geometry_transformation_data_t * geometry_data,
                                                          const double reference_vertex[3], double image_vertex[3])
{
  double              coarse_vertex[3];
  /* transform into coarse reference element */
  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (geometry_data, reference_vertex, coarse_vertex);

  /* tree vertices are application data for linear geometry */
  t8dg_coarse_geometry_apply (geometry_data->coarse_geometry, geometry_data->forest, geometry_data->itree, coarse_vertex, image_vertex);
}

double
t8dg_geometry_calculate_gram_determinant (const t8dg_geometry_transformation_data_t * geometry_data, const double reference_vertex[3])
{
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;
  t8_element_t       *element;

  int                 level;
  double              scaling_factor;
  double              coarse_vertex[3];

  eclass = t8_forest_get_eclass (geometry_data->forest, geometry_data->itree);
  scheme = t8_forest_get_eclass_scheme (geometry_data->forest, eclass);
  element = t8_forest_get_element_in_tree (geometry_data->forest, geometry_data->itree, geometry_data->ielement);

  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (geometry_data, reference_vertex, coarse_vertex);

  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -level);     /*TODO: dimension dependency ? */

  return scaling_factor * t8dg_coarse_geometry_calculate_gram_determinant (geometry_data->coarse_geometry, geometry_data->forest,
                                                                           geometry_data->itree, coarse_vertex);
}

void
t8dg_geometry_calculate_transformed_gradient_tangential_vector (const t8dg_geometry_transformation_data_t * geometry_data,
                                                                const double reference_vertex[3],
                                                                const double reference_tangential_vector[3],
                                                                double transformed_gradient_tangential_vector[3])
{
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;
  t8_element_t       *element;

  int                 level;
  double              scaling_factor;
  double              coarse_vertex[3];
  double              coarse_tangential_vector[3];

  eclass = t8_forest_get_eclass (geometry_data->forest, geometry_data->itree);
  scheme = t8_forest_get_eclass_scheme (geometry_data->forest, eclass);
  element = t8_forest_get_element_in_tree (geometry_data->forest, geometry_data->itree, geometry_data->ielement);

  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (geometry_data, reference_vertex, coarse_vertex);

  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -level);     /*TODO: dimension dependency ? */

  t8_vec_axb (reference_tangential_vector, coarse_tangential_vector, 1. / scaling_factor, 0);

  t8dg_coarse_geometry_calculate_gradient_tangential_vector (geometry_data->coarse_geometry, geometry_data->forest,
                                                             geometry_data->itree, coarse_vertex, coarse_tangential_vector,
                                                             transformed_gradient_tangential_vector);
}

void
t8dg_geometry_calculate_normal_vector (const t8dg_geometry_transformation_data_t * geometry_data, const int iface,
                                       const double reference_vertex[3], double image_normal_vector[3])
{
  t8_eclass_t         eclass;

  double              coarse_normal_vector[3] = { 0, 0, 0 };
  int                 side;
  int                 normal_direction;

  eclass = t8_forest_get_eclass (geometry_data->forest, geometry_data->itree);
  if (eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX) {
    /* For Hypercubes, there is no change of the normal vector from reference element to coarse element */
    side = iface % 2;
    normal_direction = iface / 2;
    coarse_normal_vector[normal_direction] = (side) ? 1 : -1;
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
  t8dg_coarse_geometry_transform_normal_vector (geometry_data->coarse_geometry, geometry_data->forest, geometry_data->itree,
                                                reference_vertex, coarse_normal_vector, image_normal_vector);
  t8_vec_ax (image_normal_vector, 1. / t8_vec_norm (image_normal_vector));
}

double
t8dg_geometry_calculate_face_gram_determinant (const t8dg_geometry_transformation_data_t * geometry_data,
                                               const int iface, const double reference_vertex[3])
{
  /*TODO: call coarse function for higher dimension */
  return 1;
}
