/*
 * t8dg_local_values.c
 *
 *  Created on: Apr 5, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_local_values.h"
#include "t8dg_sc_array.h"
#include "t8dg_geometry.h"
#include "t8dg_flux.h"
#include "t8dg_global_values.h"
#include "t8dg_dof.h"
#include "t8dg_quad.h"

#include <t8_vec.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>

struct t8dg_local_values
{
  t8dg_global_values_t **global_values;

  t8_forest_t         forest;

  t8dg_coarse_geometry_t *coarse_geometry;

  int                 max_num_faces, max_num_elem_values, max_num_face_values, dim;

  t8dg_quad_values_t *element_trafo_quad_weight;                        /**< For each element, an array of length num_quad_points with the combined weight
								of quadrature and transformation */
  t8dg_quad_values_t *element_transformed_gradient_tangential_vectors[DIM3][DIM3];    /**< For each element, (d_x F^-1)(tau_d), needed to calculate the gradient */

  /* To avoid another dimension when flattening the array, and since the number of Faces is bound,
   *  have an sc_array for each faceindex of length
   * of the processorlocal elements*/
  t8dg_quad_values_t *face_trafo_quad_weight[MAX_FACES];                /**< for each element, an array of length num_quad_points with the combined weight
								of quadrature and transformation */

  t8dg_dof_values_t  *face_normal_vectors[MAX_FACES][DIM3];                  /**< For each face and nodal basis point the 3-dim normal vector in image space*/
};

/*
void
t8dg_local_values_print_debug (const t8dg_local_values_t * values, const t8_locidx_t idata)
{
  int                 iface;
  sc_array_t         *tmp;
  t8dg_debugf ("Element trafo quad weight: \n");
  tmp = t8dg_sc_array_block_double_new_view (values->element_trafo_quad_weight, idata);
  t8dg_sc_array_block_double_print (tmp);
  sc_array_destroy (tmp);
  t8dg_debugf ("element_transformed_gradient_tangential_vectors: \n");
  tmp = t8dg_sc_array_block_double_new_view (values->element_transformed_gradient_tangential_vectors, idata);
  t8dg_sc_array_block_double_print (tmp);
  sc_array_destroy (tmp);
  for (iface = 0; iface < values->max_num_faces; iface++) {
    t8dg_debugf ("face trafo quad weight: \n");
    tmp = t8dg_sc_array_block_double_new_view (values->face_trafo_quad_weight[iface], idata);
    t8dg_sc_array_block_double_print (tmp);
    sc_array_destroy (tmp);
    t8dg_debugf ("face_normal_vector: \n");
    tmp = t8dg_sc_array_block_double_new_view (values->face_normal_vectors[iface], idata);
    t8dg_sc_array_block_double_print (tmp);
    sc_array_destroy (tmp);
  }
}
*/

static t8dg_face_dof_values_t *
t8dg_local_values_get_face_trafo_quad_weights_view (t8dg_local_values_t * values, const t8_locidx_t itree, const t8_locidx_t ielement,
                                                    const int faceindex)
{
  return t8dg_quad_values_new_face_quad_values_view (values->face_trafo_quad_weight[faceindex], faceindex, itree, ielement);
}

static t8dg_face_dof_values_t *
t8dg_local_values_get_face_trafo_quad_weights_view_idata_eclass (t8dg_local_values_t * values, const t8_locidx_t idata,
                                                                 const t8_eclass_t element_eclass, const int faceindex)
{
  return t8dg_quad_values_new_face_quad_values_view_idata_eclass (values->face_trafo_quad_weight[faceindex], faceindex, idata,
                                                                  element_eclass);
}

static t8dg_element_dof_values_t *
t8dg_local_values_get_element_trafo_quad_weights_view (const t8dg_local_values_t * values, const t8_locidx_t itree,
                                                       const t8_locidx_t ielement)
{
  return t8dg_quad_values_new_element_quad_values_view (values->element_trafo_quad_weight, itree, ielement);
}

void
t8dg_local_values_set_transformed_gradient_tangential_vector (const t8dg_local_values_t * values,
                                                              const t8_locidx_t itree, const t8_locidx_t ielement, const int iquad,
                                                              const int idim, double vector[3])
{
  int                 icomp;
  t8dg_element_quad_values_t *element_transformed_gradient_tangential_component;
  for (icomp = 0; icomp < DIM3; icomp++) {
    element_transformed_gradient_tangential_component =
      t8dg_quad_values_new_element_quad_values_view (values->element_transformed_gradient_tangential_vectors[idim][icomp], itree, ielement);
    t8dg_element_quad_values_set_value (element_transformed_gradient_tangential_component, iquad, vector[icomp]);
    t8dg_element_quad_values_destroy (&element_transformed_gradient_tangential_component);
  }
}

void
t8dg_local_values_get_transformed_gradient_tangential_vector (const t8dg_local_values_t * values,
                                                              const t8_locidx_t itree, const t8_locidx_t ielement, const int iquad,
                                                              const int idim, double vector[3])
{
  int                 icomp;
  t8dg_element_quad_values_t *element_transformed_gradient_tangential_component;
  for (icomp = 0; icomp < DIM3; icomp++) {
    element_transformed_gradient_tangential_component =
      t8dg_quad_values_new_element_quad_values_view (values->element_transformed_gradient_tangential_vectors[idim][icomp], itree, ielement);
    vector[icomp] = t8dg_element_quad_values_get_value (element_transformed_gradient_tangential_component, iquad);
    t8dg_element_quad_values_destroy (&element_transformed_gradient_tangential_component);
  }
}

void
t8dg_local_values_set_face_normal_vector (const t8dg_local_values_t * values,
                                          const t8_locidx_t itree, const t8_locidx_t ielement, const int iface, const int idof,
                                          double vector[3])
{
  int                 icomp;
  t8dg_face_dof_values_t *face_normal_component;
  for (icomp = 0; icomp < DIM3; icomp++) {
    face_normal_component = t8dg_dof_values_new_face_dof_values_view (values->face_normal_vectors[iface][icomp], iface, itree, ielement);
    t8dg_face_dof_values_set_value (face_normal_component, idof, vector[icomp]);
    t8dg_face_dof_values_destroy (&face_normal_component);
  }
}

void
t8dg_local_values_set_face_normal_vector_idata_eclass (const t8dg_local_values_t * values,
                                                       const t8_locidx_t idata, const t8_eclass_t element_eclass, const int iface,
                                                       const int idof, double vector[3])
{
  int                 icomp;
  t8dg_face_dof_values_t *face_normal_component;
  for (icomp = 0; icomp < DIM3; icomp++) {
    face_normal_component =
      t8dg_dof_values_new_face_dof_values_view_idata_eclass (values->face_normal_vectors[iface][icomp], iface, idata, element_eclass);
    t8dg_element_dof_values_set_value (face_normal_component, idof, vector[icomp]);
    t8dg_element_dof_values_destroy (&face_normal_component);
  }
}

void
t8dg_local_values_get_face_normal_vector (const t8dg_local_values_t * values,
                                          const t8_locidx_t itree, const t8_locidx_t ielement, const int iface, const int idof,
                                          double vector[3])
{
  int                 icomp;
  t8dg_face_dof_values_t *face_normal_component;
  for (icomp = 0; icomp < DIM3; icomp++) {
    face_normal_component = t8dg_dof_values_new_face_dof_values_view (values->face_normal_vectors[iface][icomp], iface, itree, ielement);
    vector[icomp] = t8dg_face_dof_values_get_value (face_normal_component, idof);
    t8dg_face_dof_values_destroy (&face_normal_component);
  }
}

void
t8dg_local_values_get_face_normal_vector_idata_eclass (const t8dg_local_values_t * values,
                                                       const t8_locidx_t idata, const t8_eclass_t element_eclass, const int iface,
                                                       const int idof, double vector[3])
{
  int                 icomp;
  t8dg_face_dof_values_t *face_normal_component;
  for (icomp = 0; icomp < DIM3; icomp++) {
    face_normal_component =
      t8dg_dof_values_new_face_dof_values_view_idata_eclass (values->face_normal_vectors[iface][icomp], iface, idata, element_eclass);
    vector[icomp] = t8dg_face_dof_values_get_value (face_normal_component, idof);
    t8dg_face_dof_values_destroy (&face_normal_component);
  }
}

void
t8dg_local_values_copy_element_values (t8dg_local_values_t * src_values, t8_locidx_t src_idata,
                                       t8dg_local_values_t * dest_values, t8_locidx_t dest_idata)
{
  int                 iface, idim, icomp;

  t8dg_quad_values_copy_from_index_to_index (src_values->element_trafo_quad_weight, src_idata,
                                             dest_values->element_trafo_quad_weight, dest_idata);

  for (idim = 0; idim < src_values->dim; idim++) {
    for (icomp = 0; icomp < DIM3; icomp++) {
      t8dg_quad_values_copy_from_index_to_index (src_values->element_transformed_gradient_tangential_vectors[idim][icomp], src_idata,
                                                 dest_values->element_transformed_gradient_tangential_vectors[idim][icomp], dest_idata);

    }
  }

  for (iface = 0; iface < src_values->max_num_faces; iface++) {
    t8dg_quad_values_copy_from_index_to_index (src_values->face_trafo_quad_weight[iface], src_idata,
                                               dest_values->face_trafo_quad_weight[iface], dest_idata);

    for (icomp = 0; icomp < DIM3; icomp++) {

      t8dg_dof_values_copy_from_index_to_index (src_values->face_normal_vectors[iface][icomp], src_idata,
                                                dest_values->face_normal_vectors[iface][icomp], dest_idata);
    }
  }
}

void
t8dg_local_values_set_all_ghost_elements (t8dg_local_values_t * local_values)
{
  t8_locidx_t         ighosttree, ielement, num_ghost_trees, num_elements_in_ghost_tree;

  num_ghost_trees = t8_forest_get_num_ghost_trees (local_values->forest);
  for (ighosttree = 0; ighosttree < num_ghost_trees; ighosttree++) {
    num_elements_in_ghost_tree = t8_forest_ghost_tree_num_elements (local_values->forest, ighosttree);
    for (ielement = 0; ielement < num_elements_in_ghost_tree; ielement++) {
      t8dg_local_values_set_ghost_element (local_values, ighosttree, ielement);
    }
  }
}

void
t8dg_local_values_set_ghost_element (t8dg_local_values_t * local_values, t8_locidx_t ighosttree, t8_locidx_t ielement)
{
  double              sqrt_gram_det, face_sqrt_gram_det;
  int                 iquad, iface, idim, idof;

  /*pointer on the values to fill */
  t8dg_face_quad_values_t *face_trafo_quad;     /*size: num_face_quad */
  double              image_normal_vector[3];   /*size: 3, gets evaluated for each face dof point */

  double              reference_vertex[3];
  double              reference_tangential_vector[3];

  int                 num_face_quad;
  int                 num_faces, num_face_dof;

  t8dg_global_values_t *global_values;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;
  t8_element_t       *element;
  element = t8_forest_ghost_get_element (local_values->forest, ighosttree, ielement);

  t8_gloidx_t         iglobaltree = t8_forest_ghost_get_global_treeid (local_values->forest, ighosttree);
  t8_eclass_t         element_eclass = t8dg_eclass_from_gloidx_element (local_values->forest, iglobaltree, element);

  t8_locidx_t         idata = t8dg_ighosttree_ielement_to_idata (local_values->forest, ighosttree, ielement);

  global_values = local_values->global_values[element_eclass];
  quadrature = t8dg_global_values_get_quadrature (global_values);
  functionbasis = t8dg_global_values_get_functionbasis (global_values);

  num_faces = t8dg_global_values_get_num_faces (global_values);

  for (iface = 0; iface < num_faces; iface++) {
    face_trafo_quad = t8dg_local_values_get_face_trafo_quad_weights_view_idata_eclass (local_values, idata, element_eclass, iface);
    num_face_quad = t8dg_quadrature_get_num_face_vertices (quadrature, iface);
    for (iquad = 0; iquad < num_face_quad; iquad++) {
      t8dg_quadrature_get_face_vertex (quadrature, iface, iquad, reference_vertex);

      face_sqrt_gram_det =
        t8dg_geometry_calculate_face_sqrt_gram_determinant (local_values->coarse_geometry, local_values->forest, iglobaltree, element,
                                                            iface, reference_vertex);
      t8dg_face_quad_values_set_value (face_trafo_quad, iquad,
                                       face_sqrt_gram_det * t8dg_quadrature_get_face_weight (quadrature, iface, iquad));
    }

    num_face_dof = t8dg_functionbasis_get_num_face_dof (functionbasis, iface);
    for (idof = 0; idof < num_face_dof; idof++) {
      t8dg_functionbasis_get_lagrange_face_vertex (functionbasis, iface, idof, reference_vertex);
      t8dg_geometry_calculate_normal_vector (local_values->coarse_geometry, local_values->forest, iglobaltree, element, iface,
                                             reference_vertex, image_normal_vector);
      t8dg_local_values_set_face_normal_vector_idata_eclass (local_values, idata, element_eclass, iface, idof, image_normal_vector);
    }
    t8dg_face_quad_values_destroy (&face_trafo_quad);
  }
}

void
t8dg_local_values_set_all_local_elements (t8dg_local_values_t * local_values)
{
  t8_locidx_t         num_trees, itree, num_elems_in_tree, ielement;
  num_trees = t8_forest_get_num_local_trees (local_values->forest);
  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (local_values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      t8dg_local_values_set_element (local_values, itree, ielement);
    }
  }
}

void
t8dg_local_values_set_all_elements (t8dg_local_values_t * local_values)
{
  t8dg_local_values_set_all_ghost_elements (local_values);
  t8dg_local_values_set_all_local_elements (local_values);
}

void
t8dg_local_values_set_element (t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement)
{
  double              sqrt_gram_det, face_sqrt_gram_det;
  int                 iquad, iface, idim, idof;

  /*pointer on the values to fill */
  t8dg_element_quad_values_t *element_trafo_quad;       /*size: num_elem_quad */
  t8dg_face_quad_values_t *face_trafo_quad;     /*size: num_face_quad */
  double              transformed_gradient_tangential_vector[3];        /*size: 3, gets evaluated for each element quad point */
  double              image_normal_vector[3];   /*size: 3, gets evaluated for each face dof point */

  double              reference_vertex[3];
  double              reference_tangential_vector[3];

  int                 num_elem_quad, num_face_quad;
  int                 num_faces, num_face_dof;

  t8dg_global_values_t *global_values;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;

  t8_gloidx_t         iglobaltree = t8_forest_global_tree_id (local_values->forest, itree);
  t8_element_t       *element = t8_forest_get_element_in_tree (local_values->forest, itree, ielement);

  global_values = t8dg_global_values_array_get_global_values (local_values->global_values, local_values->forest, itree, ielement);
  quadrature = t8dg_global_values_get_quadrature (global_values);
  functionbasis = t8dg_global_values_get_functionbasis (global_values);

  num_elem_quad = t8dg_global_values_get_num_elem_quad (global_values);
  num_faces = t8dg_global_values_get_num_faces (global_values);

  element_trafo_quad = t8dg_local_values_get_element_trafo_quad_weights_view (local_values, itree, ielement);

  for (iquad = 0; iquad < num_elem_quad; iquad++) {
    t8dg_quadrature_get_element_vertex (quadrature, iquad, reference_vertex);
    sqrt_gram_det =
      t8dg_geometry_calculate_sqrt_gram_determinant (local_values->coarse_geometry, local_values->forest, iglobaltree, element,
                                                     reference_vertex);

    t8dg_element_quad_values_set_value (element_trafo_quad, iquad, sqrt_gram_det * t8dg_quadrature_get_element_weight (quadrature, iquad));

    for (idim = 0; idim < local_values->dim; idim++) {
      reference_tangential_vector[0] = reference_tangential_vector[1] = reference_tangential_vector[2] = 0;
      reference_tangential_vector[idim] = 1;

      t8dg_geometry_calculate_transformed_gradient_tangential_vector (local_values->coarse_geometry, local_values->forest, iglobaltree,
                                                                      element, reference_vertex, reference_tangential_vector,
                                                                      transformed_gradient_tangential_vector);
      t8dg_local_values_set_transformed_gradient_tangential_vector (local_values, itree, ielement, iquad, idim,
                                                                    transformed_gradient_tangential_vector);
    }
  }
  t8dg_element_quad_values_destroy (&element_trafo_quad);

  for (iface = 0; iface < num_faces; iface++) {
    face_trafo_quad = t8dg_local_values_get_face_trafo_quad_weights_view (local_values, itree, ielement, iface);
    num_face_quad = t8dg_quadrature_get_num_face_vertices (quadrature, iface);
    for (iquad = 0; iquad < num_face_quad; iquad++) {
      t8dg_quadrature_get_face_vertex (quadrature, iface, iquad, reference_vertex);
      face_sqrt_gram_det =
        t8dg_geometry_calculate_face_sqrt_gram_determinant (local_values->coarse_geometry, local_values->forest, iglobaltree, element,
                                                            iface, reference_vertex);
      t8dg_face_quad_values_set_value (face_trafo_quad, iquad,
                                       face_sqrt_gram_det * t8dg_quadrature_get_face_weight (quadrature, iface, iquad));
    }

    num_face_dof = t8dg_functionbasis_get_num_face_dof (functionbasis, iface);
    for (idof = 0; idof < num_face_dof; idof++) {
      t8dg_functionbasis_get_lagrange_face_vertex (functionbasis, iface, idof, reference_vertex);
      t8dg_geometry_calculate_normal_vector (local_values->coarse_geometry, local_values->forest, iglobaltree, element, iface,
                                             reference_vertex, image_normal_vector);
      t8dg_local_values_set_face_normal_vector (local_values, itree, ielement, iface, idof, image_normal_vector);
    }
    t8dg_face_quad_values_destroy (&face_trafo_quad);
  }
}

t8dg_local_values_t *
t8dg_local_values_new (t8_forest_t forest, t8dg_global_values_t ** global_values_array, t8dg_coarse_geometry_t * coarse_geometry)
{
  t8dg_local_values_t *local_values = T8DG_ALLOC_ZERO (t8dg_local_values_t, 1);

  int                 iface, idim, icomp;
  int                 eclass;

  local_values->forest = forest;
  local_values->coarse_geometry = coarse_geometry;

  /*TODO: Move to global_values_array */
  for (eclass = 0; eclass < (int) T8_ECLASS_COUNT; eclass++) {
    if (global_values_array[eclass] != NULL) {
      local_values->max_num_faces = SC_MAX (t8dg_global_values_get_num_faces (global_values_array[eclass]), local_values->max_num_faces);
      local_values->max_num_elem_values =
        SC_MAX (t8dg_global_values_get_num_elem_quad (global_values_array[eclass]), local_values->max_num_elem_values);
      local_values->max_num_face_values =
        SC_MAX (t8dg_global_values_get_max_num_facevalues (global_values_array[eclass]), local_values->max_num_face_values);

      /* Only elements of the same dimension ? */
      local_values->dim = t8dg_global_values_get_dim (global_values_array[eclass]);
    }
  }
  local_values->global_values = global_values_array;

  /*for each element an array of double values */
  local_values->element_trafo_quad_weight = t8dg_quad_values_new_local (forest, global_values_array);   /*? */

  for (idim = 0; idim < local_values->dim; idim++) {
    for (icomp = 0; icomp < DIM3; icomp++) {
      local_values->element_transformed_gradient_tangential_vectors[idim][icomp] = t8dg_quad_values_new_local (forest, global_values_array);
    }
  }

  for (iface = 0; iface < local_values->max_num_faces; iface++) {
    local_values->face_trafo_quad_weight[iface] = t8dg_quad_values_new (forest, global_values_array);

    for (icomp = 0; icomp < DIM3; icomp++) {
      local_values->face_normal_vectors[iface][icomp] = t8dg_dof_values_new (forest, global_values_array);
    }
  }
  t8_locidx_t         itree, ielement, num_trees, num_elems_in_tree;

  num_trees = t8_forest_get_num_local_trees (local_values->forest);

  //set_elements selber

  return local_values;
}

void
t8dg_local_values_destroy (t8dg_local_values_t ** pvalues)
{
  t8dg_local_values_t *values = *pvalues;
  int                 iface, idim, icomp;
  t8dg_quad_values_destroy (&values->element_trafo_quad_weight);
  for (idim = 0; idim < values->dim; idim++) {
    for (icomp = 0; icomp < DIM3; icomp++) {
      t8dg_quad_values_destroy (&values->element_transformed_gradient_tangential_vectors[idim][icomp]);
    }
  }
  for (iface = 0; iface < values->max_num_faces; iface++) {
    t8dg_quad_values_destroy (&values->face_trafo_quad_weight[iface]);
    for (icomp = 0; icomp < DIM3; icomp++) {
      t8dg_dof_values_destroy (&values->face_normal_vectors[iface][icomp]);
    }
  }
  T8DG_FREE (values);
  *pvalues = NULL;
}

void
t8dg_local_values_element_multiply_trafo_quad_weight (const t8dg_local_values_t * local_values,
                                                      t8_locidx_t itree, t8_locidx_t ielement,
                                                      t8dg_element_quad_values_t * src_element_quad,
                                                      t8dg_element_quad_values_t * dest_element_quad)
{
  int                 iquad, num_element_quad;
  t8dg_element_quad_values_t *element_trafo_quad_weights;

  double              trafo_quad_weight;
  double              old_value;

  element_trafo_quad_weights = t8dg_quad_values_new_element_quad_values_view (local_values->element_trafo_quad_weight, itree, ielement);

  num_element_quad = t8dg_element_quad_values_get_num_element_quad_points (element_trafo_quad_weights);
  T8DG_ASSERT (num_element_quad <= local_values->max_num_elem_values);

  for (iquad = 0; iquad < num_element_quad; iquad++) {
    trafo_quad_weight = t8dg_element_quad_values_get_value (element_trafo_quad_weights, iquad);
    old_value = t8dg_element_quad_values_get_value (src_element_quad, iquad);

    t8dg_element_quad_values_set_value (dest_element_quad, iquad, old_value * trafo_quad_weight);
  }
  t8dg_element_quad_values_destroy (&element_trafo_quad_weights);
}

void
t8dg_local_values_element_divide_trafo_quad_weight (const t8dg_local_values_t * local_values,
                                                    t8_locidx_t itree, t8_locidx_t ielement, t8dg_element_quad_values_t * src_element_quad,
                                                    t8dg_element_quad_values_t * dest_element_quad)
{
  int                 iquad, num_element_quad;
  t8dg_element_quad_values_t *element_trafo_quad_weights;

  double              trafo_quad_weight;
  double              old_value;

  element_trafo_quad_weights = t8dg_quad_values_new_element_quad_values_view (local_values->element_trafo_quad_weight, itree, ielement);

  num_element_quad = t8dg_element_quad_values_get_num_element_quad_points (element_trafo_quad_weights);
  T8DG_ASSERT (num_element_quad <= local_values->max_num_elem_values);

  for (iquad = 0; iquad < num_element_quad; iquad++) {
    trafo_quad_weight = t8dg_element_quad_values_get_value (element_trafo_quad_weights, iquad);
    old_value = t8dg_element_quad_values_get_value (src_element_quad, iquad);

    t8dg_element_quad_values_set_value (dest_element_quad, iquad, old_value / trafo_quad_weight);
  }
  t8dg_element_quad_values_destroy (&element_trafo_quad_weights);
}

void
t8dg_local_values_face_multiply_trafo_quad_weight (const t8dg_local_values_t * local_values,
                                                   t8_locidx_t itree, t8_locidx_t ielement, const int iface,
                                                   t8dg_face_quad_values_t * src_face_quad, t8dg_face_quad_values_t * dest_face_quad)
{
  int                 ifacequad, num_face_quad;
  t8dg_face_quad_values_t *face_trafo_quad_weights;
  double              trafo_quad_weight, old_value;

  face_trafo_quad_weights =
    t8dg_quad_values_new_face_quad_values_view (local_values->face_trafo_quad_weight[iface], iface, itree, ielement);

  num_face_quad = t8dg_face_quad_values_get_num_face_quad_points (face_trafo_quad_weights);

  T8DG_ASSERT (num_face_quad <= local_values->max_num_face_values);

  for (ifacequad = 0; ifacequad < num_face_quad; ifacequad++) {
    trafo_quad_weight = t8dg_face_quad_values_get_value (face_trafo_quad_weights, ifacequad);
    old_value = t8dg_face_quad_values_get_value (src_face_quad, ifacequad);
    t8dg_face_quad_values_set_value (dest_face_quad, ifacequad, old_value * trafo_quad_weight);
  }
  t8dg_face_quad_values_destroy (&face_trafo_quad_weights);
}

void
t8dg_local_values_face_divide_trafo_quad_weight (const t8dg_local_values_t * local_values,
                                                 t8_locidx_t itree, t8_locidx_t ielement, const int iface,
                                                 t8dg_face_quad_values_t * src_face_quad, t8dg_face_quad_values_t * dest_face_quad)
{
  int                 ifacequad, num_face_quad;
  t8dg_face_quad_values_t *face_trafo_quad_weights;
  double              trafo_quad_weight, old_value;

  face_trafo_quad_weights =
    t8dg_quad_values_new_face_quad_values_view (local_values->face_trafo_quad_weight[iface], iface, itree, ielement);

  num_face_quad = t8dg_face_quad_values_get_num_face_quad_points (face_trafo_quad_weights);

  T8DG_ASSERT (num_face_quad <= local_values->max_num_face_values);

  for (ifacequad = 0; ifacequad < num_face_quad; ifacequad++) {
    trafo_quad_weight = t8dg_face_quad_values_get_value (face_trafo_quad_weights, ifacequad);
    old_value = t8dg_face_quad_values_get_value (src_face_quad, ifacequad);
    t8dg_face_quad_values_set_value (dest_face_quad, ifacequad, old_value / trafo_quad_weight);
  }
  t8dg_face_quad_values_destroy (&face_trafo_quad_weights);
}

void
t8dg_local_values_face_divide_trafo_quad_weight_idata_eclass (const t8dg_local_values_t * local_values,
                                                              t8_locidx_t idata, t8_eclass_t element_eclass, const int iface,
                                                              t8dg_face_quad_values_t * src_face_quad,
                                                              t8dg_face_quad_values_t * dest_face_quad)
{
  int                 ifacequad, num_face_quad;
  t8dg_face_quad_values_t *face_trafo_quad_weights;
  double              trafo_quad_weight, old_value;

  face_trafo_quad_weights =
    t8dg_quad_values_new_face_quad_values_view_idata_eclass (local_values->face_trafo_quad_weight[iface], iface, idata, element_eclass);

  num_face_quad = t8dg_face_quad_values_get_num_face_quad_points (face_trafo_quad_weights);

  T8DG_ASSERT (num_face_quad <= local_values->max_num_face_values);

  for (ifacequad = 0; ifacequad < num_face_quad; ifacequad++) {
    trafo_quad_weight = t8dg_face_quad_values_get_value (face_trafo_quad_weights, ifacequad);
    old_value = t8dg_face_quad_values_get_value (src_face_quad, ifacequad);
    t8dg_face_quad_values_set_value (dest_face_quad, ifacequad, old_value / trafo_quad_weight);
  }
  t8dg_face_quad_values_destroy (&face_trafo_quad_weights);
}

void
t8dg_local_values_face_multiply_trafo_quad_weight_idata_eclass (const t8dg_local_values_t * local_values,
                                                                t8_locidx_t idata, t8_eclass_t element_eclass, const int iface,
                                                                t8dg_face_quad_values_t * src_face_quad,
                                                                t8dg_face_quad_values_t * dest_face_quad)
{
  int                 ifacequad, num_face_quad;
  t8dg_face_quad_values_t *face_trafo_quad_weights;
  double              trafo_quad_weight, old_value;

  face_trafo_quad_weights =
    t8dg_quad_values_new_face_quad_values_view_idata_eclass (local_values->face_trafo_quad_weight[iface], iface, idata, element_eclass);

  num_face_quad = t8dg_face_quad_values_get_num_face_quad_points (face_trafo_quad_weights);

  T8DG_ASSERT (num_face_quad <= local_values->max_num_face_values);

  for (ifacequad = 0; ifacequad < num_face_quad; ifacequad++) {
    trafo_quad_weight = t8dg_face_quad_values_get_value (face_trafo_quad_weights, ifacequad);
    old_value = t8dg_face_quad_values_get_value (src_face_quad, ifacequad);
    t8dg_face_quad_values_set_value (dest_face_quad, ifacequad, old_value * trafo_quad_weight);
  }
  t8dg_face_quad_values_destroy (&face_trafo_quad_weights);
}

void
t8dg_local_values_partition (t8dg_local_values_t * local_values_old, t8dg_local_values_t * local_values_partition)
{
  int                 iface, idim, icomp;

  t8dg_quad_values_partition (local_values_old->element_trafo_quad_weight, local_values_partition->element_trafo_quad_weight);

  for (idim = 0; idim < local_values_old->dim; idim++) {
    for (icomp = 0; icomp < DIM3; icomp++) {
      t8dg_quad_values_partition (local_values_old->element_transformed_gradient_tangential_vectors[idim][icomp],
                                  local_values_partition->element_transformed_gradient_tangential_vectors[idim][icomp]);

    }
  }

  for (iface = 0; iface < local_values_old->max_num_faces; iface++) {
    t8dg_quad_values_partition (local_values_old->face_trafo_quad_weight[iface], local_values_partition->face_trafo_quad_weight[iface]);

    for (icomp = 0; icomp < DIM3; icomp++) {
      t8dg_dof_values_partition (local_values_old->face_normal_vectors[iface][icomp],
                                 local_values_partition->face_normal_vectors[iface][icomp]);
    }
  }
}

void
t8dg_local_values_ghost_exchange (t8dg_local_values_t * local_values)
{
  int                 iface, idim;
  for (iface = 0; iface < local_values->max_num_faces; iface++) {
    t8dg_quad_values_ghost_exchange (local_values->face_trafo_quad_weight[iface]);
    for (idim = 0; idim < DIM3; idim++) {
      t8dg_dof_values_ghost_exchange (local_values->face_normal_vectors[iface][idim]);
    }
  }
}

void
t8dg_local_values_element_multiply_directional_transformed_gradient_tangential_vector_component (const t8dg_local_values_t * local_values,
                                                                                                 t8_locidx_t itree, t8_locidx_t ielement,
                                                                                                 int derivative_direction, int icomp,
                                                                                                 t8dg_element_quad_values_t *
                                                                                                 src_element_quad_values,
                                                                                                 t8dg_element_quad_values_t *
                                                                                                 dest_element_quad_values)
{
  int                 iquad, num_quad_vertices;
  double              transformed_gradient_tangential_vector[3];
  double              flux_value;

  t8dg_global_values_t *global_values;
  global_values = t8dg_global_values_array_get_global_values (local_values->global_values, local_values->forest, itree, ielement);
  num_quad_vertices = t8dg_global_values_get_num_elem_quad (global_values);

  for (iquad = 0; iquad < num_quad_vertices; iquad++) {
    t8dg_local_values_get_transformed_gradient_tangential_vector (local_values, itree, ielement, iquad, derivative_direction,
                                                                  transformed_gradient_tangential_vector);

    flux_value = t8dg_element_quad_values_get_value (src_element_quad_values, iquad) * transformed_gradient_tangential_vector[icomp];

    t8dg_element_quad_values_set_value (dest_element_quad_values, iquad, flux_value);
  }
}

void
t8dg_local_values_apply_element_component_stiffness_matrix_dof (t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement,
                                                                int icomp, t8dg_element_dof_values_t * src_element_dof,
                                                                t8dg_element_dof_values_t * dest_element_dof)
{
  t8dg_global_values_t *global_values;
  int                 idirection;
  t8dg_quadidx_t      num_quad_vertices;

  t8dg_element_dof_values_t *element_res_dof_values;
  t8dg_element_dof_values_t *element_derivative_dof_values;
  t8dg_element_dof_values_t *summand;

  t8dg_element_quad_values_t *element_quad_values;
  t8dg_element_quad_values_t *element_derivative_quad_values;

  global_values = t8dg_local_values_get_global_values (local_values, itree, ielement);
  num_quad_vertices = t8dg_global_values_get_num_elem_quad (global_values);

  element_quad_values = t8dg_element_quad_values_new (num_quad_vertices);
  element_derivative_quad_values = t8dg_element_quad_values_new (num_quad_vertices);

  element_res_dof_values = t8dg_element_dof_values_duplicate (dest_element_dof);
  element_derivative_dof_values = t8dg_element_dof_values_duplicate (dest_element_dof);
  summand = t8dg_element_dof_values_duplicate (dest_element_dof);

  t8dg_element_dof_values_set_zero (element_res_dof_values);

//      t8_debugf ("Calculate stiffness matrix for itree: %i, ielement: %i\n", itree, ielement);

//      t8_debugf ("src_element_dof\n");
//      t8dg_element_dof_values_debug_print (src_element_dof);

  t8dg_global_values_transform_element_dof_to_element_quad (global_values, src_element_dof, element_quad_values);

  t8dg_local_values_element_multiply_trafo_quad_weight (local_values, itree, ielement, element_quad_values, element_quad_values);

  //     t8_debugf ("element_quad_values*qtw\n");
//      t8dg_element_dof_values_debug_print (element_quad_values);

  for (idirection = 0; idirection < local_values->dim; idirection++) {
    t8dg_local_values_element_multiply_directional_transformed_gradient_tangential_vector_component (local_values, itree, ielement,
                                                                                                     idirection, icomp, element_quad_values,
                                                                                                     element_derivative_quad_values);

    t8dg_global_values_transform_element_quad_to_element_dof (global_values, element_derivative_quad_values, element_derivative_dof_values);

//          t8_debugf ("element_dof_derivative_values\n");
//          t8dg_element_dof_values_debug_print (element_derivative_dof_values);

    t8dg_global_values_element_apply_derivative_matrix_transpose (global_values, idirection, element_derivative_dof_values, summand);
//          t8_debugf ("summand, direction: %i\n", idirection);
//          t8dg_element_dof_values_debug_print (summand);

    t8dg_element_dof_values_axpy (1, summand, element_res_dof_values);
//          t8_debugf ("element_res_dof_values, direction: %i\n", idirection);
//          t8dg_element_dof_values_debug_print (element_res_dof_values);
  }

  t8dg_element_dof_values_copy (element_res_dof_values, dest_element_dof);

  t8dg_element_dof_values_destroy (&element_res_dof_values);
  t8dg_element_dof_values_destroy (&element_derivative_dof_values);
  t8dg_element_dof_values_destroy (&summand);
  t8dg_element_quad_values_destroy (&element_quad_values);
  t8dg_element_quad_values_destroy (&element_derivative_quad_values);
}

void
t8dg_local_values_apply_face_inverse_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement, int iface,
                                                  t8dg_face_dof_values_t * src_face_dof, t8dg_face_dof_values_t * dest_face_dof)
{
  int                 eclass = t8dg_forest_get_eclass (values->forest, itree, ielement);
  if (t8dg_global_values_simplifies (values->global_values[eclass])) {
    t8dg_local_values_face_divide_trafo_quad_weight (values, itree, ielement, iface, (t8dg_face_quad_values_t *) src_face_dof,
                                                     (t8dg_face_quad_values_t *) dest_face_dof);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_local_values_apply_face_inverse_mass_matrix_idata_eclass (t8dg_local_values_t * values, t8_locidx_t idata, t8_eclass_t element_eclass,
                                                               int iface, t8dg_face_dof_values_t * src_face_dof,
                                                               t8dg_face_dof_values_t * dest_face_dof)
{
  if (t8dg_global_values_simplifies (values->global_values[element_eclass])) {
    t8dg_local_values_face_divide_trafo_quad_weight_idata_eclass (values, idata, element_eclass, iface,
                                                                  (t8dg_face_quad_values_t *) src_face_dof,
                                                                  (t8dg_face_quad_values_t *) dest_face_dof);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_local_values_apply_face_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement, int iface,
                                          t8dg_face_dof_values_t * src_face_dof, t8dg_face_dof_values_t * dest_face_dof)
{
  int                 eclass = t8dg_forest_get_eclass (values->forest, itree, ielement);
  if (t8dg_global_values_simplifies (values->global_values[eclass])) {
    t8dg_local_values_face_multiply_trafo_quad_weight (values, itree, ielement, iface, (t8dg_face_quad_values_t *) src_face_dof,
                                                       (t8dg_face_quad_values_t *) dest_face_dof);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_local_values_apply_face_mass_matrix_idata_eclass (t8dg_local_values_t * values, t8_locidx_t idata, t8_eclass_t element_eclass,
                                                       int iface, t8dg_face_dof_values_t * src_face_dof,
                                                       t8dg_face_dof_values_t * dest_face_dof)
{
  if (t8dg_global_values_simplifies (values->global_values[element_eclass])) {
    t8dg_local_values_face_multiply_trafo_quad_weight_idata_eclass (values, idata, element_eclass, iface,
                                                                    (t8dg_face_quad_values_t *) src_face_dof,
                                                                    (t8dg_face_quad_values_t *) dest_face_dof);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_local_values_apply_element_inverse_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement,
                                                     t8dg_element_dof_values_t * src_element_dof,
                                                     t8dg_element_dof_values_t * dest_element_dof)
{
  int                 eclass = t8dg_forest_get_eclass (values->forest, itree, ielement);
  if (t8dg_global_values_simplifies (values->global_values[eclass])) {
    t8dg_local_values_element_divide_trafo_quad_weight (values, itree, ielement, (t8dg_element_quad_values_t *) src_element_dof,
                                                        (t8dg_element_quad_values_t *) dest_element_dof);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_local_values_apply_element_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement,
                                             t8dg_element_dof_values_t * src_element_dof, t8dg_element_dof_values_t * dest_element_dof)
{
  int                 eclass = t8dg_forest_get_eclass (values->forest, itree, ielement);
  if (t8dg_global_values_simplifies (values->global_values[eclass])) {
    t8dg_local_values_element_multiply_trafo_quad_weight (values, itree, ielement, (t8dg_element_quad_values_t *) src_element_dof,
                                                          (t8dg_element_quad_values_t *) dest_element_dof);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

t8dg_global_values_t *
t8dg_local_values_get_global_values (t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement)
{
  return t8dg_global_values_array_get_global_values (local_values->global_values, local_values->forest, itree, ielement);
}

int
t8dg_local_values_get_max_num_faces (const t8dg_local_values_t * local_values)
{
  return local_values->max_num_faces;
}

void
t8dg_local_values_transform_child_dof_to_parent_dof (t8dg_local_values_t * local_values_old,
                                                     t8dg_local_values_t * local_values_new,
                                                     t8dg_element_dof_values_t * child_dof[MAX_SUBELEMENTS],
                                                     t8dg_element_dof_values_t * parent_dof, int num_children,
                                                     t8_locidx_t itree, t8_locidx_t ielem_first_child, t8_locidx_t ielem_parent)
{
  t8dg_element_dof_values_t *summand;
  t8dg_element_dof_values_t *mass_times_child_dof;
  int                 ichild;

  summand = t8dg_element_dof_values_duplicate (parent_dof);
  mass_times_child_dof = t8dg_element_dof_values_duplicate (child_dof[0]);

  t8dg_element_dof_values_set_zero (parent_dof);

  t8dg_global_values_t *global_values;
  t8dg_functionbasis_t *functionbasis;

  global_values =
    t8dg_global_values_array_get_global_values (local_values_new->global_values, local_values_new->forest, itree, ielem_parent);
  functionbasis = t8dg_global_values_get_functionbasis (global_values);

  for (ichild = 0; ichild < num_children; ichild++) {
    t8dg_local_values_apply_element_mass_matrix (local_values_old, itree, ielem_first_child + ichild, child_dof[ichild],
                                                 mass_times_child_dof);

    t8dg_functionbasis_apply_child_interpolation_matrix_transpose (functionbasis, ichild, mass_times_child_dof, summand);

    t8dg_element_dof_values_axpy (1, summand, parent_dof);
  }

  t8dg_local_values_apply_element_inverse_mass_matrix (local_values_new, itree, ielem_parent, parent_dof, parent_dof);
  t8dg_element_dof_values_destroy (&mass_times_child_dof);
  t8dg_element_dof_values_destroy (&summand);
}

void
t8dg_local_values_transform_orient_face_child_dof_to_parent_dof_hanging_nodes (t8dg_local_values_t * local_values,
                                                                               t8dg_face_dof_values_t * child_face_dof[MAX_SUBFACES],
                                                                               t8dg_face_dof_values_t * parent_face_dof,
                                                                               const int num_face_children,
                                                                               t8_locidx_t idata_child_neighbour[MAX_SUBFACES],
                                                                               t8_eclass_t element_eclass_children[MAX_SUBFACES],
                                                                               int iface_child_neighbour[MAX_SUBFACES],
                                                                               t8_locidx_t idata_parent, t8_eclass_t element_eclass_parent,
                                                                               int iface_parent)
{
  t8dg_face_dof_values_t *summand;
  t8dg_face_dof_values_t *mass_times_child_dof;
  int                 ichild;

  t8dg_global_values_t *global_values;
  t8dg_functionbasis_t *face_functionbasis;

  summand = t8dg_face_dof_values_duplicate (parent_face_dof);
  mass_times_child_dof = t8dg_face_dof_values_duplicate (child_face_dof[0]);

  t8dg_face_dof_values_set_zero (parent_face_dof);

  global_values = local_values->global_values[element_eclass_parent];
  face_functionbasis = t8dg_functionbasis_get_face_functionbasis (t8dg_global_values_get_functionbasis (global_values), iface_parent);

  for (ichild = 0; ichild < num_face_children; ichild++) {
    t8dg_local_values_apply_face_mass_matrix_idata_eclass (local_values, idata_child_neighbour[ichild], element_eclass_children[ichild],
                                                           iface_child_neighbour[ichild], child_face_dof[ichild], mass_times_child_dof);

    /*TODO: orient */

    t8dg_functionbasis_apply_child_interpolation_matrix_transpose (face_functionbasis, ichild, mass_times_child_dof, summand);

    t8dg_face_dof_values_axpy (1, summand, parent_face_dof);
  }
  t8dg_local_values_apply_face_inverse_mass_matrix_idata_eclass (local_values, idata_parent, element_eclass_parent, iface_parent,
                                                                 parent_face_dof, parent_face_dof);
  t8dg_face_dof_values_destroy (&mass_times_child_dof);
  t8dg_face_dof_values_destroy (&summand);
}

double
t8dg_local_values_element_norm_l2_squared (t8dg_local_values_t * local_values, t8dg_element_dof_values_t * element_dof_values,
                                           t8_locidx_t itree, t8_locidx_t ielement)
{
  double              norm = 0;
  int                 idof, num_dof;
  t8dg_element_dof_values_t *mass_times_square_dof;
  t8dg_element_dof_values_t *element_dof_square_values;

  num_dof = t8dg_element_dof_values_get_num_dof (element_dof_values);

  element_dof_square_values = t8dg_element_dof_values_duplicate (element_dof_values);
  mass_times_square_dof = t8dg_element_dof_values_duplicate (element_dof_values);

  t8dg_element_dof_values_square_values (element_dof_values, element_dof_square_values);
  t8dg_local_values_apply_element_mass_matrix (local_values, itree, ielement, element_dof_square_values, mass_times_square_dof);

  for (idof = 0; idof < num_dof; idof++) {
    norm += t8dg_element_dof_values_get_value (mass_times_square_dof, idof);
  }
  t8dg_element_dof_values_destroy (&element_dof_square_values);
  t8dg_element_dof_values_destroy (&mass_times_square_dof);
  return norm;
}

void
t8dg_local_values_element_error_ana_l2_squared (t8dg_local_values_t * local_values, t8dg_dof_values_t * dof_values,
                                                t8dg_dof_values_t * analytical_sol_dof, t8_locidx_t itree, t8_locidx_t ielement,
                                                double time, double *error_squared_summand, double *ana_norm_squared_summand)
{
  t8dg_element_dof_values_t *elem_dof_val;
  t8dg_element_dof_values_t *elem_ana_sol;
  t8dg_element_dof_values_t *elem_error;

  elem_dof_val = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielement);
  elem_ana_sol = t8dg_dof_values_new_element_dof_values_view (analytical_sol_dof, itree, ielement);
  elem_error = t8dg_element_dof_values_duplicate (elem_dof_val);

  t8dg_element_dof_values_axpyz (-1, elem_ana_sol, elem_dof_val, elem_error);
  *error_squared_summand = t8dg_local_values_element_norm_l2_squared (local_values, elem_error, itree, ielement);
  *ana_norm_squared_summand = t8dg_local_values_element_norm_l2_squared (local_values, elem_ana_sol, itree, ielement);

  t8dg_element_dof_values_destroy (&elem_ana_sol);
  t8dg_element_dof_values_destroy (&elem_dof_val);
  t8dg_element_dof_values_destroy (&elem_error);
}

t8dg_coarse_geometry_t *
t8dg_local_values_get_coarse_geometry (t8dg_local_values_t * local_values)
{
  return local_values->coarse_geometry;
}
