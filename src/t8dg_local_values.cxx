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

typedef struct t8dg_local_values_fn_evaluation_data
{
  t8dg_coarse_geometry_t *coarse_geometry;
  t8_forest_t         forest;
  t8dg_scalar_function_3d_time_fn function;
  double              time;
} t8dg_local_values_fn_evaluation_data_t;

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

/*TODO: change getter functions to t8dg_dof/quad equivalent*/

static t8dg_face_dof_values_t *
t8dg_local_values_get_face_trafo_quad_weights_view (t8dg_local_values_t * values, const t8_locidx_t itree, const t8_locidx_t ielement,
                                                    const int faceindex)
{
  return t8dg_quad_values_new_face_quad_values_view (values->face_trafo_quad_weight[faceindex], faceindex, itree, ielement);
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
  T8DG_ABORT ("Not implemented \n ");
  /*
     T8_ASSERT (idim >= 0 && idim < values->dim);
     T8_ASSERT (iquad >= 0 && (size_t) iquad < values->element_transformed_gradient_tangential_vectors->elem_size / (DIM3 * sizeof (double)));
     return ((double *) t8_sc_array_index_locidx (values->element_transformed_gradient_tangential_vectors, idata)) +
     DIM3 * (values->dim * iquad + idim);
   */
}

void
t8dg_local_values_set_face_normal_vector (const t8dg_local_values_t * values,
                                          const t8_locidx_t itree, const t8_locidx_t ielement, const int iface, const int idof,
                                          double vector[3])
{
  /*
     T8_ASSERT (iface >= 0 && iface < MAX_FACES);
     T8_ASSERT (iquad >= 0 && (size_t) iquad < values->face_normal_vectors[iface]->elem_size / (DIM3 * sizeof (double)));
     return ((double *) t8_sc_array_index_locidx (values->face_normal_vectors[iface], idata)) + DIM3 * iquad;
   */
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

  global_values = t8dg_global_values_array_get_global_values (local_values->global_values, local_values->forest, itree, ielement);
  quadrature = t8dg_global_values_get_quadrature (global_values);
  functionbasis = t8dg_global_values_get_functionbasis (global_values);

  num_elem_quad = t8dg_global_values_get_num_elem_quad (global_values);
  num_faces = t8dg_global_values_get_num_faces (global_values);

  element_trafo_quad = t8dg_local_values_get_element_trafo_quad_weights_view (local_values, itree, ielement);

  for (iquad = 0; iquad < num_elem_quad; iquad++) {
    t8dg_quadrature_get_element_vertex (quadrature, iquad, reference_vertex);
    sqrt_gram_det =
      t8dg_geometry_calculate_sqrt_gram_determinant (local_values->coarse_geometry, local_values->forest, itree, ielement,
                                                     reference_vertex);

    t8dg_element_quad_values_set_value (element_trafo_quad, iquad, sqrt_gram_det * t8dg_quadrature_get_element_weight (quadrature, iquad));

    for (idim = 0; idim < local_values->dim; idim++) {
      reference_tangential_vector[0] = reference_tangential_vector[1] = reference_tangential_vector[2] = 0;
      reference_tangential_vector[idim] = 1;

      t8dg_geometry_calculate_transformed_gradient_tangential_vector (local_values->coarse_geometry, local_values->forest, itree, ielement,
                                                                      reference_vertex, reference_tangential_vector,
                                                                      transformed_gradient_tangential_vector);
      t8dg_local_values_set_transformed_gradient_tangential_vector (local_values, itree, ielement, iquad, idim,
                                                                    transformed_gradient_tangential_vector);
    }
  }

  for (iface = 0; iface < num_faces; iface++) {
    face_trafo_quad = t8dg_local_values_get_face_trafo_quad_weights_view (local_values, itree, ielement, iface);
    num_face_quad = t8dg_quadrature_get_num_face_vertices (quadrature, iface);
    for (iquad = 0; iquad < num_face_quad; iquad++) {
      t8dg_quadrature_get_face_vertex (quadrature, iface, iquad, reference_vertex);
      face_sqrt_gram_det =
        t8dg_geometry_calculate_face_sqrt_gram_determinant (local_values->coarse_geometry, local_values->forest, itree, ielement, iface,
                                                            reference_vertex);
      t8dg_face_quad_values_set_value (face_trafo_quad, iquad,
                                       face_sqrt_gram_det * t8dg_quadrature_get_face_weight (quadrature, iface, iquad));
    }

    num_face_dof = t8dg_functionbasis_get_num_face_dof (functionbasis, iface);
    for (idof = 0; idof < num_face_dof; idof++) {
      t8dg_functionbasis_get_lagrange_face_vertex (functionbasis, iface, idof, reference_vertex);
      t8dg_geometry_calculate_normal_vector (local_values->coarse_geometry, local_values->forest, itree, ielement, iface, reference_vertex,
                                             image_normal_vector);
      t8dg_local_values_set_face_normal_vector (local_values, itree, ielement, iface, idof, image_normal_vector);
    }
  }
}

t8dg_local_values_t *
t8dg_local_values_new (t8_forest_t forest, t8dg_global_values_t ** global_values_array, t8dg_coarse_geometry_t * coarse_geometry)
{
  t8dg_local_values_t *local_values = T8DG_ALLOC_ZERO (t8dg_local_values_t, 1);

  t8_locidx_t         num_local_elems;

  int                 iface, idim, icomp;
  int                 eclass;

  local_values->forest = forest;

  num_local_elems = t8_forest_get_num_element (forest);

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
    else {
      t8dg_debugf ("global_values has no entry for eclass = %i\n", eclass);
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
  double              quad_trafo_weight;
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
t8dg_local_values_partition (t8dg_local_values_t * local_values_old, t8dg_local_values_t * local_values_partition)
{
  int                 iface, idim, icomp;
  /*TODO!!!! */
  sc_array_t         *local_view_normal_old;
  sc_array_t         *local_view_normal_new;
  sc_array_t         *local_view_face_weight_old;
  sc_array_t         *local_view_face_weight_new;

  t8dg_quad_values_partition (local_values_old->element_trafo_quad_weight, local_values_partition->element_trafo_quad_weight);

  for (idim = 0; idim < local_values_old->dim; idim++) {
    for (icomp = 0; icomp < DIM3; icomp++) {
      t8dg_quad_values_partition (local_values_old->element_transformed_gradient_tangential_vectors[idim][icomp],
                                  local_values_partition->element_transformed_gradient_tangential_vectors[idim][icomp]);

    }
  }

  for (iface = 0; iface < local_values_old->max_num_faces; iface++) {
    t8dg_quad_values_partition (local_values_old->face_trafo_quad_weight[iface], local_values_partition->face_trafo_quad_weight[iface]);

    for (idim = 0; icomp < DIM3; icomp++) {
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
    t8_forest_ghost_exchange_data (local_values->forest, (sc_array_t *) local_values->face_trafo_quad_weight[iface]);
    for (idim = 0; idim < DIM3; idim++) {
      t8_forest_ghost_exchange_data (local_values->forest, (sc_array_t *) local_values->face_normal_vectors[iface][idim]);
    }
  }
}

#if 0
void
t8dg_local_values_element_multiply_flux_value (const t8dg_local_values_t * local_values, const t8dg_flux_t * flux,
                                               const t8dg_geometry_transformation_data_t * geometry_data,
                                               t8dg_quadrature_t * quadrature, double current_time,
                                               int idim, sc_array_t * element_quad_values, sc_array_t * element_flux_quad_values)
{
  int                 iquad, num_quad_vertices;
  double             *transformed_gradient_tangential_vector;
  double              reference_vertex[3];
  double              image_vertex[3];
  double              flux_vec[3];
  double              flux_value;
  t8_locidx_t         idata;
  num_quad_vertices = t8dg_quadrature_get_num_element_vertices (quadrature);

  idata = t8dg_itree_ielement_to_idata (geometry_data->forest, geometry_data->itree, geometry_data->ielement);
  for (iquad = 0; iquad < num_quad_vertices; iquad++) {
    t8dg_quadrature_get_element_vertex (quadrature, iquad, reference_vertex);

    t8dg_geometry_transform_reference_vertex_to_image_vertex (geometry_data, reference_vertex, image_vertex);

    t8dg_flux_calulate_flux (flux, image_vertex, flux_vec, current_time);

    transformed_gradient_tangential_vector =
      t8dg_local_values_get_transformed_gradient_tangential_vector (local_values, idata, iquad, idim);

    flux_value = t8_vec_dot (flux_vec, transformed_gradient_tangential_vector);

    *(double *) sc_array_index_int (element_flux_quad_values, iquad) = flux_value *
      *(double *) sc_array_index_int (element_quad_values, iquad);
  }
}

#endif

void
t8dg_local_values_apply_face_inverse_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement, int iface,
                                                  t8dg_face_dof_values_t * src_face_dof, t8dg_face_dof_values_t * dest_face_dof)
{
  int                 eclass = t8_forest_get_eclass (values->forest, itree);
  if (t8dg_global_values_simplifies (values->global_values[eclass])) {
    t8dg_local_values_face_divide_trafo_quad_weight (values, itree, ielement, iface, (t8dg_face_quad_values_t *) src_face_dof,
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
  int                 eclass = t8_forest_get_eclass (values->forest, itree);
  if (t8dg_global_values_simplifies (values->global_values[eclass])) {
    t8dg_local_values_face_multiply_trafo_quad_weight (values, itree, ielement, iface, (t8dg_face_quad_values_t *) src_face_dof,
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
  int                 eclass = t8_forest_get_eclass (values->forest, itree);
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
  int                 eclass = t8_forest_get_eclass (values->forest, itree);
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

double             *
t8dg_local_values_get_face_normal_vector (const t8dg_local_values_t * values, const t8_locidx_t idata, const int iface, const int iquad)
{
  T8DG_ABORT ("Not yet implemented");

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

    t8dg_functionbasis_apply_child_interpolation_matrix_transpose (t8dg_global_values_get_functionbasis (global_values), ichild,
                                                                   mass_times_child_dof, summand);

    t8dg_element_dof_values_axpy (1, summand, parent_dof);
  }

  t8dg_local_values_apply_element_inverse_mass_matrix (local_values_new, itree, ielem_parent, parent_dof, parent_dof);
  t8dg_element_dof_values_destroy (&mass_times_child_dof);
  t8dg_element_dof_values_destroy (&summand);
}

void
t8dg_local_values_transform_face_child_dof_to_parent_dof (t8dg_local_values_t * local_values_old,
                                                          t8dg_local_values_t * local_values_new,
                                                          t8dg_face_dof_values_t * child_face_dof[MAX_SUBFACES],
                                                          t8dg_face_dof_values_t * parent_face_dof, const int num_face_children,
                                                          t8_locidx_t itree,
                                                          t8_locidx_t ielem_child_old[MAX_SUBFACES], t8_locidx_t ielem_parent_new,
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

  global_values =
    t8dg_global_values_array_get_global_values (local_values_new->global_values, local_values_new->forest, itree, ielem_parent_new);
  face_functionbasis = t8dg_functionbasis_get_face_functionbasis (t8dg_global_values_get_functionbasis (global_values), iface_parent);
  for (ichild = 0; ichild < num_face_children; ichild++) {
    t8dg_local_values_apply_face_mass_matrix (local_values_old, itree, ielem_child_old[ichild], iface_parent, child_face_dof[ichild],
                                              mass_times_child_dof);

    t8dg_functionbasis_apply_child_interpolation_matrix_transpose (face_functionbasis, ichild, mass_times_child_dof, summand);

    t8dg_face_dof_values_axpy (1, summand, parent_face_dof);
  }
  t8dg_local_values_apply_face_inverse_mass_matrix (local_values_new, itree, ielem_parent_new, iface_parent, parent_face_dof,
                                                    parent_face_dof);
  t8dg_face_dof_values_destroy (&mass_times_child_dof);
  t8dg_face_dof_values_destroy (&summand);
}

#if 0

double
t8dg_local_values_element_norm_infty (sc_array_t * element_dof_values)
{
  double              norm = 0;
  size_t              idof;
  for (idof = 0; idof < element_dof_values->elem_count; idof++) {
    norm = SC_MAX (norm, fabs (*(double *) sc_array_index (element_dof_values, idof)));
  }
  return norm;
}

double
t8dg_local_values_element_norm_l2_squared (sc_array_t * element_dof_values, t8dg_global_values_t * global_values,
                                           t8dg_local_values_t * local_values, t8_locidx_t idata)
{
  T8DG_ASSERT (element_dof_values->elem_count == (size_t) t8dg_global_values_get_num_dof (global_values));
  double              norm = 0;
  int                 idof, num_dof;
  sc_array_t         *mass_times_square_dof;
  sc_array_t         *element_dof_square_values;

  num_dof = t8dg_global_values_get_num_dof (global_values);

  element_dof_square_values = t8dg_sc_array_duplicate (element_dof_values);
  mass_times_square_dof = t8dg_sc_array_duplicate (element_dof_values);

  t8dg_sc_array_block_square_values (element_dof_values, element_dof_square_values);
  t8dg_local_values_apply_element_mass_matrix (global_values, local_values, idata, element_dof_square_values, mass_times_square_dof);

  for (idof = 0; idof < num_dof; idof++) {
    norm += *(double *) sc_array_index_int (mass_times_square_dof, idof);
  }
  sc_array_destroy (element_dof_square_values);
  sc_array_destroy (mass_times_square_dof);
  return norm;
}

/*Testing purpose*/
double
t8dg_local_values_transform_reference_vertex_and_evaluate (const double reference_vertex[3], void *scalar_fn_data)
{
  double              image_vertex[DIM3];
  t8dg_local_values_fn_evaluation_data_t *data;

  data = (t8dg_values_fn_evaluation_data_t *) scalar_fn_data;

  t8dg_geometry_transform_reference_vertex_to_image_vertex (data->geometry_data, reference_vertex, image_vertex);

  /* apply initial condition function at image vertex and start time */
  return data->function (image_vertex, data->time);

}

#endif
