/*
 * t8dg_local_precomputed_values.c
 *
 *  Created on: Apr 5, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_local_precomputed_values.h"
#include "t8dg_sc_array.h"
#include "t8dg_geometry.h"
#include "t8dg_flux.h"

#include <t8_vec.h>
#include <t8_forest/t8_forest_partition.h>

struct t8dg_local_precomputed_values
{
  int                 num_faces, num_elem_quad, num_face_quad[MAX_FACES], dim;

  sc_array_t         *element_trafo_quad_weight;                /**< For each element, an array of length num_quad_points with the combined weight
								of quadrature and transformation */
  sc_array_t         *element_transformed_gradient_tangential_vectors;  /**< For each element, (d_x F^-1)(tau_d), needed to calculate the gradient */

  /* To avoid another dimension when flattening the array, and since the number of Faces is bound,
   *  have an sc_array for each faceindex of length
   * of the processorlocal elements*/
  sc_array_t         *face_trafo_quad_weight[MAX_FACES];        /**< for each element, an array of length num_quad_points with the combined weight
								of quadrature and transformation */
  sc_array_t         *face_normal_vectors[MAX_FACES];           /**< For each face and quadrature point the 3-dim normal vector in image space*/

};

static double      *
t8dg_local_precomputed_values_get_face_quad_trafo_weights (const t8dg_local_precomputed_values_t * values, const t8_locidx_t idata,
                                                           const int faceindex)
{
  return ((double *) t8_sc_array_index_locidx (values->face_trafo_quad_weight[faceindex], idata));
}

static double      *
t8dg_local_precomputed_values_get_element_quad_trafo_weights (const t8dg_local_precomputed_values_t * values, const t8_locidx_t idata)
{
  return ((double *) t8_sc_array_index_locidx (values->element_trafo_quad_weight, idata));
}

double             *
t8dg_local_precomputed_values_get_transformed_gradient_tangential_vector (const t8dg_local_precomputed_values_t * values,
                                                                          const t8_locidx_t idata, const t8dg_quad_idx_t iquad,
                                                                          const int idim)
{
  T8_ASSERT (idim >= 0 && idim < values->dim);
  T8_ASSERT (iquad >= 0 && (size_t) iquad < values->element_transformed_gradient_tangential_vectors->elem_size / (DIM3 * sizeof (double)));
  return ((double *) t8_sc_array_index_locidx (values->element_transformed_gradient_tangential_vectors, idata)) +
    DIM3 * (values->dim * iquad + idim);
}

double             *
t8dg_local_precomputed_values_get_face_normal_vector (const t8dg_local_precomputed_values_t * values, const t8_locidx_t idata,
                                                      const int iface, const t8dg_quad_idx_t iquad)
{
  T8_ASSERT (iface >= 0 && iface < MAX_FACES);
  T8_ASSERT (iquad >= 0 && (size_t) iquad < values->face_normal_vectors[iface]->elem_size / (DIM3 * sizeof (double)));
  return ((double *) t8_sc_array_index_locidx (values->face_normal_vectors[iface], idata)) + DIM3 * iquad;
}

void
t8dg_local_precomputed_values_copy_element_values (t8dg_local_precomputed_values_t * incoming_values, t8_locidx_t incoming_idata,
                                                   t8dg_local_precomputed_values_t * outgoing_values, t8_locidx_t outgoing_idata)
{
  int                 iface;
  t8dg_sc_array_copy_only_at_indices (incoming_values->element_trafo_quad_weight, incoming_idata,
                                      outgoing_values->element_trafo_quad_weight, outgoing_idata);

  t8dg_sc_array_copy_only_at_indices (incoming_values->element_transformed_gradient_tangential_vectors, incoming_idata,
                                      outgoing_values->element_transformed_gradient_tangential_vectors, outgoing_idata);

  for (iface = 0; iface < incoming_values->num_faces; iface++) {
    t8dg_sc_array_copy_only_at_indices (incoming_values->face_trafo_quad_weight[iface], incoming_idata,
                                        outgoing_values->face_trafo_quad_weight[iface], outgoing_idata);

    t8dg_sc_array_copy_only_at_indices (incoming_values->face_normal_vectors[iface], incoming_idata,
                                        outgoing_values->face_normal_vectors[iface], outgoing_idata);
  }
}

void
t8dg_local_precomputed_values_set_element (t8dg_local_precomputed_values_t * values,
                                           const t8dg_geometry_transformation_data_t * geometry_data, const t8dg_quadrature_t * quadrature)
{
  double              sqrt_gram_det, face_sqrt_gram_det;
  int                 iquad, iface, idim;
  t8_locidx_t         idata;

  /*pointer on the values to fill */
  double             *element_quad_trafo;       /*size: num_elem_quad */
  double             *face_quad_trafo;  /*size: num_face_quad */
  double             *transformed_gradient_tangential_vector;   /*size: 3, gets evaluated for each element quad point */
  double             *image_normal_vector;      /*size: 3, gets evaluated for each element quad point */

  double              reference_vertex[3];
  double              reference_tangential_vector[3] = { 1, 0, 0 };

  idata = t8dg_itree_ielement_to_idata (geometry_data->forest, geometry_data->itree, geometry_data->ielement);

  element_quad_trafo = t8dg_local_precomputed_values_get_element_quad_trafo_weights (values, idata);
  for (iquad = 0; iquad < values->num_elem_quad; iquad++) {
    t8dg_quadrature_get_element_vertex (quadrature, iquad, reference_vertex);
    sqrt_gram_det = t8dg_geometry_calculate_sqrt_gram_determinant (geometry_data, reference_vertex);
    element_quad_trafo[iquad] = sqrt_gram_det * t8dg_quadrature_get_element_weight (quadrature, iquad);
    for (idim = 0; idim < t8dg_quadrature_get_dim (quadrature); idim++) {
      transformed_gradient_tangential_vector =
        t8dg_local_precomputed_values_get_transformed_gradient_tangential_vector (values, idata, iquad, idim);

      reference_tangential_vector[0] = reference_tangential_vector[1] = reference_tangential_vector[2] = 0;
      reference_tangential_vector[idim] = 1;

      t8dg_geometry_calculate_transformed_gradient_tangential_vector (geometry_data, reference_vertex, reference_tangential_vector,
                                                                      transformed_gradient_tangential_vector);
    }
  }

  for (iface = 0; iface < values->num_faces; iface++) {
    face_quad_trafo = t8dg_local_precomputed_values_get_face_quad_trafo_weights (values, idata, iface);
    for (iquad = 0; iquad < values->num_face_quad[iface]; iquad++) {
      /*for 1D elements the faceintegrals are just the value at the facequadrature point */
      t8dg_quadrature_get_face_vertex (quadrature, iface, iquad, reference_vertex);
      face_sqrt_gram_det = t8dg_geometry_calculate_face_sqrt_gram_determinant (geometry_data, iface, reference_vertex);
      face_quad_trafo[iquad] = face_sqrt_gram_det * t8dg_quadrature_get_face_weight (quadrature, iface, iquad);
      image_normal_vector = t8dg_local_precomputed_values_get_face_normal_vector (values, idata, iface, iquad);
      t8dg_geometry_calculate_normal_vector (geometry_data, iface, reference_vertex, image_normal_vector);
    }
  }
}

t8dg_local_precomputed_values_t *
t8dg_local_precomputed_values_new (const t8dg_quadrature_t * quadrature, const t8_locidx_t num_local_elems)
{
  int                 iface;
  t8dg_local_precomputed_values_t *values = T8DG_ALLOC (t8dg_local_precomputed_values_t, 1);

  values->num_faces = t8dg_quadrature_get_num_faces (quadrature);
  values->num_elem_quad = t8dg_quadrature_get_num_element_vertices (quadrature);
  values->dim = t8dg_quadrature_get_dim (quadrature);

  /*for each element an array of double values */
  values->element_trafo_quad_weight =
    sc_array_new_count (sizeof (double) * t8dg_quadrature_get_num_element_vertices (quadrature), num_local_elems);

  for (iface = 0; iface < values->num_faces; iface++) {
    values->num_face_quad[iface] = t8dg_quadrature_get_num_face_vertices (quadrature, iface);

    values->face_trafo_quad_weight[iface] = sc_array_new_count (sizeof (double) * values->num_face_quad[iface], num_local_elems);

    values->face_normal_vectors[iface] = sc_array_new_count (sizeof (double) * DIM3 * values->num_face_quad[iface], num_local_elems);

  }

  values->element_transformed_gradient_tangential_vectors =
    sc_array_new_count (sizeof (double) * DIM3 * values->dim * values->num_elem_quad, num_local_elems);

  return values;
}

void
t8dg_local_precomputed_values_destroy (t8dg_local_precomputed_values_t ** pvalues)
{
  t8dg_local_precomputed_values_t *values = *pvalues;
  int                 iface;
  sc_array_destroy (values->element_trafo_quad_weight);
  sc_array_destroy (values->element_transformed_gradient_tangential_vectors);
  for (iface = 0; iface < values->num_faces; iface++) {
    sc_array_destroy (values->face_trafo_quad_weight[iface]);
    sc_array_destroy (values->face_normal_vectors[iface]);
  }
  T8DG_FREE (values);
  *pvalues = NULL;
}

void
t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (const t8dg_local_precomputed_values_t * local_values,
                                                                  sc_array_t * element_quad, t8_locidx_t idata)
{
  t8dg_quad_idx_t     iquad;
  sc_array_t         *element_trafo_quad_weights;
  double              quad_trafo_weight;
  element_trafo_quad_weights = t8dg_sc_array_block_double_new_view (local_values->element_trafo_quad_weight, idata);

  for (iquad = 0; iquad < local_values->num_elem_quad; iquad++) {
    quad_trafo_weight = *(double *) t8dg_sc_array_index_quadidx (element_trafo_quad_weights, iquad);
    *(double *) t8dg_sc_array_index_quadidx (element_quad, iquad) *= quad_trafo_weight;
  }
  sc_array_destroy (element_trafo_quad_weights);
}

void
t8dg_local_precomputed_values_element_divide_trafo_quad_weight (const t8dg_local_precomputed_values_t * local_values,
                                                                sc_array_t * element_quad, t8_locidx_t idata)
{
  t8dg_quad_idx_t     iquad;
  sc_array_t         *element_trafo_quad_weights;
  double              quad_trafo_weight;
  element_trafo_quad_weights = t8dg_sc_array_block_double_new_view (local_values->element_trafo_quad_weight, idata);

  for (iquad = 0; iquad < local_values->num_elem_quad; iquad++) {
    quad_trafo_weight = *(double *) t8dg_sc_array_index_quadidx (element_trafo_quad_weights, iquad);
    *(double *) t8dg_sc_array_index_quadidx (element_quad, iquad) /= quad_trafo_weight;
  }
  sc_array_destroy (element_trafo_quad_weights);
}

void
t8dg_local_precomputed_values_partition (t8_forest_t forest_old, t8_forest_t forest_partition,
                                         t8dg_local_precomputed_values_t * local_values_old,
                                         t8dg_local_precomputed_values_t * local_values_partition)
{
  int                 iface;
  t8_forest_partition_data (forest_old, forest_partition,
                            local_values_old->element_trafo_quad_weight, local_values_partition->element_trafo_quad_weight);

  t8_forest_partition_data (forest_old, forest_partition,
                            local_values_old->element_transformed_gradient_tangential_vectors,
                            local_values_partition->element_transformed_gradient_tangential_vectors);

  for (iface = 0; iface < local_values_old->num_faces; iface++) {
    t8_forest_partition_data (forest_old, forest_partition,
                              local_values_old->face_normal_vectors[iface], local_values_partition->face_normal_vectors[iface]);

    t8_forest_partition_data (forest_old, forest_partition,
                              local_values_old->face_trafo_quad_weight[iface], local_values_partition->face_trafo_quad_weight[iface]);
  }
}

void
t8dg_local_precomputed_values_element_multiply_flux_value (const t8dg_local_precomputed_values_t * local_values, const t8dg_flux_t * flux,
                                                           const t8dg_geometry_transformation_data_t * geometry_data,
                                                           t8dg_quadrature_t * quadrature, double current_time,
                                                           sc_array_t * element_quad_values)
{
  t8dg_quad_idx_t     iquad, num_quad_vertices;
  double             *transformed_gradient_tangential_vector;
  double              reference_vertex[3];
  double              image_vertex[3];
  double              flux_vec[3];
  double              flux_value;
  t8_locidx_t         idata;
  int                 idim = 0; /*TODO: loop over dim */
  num_quad_vertices = t8dg_quadrature_get_num_element_vertices (quadrature);

  idata = t8dg_itree_ielement_to_idata (geometry_data->forest, geometry_data->itree, geometry_data->ielement);
  for (iquad = 0; iquad < num_quad_vertices; iquad++) {
    transformed_gradient_tangential_vector =
      t8dg_local_precomputed_values_get_transformed_gradient_tangential_vector (local_values, idata, iquad, idim);

    t8dg_quadrature_get_element_vertex (quadrature, iquad, reference_vertex);

    t8dg_geometry_transform_reference_vertex_to_image_vertex (geometry_data, reference_vertex, image_vertex);

    t8dg_flux_calulate_flux (flux, image_vertex, flux_vec, current_time);
    flux_value = t8_vec_dot (flux_vec, transformed_gradient_tangential_vector);

    *(double *) t8dg_sc_array_index_quadidx (element_quad_values, iquad) *= flux_value;
  }
}
