/*
 * t8dg_mortar.c
 *
 *  Created on: Apr 3, 2020
 *      Author: lukas
 */

#include <t8_forest.h>
#include <sc_containers.h>

#include "t8dg_mortar.h"
#include "t8dg.h"
#include "t8dg_timestepping.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_flux.h"
#include "t8dg_sc_array.h"

/** struct used to save the calculated numerical fluxes at quadrature points*/
struct t8dg_mortar
{
  int                 number_face_quadrature_points;      /**< The number of face quadrature points*/
  t8_locidx_t         elem_idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_idata_plus;                    /**< Local index of the element corresponding to u_plus */

  int                 iface_plus, iface_minus;

  /*one value for each quadrature point */
  sc_array_t         *u_minus;                          /**< value of u on elem_minus at face quadrature points */
  sc_array_t         *u_plus;                           /**< value of u on elem_plus at face quadrature points */

  sc_array_t         *fluxes;                           /**< value of (cu)*.n at face quadrature points */
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/

};                              /*maybe change to opaque handle */

sc_array_t         *
t8dg_mortar_get_flux (t8dg_mortar_t * mortar)
{
  return mortar->fluxes;
}

void
t8dg_mortar_get_idata_iface (t8dg_mortar_t * mortar, t8_locidx_t * pidata, int *piface, int side)
{
  T8DG_ASSERT (side == 0 || side == 1);
  T8DG_ASSERT (pidata != NULL && piface != NULL);
  T8DG_ASSERT (mortar != NULL);
  if (side == 0) {
    *piface = mortar->iface_minus;
    *pidata = mortar->elem_idata_minus;
  }
  else if (side == 1) {
    *piface = mortar->iface_plus;
    *pidata = mortar->elem_idata_plus;
  }
}

int
t8dg_mortar_is_valid (t8dg_mortar_t * mortar)
{
  return mortar->valid;
}

int
t8dg_mortar_get_side (t8dg_mortar_t * mortar, t8_locidx_t idata)
{
  if (idata == mortar->elem_idata_minus) {
    return -1;
  }
  else if (idata == mortar->elem_idata_plus) {
    return 1;
  }
  return 0;
}

void
t8dg_mortar_destroy (t8dg_mortar_t ** pmortar)
{
  sc_array_destroy ((*pmortar)->fluxes);
  sc_array_destroy ((*pmortar)->u_minus);
  sc_array_destroy ((*pmortar)->u_plus);
  T8DG_FREE (*pmortar);
  *pmortar = NULL;
}

void
t8dg_mortar_invalidate (t8dg_mortar_t * mortar)
{
  mortar->valid = 0;
}

t8dg_mortar_t      *
t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface)
{
  t8dg_mortar_t      *mortar = T8_ALLOC (t8dg_mortar_t, 1);
  t8_locidx_t         idata;
  int                *neigh_ifaces;
  int                 num_neighs;

  t8_element_t       *element, **neigh_elems;
  t8_eclass_scheme_c *neigh_scheme;
  t8_locidx_t        *elem_indices;

  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  /* use this function to also get the data-indices of the neighbouring element */
  t8_forest_leaf_face_neighbors (forest, itree, element, &neigh_elems, iface, &neigh_ifaces, &num_neighs, &elem_indices, &neigh_scheme, 1);

  SC_CHECK_ABORT (num_neighs == 1, "only one faceneighbour currently implemented !");

  idata = t8dg_itree_ielement_to_idata (forest, itree, ielement);

  mortar->elem_idata_minus = idata;
  mortar->elem_idata_plus = elem_indices[0];    /*could be greater than number of local elements -> ghost */
  mortar->iface_minus = iface;
  mortar->iface_plus = neigh_ifaces[0]; /* get neighbouring face index */
  mortar->number_face_quadrature_points = 1;    /* Only viable for 1D case, getter function on problem? */

  /* allocate memory for sc_arrays */
  mortar->u_minus = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->u_plus = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->fluxes = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->valid = 0;

  T8_FREE (neigh_elems);
  T8_FREE (elem_indices);
  T8_FREE (neigh_ifaces);
  return mortar;
}

void
t8dg_mortar_fill (t8dg_mortar_t * mortar,
                  sc_array_t * element_dof_values,
                  t8dg_timestepping_data_t * time_data,
                  t8dg_global_precomputed_values_t * global_values,
                  t8dg_local_precomputed_values_t * local_values,
                  t8dg_linear_flux_t * flux, t8dg_linear_numerical_flux_fn numerical_flux_fn)
{
  sc_array_t         *elem_dof_values_minus;
  sc_array_t         *elem_dof_values_plus;

  t8dg_quad_idx_t     iquad;
  double              fluxvalue;
  double              u_minus_quad;
  double              u_plus_quad;
  double             *normal_vector;
  double              flux_vec[3];
  double              time = t8dg_timestepping_data_get_current_time (time_data);

  elem_dof_values_minus = t8dg_sc_array_block_double_new_view (element_dof_values, mortar->elem_idata_minus);
  elem_dof_values_plus = t8dg_sc_array_block_double_new_view (element_dof_values, mortar->elem_idata_plus);

  t8dg_global_precomputed_values_transform_element_dof_to_face_quad (global_values, mortar->iface_minus,
                                                                     elem_dof_values_minus, mortar->u_minus);
  t8dg_global_precomputed_values_transform_element_dof_to_face_quad (global_values, mortar->iface_plus,
                                                                     elem_dof_values_plus, mortar->u_plus);

  /*normal vector from element_minus to element_plus */
  for (iquad = 0; iquad < mortar->number_face_quadrature_points; iquad++) {
    normal_vector =
      t8dg_local_precomputed_values_get_face_normal_vector (local_values, mortar->elem_idata_minus, mortar->iface_minus, iquad);

    /*TODO: t8dg_geometry, current Flux is independent of x and time */
    double              x_vec[3] = { 0, 0, 0 };
    t8dg_linear_flux_calulate_flux (flux, x_vec, flux_vec, time);

    /*TODO: Orientation!! */
    u_minus_quad = *(double *) t8dg_sc_array_index_quadidx (mortar->u_minus, iquad);
    u_plus_quad = *(double *) t8dg_sc_array_index_quadidx (mortar->u_plus, iquad);

    fluxvalue = numerical_flux_fn (u_minus_quad, u_plus_quad, flux_vec, normal_vector);
    *(double *) t8dg_sc_array_index_quadidx (mortar->fluxes, iquad) = fluxvalue;
  }
  sc_array_destroy (elem_dof_values_minus);
  sc_array_destroy (elem_dof_values_plus);
  mortar->valid = 1;
}
