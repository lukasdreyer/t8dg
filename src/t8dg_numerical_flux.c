/*
 * t8dg_numerical_flux.c
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_numerical_flux.h"
#include <t8_vec.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest.h>

double
t8dg_upwind_flux_1D (const double u_minus, const double u_plus, const double flow_vector[3], const double normal_vector[3])
{
  /*TODO: ASSERT linear dependence of flow_vector and normal_vector */
  if (t8_vec_dot (flow_vector, normal_vector) > 0)
    return u_minus * t8_vec_dot (flow_vector, normal_vector);
  else
    return u_plus * t8_vec_dot (flow_vector, normal_vector);
}

sc_array_t         *
t8dg_mortar_get_flux (t8dg_mortar_t * mortar)
{
  return mortar->fluxes;
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

t8dg_mortar_t      *
t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface)
{
  t8dg_mortar_t      *mortar = T8_ALLOC (t8dg_mortar_t, 1);
  t8_locidx_t         idata;
  int                *neigh_ifaces;
  int                 num_neighs;
  t8_locidx_t         offset;

  t8_element_t       *element, **neigh_elems;
  t8_eclass_scheme_c *neigh_scheme;
  t8_locidx_t        *elem_indices;

  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  /* use this function to also get the data-indices of the neighbouring element */
  t8_forest_leaf_face_neighbors (forest, itree, element, &neigh_elems, iface, &neigh_ifaces, &num_neighs, &elem_indices, &neigh_scheme, 1);

  SC_CHECK_ABORT (num_neighs == 1, "only one faceneighbour currently implemented !");

  /* TODO: outsource as function */
  offset = t8_forest_get_tree_element_offset (forest, itree);
  idata = offset + ielement;

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
