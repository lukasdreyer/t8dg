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
#include "t8dg_global_precomputed_values.h"
#include "t8dg_flux.h"
#include "t8dg_sc_array.h"

/** struct used to save the calculated numerical fluxes at quadrature points*/
struct t8dg_mortar
{
  int                 number_face_quadrature_points;      /**< The number of (sub-)face quadrature points*/
  t8_locidx_t         elem_idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_idata_plus;                    /**< Local index of the element corresponding to u_plus */

  int                 iface_minus;
  int                 iface_plus;
  int                 num_subfaces;

  t8_eclass_t         eclass;
  int                 orientation;
  /*one value for each quadrature point, orientation of u_ */
  sc_array_t         *u_minus;                          /**< value of u on elem_minus at face quadrature points */
  sc_array_t         *u_plus;                           /**< value of u on elem_plus at face quadrature points */
  sc_array_t         *fluxvalue_minus;                           /**< value of (cu)*.n at face quadrature points */
  /*orientation and sign of right side */
  sc_array_t         *fluxvalue_plus;
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/

};                              /*maybe change to opaque handle */

struct t8dg_mortar_array
{
  t8_forest_t         forest;
  sc_array_t         *mortars[MAX_FACES];
  int                 num_local_elements;
  int                 num_total_elements;
  int                 num_faces;
  /*for hybrid the number of faces for each tree has to be known */
};

void
t8dg_mortar_destroy (t8dg_mortar_t ** pmortar)
{
  sc_array_destroy ((*pmortar)->fluxvalue_plus);
  sc_array_destroy ((*pmortar)->fluxvalue_minus);
  sc_array_destroy ((*pmortar)->u_minus);
  sc_array_destroy ((*pmortar)->u_plus);
  T8DG_FREE (*pmortar);
  *pmortar = NULL;
}

t8dg_mortar_t      *
t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface, t8dg_quadrature_t * quadrature)
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
  mortar->number_face_quadrature_points = t8dg_quadrature_get_num_face_vertices (quadrature, iface);

  mortar->eclass = t8dg_quadrature_get_face_eclass (quadrature, iface);
  mortar->orientation = 0;      //t8dg_mortar_calculate_orientation(mortar);

  /* allocate memory for sc_arrays */
  mortar->u_minus = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->u_plus = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->fluxvalue_plus = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->fluxvalue_minus = sc_array_new_count (sizeof (double), mortar->number_face_quadrature_points);
  mortar->valid = 0;

  T8_FREE (neigh_elems);
  T8_FREE (elem_indices);
  T8_FREE (neigh_ifaces);
  return mortar;
}

void
t8dg_mortar_sc_array_orient (sc_array_t * array, t8_eclass_t eclass, int orientation)
{
  switch (orientation) {
  case 0:
    return;
  default:
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_mortar_fill (t8dg_mortar_t * mortar, t8dg_mortar_fill_data_t * mortar_fill_data)
{
  sc_array_t         *elem_dof_values_minus;
  sc_array_t         *elem_dof_values_plus;

  t8dg_quad_idx_t     iquad;
  double              fluxvalue;
  double              u_minus_quad;
  double              u_plus_quad;
  double             *normal_vector;
  double              flux_vec[3];
  double              reference_vertex[3] = { 0, 0, 0 };
  double              image_vertex[3];

  elem_dof_values_minus = t8dg_sc_array_block_double_new_view (mortar_fill_data->dof_values, mortar->elem_idata_minus);
  elem_dof_values_plus = t8dg_sc_array_block_double_new_view (mortar_fill_data->dof_values, mortar->elem_idata_plus);

  t8dg_global_precomputed_values_transform_element_dof_to_face_quad (mortar_fill_data->global_values, mortar->iface_minus,
                                                                     elem_dof_values_minus, mortar->u_minus);
  t8dg_global_precomputed_values_transform_element_dof_to_face_quad (mortar_fill_data->global_values, mortar->iface_plus,
                                                                     elem_dof_values_plus, mortar->u_plus);
  /*Orient u_plus */
  t8dg_mortar_sc_array_orient (mortar->u_plus, mortar->eclass, mortar->orientation);

  /*normal vector from element_minus to element_plus */
  for (iquad = 0; iquad < mortar->number_face_quadrature_points; iquad++) {
    normal_vector =
      t8dg_local_precomputed_values_get_face_normal_vector (mortar_fill_data->local_values, mortar->elem_idata_minus, mortar->iface_minus,
                                                            iquad);

    t8dg_quadrature_get_face_vertex (t8dg_global_precomputed_values_get_quadrature (mortar_fill_data->global_values),
                                     mortar->iface_minus, iquad, reference_vertex);
    t8dg_geometry_transform_reference_vertex_to_image_vertex (mortar_fill_data->geometry_data, reference_vertex, image_vertex);

    t8dg_flux_calulate_flux (mortar_fill_data->flux, image_vertex, flux_vec, mortar_fill_data->time);

    u_minus_quad = *(double *) t8dg_sc_array_index_quadidx (mortar->u_minus, iquad);
    u_plus_quad = *(double *) t8dg_sc_array_index_quadidx (mortar->u_plus, iquad);

    fluxvalue = t8dg_flux_calculate_numerical_flux_value (mortar_fill_data->flux, u_minus_quad, u_plus_quad, flux_vec, normal_vector);
    *(double *) t8dg_sc_array_index_quadidx (mortar->fluxvalue_minus, iquad) = fluxvalue;
    *(double *) t8dg_sc_array_index_quadidx (mortar->fluxvalue_plus, iquad) = -fluxvalue;
  }
  /*Orient fluxvalue_plus */
  t8dg_mortar_sc_array_orient (mortar->fluxvalue_plus, mortar->eclass, mortar->orientation);    /*TODO: different orientation for higher dimension */

  sc_array_destroy (elem_dof_values_minus);
  sc_array_destroy (elem_dof_values_plus);
  mortar->valid = 1;
}

t8dg_mortar_t      *
t8dg_mortar_array_get_mortar (const t8dg_mortar_array_t * mortar_array, const t8_locidx_t idata, const int iface)
{
  return *(t8dg_mortar_t **) sc_array_index (mortar_array->mortars[iface], idata);
}

static void
t8dg_mortar_array_set_mortar (t8dg_mortar_array_t * mortar_array, const t8_locidx_t idata, const int iface, t8dg_mortar_t * mortar)
{
  T8DG_ASSERT (idata < mortar_array->num_total_elements);
  T8DG_ASSERT (iface < mortar_array->num_faces);
  if (idata < mortar_array->num_local_elements) {
    *(t8dg_mortar_t **) sc_array_index (mortar_array->mortars[iface], idata) = mortar;
  }
}

static void
t8dg_mortar_array_set_all_pointers (t8dg_mortar_array_t * mortar_array, t8dg_mortar_t * mortar)
{
  t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_minus, mortar->iface_minus, mortar);
  t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_plus, mortar->iface_plus, mortar);
}

static void
t8dg_mortar_array_set_all_pointers_to_NULL (t8dg_mortar_array_t * mortar_array, t8dg_mortar_t * mortar)
{
  t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_minus, mortar->iface_minus, NULL);
  t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_plus, mortar->iface_plus, NULL);
}

sc_array_t         *
t8dg_mortar_array_get_oriented_flux (t8dg_mortar_array_t * mortar_array, t8_locidx_t idata, int iface)
{
  t8dg_mortar_t      *mortar;
  mortar = *(t8dg_mortar_t **) t8_sc_array_index_locidx (mortar_array->mortars[iface], idata);
  if (idata == mortar->elem_idata_minus) {
    return mortar->fluxvalue_minus;
  }
  else {
    return mortar->fluxvalue_plus;
  }
}

void
t8dg_mortar_array_fill (t8dg_mortar_array_t * mortar_array, t8dg_mortar_fill_data_t * mortar_fill_data)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface;
  t8dg_mortar_t      *mortar;

  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    mortar_fill_data->geometry_data->itree = itree;

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      mortar_fill_data->geometry_data->ielement = ielement;

      for (iface = 0; iface < mortar_array->num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);

        /*for all faces check wether mortar is already allocated/computed */
        if (mortar == NULL) {
          mortar = t8dg_mortar_new (mortar_array->forest, itree, ielement, iface,
                                    t8dg_global_precomputed_values_get_quadrature (mortar_fill_data->global_values));
          t8dg_mortar_array_set_all_pointers (mortar_array, mortar);
        }
        if (!mortar->valid) {
          t8dg_mortar_fill (mortar, mortar_fill_data);
        }
      }
    }
  }
}

t8dg_mortar_array_t *
t8dg_mortar_array_new_empty (t8_forest_t forest, int num_faces)
{
  int                 iface;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_locidx_t         itree, ielement, idata;
  t8dg_mortar_array_t *mortar_array = T8DG_ALLOC_ZERO (t8dg_mortar_array_t, 1);

  mortar_array->forest = forest;
  mortar_array->num_local_elements = t8_forest_get_num_element (forest);
  mortar_array->num_total_elements = mortar_array->num_local_elements + t8_forest_get_num_ghosts (forest);
  mortar_array->num_faces = num_faces;

  for (iface = 0; iface < num_faces; iface++) {
    mortar_array->mortars[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), mortar_array->num_local_elements);
    num_trees = t8_forest_get_num_local_trees (forest);

    for (itree = 0, idata = 0; itree < num_trees; itree++) {
      num_elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
      for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
        t8dg_mortar_array_set_mortar (mortar_array, idata, iface, NULL);
      }
    }
  }
  return mortar_array;
}

void
t8dg_mortar_array_invalidate_all (t8dg_mortar_array_t * mortar_array)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface;
  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);

  t8dg_mortar_t      *mortar;

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      for (iface = 0; iface < mortar_array->num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);
        if (mortar != NULL) {
          mortar->valid = 0;
        }
      }
    }
  }
}

void
t8dg_mortar_array_destroy (t8dg_mortar_array_t ** pmortar_array)
{

  t8_locidx_t         itree, ielement, idata, iface;
  t8_locidx_t         num_trees, num_elems_in_tree;

  t8dg_mortar_t      *mortar;
  t8dg_mortar_array_t *mortar_array;
  mortar_array = *pmortar_array;

  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      for (iface = 0; iface < mortar_array->num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);
        if (mortar != NULL) {
          t8dg_mortar_array_set_all_pointers_to_NULL (mortar_array, mortar);
          t8dg_mortar_destroy (&mortar);
        }
      }

    }
  }
  for (iface = 0; iface < mortar_array->num_faces; iface++) {
    sc_array_destroy (mortar_array->mortars[iface]);
    mortar_array->mortars[iface] = 0;
  }
  T8DG_FREE (mortar_array);
  *pmortar_array = NULL;
}
