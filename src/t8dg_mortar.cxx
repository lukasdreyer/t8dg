/*
 * t8dg_mortar.c
 *
 *  Created on: Apr 3, 2020
 *      Author: lukas
 */

#include <t8_forest.h>
#include <sc_containers.h>
#include <t8_element_cxx.hxx>

#include "t8dg_mortar.h"
#include "t8dg.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_precomputed_values.h"
#include "t8dg_flux.h"
#include "t8dg_sc_array.h"

/** struct used to save the calculated numerical fluxes at quadrature points*/
struct t8dg_mortar
{
  t8_locidx_t         elem_idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_idata_plus[MAX_SUBFACES];                    /**< Local index of the element corresponding to u_plus */

  int                 iface_minus;
  int                 iface_plus[MAX_SUBFACES];
  int                 num_subfaces;

  t8dg_functionbasis_t *functionbasis;
  int                 number_face_dof;

  int                 bigface_is_ghost;
  int                 subface_is_local[MAX_SUBFACES];

  t8_eclass_t         eclass;
  int                 orientation;
  /*one value for each quadrature point, orientation of u_ */
  sc_array_t         *fluxvalue_minus;                           /**< value of (cu)*.n at face dof */
  /*orientation and sign of right side */
  sc_array_t         *fluxvalue_plus[MAX_SUBFACES];
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/

};

struct t8dg_mortar_array
{
  t8_forest_t         forest;
  sc_array_t         *mortars[MAX_FACES];
  int                 num_local_elements;
  int                 num_total_elements;       /*including ghosts */
  int                 max_num_faces;
  /*for hybrid the number of faces for each tree has to be known */
};

void
t8dg_mortar_destroy (t8dg_mortar_t ** pmortar)
{
  t8dg_mortar        *mortar = *pmortar;
  int                 isubface;
  sc_array_destroy (mortar->fluxvalue_minus);
  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    sc_array_destroy (mortar->fluxvalue_plus[isubface]);
  }
  t8dg_functionbasis_unref (&mortar->functionbasis);
  T8DG_FREE (*pmortar);
  *pmortar = NULL;
}

/*make accessible to ghost elements*/
static int
t8dg_mortar_calculate_face_orientation (const t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface)
{
  t8_element_t       *element;
  t8_eclass_scheme_c *scheme;
  int                 orientation;
  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  scheme = t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest, itree));
  if (scheme->t8_element_is_root_boundary (element, iface)) {
    t8_cmesh_get_face_neighbor (t8_forest_get_cmesh (forest), itree, scheme->t8_element_tree_face (element, iface), NULL, &orientation);
    return orientation;
  }
  else {
    return 0;
  }
}

/*Is only called by local elements*/
t8dg_mortar_t      *
t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface, t8dg_functionbasis_t * element_functionbasis)
{
  t8dg_mortar_t      *mortar;
  t8_locidx_t         idata;
  int                 neigh_itree, neigh_ielement, isubface;
  int                *neigh_ifaces;
  int                 num_neighs;
  int                 own_level;
  t8_element_t       *element, **neigh_elems;
  t8_eclass_scheme_c *neigh_scheme;
  t8_locidx_t        *neigh_idatas;

  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  own_level = t8_forest_get_eclass_scheme (forest, t8_forest_get_eclass (forest, itree))->t8_element_level (element);

  /* use this function to also get the data-indices of the neighbouring element */
  t8_forest_leaf_face_neighbors (forest, itree, element, &neigh_elems, iface, &neigh_ifaces, &num_neighs, &neigh_idatas, &neigh_scheme, 1);

  if (num_neighs == 1 && (neigh_scheme->t8_element_level (neigh_elems[0]) < own_level)) {
    /*the neighbour element is the bigger one */
    if (neigh_idatas[0] < t8_forest_get_num_element (forest)) {
      /*The neighbour element is local */
      t8_forest_get_element (forest, neigh_idatas[0], &neigh_itree);
      neigh_ielement = neigh_idatas[0] - t8_forest_get_tree_element_offset (forest, neigh_itree);
      mortar = t8dg_mortar_new (forest, neigh_itree, neigh_ielement, neigh_ifaces[0], element_functionbasis);
    }
    else {
      T8DG_ABORT ("Not yet implemented");
      //TODO:
      /*The bigger neighbour element is a ghost element */
      /*Check if there is already an allocated mortar for the ghost element */

      /*Add this local subface to the mortar */

      /*Create a new mortar */
    }
  }
  else {
    mortar = T8DG_ALLOC_ZERO (t8dg_mortar_t, 1);
    mortar->bigface_is_ghost = 0;
    idata = t8dg_itree_ielement_to_idata (forest, itree, ielement);

    mortar->elem_idata_minus = idata;
    mortar->num_subfaces = num_neighs;
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      mortar->subface_is_local[isubface] = 1;   /*Since the big face is local, all subfaces are local */
      mortar->elem_idata_plus[isubface] = neigh_idatas[isubface];       /*could be greater than number of local elements -> ghost */
      mortar->iface_plus[isubface] = neigh_ifaces[isubface];    /* get neighbouring face index */
    }
    mortar->iface_minus = iface;
    mortar->functionbasis = t8dg_functionbasis_get_face_functionbasis (element_functionbasis, iface);
    t8dg_functionbasis_ref (mortar->functionbasis);
    mortar->number_face_dof = t8dg_functionbasis_get_num_dof (mortar->functionbasis);

    mortar->eclass = t8dg_functionbasis_get_eclass (mortar->functionbasis);
    mortar->orientation = t8dg_mortar_calculate_face_orientation (forest, itree, ielement, iface);
    /* allocate memory for sc_arrays */
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      mortar->fluxvalue_plus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
    }
    mortar->fluxvalue_minus = sc_array_new_count (sizeof (double), mortar->number_face_dof);
    mortar->valid = 0;
  }

  T8_FREE (neigh_elems);
  T8_FREE (neigh_idatas);
  T8_FREE (neigh_ifaces);
  return mortar;
}

void
t8dg_mortar_sc_array_orient (sc_array_t * array, t8_eclass_t eclass, int orientation)
{
  T8DG_CHECK_ABORT (eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_VERTEX, "Currently only implemented for line faces");
  switch (orientation) {
  case 0:
    return;
  case 1:
    size_t idx;
    double              tmp;
    for (idx = 0; idx < array->elem_count / 2; idx++) {
      tmp = *(double *) sc_array_index (array, idx);
      *(double *) sc_array_index (array, idx) = *(double *) sc_array_index (array, array->elem_count - idx - 1);
      *(double *) sc_array_index (array, array->elem_count - idx - 1) = tmp;
    }
    return;
  default:
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_mortar_fill (t8dg_mortar_t * mortar, t8dg_mortar_fill_data_t * mortar_fill_data)
{
  sc_array_t         *elem_dof_values_minus;
  sc_array_t         *elem_dof_values_plus[MAX_SUBFACES];

  sc_array_t         *face_dof_values_minus_big;
  sc_array_t         *face_dof_values_minus[MAX_SUBFACES];
  sc_array_t         *face_dof_values_plus[MAX_SUBFACES];
  sc_array_t         *face_flux_values_minus[MAX_SUBFACES];

  int                 idof;
  int                 isubface;
  double              fluxvalue;
  double              u_minus_val;
  double              u_plus_val;
  double             *normal_vector;
  double              flux_vec[3];
  double              reference_vertex[3] = { 0, 0, 0 };
  double              image_vertex[3];

  /*TODO: Where to take fb from */
  t8dg_functionbasis_t *functionbasis = t8dg_global_precomputed_values_get_functionbasis (mortar_fill_data->global_values);

  face_dof_values_minus_big = sc_array_new_count (sizeof (double), mortar->number_face_dof);
  elem_dof_values_minus = t8dg_sc_array_block_double_new_view (mortar_fill_data->dof_values, mortar->elem_idata_minus);
  t8dg_functionbasis_transform_element_dof_to_face_dof (functionbasis, mortar->iface_minus, elem_dof_values_minus,
                                                        face_dof_values_minus_big);

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      elem_dof_values_plus[isubface] =
        t8dg_sc_array_block_double_new_view (mortar_fill_data->dof_values, mortar->elem_idata_plus[isubface]);
      face_dof_values_plus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      t8dg_functionbasis_transform_element_dof_to_face_dof (functionbasis, mortar->iface_plus[isubface], elem_dof_values_plus[isubface],
                                                            face_dof_values_plus[isubface]);

      t8dg_mortar_sc_array_orient (face_dof_values_plus[isubface], mortar->eclass, mortar->orientation);
    }
  }

  if (mortar->num_subfaces > 1) {
    /*Interpolate to children that are local */
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (mortar->subface_is_local[isubface]) {
        face_dof_values_minus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
        t8dg_functionbasis_apply_child_interpolation_matrix (mortar->functionbasis, isubface, face_dof_values_minus_big,
                                                             face_dof_values_minus[isubface]);
      }
    }
  }
  else {
    T8DG_ASSERT (mortar->num_subfaces == 1);
    face_dof_values_minus[0] = face_dof_values_minus_big;
  }

  /*Calculate Flux values */
  if (t8dg_functionbasis_is_lagrange (mortar->functionbasis)) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      face_flux_values_minus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      for (idof = 0; idof < mortar->number_face_dof; idof++) {
        double              outward_normal[3];
        normal_vector = t8dg_local_precomputed_values_get_face_normal_vector (mortar_fill_data->local_values, mortar->elem_idata_plus[isubface], mortar->iface_plus[isubface], idof);   /*WRONG!!!! idof must be oriented */
        for (int idim = 0; idim < 3; idim++) {
          outward_normal[idim] = -normal_vector[idim];
        }
        /*TODO: Wrong, use smaller functionbasis and geometry */
        t8dg_functionbasis_get_lagrange_vertex (mortar->functionbasis, idof, reference_vertex);

        t8dg_geometry_transform_reference_vertex_to_image_vertex (mortar_fill_data->geometry_data, reference_vertex, image_vertex);

        t8dg_flux_calulate_flux (mortar_fill_data->flux, image_vertex, flux_vec, mortar_fill_data->time);

        u_minus_val = *(double *) sc_array_index_int (face_dof_values_minus[isubface], idof);
        u_plus_val = *(double *) sc_array_index_int (face_dof_values_plus[isubface], idof);

        fluxvalue = t8dg_flux_calculate_numerical_flux_value (mortar_fill_data->flux, u_minus_val, u_plus_val, flux_vec, outward_normal);
        *(double *) sc_array_index_int (face_flux_values_minus[isubface], idof) = +fluxvalue;
        *(double *) sc_array_index_int (mortar->fluxvalue_plus[isubface], idof) = -fluxvalue;
      }
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  /*Project on big edge!! */
  if (mortar->num_subfaces > 1) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {

      t8dg_precomputed_values_transform_face_child_dof_to_parent_dof (mortar_fill_data->global_values, face_flux_values_minus,
                                                                      mortar->fluxvalue_minus, mortar->num_subfaces,
                                                                      mortar_fill_data->local_values, mortar->elem_idata_plus,
                                                                      mortar->elem_idata_minus, mortar->iface_minus);
    }
  }
  else {
    sc_array_copy (mortar->fluxvalue_minus, face_flux_values_minus[0]);
  }

#if 0
  /*Orient fluxvalue_plus */
  t8dg_mortar_sc_array_orient (mortar->fluxvalue_plus, mortar->eclass, mortar->orientation);    /*TODO: different orientation for higher dimension */
#endif

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    sc_array_destroy (elem_dof_values_plus[isubface]);
    sc_array_destroy (face_dof_values_minus[isubface]);
    sc_array_destroy (face_dof_values_plus[isubface]);
    sc_array_destroy (face_flux_values_minus[isubface]);
  }
  if (mortar->num_subfaces > 1) {
    sc_array_destroy (face_dof_values_minus_big);

  }
  sc_array_destroy (elem_dof_values_minus);

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
  T8DG_ASSERT (iface < mortar_array->max_num_faces);
  *(t8dg_mortar_t **) sc_array_index (mortar_array->mortars[iface], idata) = mortar;
}

static void
t8dg_mortar_array_set_all_pointers (t8dg_mortar_array_t * mortar_array, t8dg_mortar_t * mortar)
{
  int                 isubface;
  if (mortar->bigface_is_ghost) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (mortar->subface_is_local[isubface]) {
        t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_plus[isubface], mortar->iface_plus[isubface], mortar);
      }
    }
  }
  else {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_plus[isubface], mortar->iface_plus[isubface], mortar);
    }
  }
  t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_minus, mortar->iface_minus, mortar);
}

static void
t8dg_mortar_array_set_all_pointers_to_NULL (t8dg_mortar_array_t * mortar_array, t8dg_mortar_t * mortar)
{
  int                 isubface;
  if (mortar->bigface_is_ghost) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (mortar->subface_is_local[isubface]) {
        t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_plus[isubface], mortar->iface_plus[isubface], NULL);
      }
    }
  }
  else {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_plus[isubface], mortar->iface_plus[isubface], NULL);
    }
  }
  t8dg_mortar_array_set_mortar (mortar_array, mortar->elem_idata_minus, mortar->iface_minus, NULL);
}

sc_array_t         *
t8dg_mortar_array_get_oriented_flux (t8dg_mortar_array_t * mortar_array, t8_locidx_t idata, int iface)
{
  t8dg_mortar_t      *mortar;
  int                 isubface;
  mortar = *(t8dg_mortar_t **) t8_sc_array_index_locidx (mortar_array->mortars[iface], idata);
  if (idata == mortar->elem_idata_minus && iface == mortar->iface_minus) {
    return mortar->fluxvalue_minus;
  }
  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface] && mortar->elem_idata_plus[isubface] == idata && mortar->iface_plus[isubface] == iface) {
      return mortar->fluxvalue_plus[isubface];
    }
  }
  T8DG_ABORT ("The element does not border the mortar");
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

      /*TODO: num_faces eclass dependent */
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);

        /*for all faces check wether mortar is already allocated/computed */
        if (mortar == NULL) {
          mortar = t8dg_mortar_new (mortar_array->forest, itree, ielement, iface,
                                    t8dg_global_precomputed_values_get_functionbasis (mortar_fill_data->global_values));
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
  mortar_array->max_num_faces = num_faces;

  /*TODO: split in max_num_faces and tree dependent actual num_faces */
  for (iface = 0; iface < num_faces; iface++) {
    mortar_array->mortars[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), mortar_array->num_total_elements);
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
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
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
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);
        if (mortar != NULL) {
          t8dg_mortar_array_set_all_pointers_to_NULL (mortar_array, mortar);
          t8dg_mortar_destroy (&mortar);
        }
      }

    }
  }
  for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
    sc_array_destroy (mortar_array->mortars[iface]);
    mortar_array->mortars[iface] = 0;
  }
  T8DG_FREE (mortar_array);
  *pmortar_array = NULL;
}
