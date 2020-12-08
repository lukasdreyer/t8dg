/*
 * t8dg_mortar.c
 *
 *  Created on: Apr 3, 2020
 *      Author: lukas
 */

#include <t8_forest.h>
#include <sc_containers.h>
#include <t8_element_cxx.hxx>
#include <t8_vec.h>

#include "t8dg_mortar.h"
#include "t8dg.h"
#include "t8dg_global_values.h"
#include "t8dg_values.h"
#include "t8dg_flux.h"
#include "t8dg_sc_array.h"
#include "t8dg_geometry.h"

/** struct used to save the calculated numerical fluxes at quadrature points*/
struct t8dg_mortar
{
  t8_locidx_t         elem_itree_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_itree_plus;                    /**< All faceneighbours are on the same tree */
  t8_locidx_t         elem_ielement_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_ielement_plus[MAX_SUBFACES];                    /**< Local index of the element corresponding to u_plus */

  t8_locidx_t         elem_idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_idata_plus[MAX_SUBFACES];                    /**< Local index of the element corresponding to u_plus */

  int                 iface_minus;
  int                 iface_plus[MAX_SUBFACES];
  int                 num_subfaces;

  t8dg_functionbasis_t *face_functionbasis;    /**< face functionbasis with orientation from elem_minus*/
  t8dg_functionbasis_t *functionbasis_minus;

  int                 number_face_dof;

  int                 bigface_is_ghost;
  int                 subface_is_local[MAX_SUBFACES];

  t8_eclass_t         eclass;
  int                 orientation;
  /*one value for each quadrature point, orientation of u_ */
  t8dg_face_quad_values_t *fluxvalue_minus;                                   /**< value of (cu)*.n at face dof */
  /*orientation and sign of right side */
  t8dg_face_quad_values_t *fluxvalue_plus[MAX_SUBFACES];
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/

  t8dg_mortar_array_t *mortar_array;
};

struct t8dg_mortar_array
{
  t8_forest_t         forest;
  sc_array_t         *mortars[MAX_FACES];

  int                 num_local_elements;
  int                 num_total_elements;       /*including ghosts */
  int                 max_num_faces;
  /*for hybrid the number of faces for each tree has to be known */
  t8dg_local_values_t *local_values;
};

void
t8dg_mortar_destroy (t8dg_mortar_t ** pmortar)
{
  t8dg_mortar        *mortar = *pmortar;
  int                 isubface;
  if (!mortar->bigface_is_ghost) {
    t8dg_debugf ("destroy fluxvalue_minus %p, idata: %i, iface: %i\n", (void *) mortar->fluxvalue_minus, mortar->elem_idata_minus,
                 mortar->iface_minus);
    sc_array_destroy (mortar->fluxvalue_minus);
    t8dg_debugf ("destroyed fluxvalue_minus, idata: %i, iface: %i\n", mortar->elem_idata_minus, mortar->iface_minus);
  }
  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      sc_array_destroy (mortar->fluxvalue_plus[isubface]);
    }
  }
  t8dg_functionbasis_unref (&mortar->functionbasis_minus);
  t8dg_functionbasis_unref (&mortar->face_functionbasis);

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

static t8dg_mortar_t *
t8dg_mortar_array_get_mortar (const t8dg_mortar_array_t * mortar_array, const t8_locidx_t idata, const int iface)
{
  return *(t8dg_mortar_t **) sc_array_index (mortar_array->mortars[iface], idata);
}

/*Is only called by local elements*/
static t8dg_mortar_t *
t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface, t8dg_mortar_array_t * mortar_array)
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
  t8dg_global_values_t *global_values;
  t8dg_functionbasis_t *element_functionbasis;

  global_values = t8dg_local_values_get_global_values (mortar_array->local_values, itree, ielement);
  element_functionbasis = t8dg_global_values_get_functionbasis (global_values);

  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  own_level = t8_forest_get_eclass_scheme (forest, t8_forest_get_eclass (forest, itree))->t8_element_level (element);

  /* use this function to also get the data-indices of the neighbouring element */
  t8_forest_leaf_face_neighbors (forest, itree, element, &neigh_elems, iface, &neigh_ifaces, &num_neighs, &neigh_idatas, &neigh_scheme, 1);

  if (num_neighs == 1 && (neigh_scheme->t8_element_level (neigh_elems[0]) < own_level)) {
    /*the neighbour element is the bigger one */
    if (neigh_idatas[0] < t8_forest_get_num_element (forest)) {
      /*The neighbour element is local */
      t8_forest_get_element (forest, neigh_idatas[0], &neigh_itree);    /*Only needed for neighbour itree */
      neigh_ielement = neigh_idatas[0] - t8_forest_get_tree_element_offset (forest, neigh_itree);
      mortar = t8dg_mortar_new (forest, neigh_itree, neigh_ielement, neigh_ifaces[0], mortar_array);
    }
    else {
      mortar = t8dg_mortar_array_get_mortar (mortar_array, neigh_idatas[0], neigh_ifaces[0]);
      if (mortar == NULL) {
        /*Create new mortar, since neighbour element is not local */
        mortar = T8DG_ALLOC_ZERO (t8dg_mortar_t, 1);
        mortar->mortar_array = mortar_array;
        mortar->bigface_is_ghost = 1;
        mortar->elem_itree_minus = neigh_itree;

        mortar->elem_idata_minus = neigh_idatas[0];
        mortar->elem_itree_minus = -1;  /*since elem_minus is not local */
        mortar->elem_ielement_minus = mortar->elem_idata_minus;

        mortar->iface_minus = neigh_ifaces[0];
        mortar->face_functionbasis = t8dg_functionbasis_get_face_functionbasis (element_functionbasis, iface);
        mortar->functionbasis_minus = element_functionbasis;
        mortar->num_subfaces = t8dg_functionbasis_get_num_children (mortar->face_functionbasis);
        t8dg_functionbasis_ref (mortar->face_functionbasis);
        t8dg_functionbasis_ref (mortar->functionbasis_minus);
        mortar->number_face_dof = t8dg_functionbasis_get_num_dof (mortar->face_functionbasis);
        mortar->eclass = t8dg_functionbasis_get_eclass (mortar->face_functionbasis);
        mortar->orientation = t8dg_mortar_calculate_face_orientation (forest, itree, ielement, iface);
      }
      t8dg_debugf ("Adress: %p, idata: %i, num_local: %i\n", (void *) mortar, neigh_idatas[0], t8_forest_get_num_element (forest));
      t8dg_debugf ("mortar: num_subfaces: %i\n", mortar->num_subfaces);
      t8dg_debugf ("Adress: %p\n", (void *) mortar);

      t8_eclass_scheme_c *boundary_scheme, *element_scheme;
      t8_element_t       *boundary_element;

      boundary_scheme = t8_forest_get_eclass_scheme (forest, mortar->eclass);
      boundary_scheme->t8_element_new (1, &boundary_element);
      element = t8_forest_get_element_in_tree (forest, itree, ielement);
      element_scheme = t8_forest_get_eclass_scheme (forest, t8dg_functionbasis_get_eclass (element_functionbasis));
      element_scheme->t8_element_boundary_face (element, iface, boundary_element, boundary_scheme);

      isubface = boundary_scheme->t8_element_child_id (boundary_element);
      boundary_scheme->t8_element_destroy (1, &boundary_element);

      idata = t8dg_itree_ielement_to_idata (forest, itree, ielement);

      t8dg_debugf ("Halfmortar set: isubface: %i, idata: %i, iface: %i\n", isubface, idata, iface);

      mortar->elem_idata_plus[isubface] = idata;
      mortar->elem_itree_plus = itree;
      mortar->elem_ielement_plus[isubface] = ielement;
      mortar->iface_plus[isubface] = iface;
      T8DG_ASSERT (mortar->subface_is_local[isubface] == 0);
      mortar->subface_is_local[isubface] = 1;
      mortar->fluxvalue_plus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      t8dg_debugf ("add fluxvalue_plus[isubface] = %p", (void *) mortar->fluxvalue_plus[isubface]);
      mortar->valid = 0;
    }
  }
  else {
    /*own element is bigger or same size */
    mortar = T8DG_ALLOC_ZERO (t8dg_mortar_t, 1);
    mortar->mortar_array = mortar_array;
    mortar->bigface_is_ghost = 0;
    idata = t8dg_itree_ielement_to_idata (forest, itree, ielement);

    mortar->elem_idata_minus = idata;
    mortar->elem_ielement_minus = ielement;
    mortar->elem_itree_minus = itree;
    mortar->num_subfaces = num_neighs;
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      mortar->subface_is_local[isubface] = 1;   /*Since the big face is local, all subfaces are local */
      mortar->elem_idata_plus[isubface] = neigh_idatas[isubface];       /*could be greater than number of local elements -> ghost */

      if (neigh_idatas[isubface] >= mortar_array->num_local_elements) {
        mortar->elem_itree_plus = -1;
        mortar->elem_ielement_plus[isubface] = mortar->elem_idata_plus[isubface];       /*neighbouring element is not local */
      }
      else {
        t8_forest_get_element (forest, neigh_idatas[isubface], &neigh_itree);   /*Only needed for neighbour itree */
        neigh_ielement = neigh_idatas[isubface] - t8_forest_get_tree_element_offset (forest, neigh_itree);

        mortar->elem_ielement_plus[isubface] = neigh_ielement;
        mortar->elem_itree_plus = neigh_itree;
      }
      mortar->iface_plus[isubface] = neigh_ifaces[isubface];    /* get neighbouring face index */
    }
    mortar->iface_minus = iface;
    mortar->face_functionbasis = t8dg_functionbasis_get_face_functionbasis (element_functionbasis, iface);
    mortar->functionbasis_minus = element_functionbasis;
    t8dg_functionbasis_ref (mortar->face_functionbasis);
    t8dg_functionbasis_ref (mortar->functionbasis_minus);
    mortar->number_face_dof = t8dg_functionbasis_get_num_dof (mortar->face_functionbasis);

    mortar->eclass = t8dg_functionbasis_get_eclass (mortar->face_functionbasis);
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

static void
t8dg_mortar_calculate_linear_flux3D (t8dg_mortar_t * mortar, t8dg_dof_values_t * dof_values, t8dg_linear_flux3D_fn linear_flux,
                                     void *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux, void *numerical_flux_data,
                                     double time)
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
  double              flux_vec[3];
  double              reference_vertex[3] = { 0, 0, 0 };
  double              image_vertex[3];

  t8dg_mortar_array_t *mortar_array;
  t8dg_functionbasis_t *element_functionbasis = mortar->functionbasis_minus;

  mortar_array = mortar->mortar_array;

  face_dof_values_minus_big = sc_array_new_count (sizeof (double), mortar->number_face_dof);
  elem_dof_values_minus = t8dg_dof_values_new_element_dof_values_view (dof_values, mortar->elem_itree_minus, mortar->elem_ielement_minus);
  t8dg_functionbasis_transform_element_dof_to_face_dof (element_functionbasis, mortar->iface_minus, elem_dof_values_minus,
                                                        face_dof_values_minus_big);

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      elem_dof_values_plus[isubface] =
        t8dg_dof_values_new_element_dof_values_view (dof_values, mortar->elem_itree_plus, mortar->elem_ielement_plus[isubface]);
      face_dof_values_plus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      t8dg_functionbasis_transform_element_dof_to_face_dof (element_functionbasis, mortar->iface_plus[isubface],
                                                            elem_dof_values_plus[isubface], face_dof_values_plus[isubface]);
    }
  }

  if (mortar->num_subfaces > 1) {
    /*Interpolate to children that are local */
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (mortar->subface_is_local[isubface]) {
        face_dof_values_minus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
        t8dg_functionbasis_apply_child_interpolation_matrix (mortar->face_functionbasis, isubface, face_dof_values_minus_big,
                                                             face_dof_values_minus[isubface]);
        t8dg_mortar_sc_array_orient (face_dof_values_minus[isubface], mortar->eclass, mortar->orientation);
      }
    }
  }
  else {
    T8DG_ASSERT (mortar->num_subfaces == 1);
    face_dof_values_minus[0] = face_dof_values_minus_big;
    t8dg_mortar_sc_array_orient (face_dof_values_minus[0], mortar->eclass, mortar->orientation);
  }

  /*Calculate Flux values */
  t8dg_debugf ("calculate flux for mortar %p\n", (void *) mortar);
  if (t8dg_functionbasis_is_lagrange (mortar->face_functionbasis)) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (!mortar->subface_is_local[isubface]) {
        continue;
      }
      t8dg_debugf ("calculate flux for idata %i, iface %i ,isubface %i\n", mortar->elem_idata_plus[isubface], mortar->iface_plus[isubface],
                   isubface);
      face_flux_values_minus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      T8DG_ASSERT (mortar->fluxvalue_plus[isubface]->elem_count == (size_t) mortar->number_face_dof);
      for (idof = 0; idof < mortar->number_face_dof; idof++) {
        double              outward_normal[3] = { 0, 0, 0 };

        t8dg_local_values_get_face_normal_vector (mortar_array->local_values, mortar->elem_itree_plus, mortar->elem_ielement_plus[isubface],
                                                  mortar->iface_plus[isubface], idof, outward_normal);
        t8_vec_ax (outward_normal, -1); //need outward normal to elem_minus at the face nodal_basis_vertices of elem_plus,

        t8dg_functionbasis_t *functionbasis_plus;
        //t8dg_global_values_t *global_values_plus;
        /*TODO: Wrong, use smaller functionbasis and geometry */
        functionbasis_plus = mortar->functionbasis_minus;

        t8dg_functionbasis_get_lagrange_vertex (functionbasis_plus, idof, reference_vertex);

        t8dg_geometry_transform_reference_vertex_to_image_vertex (t8dg_local_values_get_coarse_geometry (mortar_array->local_values),
                                                                  mortar_array->forest, mortar->elem_itree_plus,
                                                                  mortar->elem_ielement_plus[isubface], reference_vertex, image_vertex);
        linear_flux (image_vertex, flux_vec, time, flux_data);

        u_minus_val = t8dg_face_dof_values_get_value (face_dof_values_minus[isubface], idof);
        u_plus_val = t8dg_face_dof_values_get_value (face_dof_values_plus[isubface], idof);

        fluxvalue = numerical_flux (u_minus_val, u_plus_val, flux_vec, outward_normal, numerical_flux_data);
        t8dg_debugf ("fluxvalue: %f\n", fluxvalue);
//        T8DG_ASSERT (fluxvalue == fluxvalue && fabs (fluxvalue) < 1e200);
        *(double *) sc_array_index_int (face_flux_values_minus[isubface], idof) = +fluxvalue;   /*TODO: check! */
        *(double *) sc_array_index_int (mortar->fluxvalue_plus[isubface], idof) = -fluxvalue;
      }
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (face_flux_values_minus[isubface]));
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_plus[isubface]));
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  /*Project on big face if it is local */
  if (!(mortar->bigface_is_ghost)) {
    if (mortar->num_subfaces > 1) {
      for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
        t8dg_mortar_sc_array_orient (face_flux_values_minus[isubface], mortar->eclass, mortar->orientation);    /*TODO: Move to transform_orient! */
      }

      t8dg_local_values_transform_orient_face_child_dof_to_parent_dof_hanging_nodes (mortar_array->local_values, face_flux_values_minus,
                                                                                     mortar->fluxvalue_minus, mortar->num_subfaces,
                                                                                     mortar->elem_itree_plus, mortar->elem_ielement_plus,
                                                                                     mortar->iface_plus, mortar->elem_itree_minus,
                                                                                     mortar->elem_ielement_minus, mortar->iface_minus);
    }
    else {
      t8dg_mortar_sc_array_orient (face_flux_values_minus[0], mortar->eclass, mortar->orientation);
      sc_array_copy (mortar->fluxvalue_minus, face_flux_values_minus[0]);

    }
  }

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      sc_array_destroy (elem_dof_values_plus[isubface]);
      sc_array_destroy (face_dof_values_minus[isubface]);
      sc_array_destroy (face_flux_values_minus[isubface]);
      sc_array_destroy (face_dof_values_plus[isubface]);
    }
  }
  if (mortar->num_subfaces > 1) {
    sc_array_destroy (face_dof_values_minus_big);

  }
  sc_array_destroy (elem_dof_values_minus);

  mortar->valid = 1;
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
        t8dg_debugf ("set mortar %p at: idata:%i, iface:%i\n", (void *) mortar, mortar->elem_idata_plus[isubface],
                     mortar->iface_plus[isubface]);
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
  t8dg_debugf ("idata: %i, iface: %i \n", idata, iface);
  mortar = *(t8dg_mortar_t **) t8_sc_array_index_locidx (mortar_array->mortars[iface], idata);
  T8DG_ASSERT (mortar != NULL);
  t8dg_debugf ("adress: %p \n", (void *) mortar);
  if (idata == mortar->elem_idata_minus && iface == mortar->iface_minus) {
    T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_minus));
    return mortar->fluxvalue_minus;
  }
  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface] && mortar->elem_idata_plus[isubface] == idata && mortar->iface_plus[isubface] == iface) {
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_plus[isubface]));
      return mortar->fluxvalue_plus[isubface];
    }
  }
  T8DG_ABORT ("The element does not border the mortar");
}

void
t8dg_mortar_array_calculate_linear_flux3D (t8dg_mortar_array_t * mortar_array, t8dg_dof_values_t * dof_values,
                                           t8dg_linear_flux3D_fn linear_flux, void *flux_data,
                                           t8dg_numerical_linear_flux3D_fn numerical_flux, void *numerical_flux_data, double time)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface;
  t8dg_mortar_t      *mortar;

  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {

      /*TODO: num_faces eclass dependent */
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);

        /*for all faces check wether mortar is already allocated/computed */
        if (mortar == NULL) {
          mortar = t8dg_mortar_new (mortar_array->forest, itree, ielement, iface, mortar_array);
          t8dg_mortar_array_set_all_pointers (mortar_array, mortar);
        }
        if (!mortar->valid) {
          t8dg_mortar_calculate_linear_flux3D (mortar, dof_values, linear_flux, flux_data, numerical_flux, numerical_flux_data, time);
        }
        T8DG_ASSERT (t8dg_mortar_array_get_mortar (mortar_array, idata, iface) != NULL);
      }
    }
  }
}

t8dg_mortar_array_t *
t8dg_mortar_array_new_empty (t8_forest_t forest, t8dg_local_values_t * local_values)
{
  int                 iface;
  t8_locidx_t         idata;
  t8dg_mortar_array_t *mortar_array = T8DG_ALLOC_ZERO (t8dg_mortar_array_t, 1);

  t8_forest_ref (forest);
  mortar_array->forest = forest;
  mortar_array->num_local_elements = t8_forest_get_num_element (forest);
  mortar_array->num_total_elements = mortar_array->num_local_elements + t8_forest_get_num_ghosts (forest);
  mortar_array->local_values = local_values;
  mortar_array->max_num_faces = t8dg_local_values_get_max_num_faces (local_values);

  /*TODO: split in max_num_faces and tree dependent actual num_faces */
  for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
    mortar_array->mortars[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), mortar_array->num_total_elements);

    for (idata = 0; idata < mortar_array->num_total_elements; idata++) {
      t8dg_mortar_array_set_mortar (mortar_array, idata, iface, NULL);
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
  t8_forest_unref (&mortar_array->forest);
  T8DG_FREE (mortar_array);
  *pmortar_array = NULL;
}

/*TODO: Who does ghost exchange?*/
void
t8dg_mortar_array_apply_element_boundary_integral (t8dg_mortar_array_t * mortar_array,
                                                   t8_locidx_t itree, t8_locidx_t ielement, sc_array_t * element_result_dof)
{
  int                 iface, num_faces;

  sc_array_t         *face_flux_dof;
  sc_array_t         *face_dof;
  sc_array_t         *summand;

  t8_locidx_t         idata = t8dg_itree_ielement_to_idata (mortar_array->forest, itree, ielement);

  t8_eclass_t         eclass;
  eclass = t8_forest_get_eclass (mortar_array->forest, itree);

  num_faces = t8_eclass_num_faces[eclass];

  summand = t8dg_element_dof_values_duplicate (element_result_dof);
  t8dg_element_dof_values_set_zero (element_result_dof);

  for (iface = 0; iface < num_faces; iface++) {
    face_flux_dof = t8dg_mortar_array_get_oriented_flux (mortar_array, idata, iface);
    face_dof = t8dg_face_dof_values_duplicate (face_flux_dof);

    t8dg_local_values_apply_face_mass_matrix (mortar_array->local_values, itree, ielement, iface, face_flux_dof, face_dof);

    t8dg_functionbasis_transform_face_dof_to_element_dof (t8dg_global_values_get_functionbasis
                                                          (t8dg_local_values_get_global_values
                                                           (mortar_array->local_values, itree, ielement)), iface, face_dof, summand);

    T8DG_ASSERT (t8dg_element_dof_values_is_valid (summand));
    t8dg_element_dof_values_axpy (1, summand, element_result_dof);
    sc_array_destroy (face_dof);
  }
  sc_array_destroy (summand);
}
