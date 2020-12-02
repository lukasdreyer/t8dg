#include "t8dg.h"
#include "t8dg_global_values.h"
#include "t8dg_local_values.h"
#include "t8dg_values.h"
#include "t8dg_sc_array.h"
#include "t8dg_mortar.h"
#include "t8dg_dof.h"
#include "t8dg_flux.h"

typedef struct t8dg_values t8dg_values_t;

struct t8dg_values
{
  int                 dim;
  t8dg_local_values_t *local_values;
  t8dg_local_values_t *local_values_adapt;

  t8dg_mortar_array_t *mortar_array;
  t8dg_global_values_t **global_values_array;
  t8dg_coarse_geometry_t *coarse_geometry;

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
};

t8dg_values_t      *
t8dg_values_new_LGL_hypercube (int dim, int num_LGL_vertices, t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest)
{
  t8dg_values_t      *return_values;
  int                 eclass;
  switch (dim) {
  case 1:
    eclass = T8_ECLASS_LINE;
    break;
  case 2:
    eclass = T8_ECLASS_QUAD;
    break;
  case 3:
    eclass = T8_ECLASS_HEX;
    break;

  default:
    eclass = -1;
    break;
  }

  return_values = T8DG_ALLOC_ZERO (t8dg_values_t, 1);
  return_values->dim = dim;
  return_values->global_values_array = T8DG_ALLOC_ZERO (t8dg_global_values_t *, T8_ECLASS_COUNT);
  return_values->global_values_array[eclass] = t8dg_global_values_new_hypercube_LGL (dim, num_LGL_vertices);
  return_values->local_values = t8dg_local_values_new (forest, return_values->global_values_array, coarse_geometry);
  return_values->coarse_geometry = coarse_geometry;
  return_values->mortar_array = t8dg_mortar_array_new_empty (forest, return_values->local_values);
  return_values->forest = forest;
  t8_forest_ref (forest);
  return return_values;
}

void
t8dg_values_destroy (t8dg_values_t ** p_values)
{
  t8dg_values_t      *values;
  values = *p_values;
  t8_eclass_t         eclass;

  t8dg_mortar_array_destroy (&values->mortar_array);

  t8dg_local_values_destroy (&values->local_values);

  for (eclass = 0; eclass < T8_ECLASS_COUNT; eclass++) {
    if (values->global_values_array[eclass] != NULL) {
      t8dg_global_values_destroy (&values->global_values_array[eclass]);
    }
  }
  t8_forest_unref (&values->forest);
  T8DG_FREE (values->global_values_array);
  T8DG_FREE (values);
  *p_values = NULL;
}

void
t8dg_values_apply_mass_matrix (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8dg_element_dof_values_t *src_element_dofs;
  t8dg_element_dof_values_t *dest_element_dofs;

  for (itree = 0, idata = 0; itree < t8_forest_get_num_local_trees (values->forest); itree++) {
    for (ielement = 0; ielement < t8_forest_get_tree_num_elements (values->forest, itree); ielement++, idata++) {
      src_element_dofs = t8dg_dof_values_new_element_dof_values_view (src_dof, itree, ielement);
      dest_element_dofs = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);

      t8dg_local_values_apply_element_mass_matrix (values->local_values, itree, ielement, src_element_dofs, dest_element_dofs);

      t8dg_element_dof_values_destroy (&src_element_dofs);
      t8dg_element_dof_values_destroy (&dest_element_dofs);
    }
  }
}

/*same for inverse*/
void
t8dg_values_apply_inverse_mass_matrix (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_elements;
  t8dg_element_dof_values_t *src_element_dofs;
  t8dg_element_dof_values_t *dest_element_dofs;

  for (itree = 0, idata = 0; itree < t8_forest_get_num_local_trees (values->forest); itree++) {
    for (ielement = 0; ielement < t8_forest_get_tree_num_elements (values->forest, itree); ielement++, idata++) {
      src_element_dofs = t8dg_dof_values_new_element_dof_values_view (src_dof, itree, ielement);
      dest_element_dofs = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);

      t8dg_local_values_apply_element_inverse_mass_matrix (values->local_values, itree, ielement, src_element_dofs, dest_element_dofs);

      t8dg_element_dof_values_destroy (&src_element_dofs);
      t8dg_element_dof_values_destroy (&dest_element_dofs);
    }
  }
}

void
t8dg_values_apply_directional_stiffness_matrix_linear_flux_fn (t8dg_values_t * values, int idirection, t8dg_linear_flux1D_fn * flux_fn,
                                                               void *flux_data, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof)
{

}

void
t8dg_values_apply_stiffness_matrix_linear_flux_fn3D (t8dg_values_t * values, t8dg_linear_flux3D_fn flux_fn, void *flux_data,
                                                     t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_dofidx_t       num_dof;
  t8dg_quadidx_t      num_quad_vertices;
  int                 idim;

  sc_array_t         *element_quad_values;
  sc_array_t         *element_flux_quad_values;
  sc_array_t         *element_dof_values;
  sc_array_t         *element_res_dof_values;
  sc_array_t         *element_res_summand_dof_values;
  t8_eclass_t         eclass;

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    eclass = t8_forest_get_eclass (values->forest, itree);
    num_quad_vertices = t8dg_global_values_get_num_elem_quad (values->global_values_array[eclass]);
    num_dof = t8dg_global_values_get_num_dof (values->global_values_array[eclass]);

    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dof_values = t8dg_dof_values_new_element_dof_values_view (src_dof, itree, ielement);
//      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dof_values));
      element_quad_values = sc_array_new_count (sizeof (double), num_quad_vertices);
      element_flux_quad_values = sc_array_new_count (sizeof (double), num_quad_vertices);
      element_res_summand_dof_values = sc_array_new_count (sizeof (double), num_dof);
      element_res_dof_values = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);

      t8dg_local_values_element_multiply_trafo_quad_weight (values->local_values, itree, ielement, element_dof_values, element_quad_values);    /*TODO: Vandermonde? */
//      t8_debugf ("element_quad_values*qtw\n");
//      t8dg_sc_array_block_double_debug_print (element_quad_values);

      t8dg_element_dof_values_set_zero (element_res_dof_values);
      for (idim = 0; idim < values->dim; idim++) {
/* TODO: 
        t8dg_local_values_element_multiply_flux_value (problem->local_values, problem->description.flux, &geometry_data,
                                                                   t8dg_global_values_get_quadrature (problem->global_values),
                                                                   t8dg_timestepping_data_get_current_time (problem->time_data),
                                                                   idim, element_quad_values, element_flux_quad_values);

*/

//        t8_debugf ("element_dof_derivative_values\n");
//        t8dg_sc_array_block_double_debug_print (element_dof_derivative_values);
        t8dg_global_values_element_apply_derivative_matrix_transpose (values->global_values_array[eclass], idim,
                                                                      element_flux_quad_values, element_res_summand_dof_values);
//        t8_debugf ("element_res_summand_dof_values\n");
//        t8dg_sc_array_block_double_debug_print (element_res_summand_dof_values);
        t8dg_element_dof_values_axpy (1, element_res_summand_dof_values, element_res_dof_values);
      }

//      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_res_dof_values));

      sc_array_destroy (element_dof_values);
      sc_array_destroy (element_quad_values);
      sc_array_destroy (element_flux_quad_values);
      sc_array_destroy (element_res_dof_values);
      sc_array_destroy (element_res_summand_dof_values);
    }
  }
}

void
t8dg_values_apply_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof,
                                      t8dg_linear_flux3D_fn linear_flux, void *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux,
                                      void *numerical_flux_data, double time)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_element_dof_values_t *element_dest_dof;

  t8dg_mortar_array_calculate_linear_flux3D (values->mortar_array, src_dof, linear_flux, numerical_flux, flux_data, numerical_flux_data,
                                             time);

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dest_dof = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);
      t8dg_mortar_array_apply_element_boundary_integral (values->mortar_array, itree, ielement, element_dest_dof);

      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dest_dof));
      sc_array_destroy (element_dest_dof);
    }
  }

  t8dg_mortar_array_invalidate_all (values->mortar_array);

}

void
t8dg_values_ghost_exchange (t8dg_values_t * values)
{
  t8dg_local_values_ghost_exchange (values->local_values);
}

t8dg_global_values_t *
t8dg_values_get_global_values (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement)
{
  t8_eclass_t         eclass;
  eclass = t8_forest_get_eclass (values->forest, itree);
  return values->global_values_array[eclass];
}

t8dg_global_values_t *
t8dg_values_get_global_values_adapt (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement)
{
  T8DG_ASSERT (values->forest_adapt != NULL);
  t8_eclass_t         eclass;
  eclass = t8_forest_get_eclass (values->forest, itree);
  return values->global_values_array[eclass];
}

t8dg_global_values_t **
t8dg_values_get_global_values_array (t8dg_values_t * values)
{
  return values->global_values_array;
}

void
t8dg_values_copy_element_values (t8dg_values_t * values, t8_locidx_t idata_old, t8_locidx_t idata_new)
{
  T8DG_ASSERT (values->local_values_adapt != NULL);
  t8dg_local_values_copy_element_values (values->local_values, idata_old, values->local_values_adapt, idata_new);
}

void
t8dg_values_set_element_adapt (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement_new)
{
  t8dg_local_values_set_element (values->local_values_adapt, itree, ielement_new);
}

/*Interpolation only dependent on functionbasis*/
void
t8dg_values_transform_parent_dof_to_child_dof (t8dg_values_t * values, t8dg_dof_values_t * dof_values, t8dg_dof_values_t * dof_values_adapt,
                                               t8_locidx_t itree, t8_locidx_t ielem_parent_old, t8_locidx_t ielem_child_new, int ichild)
{
  t8dg_global_values_t *global_values;
  t8dg_element_dof_values_t *parent_dof;
  t8dg_element_dof_values_t *child_dof;

  parent_dof = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielem_parent_old);
  child_dof = t8dg_dof_values_new_element_dof_values_view (dof_values_adapt, itree, ielem_child_new);

  t8dg_global_values_transform_element_dof_to_child_dof (global_values, parent_dof, child_dof, ichild);
  t8dg_element_dof_values_destroy (&parent_dof);
  t8dg_element_dof_values_destroy (&child_dof);
}

/*L2 projection needs geometry information */
void
t8dg_values_transform_child_dof_to_parent_dof (t8dg_values_t * values, t8dg_dof_values_t * dof_values, t8dg_dof_values_t * dof_values_adapt,
                                               t8_locidx_t itree, t8_locidx_t num_children, t8_locidx_t ielem_first_child_old,
                                               t8_locidx_t ielem_parent_new)
{
  t8dg_element_dof_values_t *child_dof[MAX_SUBELEMENTS];
  t8dg_element_dof_values_t *parent_dof;
  int                 ichild;
  for (ichild = 0; ichild < num_children; ichild++) {
    child_dof[ichild] = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielem_first_child_old + ichild);
  }
  parent_dof = t8dg_dof_values_new_element_dof_values_view (dof_values_adapt, itree, ielem_parent_new);

  t8dg_local_values_transform_child_dof_to_parent_dof (values->local_values, values->local_values_adapt, child_dof, parent_dof,
                                                       num_children, itree, ielem_first_child_old, ielem_parent_new);
  for (ichild = 0; ichild < num_children; ichild++) {
    t8dg_element_dof_values_destroy (&child_dof[ichild]);
  }
  t8dg_element_dof_values_destroy (&parent_dof);
}

void
t8dg_values_allocate_adapt (t8dg_values_t * values, t8_forest_t forest_adapt)
{
  values->local_values_adapt = t8dg_local_values_new (forest_adapt, values->global_values_array, values->coarse_geometry);
  values->forest_adapt = forest_adapt;
  t8_forest_ref (values->forest_adapt);
}

void
t8dg_values_cleanup_adapt (t8dg_values_t * values)
{
  /* clean the old element data */
  t8_forest_unref (&values->forest);
  values->forest = values->forest_adapt;
  values->forest_adapt = NULL;

  /* Set the elem data to the adapted elem data */
  t8dg_local_values_destroy (&values->local_values);
  values->local_values = values->local_values_adapt;
  values->local_values_adapt = NULL;

/*TODO: Is the order important?*/

  /*Create new mortar arrays */
  t8dg_mortar_array_destroy (&values->mortar_array);
  values->mortar_array = t8dg_mortar_array_new_empty (values->forest, values->local_values_adapt);
}

void
t8dg_values_partition (t8dg_values_t * values, t8_forest_t forest_partition)
{
  /* Partition local precomputed values */
  t8dg_local_values_t *local_values_partition;
  local_values_partition = t8dg_local_values_new (forest_partition, values->global_values_array, values->coarse_geometry);
  t8dg_local_values_partition (values->local_values, local_values_partition);

  t8dg_local_values_destroy (&values->local_values);
  values->local_values = local_values_partition;

  t8_debugf ("[VALUES] Done partition local_data\n");

  t8_debugf (" [ADVECT] begin mortars destroy\n");
  t8dg_mortar_array_destroy (&values->mortar_array);

  t8_forest_unref (&values->forest);
  values->forest = forest_partition;
  t8_forest_ref (values->forest);

  values->mortar_array = t8dg_mortar_array_new_empty (forest_partition, values->local_values);
  t8_debugf (" [ADVECT] Done partition\n");

}

double
t8dg_values_element_norm_l2_squared (t8dg_values_t * values, t8dg_element_dof_values_t * element_dof, t8_locidx_t itree,
                                     t8_locidx_t ielement)
{
  T8DG_ABORT ("Not implemented\n");
}
