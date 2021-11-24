#include "t8dg.h"
#include "t8dg_global_values.h"
#include "t8dg_local_values.h"
#include "t8dg_values.h"
#include "t8dg_sc_array.h"
#include "t8dg_mortar.h"
#include "t8dg_dof.h"
#include "t8dg_flux.h"
#include "t8dg_flux_implementation.h"

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

  double              ghost_exchange_data_time;
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
  t8dg_local_values_set_all_elements (return_values->local_values);
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
  int         eclass;

  t8dg_mortar_array_destroy (&values->mortar_array);

  t8dg_local_values_destroy (&values->local_values);

  t8dg_coarse_geometry_destroy (&values->coarse_geometry);

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
t8dg_values_apply_stiffness_matrix_linear_flux_fn3D (t8dg_values_t * values, t8dg_linear_flux3D_fn flux_fn, t8dg_flux_data_base *flux_data, double time,
                                                     t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof)
{
  int                 direction;
  t8dg_dofidx_t       idof;
  t8_locidx_t         ielement, itree, num_trees, num_elems_in_tree;
  t8dg_dof_values_t  *result_dof, *summand;
  t8dg_dof_values_t **flux_dof;
  t8dg_element_dof_values_t *element_flux_dof;
  t8dg_element_dof_values_t *src_element_dof;
  double              flux_vec[3];
  double              reference_vertex[3];
  double              image_vertex[3];
  t8dg_global_values_t *global_values;
  t8dg_functionbasis_t *functionbasis;

  t8_gloidx_t         iglobaltree;
  t8_element_t       *element;

  flux_dof = T8DG_ALLOC (t8dg_dof_values_t *, values->dim);
  for (direction = 0; direction < values->dim; direction++) {
    flux_dof[direction] = t8dg_dof_values_new (values->forest, values->global_values_array);
  }

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);
    iglobaltree = t8_forest_global_tree_id (values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      flux_data->before_first_call_on_element (values->forest, itree, ielement);
      for (direction = 0; direction < values->dim; direction++) {
        element = t8_forest_get_element_in_tree (values->forest, itree, ielement);
        element_flux_dof = t8dg_dof_values_new_element_dof_values_view (flux_dof[direction], itree, ielement);
        src_element_dof = t8dg_dof_values_new_element_dof_values_view (src_dof, itree, ielement);

        global_values = t8dg_values_get_global_values (values, itree, ielement);
        functionbasis = t8dg_global_values_get_functionbasis (global_values);
        for (idof = 0; idof < t8dg_element_dof_values_get_num_dof (element_flux_dof); idof++) {
          t8dg_functionbasis_get_lagrange_vertex (functionbasis, idof, reference_vertex);
          t8dg_geometry_transform_reference_vertex_to_image_vertex (values->coarse_geometry, values->forest, iglobaltree, element,
                                                                    reference_vertex, image_vertex);
          flux_fn (image_vertex, flux_vec, time, flux_data);
          t8dg_element_dof_values_set_value (element_flux_dof, idof,
                                             flux_vec[direction] * t8dg_element_dof_values_get_value (src_element_dof, idof));
        }
        t8dg_element_dof_values_destroy (&element_flux_dof);
        t8dg_element_dof_values_destroy (&src_element_dof);
      }
      flux_data->after_last_call_on_element (values->forest, itree, ielement);
    }
  }

  result_dof = t8dg_dof_values_duplicate (dest_dof);
  summand = t8dg_dof_values_duplicate (dest_dof);
  t8dg_dof_values_set_zero (result_dof);
  for (direction = 0; direction < values->dim; direction++) {
    t8dg_values_apply_component_stiffness_matrix_dof (values, direction, flux_dof[direction], summand);
    t8dg_dof_values_axpy (1, summand, result_dof);
  }
  t8dg_dof_values_copy (result_dof, dest_dof);
  t8dg_dof_values_destroy (&result_dof);
  t8dg_dof_values_destroy (&summand);
  for (direction = 0; direction < values->dim; direction++) {
    t8dg_dof_values_destroy (&flux_dof[direction]);
  }
  T8DG_FREE (flux_dof);
}

void
t8dg_values_apply_component_stiffness_matrix_dof (t8dg_values_t * values, int icomp, t8dg_dof_values_t * flux_dof,
                                                  t8dg_dof_values_t * dest_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;

  t8dg_element_dof_values_t *element_flux_dof_values;
  t8dg_element_dof_values_t *element_dest_dof_values;

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_flux_dof_values = t8dg_dof_values_new_element_dof_values_view (flux_dof, itree, ielement);
      element_dest_dof_values = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);
      t8dg_local_values_apply_element_component_stiffness_matrix_dof (values->local_values, itree, ielement, icomp, element_flux_dof_values,
                                                                      element_dest_dof_values);
      t8dg_element_dof_values_destroy (&element_dest_dof_values);
      t8dg_element_dof_values_destroy (&element_flux_dof_values);
    }
  }
}

/*
void
t8dg_values_apply_directional_stiffness_matrix_linear_flux_fn3D (t8dg_values_t * values, int direction, t8dg_linear_flux3D_fn flux_fn, void *flux_data, double time,
                                                     t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_quadidx_t      num_quad_vertices;

  t8dg_element_quad_values_t         *element_quad_values;
  t8dg_element_dof_values_t         *element_dof_values;
  t8dg_element_dof_values_t         *element_res_dof_values;
  t8dg_global_values_t *global_values;

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      global_values = t8dg_values_get_global_values(values, itree, ielement);
      num_quad_vertices = t8dg_global_values_get_num_elem_quad (global_values);

      element_dof_values = t8dg_dof_values_new_element_dof_values_view (src_dof, itree, ielement);
//      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dof_values));
      element_quad_values = t8dg_element_quad_values_new (num_quad_vertices);
      element_res_dof_values = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);

      t8dg_global_values_transform_element_dof_to_element_quad(global_values, element_dof_values, element_quad_values);

      t8dg_local_values_element_multiply_trafo_quad_weight (values->local_values, itree, ielement, element_quad_values, element_quad_values);
//      t8_debugf ("element_quad_values*qtw\n");
//      t8dg_sc_array_block_double_debug_print (element_quad_values);

        t8dg_local_values_element_multiply_flux_value_direction(values->local_values,itree, ielement, flux_fn, flux_data, direction, time, element_quad_values, element_quad_values);

        t8dg_global_values_transform_element_quad_to_element_dof(global_values, element_quad_values, element_res_dof_values);

//        t8_debugf ("element_dof_derivative_values\n");
//        t8dg_sc_array_block_double_debug_print (element_dof_derivative_values);
        t8dg_global_values_element_apply_derivative_matrix_transpose (global_values, direction,
                                                                      element_res_dof_values, element_res_dof_values);
//        t8_debugf ("element_res_summand_dof_values\n");
//        t8dg_sc_array_block_double_debug_print (element_res_summand_dof_values);
//      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_res_dof_values));

      t8dg_element_dof_values_destroy (&element_dof_values);
      t8dg_element_dof_values_destroy (&element_quad_values);
      t8dg_element_quad_values_destroy (&element_res_dof_values);
    }
  }
}
*/

void
t8dg_values_apply_component_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof,
                                                int icomp, t8dg_numerical_flux1D_fn numerical_flux, void *numerical_flux_data, double time)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_element_dof_values_t *element_dest_dof;

  t8dg_mortar_array_calculate_flux_dof1D (values->mortar_array, src_dof, icomp, numerical_flux, numerical_flux_data, time);

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
t8dg_values_apply_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof,
                                      t8dg_linear_flux3D_fn linear_flux, t8dg_flux_data_base *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux,
                                      void *numerical_flux_data, double time)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_element_dof_values_t *element_dest_dof;

  t8dg_mortar_array_calculate_linear_flux3D (values->mortar_array, src_dof, linear_flux, flux_data, numerical_flux, numerical_flux_data,
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

#if T8_WITH_PETSC
void
t8dg_values_block_precon_apply_component_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof,
                                                             t8dg_dof_values_t * dest_dof, int icomp,
                                                             t8dg_numerical_flux1D_fn numerical_flux, void *numerical_flux_data,
                                                             double time, int selector)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_element_dof_values_t *element_dest_dof;

  t8dg_mortar_array_calculate_flux_dof1D (values->mortar_array, src_dof, icomp, numerical_flux, numerical_flux_data, time);

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dest_dof = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);
      t8dg_mortar_array_block_precon_apply_element_boundary_integral (values->mortar_array, itree, ielement, element_dest_dof, selector);

      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dest_dof));
      sc_array_destroy (element_dest_dof);
    }
  }

  t8dg_mortar_array_invalidate_all (values->mortar_array);
}
#endif

#if T8_WITH_PETSC
void
t8dg_values_block_precon_apply_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof,
                                                   t8dg_linear_flux3D_fn linear_flux, t8dg_flux_data_base *flux_data,
                                                   t8dg_numerical_linear_flux3D_fn numerical_flux, void *numerical_flux_data, double time,
                                                   int selector)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_element_dof_values_t *element_dest_dof;

  /* Calculate all flux values */
  t8dg_mortar_array_calculate_linear_flux3D (values->mortar_array, src_dof, linear_flux, flux_data, numerical_flux, numerical_flux_data,
                                             time);

  /* Iterate over all local trees */
  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);

    /* Iterate over all tree-local elements */
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      /* Get the current dofs of the current element */
      element_dest_dof = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);

      /* Apply/update the element-local degrees of freedom, but only with the block-preconditioner-conforming flux values */
      t8dg_mortar_array_block_precon_apply_element_boundary_integral (values->mortar_array, itree, ielement, element_dest_dof, selector);

      T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dest_dof));
      /* Free the element_dofs */
      sc_array_destroy (element_dest_dof);
    }
  }

  /* Make the mortars accessible for recalculating, after the block preconditioner has been applied */
  t8dg_mortar_array_invalidate_all (values->mortar_array);

}
#endif

#if T8_WITH_PETSC
/* Counts the local number of degrees of freedoms */
t8_locidx_t
t8dg_values_count_num_local_dofs (t8dg_values_t * values)
{
  t8_locidx_t         iter_tree, iter_elem;
  t8_locidx_t         number_local_dofs = 0;

  for (iter_tree = 0; iter_tree < t8_forest_get_num_local_trees (values->forest); ++iter_tree) {
    for (iter_elem = 0; iter_elem < t8_forest_get_tree_num_elements (values->forest, iter_tree); ++iter_elem) {
      number_local_dofs = number_local_dofs + t8dg_global_values_get_num_dof (t8dg_values_get_global_values (values, iter_tree, iter_elem));
    }
  }

  return number_local_dofs;
}
#endif

void
t8dg_values_ghost_exchange (t8dg_values_t * values)
{
  t8dg_local_values_ghost_exchange (values->local_values);
}

double
t8dg_values_get_ghost_exchange_time (t8dg_values_t * values)
{
  return values->ghost_exchange_data_time + t8dg_mortar_array_get_ghost_exchange_time (values->mortar_array);
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

  global_values = t8dg_values_get_global_values (values, itree, ielem_parent_old);

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
  t8dg_local_values_set_all_ghost_elements (values->local_values_adapt);        /*local values get copied or set by replace function */
  values->forest_adapt = forest_adapt;
  t8_forest_ref (values->forest_adapt);
}

#if T8_WITH_PETSC
/* Allocates local_values and mortar_arrays for each multigrid level */
void
t8dg_values_mg_lvl_allocate_properties (t8dg_values_t * values, int num_mg_lvls, t8_forest_t * forests,
                                        t8dg_local_values_t ** local_values_lvl, t8dg_mortar_array_t ** mortar_array_lvl)
{
  int                 lvl_iter;

  /* the initial members at index zero are already initialized */
  for (lvl_iter = 1; lvl_iter < num_mg_lvls; ++lvl_iter) {
    local_values_lvl[lvl_iter] = t8dg_local_values_new (forests[lvl_iter], values->global_values_array, values->coarse_geometry);
    t8dg_local_values_set_all_ghost_elements (local_values_lvl[lvl_iter]);
    mortar_array_lvl[lvl_iter] = t8dg_mortar_array_new_empty (forests[lvl_iter], local_values_lvl[lvl_iter]);
  }
}

/* Assign members of dg_values in order to perform the interpolation step during the multigrid preconditioning */
void
t8dg_values_mg_lvl_set_interpolation_step (t8dg_values_t * values, t8_forest_t forest, t8_forest_t forest_adapt,
                                           t8dg_local_values_t * local_values, t8dg_local_values_t * local_values_adapt)
{
  values->forest = forest;
  values->forest_adapt = forest_adapt;
  values->local_values = local_values;
  values->local_values_adapt = local_values_adapt;
}

/* Switch the dg_Values' members so that they are conforming with the current underlying mesh */
void
t8dg_values_mg_lvl_prepare_next_interpolation_step (t8dg_values_t * values, t8dg_mortar_array_t * mortar_array_lvl)
{
  values->local_values = values->local_values_adapt;
  /* maybe local_values set all ghosts has to be called */
  values->forest = values->forest_adapt;
  values->mortar_array = mortar_array_lvl;
}

t8dg_local_values_t **
t8dg_values_get_local_values (t8dg_values_t * values)
{
  return &(values->local_values);
}

t8dg_mortar_array_t **
t8dg_values_get_mortar_array (t8dg_values_t * values)
{
  return &(values->mortar_array);
}

t8dg_coarse_geometry_t *
t8dg_values_get_coarse_geometry (t8dg_values_t * values)
{
  return values->coarse_geometry;
}

#endif

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
  t8dg_local_values_set_all_ghost_elements (values->local_values);      /*local elements get set during replace, these could also be communicated */
  /*Create new mortar arrays */
  values->ghost_exchange_data_time += t8dg_mortar_array_get_ghost_exchange_time (values->mortar_array);
  t8dg_mortar_array_destroy (&values->mortar_array);
  values->mortar_array = t8dg_mortar_array_new_empty (values->forest, values->local_values);
}

void
t8dg_values_partition (t8dg_values_t * values, t8_forest_t forest_partition)
{
  /* Partition local precomputed values */
  t8dg_local_values_t *local_values_partition;
  local_values_partition = t8dg_local_values_new (forest_partition, values->global_values_array, values->coarse_geometry);
  t8dg_local_values_partition (values->local_values, local_values_partition);   /*partitions the local values */

  t8dg_local_values_destroy (&values->local_values);
  values->local_values = local_values_partition;

  t8dg_local_values_set_all_ghost_elements (values->local_values);      /*recalculates the ghost values */

  t8_debugf ("[VALUES] Done partition local_data\n");

  t8_debugf (" [ADVECT] begin mortars destroy\n");
  values->ghost_exchange_data_time += t8dg_mortar_array_get_ghost_exchange_time (values->mortar_array); /*Add ghost exchange time from mortars, before they are destroyed */
  t8dg_mortar_array_destroy (&values->mortar_array);

  t8_forest_unref (&values->forest);
  values->forest = forest_partition;
  t8_forest_ref (values->forest);

  values->mortar_array = t8dg_mortar_array_new_empty (forest_partition, values->local_values);
  t8_debugf (" [ADVECT] Done partition\n");
}

typedef struct t8dg_values_fn_evaluation_data
{
  t8dg_coarse_geometry_t *coarse_geometry;
  t8_forest_t         forest;
  t8_locidx_t         iglobaltree;
  t8_element_t       *element;
  t8dg_scalar_function_3d_time_fn function;
  void               *function_data;
  double              time;
} t8dg_values_fn_evaluation_data_t;

static double
t8dg_values_transform_reference_vertex_and_evaluate (const double reference_vertex[3], void *scalar_fn_data)
{
  double              image_vertex[DIM3];
  t8dg_values_fn_evaluation_data_t *data;

  data = (t8dg_values_fn_evaluation_data_t *) scalar_fn_data;

  t8dg_geometry_transform_reference_vertex_to_image_vertex (data->coarse_geometry, data->forest, data->iglobaltree, data->element,
                                                            reference_vertex, image_vertex);

  /* apply initial condition function at image vertex and start time */
  return data->function (image_vertex, data->time, data->function_data);
}

void
t8dg_values_interpolate_scalar_function_3d_time (t8dg_values_t * values, t8dg_scalar_function_3d_time_fn function, double time,
                                                 void *function_data, t8dg_dof_values_t * dof_values)
{
  t8dg_values_fn_evaluation_data_t data;
  t8_locidx_t         num_trees, itree, ielement, num_elems_in_tree;

  data.coarse_geometry = values->coarse_geometry;
  data.forest = values->forest;
  data.function = function;
  data.function_data = function_data;
  data.time = time;

  t8dg_global_values_t *global_values;
  t8dg_functionbasis_t *functionbasis;
  t8dg_element_dof_values_t *element_dof_view;

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (values->forest, itree);
    data.iglobaltree = t8_forest_global_tree_id (values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      element_dof_view = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielement);
      global_values = t8dg_values_get_global_values (values, itree, ielement);
      functionbasis = t8dg_global_values_get_functionbasis (global_values);
      data.element = t8_forest_get_element_in_tree (values->forest, itree, ielement);
      t8dg_functionbasis_interpolate_scalar_fn (functionbasis, t8dg_values_transform_reference_vertex_and_evaluate, &data,
                                                element_dof_view);
      t8dg_element_dof_values_destroy (&element_dof_view);
    }
  }
}

double
t8dg_values_norm_l2 (t8dg_values_t * values, t8dg_dof_values_t * dof_values, sc_MPI_Comm comm)
{
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;

  double              norm_squared = 0, global_norm_squared;
  t8dg_element_dof_values_t *element_dof_values;

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    num_elements = t8_forest_get_tree_num_elements (values->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      element_dof_values = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielement);
      norm_squared += t8dg_values_element_norm_l2_squared (values, element_dof_values, itree, ielement);
      t8dg_element_dof_values_destroy (&element_dof_values);
    }
  }
  /* Compute the sum of the error among all processes */
  sc_MPI_Allreduce (&norm_squared, &global_norm_squared, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
  return sqrt (global_norm_squared);
}

double
t8dg_values_norm_l2_rel (t8dg_values_t * values, t8dg_dof_values_t * dof_values, t8dg_scalar_function_3d_time_fn analytical_sol_fn,
                         double time, void *analytical_sol_data, sc_MPI_Comm comm)
{
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;

  t8dg_dof_values_t  *analytical_sol_dof;

  double              error_squared = 0, global_error_squared, global_error, error_squared_summand;
  double              ana_norm_squared = 0, global_ana_norm_squared, global_ana_norm, ana_norm_squared_summand;

  analytical_sol_dof = t8dg_dof_values_duplicate (dof_values);

  t8dg_values_interpolate_scalar_function_3d_time (values, analytical_sol_fn, time, analytical_sol_data, analytical_sol_dof);

  num_trees = t8_forest_get_num_local_trees (values->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    num_elements = t8_forest_get_tree_num_elements (values->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      t8dg_local_values_element_error_ana_l2_squared (values->local_values, dof_values, analytical_sol_dof, itree, ielement, time,
                                                      &error_squared_summand, &ana_norm_squared_summand);
      error_squared += error_squared_summand;
      ana_norm_squared += ana_norm_squared_summand;
    }
  }
  /* Compute the sum of the error among all processes */
  sc_MPI_Allreduce (&ana_norm_squared, &global_ana_norm_squared, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
  sc_MPI_Allreduce (&error_squared, &global_error_squared, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);

  global_ana_norm = sqrt (global_ana_norm_squared);
  global_error = sqrt (global_error_squared);

  t8dg_dof_values_destroy (&analytical_sol_dof);

  /* Return the relative error, that is the l_2 error divided by
   * the l_2 norm of the analytical solution */
  return global_error / global_ana_norm;
}

double
t8dg_values_element_norm_l2_squared (t8dg_values_t * values, t8dg_element_dof_values_t * element_dof, t8_locidx_t itree,
                                     t8_locidx_t ielement)
{
  return t8dg_local_values_element_norm_l2_squared (values->local_values, element_dof, itree, ielement);
}

double
t8dg_values_element_area (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement)
{
  double              area;
  t8dg_element_dof_values_t *one_array;
  t8dg_dofidx_t       num_dof = t8dg_global_values_get_num_dof (t8dg_values_get_global_values (values, itree, ielement));
  one_array = t8dg_element_dof_values_new (num_dof);
  t8dg_element_dof_values_set_all (one_array, 1);
  area = t8dg_values_element_norm_l2_squared (values, one_array, itree, ielement);
  t8dg_element_dof_values_destroy (&one_array);
  return area;
}

double
t8dg_values_norm_l_infty_rel (t8dg_values_t * values, t8dg_dof_values_t * dof_values,
                              t8dg_scalar_function_3d_time_fn analytical_sol_fn, double time, void *analytical_sol_data, sc_MPI_Comm comm)
{
#if 0
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;
  sc_array_t         *elem_ana_sol;
  sc_array_t         *elem_dof_val;
  sc_array_t         *elem_error;
  double              error = 0, global_error;
  double              ana_norm = 0, global_ana_norm;

  t8dg_geometry_t     geometry = { problem->coarse_geometry, problem->forest };
  t8_eclass_t         eclass;

  elem_ana_sol = sc_array_new_count (sizeof (double), t8dg_global_values_get_num_dof (problem->global_values));
  elem_error = t8dg_sc_array_duplicate (elem_ana_sol);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    num_elements = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      elem_dof_val = t8dg_dof_values_new_element_dof_values_view (problem->dof_values, idata);
      t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_values_get_functionbasis (problem->global_values),
                                                t8dg_values_transform_reference_vertex_and_evaluate, &evaluation_data, elem_ana_sol);
      t8dg_dof_values_axpyz (-1, elem_ana_sol, elem_dof_val, elem_error);
      error = SC_MAX (error, t8dg_values_element_norm_infty (elem_error));
      ana_norm = SC_MAX (ana_norm, t8dg_values_element_norm_infty (elem_ana_sol));
      sc_array_destroy (elem_dof_val);
    }
  }
  /* Compute the maximum of the error among all processes */
  sc_MPI_Allreduce (&ana_norm, &global_ana_norm, 1, sc_MPI_DOUBLE, sc_MPI_MAX, problem->comm);
  sc_MPI_Allreduce (&error, &global_error, 1, sc_MPI_DOUBLE, sc_MPI_MAX, problem->comm);

  sc_array_destroy (elem_ana_sol);
  sc_array_destroy (elem_error);

  /* Return the relative error, that is the l_infty error divided by
   * the l_infty norm of the analytical solution */
  return global_error / global_ana_norm;
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_ERROR_INF, 0);
#endif
  return -1;
}

t8_forest_t
t8dg_values_get_forest (t8dg_values_t * values)
{
  return values->forest;
}

int
t8dg_values_get_dim (t8dg_values_t * values)
{
  return values->dim;
}
