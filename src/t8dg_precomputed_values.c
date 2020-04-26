#include "t8dg.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_local_precomputed_values.h"
#include "t8dg_precomputed_values.h"
#include "t8dg_sc_array.h"

void
t8dg_precomputed_values_apply_element_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                   const t8dg_local_precomputed_values_t * local_values,
                                                   const t8dg_locidx_t idata, const sc_array_t * dof_values, sc_array_t * result_dof_values)
{
  sc_array_t         *quad_values;

  quad_values = sc_array_new_count (sizeof (double), t8dg_global_precomputed_values_get_num_elem_quad (global_values));
  t8dg_global_precomputed_values_transform_element_dof_to_element_quad (global_values, dof_values, quad_values);
  t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (local_values, quad_values, idata);
  t8dg_global_precomputed_values_transform_element_quad_to_element_dof (global_values, quad_values, result_dof_values);
  sc_array_destroy (quad_values);
}

void
t8dg_precomputed_values_apply_element_inverse_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                           const t8dg_local_precomputed_values_t * local_values,
                                                           const t8dg_locidx_t idata, const sc_array_t * dof_values,
                                                           sc_array_t * result_dof_values)
{
  T8DG_CHECK_ABORT (t8dg_functionbasis_get_type (t8dg_global_precomputed_values_get_functionbasis (global_values)) == T8DG_LAGRANGE_LGL
                    && t8dg_quadrature_get_type (t8dg_global_precomputed_values_get_quadrature (global_values)) == T8DG_QUAD_LGL,
                    "Not yet implemented");

  T8DG_ASSERT (dof_values->elem_count == result_dof_values->elem_count);
  T8DG_ASSERT (dof_values->elem_count == (size_t) t8dg_global_precomputed_values_get_num_dof (global_values));

  t8dg_sc_array_copy (dof_values, result_dof_values);
  t8dg_local_precomputed_values_element_divide_trafo_quad_weight (local_values, result_dof_values, idata);
}

void
t8dg_precomputed_values_transform_child_dof_to_parent_dof (const t8dg_global_precomputed_values_t * global_values,
                                                           const sc_array_t * const child_dof[MAX_SUBELEMENTS],
                                                           sc_array_t * parent_dof, const int num_children,
                                                           t8dg_local_precomputed_values_t * local_values_old,
                                                           t8dg_local_precomputed_values_t * local_values_new,
                                                           t8_locidx_t idata_first_child, t8_locidx_t idata_parent)
{
  t8dg_dmatrix_t     *interpolation_matrix;
  sc_array_t         *summand;
  sc_array_t         *mass_times_child_dof;
  int                 ichild;

  summand = t8dg_sc_array_duplicate (parent_dof);
  mass_times_child_dof = t8dg_sc_array_duplicate (child_dof[0]);

  t8dg_sc_array_block_double_set_zero (parent_dof);

  for (ichild = 0; ichild < num_children; ichild++) {
    interpolation_matrix =
      t8dg_functionbasis_get_child_interpolation_matrix (t8dg_global_precomputed_values_get_functionbasis (global_values), ichild);

    t8dg_precomputed_values_apply_element_mass_matrix (global_values, local_values_old, idata_first_child + ichild,
                                                       child_dof[ichild], mass_times_child_dof);

    t8dg_dmatrix_transpose_mult_sc_array (interpolation_matrix, mass_times_child_dof, summand);
    t8dg_sc_array_block_double_axpy (1, summand, parent_dof);   /*incorporate mass matrices !!! */
  }
  t8dg_precomputed_values_apply_element_inverse_mass_matrix (global_values, local_values_new, idata_parent, parent_dof, parent_dof);
  sc_array_destroy (mass_times_child_dof);
  sc_array_destroy (summand);
}

double
t8dg_precomputed_values_element_norm_infty (sc_array_t * element_dof_values)
{
  double              norm = 0;
  size_t              idof;
  for (idof = 0; idof < element_dof_values->elem_count; idof++) {
    norm = SC_MAX (norm, fabs (*(double *) sc_array_index (element_dof_values, idof)));
  }
  return norm;
}

double
t8dg_precomputed_values_element_norm_l2_squared (sc_array_t * element_dof_values, t8dg_global_precomputed_values_t * global_values,
                                                 t8dg_local_precomputed_values_t * local_values, t8_locidx_t idata)
{
  T8DG_ASSERT (element_dof_values->elem_count == (size_t) t8dg_global_precomputed_values_get_num_dof (global_values));
  double              norm = 0;
  int                 idof, num_dof;
  sc_array_t         *mass_times_square_dof;
  sc_array_t         *element_dof_square_values;

  num_dof = t8dg_global_precomputed_values_get_num_dof (global_values);

  element_dof_square_values = t8dg_sc_array_duplicate (element_dof_values);
  mass_times_square_dof = t8dg_sc_array_duplicate (element_dof_values);

  t8dg_sc_array_block_square_values (element_dof_values, element_dof_square_values);
  t8dg_precomputed_values_apply_element_mass_matrix (global_values, local_values, idata, element_dof_square_values, mass_times_square_dof);

  for (idof = 0; idof < num_dof; idof++) {
    norm += *(double *) sc_array_index_int (mass_times_square_dof, idof);
  }
  sc_array_destroy (element_dof_square_values);
  sc_array_destroy (mass_times_square_dof);
  return norm;
}

double
t8dg_precomputed_values_transform_reference_vertex_and_evaluate (const double reference_vertex[3], void *scalar_fn_data)
{
  double              coarse_vertex[DIM3];
  double              image_vertex[DIM3];
  t8dg_precomputed_values_fn_evaluation_data_t *data;

  data = (t8dg_precomputed_values_fn_evaluation_data_t *) scalar_fn_data;
  /* transform into coarse reference element */
  t8dg_local_precomputed_values_fine_to_coarse_geometry (reference_vertex, coarse_vertex, data->scheme, data->element);

  /* tree vertices are application data for linear geometry */
  data->coarse_geometry->geometry (coarse_vertex, image_vertex, data->forest, data->itree);

  /* apply initial condition function at image vertex and start time */
  return data->function (image_vertex, data->time);

}
