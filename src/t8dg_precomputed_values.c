#include "t8dg.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_local_precomputed_values.h"
#include "t8dg_precomputed_values.h"
#include "t8dg_sc_array.h"
#include "t8dg_mortar.h"

void
t8dg_precomputed_values_apply_element_boundary_integral (t8dg_global_precomputed_values_t * global_values,
                                                         t8dg_local_precomputed_values_t * local_values, t8dg_mortar_array_t * mortar_array,
                                                         t8_locidx_t idata, sc_array_t * element_result_dof)
{
  int                 iface, num_faces;

  sc_array_t         *face_flux_dof;
  sc_array_t         *face_dof;
  sc_array_t         *summand;

  num_faces = t8dg_global_precomputed_values_get_num_faces (global_values);
  summand = t8dg_sc_array_duplicate (element_result_dof);
  t8dg_sc_array_block_double_set_zero (element_result_dof);

  for (iface = 0; iface < num_faces; iface++) {
    face_flux_dof = t8dg_mortar_array_get_oriented_flux (mortar_array, idata, iface);
    face_dof = t8dg_sc_array_duplicate (face_flux_dof);

    t8dg_precomputed_values_apply_face_mass_matrix (global_values, local_values, idata, iface, face_flux_dof, face_dof);

    t8dg_functionbasis_transform_face_dof_to_element_dof (t8dg_global_precomputed_values_get_functionbasis (global_values), iface, face_dof,
                                                          summand);
    T8DG_ASSERT (t8dg_sc_array_block_double_is_valid (summand));
    t8dg_sc_array_block_double_axpy (1, summand, element_result_dof);
    sc_array_destroy (face_dof);
  }
  sc_array_destroy (summand);
}

void
t8dg_precomputed_values_apply_face_inverse_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                        const t8dg_local_precomputed_values_t * local_values,
                                                        const t8dg_locidx_t idata, const int iface, sc_array_t * dof_values,
                                                        sc_array_t * result_dof_values)
{
  if (t8dg_global_precomputed_values_simplifies (global_values)) {
    t8dg_local_precomputed_values_face_divide_trafo_quad_weight (local_values, idata, iface, dof_values, result_dof_values);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_precomputed_values_apply_face_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                const t8dg_local_precomputed_values_t * local_values,
                                                const t8dg_locidx_t idata, const int iface, sc_array_t * dof_values,
                                                sc_array_t * result_dof_values)
{
  if (t8dg_global_precomputed_values_simplifies (global_values)) {
    t8dg_local_precomputed_values_face_multiply_trafo_quad_weight (local_values, idata, iface, dof_values, result_dof_values);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_precomputed_values_apply_element_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                   const t8dg_local_precomputed_values_t * local_values,
                                                   const t8dg_locidx_t idata, sc_array_t * dof_values, sc_array_t * result_dof_values)
{
  if (t8dg_global_precomputed_values_simplifies (global_values)) {
    t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (local_values, idata, dof_values, result_dof_values);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

void
t8dg_precomputed_values_apply_element_inverse_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                           const t8dg_local_precomputed_values_t * local_values,
                                                           const t8dg_locidx_t idata, sc_array_t * dof_values,
                                                           sc_array_t * result_dof_values)
{
  if (t8dg_global_precomputed_values_simplifies (global_values)) {
    t8dg_local_precomputed_values_element_divide_trafo_quad_weight (local_values, idata, dof_values, result_dof_values);
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

/*TODO: Add face transform child dof to parent dof*/
void
t8dg_precomputed_values_transform_face_child_dof_to_parent_dof (const t8dg_global_precomputed_values_t * global_values,
                                                                sc_array_t * child_face_dof[MAX_SUBFACES],
                                                                sc_array_t * parent_face_dof, const int num_face_children,
                                                                const t8dg_local_precomputed_values_t * local_values,
                                                                t8_locidx_t idata_child[MAX_SUBFACES], t8_locidx_t idata_parent,
                                                                int iface_parent)
{
  sc_array_t         *summand;
  sc_array_t         *mass_times_child_dof;
  int                 ichild;
  t8dg_functionbasis_t *face_functionbasis;

  summand = t8dg_sc_array_duplicate (parent_face_dof);
  mass_times_child_dof = t8dg_sc_array_duplicate (child_face_dof[0]);

  t8dg_sc_array_block_double_set_zero (parent_face_dof);

  face_functionbasis =
    t8dg_functionbasis_get_face_functionbasis (t8dg_global_precomputed_values_get_functionbasis (global_values), iface_parent);
  for (ichild = 0; ichild < num_face_children; ichild++) {
    t8dg_precomputed_values_apply_face_mass_matrix (global_values, local_values, idata_child[ichild], iface_parent,
                                                    child_face_dof[ichild], mass_times_child_dof);

    t8dg_functionbasis_apply_child_interpolation_matrix_transpose (face_functionbasis, ichild, mass_times_child_dof, summand);

    t8dg_sc_array_block_double_axpy (1, summand, parent_face_dof);
  }
  t8dg_precomputed_values_apply_face_inverse_mass_matrix (global_values, local_values, idata_parent, iface_parent, parent_face_dof,
                                                          parent_face_dof);
  sc_array_destroy (mass_times_child_dof);
  sc_array_destroy (summand);
}

void
t8dg_precomputed_values_transform_child_dof_to_parent_dof (const t8dg_global_precomputed_values_t * global_values,
                                                           sc_array_t * child_dof[MAX_SUBELEMENTS],
                                                           sc_array_t * parent_dof, const int num_children,
                                                           const t8dg_local_precomputed_values_t * local_values_old,
                                                           const t8dg_local_precomputed_values_t * local_values_new,
                                                           t8_locidx_t idata_first_child, t8_locidx_t idata_parent)
{
  sc_array_t         *summand;
  sc_array_t         *mass_times_child_dof;
  int                 ichild;

  summand = t8dg_sc_array_duplicate (parent_dof);
  mass_times_child_dof = t8dg_sc_array_duplicate (child_dof[0]);

  t8dg_sc_array_block_double_set_zero (parent_dof);

  for (ichild = 0; ichild < num_children; ichild++) {
    t8dg_precomputed_values_apply_element_mass_matrix (global_values, local_values_old, idata_first_child + ichild,
                                                       child_dof[ichild], mass_times_child_dof);

    t8dg_functionbasis_apply_child_interpolation_matrix_transpose (t8dg_global_precomputed_values_get_functionbasis (global_values), ichild,
                                                                   mass_times_child_dof, summand);

    t8dg_sc_array_block_double_axpy (1, summand, parent_dof);
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
  double              image_vertex[DIM3];
  t8dg_precomputed_values_fn_evaluation_data_t *data;

  data = (t8dg_precomputed_values_fn_evaluation_data_t *) scalar_fn_data;

  t8dg_geometry_transform_reference_vertex_to_image_vertex (data->geometry_data, reference_vertex, image_vertex);

  /* apply initial condition function at image vertex and start time */
  return data->function (image_vertex, data->time);

}
