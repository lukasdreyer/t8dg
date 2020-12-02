#include "t8dg.h"
#include "t8dg_quad.h"
#include "t8dg_global_values.h"
#include "t8_forest/t8_forest_partition.h"

struct t8dg_quad_values
{
  sc_array_t         *quads;

  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;
  t8_locidx_t         num_total_elements;

  t8dg_quadidx_t      max_num_element_quad;

  t8_forest_t         forest;   /* to convert idata into eclass */
  t8dg_global_values_t **global_values; /*to get the appropriate amount of dofs for each elementtype */
};

t8dg_element_quad_values_t *
t8dg_quad_values_new_element_quad_values_view (t8dg_quad_values_t * quad_values, t8_locidx_t itree, t8_locidx_t ielement)
{
  t8dg_dofidx_t       idata;
  idata = t8dg_itree_ielement_to_idata (quad_values->forest, itree, ielement);
  T8DG_ASSERT (idata >= 0 && (size_t) idata < quad_values->num_total_elements);

  t8dg_element_quad_values_t *element_quad_view;

  t8dg_global_values_t *global_values =
    t8dg_global_values_array_get_global_values (quad_values->global_values, quad_values->forest, itree, ielement);

  t8dg_dofidx_t       num_element_quad = t8dg_global_values_get_num_elem_quad (global_values);
  element_quad_view =
    (t8dg_element_quad_values_t *) sc_array_new_data (t8_sc_array_index_locidx (quad_values->quads, idata), sizeof (double),
                                                      num_element_quad);
  /* Create new view */
  return element_quad_view;
}

t8dg_face_quad_values_t *
t8dg_quad_values_new_face_quad_values_view (t8dg_quad_values_t * quad_values, int iface, t8_locidx_t itree, t8_locidx_t ielement)
{
  t8dg_dofidx_t       idata;
  idata = t8dg_itree_ielement_to_idata (quad_values->forest, itree, ielement);
  T8DG_ASSERT (idata >= 0 && (size_t) idata < quad_values->num_total_elements);

  t8dg_face_quad_values_t *face_quad_view;

  t8dg_global_values_t *global_values =
    t8dg_global_values_array_get_global_values (quad_values->global_values, quad_values->forest, itree, ielement);

  t8dg_dofidx_t       num_face_quad = t8dg_global_values_get_num_face_quad (global_values, iface);
  face_quad_view =
    (t8dg_element_quad_values_t *) sc_array_new_data (t8_sc_array_index_locidx (quad_values->quads, idata), sizeof (double), num_face_quad);
  /* Create new view */
  return face_quad_view;
}

void
t8dg_element_quad_values_destroy (t8dg_element_quad_values_t ** p_quad_values)
{
  sc_array_destroy_null (p_quad_values);
}

void
t8dg_face_quad_values_destroy (t8dg_face_quad_values_t ** p_quad_values)
{
  sc_array_destroy_null (p_quad_values);
}

void
t8dg_quad_values_copy_from_index_to_index (t8dg_quad_values_t * src_quad, t8_locidx_t src_idata, t8dg_quad_values_t * dest_quad,
                                           t8_locidx_t dest_idata)
{
  T8DG_ABORT ("Not implemented \n");
}

double
t8dg_element_quad_values_get_value (t8dg_element_quad_values_t * element_quad_values, t8dg_quadidx_t iquad)
{
  T8DG_ASSERT (t8dg_element_quad_values_is_valid (element_quad_values));
  T8DG_ASSERT (iquad >= 0 && iquad < element_quad_values->elem_count);
  return *(double *) sc_array_index (element_quad_values, iquad);
}

void
t8dg_element_quad_values_set_value (t8dg_element_quad_values_t * element_quad_values, t8dg_quadidx_t iquad, double value)
{
  T8DG_ASSERT (t8dg_element_quad_values_is_valid (element_quad_values));
  T8DG_ASSERT (iquad >= 0 && iquad < element_quad_values->elem_count);
  *(double *) sc_array_index (element_quad_values, iquad) = value;
}

double
t8dg_face_quad_values_get_value (t8dg_face_quad_values_t * face_quad_values, t8dg_quadidx_t iquad)
{
  return *(double *) sc_array_index (face_quad_values, iquad);
}

void
t8dg_face_quad_values_set_value (t8dg_face_quad_values_t * face_quad_values, t8dg_quadidx_t iquad, double value)
{
  *(double *) sc_array_index (face_quad_values, iquad) = value;
}

t8dg_quadidx_t
t8dg_element_quad_values_get_num_element_quad_points (t8dg_element_quad_values_t * element_quad_values)
{
  T8DG_ASSERT (t8dg_element_quad_values_is_valid (element_quad_values));
  return (t8dg_quadidx_t) element_quad_values->elem_count;
}

t8dg_quadidx_t
t8dg_face_quad_values_get_num_face_quad_points (t8dg_face_quad_values_t * face_quad_values)
{
  T8DG_ASSERT (face_quad_values->elem_size == sizeof (double));
  return face_quad_values->elem_count;
}

t8dg_quad_values_t *
t8dg_quad_values_new (t8_forest_t forest, t8dg_global_values_t * global_values_array[T8_ECLASS_COUNT])
{
  t8dg_quad_values_t *return_values;
  return_values = T8DG_ALLOC_ZERO (t8dg_quad_values_t, 1);
  return_values->forest = forest;
  return_values->global_values = global_values_array;
  return_values->max_num_element_quad = t8dg_global_values_array_get_max_num_element_quad (global_values_array);
  return_values->num_local_elements = t8_forest_get_num_element (forest);
  return_values->num_ghost_elements = t8_forest_get_num_ghosts (forest);
  return_values->num_total_elements = return_values->num_local_elements + return_values->num_ghost_elements;

  return_values->quads = sc_array_new_count (sizeof (double) * return_values->max_num_element_quad, return_values->num_total_elements);
  return return_values;
}

t8dg_quad_values_t *
t8dg_quad_values_new_local (t8_forest_t forest, t8dg_global_values_t * global_values_array[T8_ECLASS_COUNT])
{
  t8dg_quad_values_t *return_values;
  return_values = T8DG_ALLOC_ZERO (t8dg_quad_values_t, 1);
  return_values->forest = forest;
  return_values->global_values = global_values_array;
  return_values->max_num_element_quad = t8dg_global_values_array_get_max_num_element_quad (global_values_array);
  return_values->num_local_elements = t8_forest_get_num_element (forest);
  return_values->num_ghost_elements = 0;
  return_values->num_total_elements = return_values->num_local_elements;

  return_values->quads = sc_array_new_count (sizeof (double) * return_values->max_num_element_quad, return_values->num_total_elements);
  return return_values;
}

void
t8dg_quad_values_destroy (t8dg_quad_values_t ** pquad_values)
{
  t8dg_quad_values_t *quad_values = *pquad_values;
  sc_array_destroy (quad_values->quads);
  T8DG_FREE (quad_values);
  pquad_values = NULL;
}

void
t8dg_quad_values_partition (t8dg_quad_values_t * quad_values_old, t8dg_quad_values_t * quad_values_partition)
{
  sc_array_t         *quad_values_local_view;
  sc_array_t         *quad_values_partition_local_view;
  quad_values_local_view = sc_array_new_view (quad_values_old->quads, 0, quad_values_old->num_local_elements);
  quad_values_partition_local_view = sc_array_new_view (quad_values_partition->quads, 0, quad_values_partition->num_local_elements);

  t8_forest_partition_data (quad_values_old->forest, quad_values_partition->forest, quad_values_local_view,
                            quad_values_partition_local_view);

  /*destroy views */
  sc_array_destroy (quad_values_local_view);
  sc_array_destroy (quad_values_partition_local_view);
}

int
t8dg_quad_values_is_valid (t8dg_quad_values_t * values)
{
  return t8_forest_is_committed (values->forest);
}

int
t8dg_element_quad_values_is_valid (t8dg_element_quad_values_t * element_values)
{
  return (element_values->elem_size == sizeof (double));

}

int
t8dg_face_quad_values_is_valid (t8dg_element_quad_values_t * face_values)
{
  return (face_values->elem_size == sizeof (double));
}
