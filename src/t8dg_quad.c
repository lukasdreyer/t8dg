#include "t8dg.h"
#include "t8dg_quad.h"
#include "t8dg_global_values.h"

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
  T8DG_ABORT ("Not implemented \n");
}

t8dg_face_quad_values_t *
t8dg_quad_values_new_face_quad_values_view (t8dg_quad_values_t * quad_values, t8_locidx_t itree, t8_locidx_t ielement)
{
  T8DG_ABORT ("Not implemented \n");
}

void
t8dg_element_quad_values_destroy (t8dg_element_quad_values_t ** quad_values)
{
  T8DG_ABORT ("Not implemented \n");
}

void
t8dg_face_quad_values_destroy (t8dg_face_quad_values_t ** quad_values)
{
  T8DG_ABORT ("Not implemented \n");
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
  T8DG_ABORT ("Not implemented \n");
}

void
t8dg_element_quad_values_set_value (t8dg_element_quad_values_t * element_quad_values, t8dg_quadidx_t iquad, double value)
{
  T8DG_ABORT ("Not implemented \n");
}

double
t8dg_face_quad_values_get_value (t8dg_face_quad_values_t * face_quad_values, t8dg_quadidx_t iquad)
{
  T8DG_ABORT ("Not implemented \n");
}

void
t8dg_face_quad_values_set_value (t8dg_face_quad_values_t * face_quad_values, t8dg_quadidx_t iquad, double value)
{
  T8DG_ABORT ("Not implemented \n");
}

t8dg_quadidx_t
t8dg_element_quad_values_get_num_element_quad_points (t8dg_element_quad_values_t * element_quad_values)
{
  T8DG_ABORT ("Not implemented \n");
}

t8dg_quadidx_t
t8dg_face_quad_values_get_num_face_quad_points (t8dg_face_quad_values_t * face_quad_values)
{
  T8DG_ABORT ("Not implemented \n");
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
  T8DG_ABORT ("Not implemented \n");
}
