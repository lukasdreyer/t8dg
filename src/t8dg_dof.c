#include "t8dg.h"
#include "t8dg_dof.h"
#include "t8dg_global_values.h"
#include <t8_forest/t8_forest_partition.h>

struct t8dg_dof_values
{
  sc_array_t         *dofs;

  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;
  t8_locidx_t         num_total_elements;

  t8dg_dofidx_t       max_num_element_dof;

  t8_forest_t         forest;   /* to convert idata into eclass */
  t8dg_global_values_t **global_values; /*to get the appropriate amount of dofs for each elementtype */
};

/* Allocates memory for t8dg_element_dof_values_t and sc_array_view, only */
t8dg_element_dof_values_t *
t8dg_dof_values_new_element_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t idata)
{
  T8DG_ASSERT (idata >= 0 && (size_t) idata < dof_values->num_total_elements);

  t8dg_element_dof_values_t *element_dof_view;

  element_dof_view = T8DG_ALLOC (t8dg_element_dof_values_t, 1);

  t8dg_dofidx_t       num_element_dof = 0;      /*TODO!!!! */
  element_dof_view =
    (t8dg_element_dof_values_t *) sc_array_new_data (t8_sc_array_index_locidx (dof_values->dofs, idata), sizeof (double), num_element_dof);
  /* Create new view */
  return element_dof_view;
}

void
t8dg_dof_values_destroy (t8dg_dof_values_t ** dof_values)
{
  T8DG_ABORT ("Not implemented!\n");
}

void
t8dg_dof_values_ghost_exchange (t8dg_dof_values_t * dof_values)
{
  t8_forest_ghost_exchange_data (dof_values->forest, dof_values->dofs);
}

double             *
t8dg_dof_values_get_double_pointer (const t8dg_dof_values_t * dof_values, t8_locidx_t idata)
{
  return ((double *) t8_sc_array_index_locidx (dof_values->dofs, idata));
}

t8dg_dof_values_t  *
t8dg_dof_values_new (t8_forest_t forest, t8dg_global_values_t ** global_values_array)
{
  t8dg_dof_values_t  *return_values;
  return_values = T8DG_ALLOC_ZERO (t8dg_dof_values_t, 1);
  return_values->forest = forest;
  return_values->global_values = global_values_array;
  return_values->max_num_element_dof = t8dg_global_values_array_get_max_num_element_dof (global_values_array);
  return_values->num_local_elements = t8_forest_get_num_element (forest);
  return_values->num_ghost_elements = t8_forest_get_num_ghosts (forest);
  return_values->num_total_elements = return_values->num_local_elements + return_values->num_ghost_elements;

  return_values->dofs = sc_array_new_count (sizeof (double) * return_values->max_num_element_dof, return_values->num_total_elements);
  return return_values;
}

void
t8dg_dof_values_partition (t8dg_dof_values_t * dof_values_old, t8dg_dof_values_t * dof_values_partition)
{
  sc_array_t         *dof_values_local_view;
  sc_array_t         *dof_values_partition_local_view;
  dof_values_local_view = sc_array_new_view (dof_values_old->dofs, 0, dof_values_old->num_local_elements);
  dof_values_partition_local_view = sc_array_new_view (dof_values_partition->dofs, 0, dof_values_partition->num_local_elements);

  t8_forest_partition_data (dof_values_old->forest, dof_values_partition->forest, dof_values_local_view, dof_values_partition_local_view);

  /*destroy views */
  sc_array_destroy (dof_values_local_view);
  sc_array_destroy (dof_values_partition_local_view);
}

void
t8dg_element_dof_values_set_zero (t8dg_element_dof_values_t * dof_values)
{
  T8DG_ABORT ("Not implemented!\n");
}

t8dg_face_dof_values_t *
t8dg_dof_values_new_face_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t idata)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_element_dof_values_destroy (t8dg_element_dof_values_t ** dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_face_dof_values_destroy (t8dg_face_dof_values_t ** dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_copy_from_index_to_index (t8dg_dof_values_t * src_dof, t8_locidx_t src_idata, t8dg_dof_values_t * dest_dof,
                                          t8_locidx_t dest_idata)
{
  T8DG_ABORT ("Not implemented!\n");

}

t8dg_dof_values_t  *
t8dg_dof_values_duplicate (t8dg_dof_values_t * src_dof_values)
{
  T8DG_ABORT ("Not implemented!\n");
}

t8dg_dof_values_t  *
t8dg_dof_values_clone (t8dg_dof_values_t * src_dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

int
t8dg_dof_values_is_valid (t8dg_dof_values_t * dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_axpy (double a, const t8dg_dof_values_t * x, t8dg_dof_values_t * y)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_axpyz (double a, const t8dg_dof_values_t * x, const t8dg_dof_values_t * y, t8dg_dof_values_t * z)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_debug_print (t8dg_dof_values_t * array)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_print (t8dg_dof_values_t * array)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_square_values (t8dg_dof_values_t * src, t8dg_dof_values_t * dest)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_swap (t8dg_dof_values_t ** parray1, t8dg_dof_values_t ** parray2)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_copy (const t8dg_dof_values_t * src, t8dg_dof_values_t * dest)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_set_zero (t8dg_dof_values_t * array)
{
  T8DG_ABORT ("Not implemented!\n");

}

t8dg_element_dof_values_t *
t8dg_element_dof_values_duplicate (const t8dg_element_dof_values_t * element_dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

t8dg_element_dof_values_t *
t8dg_element_dof_values_clone (const t8dg_element_dof_values_t * element_dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

t8dg_face_dof_values_t *
t8dg_face_dof_values_duplicate (const t8dg_face_dof_values_t * face_dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

int
t8dg_element_dof_values_is_valid (t8dg_element_dof_values_t * element_dof_values)
{
  T8DG_ABORT ("Not implemented!\n");

}

int
t8dg_face_dof_values_is_valid (t8dg_face_dof_values_t * face_dof_values)
{
  T8DG_ABORT ("Not implemented \n ");
}

void
t8dg_element_dof_values_axpy (double a, const t8dg_element_dof_values_t * x, t8dg_element_dof_values_t * y)
{
  T8DG_ABORT ("Not implemented!\n");

}
