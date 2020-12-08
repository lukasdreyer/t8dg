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
t8dg_dof_values_new_element_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t itree, t8_locidx_t ielement)
{
//  T8DG_ASSERT(itree >= 0 && itree < t8_forest_get_num_local_trees(dof_values->forest));
  t8_locidx_t         idata;

  t8dg_element_dof_values_t *element_dof_view;
  t8dg_dofidx_t       num_element_dof;
  if (itree == -1) {
    num_element_dof = dof_values->max_num_element_dof;  /*TODO: how to get eclass and therefore right num of element_dof for ghost elements */
    idata = ielement;
  }
  else {
    t8dg_global_values_t *global_values;
    global_values = t8dg_global_values_array_get_global_values (dof_values->global_values, dof_values->forest, itree, ielement);
    num_element_dof = t8dg_global_values_get_num_dof (global_values);
    idata = t8dg_itree_ielement_to_idata (dof_values->forest, itree, ielement);
  }
  T8DG_ASSERT (idata >= 0 && (size_t) idata < dof_values->num_total_elements);
  element_dof_view =
    (t8dg_element_dof_values_t *) sc_array_new_data (t8_sc_array_index_locidx (dof_values->dofs, idata), sizeof (double), num_element_dof);
  /* Create new view */
  return element_dof_view;
}

void
t8dg_dof_values_destroy (t8dg_dof_values_t ** p_dof_values)
{
  t8dg_dof_values_t  *dof_values;
  dof_values = *p_dof_values;
  sc_array_destroy_null (&dof_values->dofs);
  T8DG_FREE (dof_values);
  *p_dof_values = NULL;
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
  T8DG_ASSERT (t8_forest_is_committed (forest));
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

t8dg_dof_values_t  *
t8dg_dof_values_new_data_local (t8_forest_t forest, t8dg_global_values_t ** global_values_array, double *array,
                                t8dg_dofidx_t num_total_values)
{
  T8DG_ASSERT (t8_forest_is_committed (forest));
  t8dg_dof_values_t  *return_values;
  return_values = T8DG_ALLOC_ZERO (t8dg_dof_values_t, 1);
  return_values->forest = forest;
  return_values->global_values = global_values_array;
  return_values->max_num_element_dof = t8dg_global_values_array_get_max_num_element_dof (global_values_array);
  return_values->num_local_elements = t8_forest_get_num_element (forest);
  return_values->num_ghost_elements = 0;
  return_values->num_total_elements = return_values->num_local_elements + return_values->num_ghost_elements;

  size_t              offset = t8_forest_get_first_local_element_id (forest);

  return_values->dofs =
    sc_array_new_data ((void *) (array + offset * return_values->max_num_element_dof), sizeof (double) * return_values->max_num_element_dof,
                       return_values->num_total_elements);
  return return_values;
}

int
t8dg_dof_values_equal (t8dg_dof_values_t * dof_values, t8dg_dof_values_t * dof_values_compare)
{
  t8dg_element_dof_values_t *element_dof_values;
  t8dg_element_dof_values_t *element_dof_values_compare;

  if (dof_values->forest != dof_values_compare->forest)
    return 0;
  if (dof_values->global_values != dof_values_compare->global_values)
    return 0;
  t8_locidx_t         itree, ielement, num_trees, num_elems_in_tree;
  num_trees = t8_forest_get_num_local_trees (dof_values->forest);
  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (dof_values->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      element_dof_values = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielement);
      element_dof_values_compare = t8dg_dof_values_new_element_dof_values_view (dof_values_compare, itree, ielement);

      if (!t8dg_element_dof_values_equal (element_dof_values, element_dof_values_compare))
        return 0;

      t8dg_element_dof_values_destroy (&element_dof_values);
      t8dg_element_dof_values_destroy (&element_dof_values_compare);
    }
  }
  return 1;
}

int
t8dg_element_dof_values_equal (t8dg_element_dof_values_t * element_dof_values, t8dg_element_dof_values_t * element_dof_values_compare)
{
  double             *double_array, *double_array_compare;
  if (element_dof_values->elem_size != sizeof (double) || element_dof_values_compare->elem_size != sizeof (double))
    return 0;
  if (element_dof_values->elem_count != element_dof_values_compare->elem_count)
    return 0;
  double_array = (double *) element_dof_values->array;
  double_array_compare = (double *) element_dof_values_compare->array;
  size_t              idx;
  for (idx = 0; idx < element_dof_values->elem_count; idx++) {
    if (double_array[idx] != double_array_compare[idx])
      return 0;
  }
  return 1;
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
t8dg_element_dof_values_set_zero (t8dg_element_dof_values_t * element_dof_values)
{
  t8dg_dofidx_t       idof;
  for (idof = 0; idof < element_dof_values->elem_count; idof++) {
    t8dg_element_dof_values_set_value (element_dof_values, idof, 0);
  }
}

double
t8dg_element_dof_values_get_value (t8dg_element_dof_values_t * element_dof_values, t8dg_dofidx_t idof)
{
  T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dof_values));
  T8DG_ASSERT (idof >= 0 && idof < element_dof_values->elem_count);
  return *(double *) sc_array_index (element_dof_values, idof);
}

void
t8dg_element_dof_values_set_value (t8dg_element_dof_values_t * element_dof_values, t8dg_dofidx_t idof, double value)
{
  T8DG_ASSERT (t8dg_element_dof_values_is_valid (element_dof_values));
  T8DG_ASSERT (idof >= 0 && idof < element_dof_values->elem_count);
  *(double *) sc_array_index (element_dof_values, idof) = value;
}

t8dg_face_dof_values_t *
t8dg_dof_values_new_face_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t idata)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_element_dof_values_destroy (t8dg_element_dof_values_t ** dof_values)
{
  sc_array_destroy_null (dof_values);
}

void
t8dg_face_dof_values_destroy (t8dg_face_dof_values_t ** dof_values)
{
  sc_array_destroy_null (dof_values);

}

void
t8dg_dof_values_copy_from_index_to_index (t8dg_dof_values_t * src_dof, t8_locidx_t src_idata, t8dg_dof_values_t * dest_dof,
                                          t8_locidx_t dest_idata)
{
  memcpy (t8_sc_array_index_locidx (dest_dof->dofs, dest_idata),
          t8_sc_array_index_locidx (src_dof->dofs, src_idata), src_dof->dofs->elem_size);
}

t8dg_dof_values_t  *
t8dg_dof_values_duplicate (t8dg_dof_values_t * src_dof_values)
{
  t8dg_dof_values_t  *dest_dof_values;
  dest_dof_values = t8dg_dof_values_new (src_dof_values->forest, src_dof_values->global_values);
  return dest_dof_values;
}

t8dg_dof_values_t  *
t8dg_dof_values_clone (t8dg_dof_values_t * src_dof_values)
{
  t8dg_dof_values_t  *dest_dof_values;
  dest_dof_values = t8dg_dof_values_duplicate (src_dof_values);
  t8dg_dof_values_copy (src_dof_values, dest_dof_values);
  return dest_dof_values;
}

int
t8dg_dof_values_is_valid (t8dg_dof_values_t * dof_values)
{
  return t8_forest_is_committed (dof_values->forest);
  /*TODO: implement more checks */
}

void
t8dg_dof_values_axpy (double a, const t8dg_dof_values_t * x, t8dg_dof_values_t * y)
{
  double             *x_double, *y_double;
  int                 double_count, i;

  /*View array as double array */
  x_double = (double *) x->dofs->array;
  y_double = (double *) y->dofs->array;
  double_count = x->max_num_element_dof * x->num_total_elements;        /*total number of doubles */

  for (i = 0; i < double_count; i++) {
    y_double[i] = a * x_double[i] + y_double[i];
  }

}

void
t8dg_dof_values_axpyz (double a, const t8dg_dof_values_t * x, const t8dg_dof_values_t * y, t8dg_dof_values_t * z)
{
  double             *x_double, *y_double, *z_double;
  int                 double_count, i;

  /*View array as double array */
  x_double = (double *) x->dofs->array;
  y_double = (double *) y->dofs->array;
  z_double = (double *) z->dofs->array;
  double_count = x->max_num_element_dof * x->num_total_elements;        /*total number of doubles */

  for (i = 0; i < double_count; i++) {
    z_double[i] = a * x_double[i] + y_double[i];
  }

}

void
t8dg_dof_values_debug_print (t8dg_dof_values_t * array)
{
#ifdef T8_ENABLE_DEBUG
  t8dg_dof_values_print (array);
#endif
}

void
t8dg_dof_values_print (t8dg_dof_values_t * array)
{
  size_t              irow, icolumn;
  double             *row_array;
  for (irow = 0; irow < array->dofs->elem_count; irow++) {
    row_array = (double *) t8_sc_array_index_locidx (array->dofs, irow);
    for (icolumn = 0; icolumn < array->dofs->elem_size / sizeof (double); icolumn++) {
      printf ("%f  ,  ", row_array[icolumn]);
    }
    printf ("\n");
  }
  printf ("\n");
}

void
t8dg_element_dof_values_debug_print (t8dg_element_dof_values_t * element_dof_values)
{
#ifdef T8_ENABLE_DEBUG
  t8dg_dofidx_t       idof;
  for (idof = 0; idof < element_dof_values->elem_count; idof++) {
    printf ("%f ", t8dg_element_dof_values_get_value (element_dof_values, idof));
  }
  printf ("\n");
#endif
}

void
t8dg_dof_values_square_values (t8dg_dof_values_t * src, t8dg_dof_values_t * dest)
{
  T8DG_ABORT ("Not implemented!\n");

}

void
t8dg_dof_values_swap (t8dg_dof_values_t ** parray1, t8dg_dof_values_t ** parray2)
{
  t8dg_dof_values_t  *temp;
  temp = *parray1;
  *parray1 = *parray2;
  *parray2 = temp;
}

void
t8dg_dof_values_copy (t8dg_dof_values_t * src, t8dg_dof_values_t * dest)
{
  T8DG_ASSERT (dest->forest == src->forest);
  T8DG_ASSERT (dest->global_values == src->global_values);
  T8DG_ASSERT (t8dg_dof_values_is_valid (src) && t8dg_dof_values_is_valid (dest));
  memcpy (dest->dofs->array, src->dofs->array, src->dofs->elem_count * src->dofs->elem_size);
}

void
t8dg_dof_values_set_zero (t8dg_dof_values_t * array)
{
  t8dg_dof_values_set_all_values (array, 0);
}

void
t8dg_dof_values_set_all_values (t8dg_dof_values_t * dof_values, double value)
{
  double             *values;
  values = (double *) dof_values->dofs->array;
  size_t              idx;
  size_t              total_dof;
  total_dof = dof_values->max_num_element_dof * dof_values->num_total_elements;
  for (idx = 0; idx < total_dof; idx++) {
    values[idx] = value;
  }
}

t8dg_element_dof_values_t *
t8dg_element_dof_values_duplicate (t8dg_element_dof_values_t * element_dof_values)
{
  T8DG_ASSERT (element_dof_values->elem_size == sizeof (double));
  t8dg_element_dof_values_t *dest_element_dof_values;
  dest_element_dof_values = sc_array_new_count (sizeof (double), element_dof_values->elem_count);
  return dest_element_dof_values;
}

t8dg_element_dof_values_t *
t8dg_element_dof_values_clone (t8dg_element_dof_values_t * element_dof_values)
{
  t8dg_element_dof_values_t *dest_element_dof_values;
  dest_element_dof_values = t8dg_element_dof_values_duplicate (element_dof_values);
  t8dg_element_dof_values_copy (element_dof_values, dest_element_dof_values);
  return dest_element_dof_values;
}

void
t8dg_element_dof_values_copy (t8dg_element_dof_values_t * src, t8dg_element_dof_values_t * dest)
{
  T8DG_ASSERT (t8dg_element_dof_values_is_valid (src) && t8dg_element_dof_values_is_valid (dest));
  memcpy (dest->array, src->array, src->elem_count * src->elem_size);
}

t8dg_face_dof_values_t *
t8dg_face_dof_values_duplicate (t8dg_face_dof_values_t * face_dof_values)
{
  T8DG_ASSERT (face_dof_values->elem_size == sizeof (double));
  t8dg_face_dof_values_t *dest_face_dof_values;
  dest_face_dof_values = sc_array_new_count (sizeof (double), face_dof_values->elem_count);
  return dest_face_dof_values;

}

void
t8dg_face_dof_values_set_zero (t8dg_face_dof_values_t * face_dof_values)
{
  t8dg_element_dof_values_set_zero (face_dof_values);
}

void
t8dg_face_dof_values_axpy (double a, t8dg_face_dof_values_t * x, t8dg_face_dof_values_t * y)
{
  t8dg_element_dof_values_axpy (a, x, y);
}

int
t8dg_element_dof_values_is_valid (t8dg_element_dof_values_t * element_dof_values)
{
  T8DG_ASSERT (element_dof_values->elem_size == sizeof (double));
  t8dg_dofidx_t       iarray;
  return 1;                     //TODO: delete!
  for (iarray = 0; iarray < element_dof_values->elem_count; iarray++) {
    if (*(double *) sc_array_index (element_dof_values, iarray) != *(double *) sc_array_index (element_dof_values, iarray)) {
      return 0;
    }
  }
  return 1;

}

double
t8dg_face_dof_values_get_value (t8dg_face_dof_values_t * face_dof_values, t8dg_dofidx_t idof)
{
  return *(double *) sc_array_index (face_dof_values, idof);
}

void
t8dg_face_dof_values_set_value (t8dg_face_dof_values_t * face_dof_values, t8dg_dofidx_t idof, double value)
{
  *(double *) sc_array_index (face_dof_values, idof) = value;
}

int
t8dg_face_dof_values_is_valid (t8dg_face_dof_values_t * face_dof_values)
{
  T8DG_ASSERT (face_dof_values->elem_size == sizeof (double));
  t8dg_dofidx_t       iarray;
  return 1;                     //TODO: delete!
  for (iarray = 0; iarray < face_dof_values->elem_count; iarray++) {
    if (*(double *) sc_array_index (face_dof_values, iarray) != *(double *) sc_array_index (face_dof_values, iarray)) {
      return 0;
    }
  }
  return 1;
}

void
t8dg_element_dof_values_axpy (double a, t8dg_element_dof_values_t * x, t8dg_element_dof_values_t * y)
{
  T8DG_ASSERT (t8dg_element_dof_values_is_valid (x));
  T8DG_ASSERT (t8dg_element_dof_values_is_valid (y));
  T8DG_ASSERT (y->elem_count == x->elem_count);
  double             *x_double, *y_double;
  int                 double_count, i;

  /*View array as double array */
  x_double = (double *) x->array;
  y_double = (double *) y->array;
  double_count = x->elem_count; /*total number of doubles */

  for (i = 0; i < double_count; i++) {
    y_double[i] = a * x_double[i] + y_double[i];
  }
}

void
t8dg_element_dof_values_axpyz (double a, t8dg_element_dof_values_t * x, t8dg_element_dof_values_t * y, t8dg_element_dof_values_t * z)
{
  T8DG_ABORT ("Not implemented \n ");
}

t8dg_dofidx_t
t8dg_element_dof_values_get_num_dof (t8dg_element_dof_values_t * element_dof_values)
{
  T8DG_ASSERT (element_dof_values != NULL);
  T8DG_ASSERT (element_dof_values->elem_size == sizeof (double));
  return element_dof_values->elem_count;
}

double
t8dg_element_dof_values_element_norm_infty (t8dg_element_dof_values_t * element_dof_values)
{
  double              norm = 0;
  size_t              idof;
  for (idof = 0; idof < element_dof_values->elem_count; idof++) {
    norm = SC_MAX (norm, fabs (t8dg_element_dof_values_get_value (element_dof_values, idof)));
  }
  return norm;
}

void
t8dg_element_dof_values_square_values (t8dg_element_dof_values_t * element_dof_values,
                                       t8dg_element_dof_values_t * element_dof_square_values)
{
  T8DG_ASSERT (element_dof_values->elem_size == sizeof (double) && element_dof_square_values->elem_size == sizeof (double));
  T8DG_ASSERT (element_dof_square_values->elem_count == element_dof_values->elem_count);
  t8dg_dofidx_t       idof, num_dof;
  double              square;
  num_dof = element_dof_values->elem_count;
  for (idof = 0; idof < num_dof; idof++) {
    square = t8dg_element_dof_values_get_value (element_dof_values, idof);
    square *= square;
    t8dg_element_dof_values_set_value (element_dof_square_values, idof, square);
  }
}

void
t8dg_element_dof_values_set_all (t8dg_element_dof_values_t * element_dof_values, double value)
{
  t8dg_dofidx_t       idof;
  for (idof = 0; idof < element_dof_values->elem_count; idof++) {
    t8dg_element_dof_values_set_value (element_dof_values, idof, value);
  }
}

t8dg_element_dof_values_t *
t8dg_element_dof_values_new (t8dg_dofidx_t num_dof)
{
  t8dg_element_dof_values_t *element_dof_values;
  element_dof_values = sc_array_new_count (sizeof (double), num_dof);
  return element_dof_values;
}
