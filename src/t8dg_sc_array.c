/*
 * t8dg_sc_array.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */
#include<t8.h>
#include <sc_containers.h>

sc_array_t         *
t8dg_sc_array_clone (sc_array_t * src)
{
  sc_array_t         *dest = sc_array_new_count (src->elem_size, src->elem_count);
  sc_array_copy (dest, src);
  return dest;
}

void
t8dg_sc_array_block_double_axpy (const double a, const sc_array_t * x, sc_array_t * y)
{
  T8_ASSERT (x != NULL && y != NULL);
  T8_ASSERT (x->array != NULL && y->array != NULL);
  T8_ASSERT (x->elem_size % sizeof (double) == 0 && y->elem_size % sizeof (double) == 0);
  T8_ASSERT (x->elem_count * x->elem_size == y->elem_count * y->elem_size);
  double             *x_double, *y_double;
  int                 double_count, i;

  /*View array as double array */
  x_double = (double *) x->array;
  y_double = (double *) y->array;
  double_count = x->elem_size / sizeof (double) * x->elem_count;        /*total number of doubles */

  for (i = 0; i < double_count; i++) {
    y_double[i] = a * x_double[i] + y_double[i];
  }
}

void
t8dg_sc_array_block_double_zaxpy (sc_array_t * z, double a, const sc_array_t * x, const sc_array_t * y)
{
  T8_ASSERT (z != NULL && x != NULL && y != NULL);
  T8_ASSERT (z->array != NULL && x->array != NULL && y->array != NULL);
  T8_ASSERT (x->elem_size % sizeof (double) == 0 && y->elem_size % sizeof (double) == 0 && z->elem_size % sizeof (double) == 0);
  T8_ASSERT (x->elem_count * x->elem_size == y->elem_count * y->elem_size);
  T8_ASSERT (x->elem_count * x->elem_size == z->elem_count * z->elem_size);

  /*View array as double array */
  double             *x_double, *y_double, *z_double;
  size_t              double_count, idouble;

  x_double = (double *) (x->array);
  y_double = (double *) (y->array);
  z_double = (double *) (z->array);

  double_count = (x->elem_size / sizeof (double)) * x->elem_count;      /*total number of doubles */

  for (idouble = 0; idouble < double_count; idouble++) {
    z_double[idouble] = a * x_double[idouble] + y_double[idouble];
  }
}

void
t8dg_sc_array_swap (sc_array_t ** parray1, sc_array_t ** parray2)
{
  sc_array_t         *temp;
  temp = *parray1;
  *parray1 = *parray2;
  *parray2 = temp;
}

void
t8dg_sc_array_block_double_print (sc_array_t * array)
{
  size_t              irow, icolumn;
  double             *row_array;
  for (irow = 0; irow < array->elem_count; irow++) {
    row_array = (double *) t8_sc_array_index_locidx (array, irow);
    for (icolumn = 0; icolumn < array->elem_size / 8; icolumn++) {
      printf ("%f  ,  ", row_array[icolumn]);
    }
    printf ("\n");
  }
  printf ("\n");
}
