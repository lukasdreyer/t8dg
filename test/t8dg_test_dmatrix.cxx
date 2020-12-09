#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_dmatrix.h"
#include <sc_containers.h>

TEST (dmatrix, mult)
{
  t8dg_dmatrix_t     *matrix;
  sc_array_t         *array;
  sc_array_t         *res_array;
  matrix = t8dg_dmatrix_new (2, 3);
  array = sc_array_new_count (sizeof (double), 3);
  res_array = sc_array_new_count (sizeof (double), 2);
  t8dg_dmatrix_set_at (matrix, 0, 0, 1);
  t8dg_dmatrix_set_at (matrix, 0, 1, 2);
  t8dg_dmatrix_set_at (matrix, 0, 2, 3);
  t8dg_dmatrix_set_at (matrix, 1, 0, 4);
  t8dg_dmatrix_set_at (matrix, 1, 1, 5);
  t8dg_dmatrix_set_at (matrix, 1, 2, 6);
  *(double *) sc_array_index (array, 0) = 1;
  *(double *) sc_array_index (array, 1) = -3;
  *(double *) sc_array_index (array, 2) = 6;
  t8dg_dmatrix_mult_sc_array (matrix, array, res_array);
#ifdef T8DG_ENABLE_MPI
  EXPECT_EQ_MPI (*(double *) sc_array_index (res_array, 0), 13);
  EXPECT_EQ_MPI (*(double *) sc_array_index (res_array, 1), 25);
#else
  EXPECT_EQ (*(double *) sc_array_index (res_array, 0), 13);
  EXPECT_EQ (*(double *) sc_array_index (res_array, 1), 25);
#endif
  t8dg_dmatrix_destroy (&matrix);
  sc_array_destroy (array);
  sc_array_destroy (res_array);

}
