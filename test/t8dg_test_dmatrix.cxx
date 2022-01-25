/*
 *  This file is part of t8dg.
 *  t8dg is a discontinuous galerkin solver on adaptive meshes
 *  that uses the adaptive meshing library t8code to distribute
 *  the workload equally among parallel processes to achieve good scaling.
 *
 *  Copyright (C) 2020 the developers
 *
 *  t8dg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  t8dg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with t8dg.  If not, see <http://www.gnu.org/licenses/>.
 */

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

  EXPECT_EQ (*(double *) sc_array_index (res_array, 0), 13);
  EXPECT_EQ (*(double *) sc_array_index (res_array, 1), 25);

  t8dg_dmatrix_destroy (&matrix);
  sc_array_destroy (array);
  sc_array_destroy (res_array);

}
