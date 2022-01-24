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
#include "../src/t8dg_quadrature.h"
#include "../src/t8dg_vertexset.h"
/* This fail just demonstrates that the test framework is working. */

static double
const_one (const double x[3], void *data)
{
  return 1;
}

static double
x_pow (const double x[3], void *data)
{
  int                 power;
  power = *(int *) data;
  return pow (x[0], power);
}

static double
x_times_one_minus_y_pow (const double x[3], void *data)
{
  int                 power;
  power = *(int *) data;
  return pow (x[0] * (1 - x[1]), power);
}

/* If this test fails, something is really wrong. */
TEST (quadrature1D, integrate_const_one)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 1; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    EXPECT_EQ (t8dg_vertexset_get_num_vertices (vertexset), num_LGL);
    quadrature = t8dg_quadrature_new_vertexset (vertexset, 0);
    EXPECT_EQ (t8dg_quadrature_get_num_element_vertices (quadrature), num_LGL);
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature, const_one, NULL), 1, 1e-10);
    t8dg_quadrature_destroy (&quadrature);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature1D, integrate_max_order)
{
  int                 num_LGL, power;
  t8dg_quadrature_t  *quadrature;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 2; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    quadrature = t8dg_quadrature_new_vertexset (vertexset, 0);
    power = 2 * num_LGL - 3;
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature, x_pow, &power), 1.0 / (power + 1), 1e-10);
    t8dg_quadrature_destroy (&quadrature);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature2D, integrate_const_one)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature_square;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 1; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    EXPECT_EQ (t8dg_vertexset_get_num_vertices (vertexset), num_LGL);
    quadrature_square = t8dg_quadrature_new_hypercube (2, vertexset, 0);
    EXPECT_EQ (t8dg_quadrature_get_num_element_vertices (quadrature_square), num_LGL * num_LGL);
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature_square, const_one, NULL), 1, 1e-10);
    t8dg_quadrature_destroy (&quadrature_square);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature2D, integrate_max_order)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature_square;
  t8dg_vertexset_t   *vertexset;
  int                 power;
  double              result;
  for (num_LGL = 2; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    EXPECT_EQ (t8dg_vertexset_get_num_vertices (vertexset), num_LGL);
    quadrature_square = t8dg_quadrature_new_hypercube (2, vertexset, 0);
    EXPECT_EQ (t8dg_quadrature_get_num_element_vertices (quadrature_square), num_LGL * num_LGL);
    power = 2 * num_LGL - 3;
    result = 1. / ((power + 1) * (power + 1));
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature_square, x_times_one_minus_y_pow, &power), result, 1e-10);
    t8dg_quadrature_destroy (&quadrature_square);
    t8dg_vertexset_destroy (&vertexset);
  }
}
