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
#include "../src/t8dg_cmesh.h"
#include <t8_cmesh_vtk.h>
#include "../src/t8dg_coarse_geometry.h"

TEST (circle_ring, construction)
{
  t8_cmesh_t          cmesh;
  cmesh = t8dg_cmesh_new_circle_ring (sc_MPI_COMM_WORLD, 1, 2, 4);
  double              expected_tree_vertices[12] = { 1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    0, 2, 0
  };
  double             *tree_vertices = t8_cmesh_get_tree_vertices (cmesh, 0);
  int                 icoordinate;
  for (icoordinate = 0; icoordinate < 12; icoordinate++) {
    EXPECT_NEAR (tree_vertices[icoordinate], expected_tree_vertices[icoordinate], 1e-12);
  }
  t8_cmesh_vtk_write_file (cmesh, "test_cmesh", 1);
  t8_cmesh_destroy (&cmesh);
}

/*
TEST (circle_ring, d_inv_trans)
{
    t8_cmesh_t cmesh = t8dg_cmesh_new_circle_ring (sc_MPI_COMM_WORLD, 1, 2, 4);
    t8_forest_t forest = t8_forest_new_uniform(cmesh,)
    t8dg_coarse_geometry_t *circle_ring_geometry;
    double coarse_vertex[3] = {0,0,0};
    double coarse_tangential_vector[3]={-1,0,0};
    double transformed_tangential_vector[3];

    circle_ring_geometry = t8dg_coarse_geometry_new_2D_circle_ring(1,2);
    t8dg_coarse_geometry_calculate_gradient_tangential_vector(circle_ring_geometry, forest, 0, coarse_vertex, coarse_tangential_vector, transformed_tangential_vector);

}
*/
