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
