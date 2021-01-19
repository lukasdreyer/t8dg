#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_cmesh.h"
#include <t8_cmesh_vtk.h>

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
