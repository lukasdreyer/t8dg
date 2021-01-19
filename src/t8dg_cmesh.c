#include "t8dg_cmesh.h"
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>

int
t8dg_cmesh_dim (int icmesh)
{
  if (icmesh >= 0 && icmesh <= 2)
    return 1;
  if (icmesh >= 3 && icmesh <= 7)
    return 2;
  if (icmesh == 8)
    return 3;
  return -1;
}

t8_cmesh_t
t8dg_cmesh_new_arg (int icmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  switch (icmesh) {
  case 0:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, comm, 0, 0, 1);
    break;
  case 1:
    cmesh = t8_cmesh_new_periodic_line_more_trees (comm);
    break;
  case 2:
    cmesh = t8dg_cmesh_new_periodic_diagonal_line_more_trees (comm);
    break;
  case 3:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 1);
    break;
  case 4:
    cmesh = t8dg_cmesh_new_square_more_trees_different_size (comm);
    break;
  case 5:
    cmesh = t8dg_cmesh_new_square_moebius (comm);
    break;
  case 6:
    cmesh = t8dg_cmesh_new_half_moebius_more_trees (comm);
    break;
  case 7:
    cmesh = t8dg_cmesh_new_square_tilted (comm);
    break;
  case 8:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 1);
    break;

  default:
    T8DG_ABORT ("Not yet implemented");
  }
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_periodic_diagonal_line_more_trees (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices[12] = {
    0, 0, 0,
    0.2, 0.2, 0.2,
    0.6, 0.6, 0.6,
    1, 1, 1
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_LINE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 2);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices + 3, 2);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0, vertices + 6, 2);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 1, 0, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_half_moebius_more_trees (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices0[3 * 4] = {
    0, 0, 0,
    0.5, 0, 0,
    0, 0.5, 0,
    0.5, 0.5, 0
  };
  double              vertices1[3 * 4] = {
    0.5, 0, 0,
    1, 0, 0,
    0.5, 0.5, 0,
    1, 0.5, 0
  };
  double              vertices2[3 * 4] = {
    0, 0.5, 0,
    0.5, 0.5, 0,
    0, 1, 0,
    0.5, 1, 0
  };
  double              vertices3[3 * 4] = {
    0.5, 0.5, 0,
    1, 0.5, 0,
    0.5, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices0, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices1, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0, vertices2, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0, vertices3, 4);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 1);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 0, 1, 0, 1);
  t8_cmesh_set_join (cmesh, 0, 2, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 3, 1, 3, 2, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_more_trees_different_size (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices0[3 * 4] = {
    0, 0, 0,
    0.4, 0, 0,
    0, 0.4, 0,
    0.4, 0.4, 0
  };
  double              vertices1[3 * 4] = {
    0.4, 0, 0,
    1, 0, 0,
    0.4, 0.4, 0,
    1, 0.4, 0
  };
  double              vertices2[3 * 4] = {
    0, 0.4, 0,
    0.4, 0.4, 0,
    0, 1, 0,
    0.4, 1, 0
  };
  double              vertices3[3 * 4] = {
    0.4, 0.4, 0,
    1, 0.4, 0,
    0.4, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices0, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices1, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0, vertices2, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0, vertices3, 4);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 0, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 2, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 2, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 3, 1, 3, 2, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_moebius (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices[3 * 4] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);
  t8_cmesh_set_join (cmesh, 0, 0, 1, 0, 1);
  t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_tilted (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices[3 * 4] = {
    0, 0, 0,
    2, 0, 0,
    -1, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);
  t8_cmesh_set_join (cmesh, 0, 0, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm, int num_trees)
{
  t8_cmesh_t          cmesh;

  double              vertices[6] = { 0 };
  int                 itree;

  t8_cmesh_init (&cmesh);

  for (itree = 0; itree < num_trees; itree++) {
    vertices[0] = 1. / num_trees * itree;
    vertices[3] = 1. / num_trees * (itree + 1);
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_LINE);
    t8_cmesh_set_tree_vertices (cmesh, itree, t8_get_package_id (), 0, vertices, 2);
    if (itree < num_trees - 1) {
      t8_cmesh_set_join (cmesh, itree, itree + 1, 1, 0, 0);
    }
  }
  t8_cmesh_set_join (cmesh, num_trees - 1, 0, 1, 0, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_circle_ring (sc_MPI_Comm comm, double inner_radius, double outer_radius, int num_trees)
{
  t8_cmesh_t          cmesh;

  double              vertices[12] = { 0 };
  int                 itree;
  double              angle = 2 * M_PI / num_trees;
  t8_cmesh_init (&cmesh);

  for (itree = 0; itree < num_trees; itree++) {
    vertices[0] = cos (itree * angle) * inner_radius;
    vertices[1] = sin (itree * angle) * inner_radius;
    vertices[3] = cos (itree * angle) * outer_radius;
    vertices[4] = sin (itree * angle) * outer_radius;
    vertices[6] = cos ((itree + 1) * angle) * inner_radius;
    vertices[7] = sin ((itree + 1) * angle) * inner_radius;
    vertices[9] = cos ((itree + 1) * angle) * outer_radius;
    vertices[10] = sin ((itree + 1) * angle) * outer_radius;

    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_vertices (cmesh, itree, t8_get_package_id (), 0, vertices, 4);
    if (itree < num_trees - 1) {
      t8_cmesh_set_join (cmesh, itree, itree + 1, 3, 2, 0);
    }
  }
  t8_cmesh_set_join (cmesh, num_trees - 1, 0, 3, 2, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
