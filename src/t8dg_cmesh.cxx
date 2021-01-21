#include "t8dg_cmesh.h"
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_helpers.h>

int
t8dg_cmesh_dim (int icmesh)
{
  if (icmesh >= 0 && icmesh <= 2)
    return 1;
  if ((icmesh >= 3 && icmesh <= 7) || icmesh == 9 || icmesh == 10)
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
  case 9:
    cmesh = t8dg_cmesh_new_circle_ring (comm, 1, 2, 4);
    break;
  case 10:
    cmesh = t8dg_cmesh_new_square_half_periodic (comm);
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
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (1);

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
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 2);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 3, 2);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 6, 2);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 1, 0, 0);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_half_moebius_more_trees (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2); /*maybe do moebius geometry */

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
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices0, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices1, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices2, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices3, 4);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 1);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 0, 1, 0, 1);
  t8_cmesh_set_join (cmesh, 0, 2, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 3, 1, 3, 2, 0);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_more_trees_different_size (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

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
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices0, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices1, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices2, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices3, 4);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 0, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 2, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 2, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 3, 1, 3, 2, 0);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_moebius (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  double              vertices[3 * 4] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_set_join (cmesh, 0, 0, 1, 0, 1);
  t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_half_periodic (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  double              vertices[3 * 4] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_set_join (cmesh, 0, 0, 1, 0, 1);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_square_tilted (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  double              vertices[3 * 4] = {
    0, 0, 0,
    2, 0, 0,
    -1, 1, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_set_join (cmesh, 0, 0, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8dg_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm, int num_trees)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (1);

  double              vertices[6] = { 0 };
  int                 itree;

  t8_cmesh_init (&cmesh);

  for (itree = 0; itree < num_trees; itree++) {
    vertices[0] = 1. / num_trees * itree;
    vertices[3] = 1. / num_trees * (itree + 1);
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_LINE);
    t8_cmesh_set_tree_vertices (cmesh, itree, vertices, 2);
    if (itree < num_trees - 1) {
      t8_cmesh_set_join (cmesh, itree, itree + 1, 1, 0, 0);
    }
  }
  t8_cmesh_set_join (cmesh, num_trees - 1, 0, 1, 0, 0);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static void
t8dg_analytic_circle_ring (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                           const double *ref_coords, double out_coords[3], const void *tree_vertices, const void *user_data)
{
/*
  const double       *tree_v = (const double *) tree_vertices;
  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  t8_eclass_t         tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  t8_geom_compute_linear_geometry (tree_class, tree_v, ref_coords,
                                   out_coords);
*/
  double              radius;
  double              angle;

  radius = 1 + ref_coords[0];
  angle = ref_coords[1] * M_PI_2 + gtreeid * M_PI_2;

  out_coords[0] = radius * cos (angle);
  out_coords[1] = radius * sin (angle);
  out_coords[2] = 0;

}

t8_cmesh_t
t8dg_cmesh_new_circle_ring (sc_MPI_Comm comm, double inner_radius, double outer_radius, int num_trees)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *geometry;

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
    t8_cmesh_set_tree_vertices (cmesh, itree, vertices, 4);
    if (itree < num_trees - 1) {
      t8_cmesh_set_join (cmesh, itree, itree + 1, 3, 2, 0);
    }
  }
  t8_cmesh_set_join (cmesh, num_trees - 1, 0, 3, 2, 0);

  geometry = new t8_geometry_analytic (2, "analytic moebius", t8dg_analytic_circle_ring, NULL, t8_geom_load_tree_data_vertices, NULL);
  t8_cmesh_register_geometry (cmesh, geometry);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
