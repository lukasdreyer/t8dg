#include <t8_cmesh.h>
#include <t8_vec.h>
#include <t8_cmesh_vtk.h>
#include "t8dg_common.h"

double
t8dg_scalar1d_hat_function (const double x[3], const double t)
{
  return 0.5 - (fabs (0.5 - x[0]));
}

double
t8dg_scalar2d_hat_function (const double x[3], const double t)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return sqrt (0.5) - t8_vec_dist (x, center);
}

double
t8dg_scalar3d_norm_function (const double x[3], const double t)
{
  return t8_vec_norm (x);
}

double
t8dg_scalar2d_step_function (const double x[3], const double t)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return t8_vec_dist (x, center) < 0.3;
}

double
t8dg_scalar2d_triangle_step_function (const double x[3], const double t)
{
  return x[0] > 0.3 && x[1] > 0.3 && x[0] + x[1] < 0.9;
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
