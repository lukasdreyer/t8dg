#include "t8dg_cmesh.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_advect_diff_problem.h"
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vec.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_helpers.h>

t8_cmesh_t
t8dg_cmesh_new_arg (int icmesh, const char *mshfile_prefix, const int mshfile_dim, int *dim, t8dg_flow_type_t *velocity_field_arg, int *geometry_arg, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  T8_ASSERT (icmesh == 12 || mshfile_prefix == NULL);
  /* Build a cmesh depending on the icmesh parameter. */
  switch (icmesh) {
  case 0:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, comm, 0, 0, 1);
    *dim = 1;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 0;
    break;
  case 1:
    cmesh = t8_cmesh_new_periodic_line_more_trees (comm);
    *dim = 1;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 0;
    break;
  case 2:
    cmesh = t8dg_cmesh_new_periodic_diagonal_line_more_trees (comm);
    *dim = 1;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;    /*flow currently to the right */
    *geometry_arg = 0;
    break;
  case 3:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 1);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 1;
    break;
  case 4:
    cmesh = t8dg_cmesh_new_square_more_trees_different_size (comm);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 1;
    break;
  case 5:
    cmesh = t8dg_cmesh_new_square_moebius (comm);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 1;
    break;
  case 6:
    cmesh = t8dg_cmesh_new_half_moebius_more_trees (comm);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 1;
    break;
  case 7:
    cmesh = t8dg_cmesh_new_square_tilted (comm);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 1;
    break;
  case 8:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 1);
    *dim = 3;
    *velocity_field_arg = T8DG_FLOW_MPTRAC_3D;
    *geometry_arg = 2;
    break;
  case 9:
    cmesh = t8dg_cmesh_new_circle_ring (comm, 1, 2, 4);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_ROTATE_2D;
    *geometry_arg = 3;
    break;
  case 10:
    cmesh = t8dg_cmesh_new_square_half_periodic (comm);
    *dim = 2;
    *velocity_field_arg = T8DG_FLOW_CONSTANT_3D;
    *geometry_arg = 1;
    break;
  case 11:
    cmesh = t8dg_cmesh_new_cylinder_ring_periodic (comm, 1, 2, 4, 1);
    *dim = 3;
    *velocity_field_arg = T8DG_FLOW_SPIRAL_3D;
    *geometry_arg = 4;
    break;
  case 12:
    /* load the cmesh from a file */
    T8DG_ASSERT (mshfile_prefix != NULL);
    cmesh = t8_cmesh_from_msh_file (mshfile_prefix, 0, comm, mshfile_dim, 0);
    *dim = mshfile_dim;
    *velocity_field_arg = T8DG_FLOW_MPTRAC_3D;
    /* TODO: What does this mean? */
    *geometry_arg = 0;
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

/*
  t8dg_cylinder_ring_data_t *geometry_data;
  geometry_data = (t8dg_cylinder_ring_data_t *) user_data;
  radius = geometry_data->inner_radius + ref_coords[0] * (geometry_data->outer_radius - geometry_data->inner_radius);
  angle = (ref_coords[1] + gtreeid) * 2 * M_PI / geometry_data->num_trees;
*/
  radius = 1 + ref_coords[0];
  angle = (ref_coords[1] + gtreeid) * 2 * M_PI / 4;

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
    t8_cmesh_set_join (cmesh, itree, itree, 1, 0, 0);   /*connect inner and outer ring */
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

static void
t8dg_analytic_cylinder_ring (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
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

//  t8dg_cylinder_ring_data_t *geometry_data = (t8dg_cylinder_ring_data_t *) user_data;

//  radius = geometry_data->inner_radius + ref_coords[0] * (geometry_data->outer_radius - geometry_data->inner_radius);
  radius = 1 + ref_coords[0];
//  angle = (ref_coords[1] + gtreeid) * 2 * M_PI / geometry_data->num_trees;
  angle = (ref_coords[1] + gtreeid) * 2 * M_PI / 4;

  out_coords[0] = radius * cos (angle);
  out_coords[1] = radius * sin (angle);
  out_coords[2] = ref_coords[2];
//  out_coords[2] = geometry_data->height * ref_coords[2];
}

t8_cmesh_t
t8dg_cmesh_new_cylinder_ring_periodic (sc_MPI_Comm comm, double inner_radius, double outer_radius, int num_trees, double height)
{
//  t8dg_cylinder_ring_data_t *geometry_data;
  t8_cmesh_t          cmesh;
  t8_geometry_c      *geometry;

  double              vertices[24] = { 0 };
  int                 itree;
  double              angle = 2 * M_PI / num_trees;

  t8_cmesh_init (&cmesh);

  for (itree = 0; itree < num_trees; itree++) {
    vertices[0] = cos (itree * angle) * inner_radius;
    vertices[1] = sin (itree * angle) * inner_radius;
    vertices[2] = 0;
    vertices[3] = cos (itree * angle) * outer_radius;
    vertices[4] = sin (itree * angle) * outer_radius;
    vertices[5] = 0;
    vertices[6] = cos ((itree + 1) * angle) * inner_radius;
    vertices[7] = sin ((itree + 1) * angle) * inner_radius;
    vertices[8] = 0;
    vertices[9] = cos ((itree + 1) * angle) * outer_radius;
    vertices[10] = sin ((itree + 1) * angle) * outer_radius;
    vertices[11] = 0;
    vertices[12] = cos (itree * angle) * inner_radius;
    vertices[13] = sin (itree * angle) * inner_radius;
    vertices[14] = height;
    vertices[15] = cos (itree * angle) * outer_radius;
    vertices[16] = sin (itree * angle) * outer_radius;
    vertices[17] = height;
    vertices[18] = cos ((itree + 1) * angle) * inner_radius;
    vertices[19] = sin ((itree + 1) * angle) * inner_radius;
    vertices[20] = height;
    vertices[21] = cos ((itree + 1) * angle) * outer_radius;
    vertices[22] = sin ((itree + 1) * angle) * outer_radius;
    vertices[23] = height;

    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_HEX);
    t8_cmesh_set_tree_vertices (cmesh, itree, vertices, 8);
    t8_cmesh_set_join (cmesh, itree, itree, 1, 0, 0);   /*connect inner and outer ring */
    t8_cmesh_set_join (cmesh, itree, itree, 4, 5, 0);   /*connect upper and lower end */
    if (itree < num_trees - 1) {
      t8_cmesh_set_join (cmesh, itree, itree + 1, 3, 2, 0);
    }
  }
  t8_cmesh_set_join (cmesh, num_trees - 1, 0, 3, 2, 0);

/*  geometry_data = T8DG_ALLOC_ZERO (t8dg_cylinder_ring_data_t, 1);
  geometry_data->inner_radius = inner_radius;
  geometry_data->outer_radius = outer_radius;
  geometry_data->num_trees = num_trees;
  geometry_data->height = height;
*/

  geometry =
    new t8_geometry_analytic (2, "analytic cylinder ring", t8dg_analytic_cylinder_ring, NULL, t8_geom_load_tree_data_vertices, NULL);
/*  geometry =
    new t8_geometry_analytic (2, "analytic cylinder ring", t8dg_analytic_cylinder_ring, NULL, t8_geom_load_tree_data_vertices,
                              geometry_data);*/
  t8_cmesh_register_geometry (cmesh, geometry);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

#if 0
static void
t8dg_analytic_sphere_extruded (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                               const double *ref_coords, double out_coords[3], const void *tree_vertices, const void *user_data)
{

  const double       *tree_v = (const double *) tree_vertices;
  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  t8_eclass_t         tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  t8_geom_compute_linear_geometry (tree_class, tree_v, ref_coords, out_coords);

  int                 direction = gtreeid / 2;
  double              radius = out_coords[direction];
  out_coords[direction] = (gtreeid % 2) * 2 - 1;        /*first of two faces in one direction maps to -1, the other to 1 */
  double              norm = t8_vec_norm (out_coords);
  int                 idim;
  for (idim = 0; idim < DIM3; idim++) {
    out_coords[idim] *= radius / norm;
  }
}

t8_cmesh_t
t8dg_cmesh_new_sphere_extruded (sc_MPI_Comm comm, void *t8dg_sphere_extruded_data)
{
  t8_cmesh_t          cmesh;
  t8_geometry_c      *geometry;

  double              vertices[24] = { 0 };
  int                 itree, num_trees = 6;
  int                 distance_coordinate, x_coordinate, y_coordinate, distance_factor, ivertex;

  double              inner_radius;
  double              outer_radius;     // get from user_data

  t8_cmesh_init (&cmesh);

  for (itree = 0; itree < num_trees; itree++) {
/* TODO:
    distance_coordinate = itree/2;
    x_coordinate = ;
    y_coordinate = ;
    distance_factor = 2* itree -1;
*/
    for (ivertex = 0; ivertex < 8; ivertex++) {
      vertices[ivertex * 3 + x_coordinate] = 2 * (ivertex % 2) - 1;
      vertices[ivertex * 3 + y_coordinate] = 2 * (ivertex % 4) - 1;
      vertices[ivertex * 3 + distance_coordinate] = (ivertex % 8) ? inner_radius : outer_radius;
      vertices[ivertex * 3 + distance_coordinate] *= distance_factor;
    }
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_HEX);
    t8_cmesh_set_tree_vertices (cmesh, itree, vertices, 8);
  }

  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);

  geometry =
    new t8_geometry_analytic (2, "analytic cylinder ring", t8dg_analytic_cylinder_ring, NULL, t8_geom_load_tree_data_vertices,
                              t8dg_sphere_extruded_data);
  t8_cmesh_register_geometry (cmesh, geometry);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
#endif
