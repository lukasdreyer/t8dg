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

/** @file t8dg_geometry.h */

#ifndef SRC_T8DG_COARSE_GEOMETRY_H_
#define SRC_T8DG_COARSE_GEOMETRY_H_

#include <sc_containers.h>
#include <t8.h>
#include <t8_forest.h>

#include "t8dg.h"
#include "t8dg_dmatrix.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_circle_ring_data
{
  double              inner_radius;
  double              outer_radius;
  int                 num_trees;
} t8dg_circle_ring_data_t;

typedef struct t8dg_cylinder_ring_data
{
  double              inner_radius;
  double              outer_radius;
  double              height;
  int                 num_trees;
} t8dg_cylinder_ring_data_t;

typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;

t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_linear (int dim);

/**Create a new coarse geometry for the linear 1D case*/
t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_1D_linear ();
/**Create a new coarse geometry for the linear 2D case*/
t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_2D_linear ();
t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_3D_linear ();
t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_2D_circle_ring (double inner_radius, double outer_radius, int num_trees);
t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_3D_cylinder_ring (double inner_radius, double outer_radius, int num_trees, double height);

t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_arg (int geometry_arg);

/**Destroy a coarse geometry*/
void                t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_t ** pgeometry);

double              t8dg_coarse_geometry_calculate_sqrt_gram_determinant
  (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest, const t8_gloidx_t iglobaltree, const double coarse_vertex[3]);

double              t8dg_coarse_geometry_calculate_sqrt_face_gram_determinant
  (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
   const t8_gloidx_t iglobaltree, const int iface, const double coarse_vertex[3]);

void                t8dg_coarse_geometry_apply (const t8dg_coarse_geometry_t * coarse_geometry,
                                                const t8_forest_t forest,
                                                const t8_gloidx_t iglobaltree, const double coarse_vertex[3], double image_vertex[3]);

void                t8dg_coarse_geometry_calculate_gradient_tangential_vector (const t8dg_coarse_geometry_t * coarse_geometry,
                                                                               const t8_forest_t forest, const t8_gloidx_t iglobaltree,
                                                                               const double coarse_vertex[3],
                                                                               const double coarse_tangential_vector[3],
                                                                               double transformed_gradient_tangential_vector[3]);

void                t8dg_coarse_geometry_transform_normal_vector
  (const t8dg_coarse_geometry_t * coarse_geometry, const t8_forest_t forest,
   const t8_gloidx_t iglobaltree, const double coarse_vertex[3], const double coarse_normal_vector[3], double image_normal_vector[3]);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_COARSE_GEOMETRY_H_ */
