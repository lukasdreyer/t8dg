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



#ifndef SRC_T8DG_GEOMETRY_H_
#define SRC_T8DG_GEOMETRY_H_

#include <t8.h>
#include <t8_forest.h>
#include "t8dg_coarse_geometry.h"

T8DG_EXTERN_C_BEGIN ();

void                t8dg_geometry_transform_reference_vertex_to_image_vertex
  (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest, const t8_gloidx_t itree, const t8_element_t * element,
   const double reference_vertex[3], double image_vertex[3]);

double              t8dg_geometry_calculate_sqrt_gram_determinant
  (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest, const t8_gloidx_t itree, const t8_element_t * element,
   const double reference_vertex[3]);

void                t8dg_geometry_calculate_transformed_gradient_tangential_vector
  (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest, const t8_gloidx_t itree, const t8_element_t * element,
   const double reference_vertex[3], const double reference_tangential_vector[3], double transformed_gradient_tangential_vector[3]);

void                t8dg_geometry_calculate_normal_vector
  (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest, const t8_gloidx_t itree, const t8_element_t * element,
   const int iface, const double reference_vertex[3], double image_normal_vector[3]);

double              t8dg_geometry_calculate_face_sqrt_gram_determinant
  (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest, const t8_gloidx_t itree, const t8_element_t * element,
   const int iface, const double reference_vertex[3]);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_GEOMETRY_H_ */
