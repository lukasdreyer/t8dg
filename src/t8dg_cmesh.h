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

#ifndef SRC_T8DG_CMESH_H_
#define SRC_T8DG_CMESH_H_

#include <t8dg.h>

T8DG_EXTERN_C_BEGIN ();

t8_cmesh_t          t8dg_cmesh_new_arg (int icmesh, int *dim, int *velocity_field_arg, int *geometry_arg, sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_periodic_diagonal_line_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_half_moebius_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_more_trees_different_size (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_moebius (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_tilted (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm, int num_trees);

t8_cmesh_t          t8dg_cmesh_new_circle_ring (sc_MPI_Comm comm, double inner_radius, double outer_radius, int num_trees);

t8_cmesh_t          t8dg_cmesh_new_square_half_periodic (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_cylinder_ring_periodic (sc_MPI_Comm comm, double inner_radius, double outer_radius, int num_trees,
                                                           double height);

T8DG_EXTERN_C_END ();

#endif
