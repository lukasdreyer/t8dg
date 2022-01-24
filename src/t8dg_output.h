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

#ifndef SRC_T8DG_OUTPUT_H_
#define SRC_T8DG_OUTPUT_H_
#include "t8dg.h"
#include "t8dg_dof.h"
#include <t8_forest.h>

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_vtk_data
{
  int                 vtk_count;
  const char         *prefix;
  int                 vtk_freq;
} t8dg_vtk_data_t;

t8dg_vtk_data_t    *t8dg_output_vtk_data_new (const char *prefix, int vtk_freq);

void                t8dg_output_vtk_data_destroy (t8dg_vtk_data_t ** p_vtk_data);

void                t8dg_output_write_vtk (t8dg_dof_values_t * dof_values, t8dg_vtk_data_t * vtk_data);

T8DG_EXTERN_C_END ();

#endif
