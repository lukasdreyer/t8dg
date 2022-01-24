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



#ifndef SRC_T8DG_MORTAR_H_
#define SRC_T8DG_MORTAR_H_

#include <t8.h>
#include <sc_containers.h>
#include <t8_forest.h>

#include "t8dg_timestepping.h"
#include "t8dg_global_values.h"
#include "t8dg_flux.h"
#include "t8dg_local_values.h"
#include "t8dg_geometry.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_mortar t8dg_mortar_t;
typedef struct t8dg_mortar_array t8dg_mortar_array_t;

void                t8dg_mortar_array_calculate_linear_flux3D (t8dg_mortar_array_t * mortar_array, t8dg_dof_values_t * dof_values,
                                                               t8dg_linear_flux3D_fn linear_flux, void *flux_data,
                                                               t8dg_numerical_linear_flux3D_fn numerical_flux, void *numerical_flux_data,
                                                               double time);

void                t8dg_mortar_array_calculate_flux_dof1D (t8dg_mortar_array_t * mortar_array, t8dg_dof_values_t * dof_values,
                                                            int icomp,
                                                            t8dg_numerical_flux1D_fn numerical_flux, void *numerical_flux_data,
                                                            double time);

t8dg_mortar_array_t *t8dg_mortar_array_new_empty (t8_forest_t forest, t8dg_local_values_t * local_values);

void                t8dg_mortar_array_invalidate_all (t8dg_mortar_array_t * mortar_array);

void                t8dg_mortar_array_destroy (t8dg_mortar_array_t ** pmortar_array);

double              t8dg_mortar_array_get_ghost_exchange_time (t8dg_mortar_array_t * mortar_array);

/*TODO!!!!!*/
sc_array_t         *t8dg_mortar_array_get_oriented_flux (t8dg_mortar_array_t * mortar_array, t8_locidx_t idata, int iface);

void
 
 
 
 
 
 
 
 t8dg_mortar_array_apply_element_boundary_integral (t8dg_mortar_array_t * mortar_array,
                                                    t8_locidx_t itree, t8_locidx_t ielement, sc_array_t * element_result_dof);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_MORTAR_H_ */
