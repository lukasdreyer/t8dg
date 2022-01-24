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

/** @file t8dg_numerical_flux.h */

#ifndef SRC_T8DG_NUMERICAL_FLUX_H_
#define SRC_T8DG_NUMERICAL_FLUX_H_

#include <t8.h>
#include <t8_forest.h>

#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

typedef double      (*t8dg_numerical_linear_flux3D_fn) (const double u_minus, const double u_plus,
                                                        const double flow_vector[3], const double normal_vector[3],
                                                        const void *numerical_flux_data);

typedef double      (*t8dg_numerical_linear_flux1D_fn) (const double u_minus, const double u_plus,
                                                        const double flow_constant, const double normal_component,
                                                        const void *numerical_flux_data);

typedef double      (*t8dg_numerical_flux1D_fn) (const double u_minus, const double u_plus,
                                                 const double normal_component, const void *numerical_flux_data);

typedef void        (*t8dg_linear_flux1D_fn) (const double x_vec[3], double *flux_velocity, const double t, const void *flux_data);

typedef void        (*t8dg_linear_flux3D_fn) (double x_vec[3], double flux_velocity[3], double t, void *flux_data);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
