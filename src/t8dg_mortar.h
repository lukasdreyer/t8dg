/*
 * t8dg_mortar.h
 *
 *  Created on: Apr 3, 2020
 *      Author: lukas
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

/*TODO!!!!!*/
sc_array_t         *t8dg_mortar_array_get_oriented_flux (t8dg_mortar_array_t * mortar_array, t8_locidx_t idata, int iface);

void
 
 
 
 
 
 
 
 t8dg_mortar_array_apply_element_boundary_integral (t8dg_mortar_array_t * mortar_array,
                                                    t8_locidx_t itree, t8_locidx_t ielement, sc_array_t * element_result_dof);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_MORTAR_H_ */
