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
#include "t8dg_global_precomputed_values.h"
#include "t8dg_flux.h"
#include "t8dg_local_precomputed_values.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_mortar t8dg_mortar_t;

sc_array_t         *t8dg_mortar_get_flux (t8dg_mortar_t * mortar);

t8dg_mortar_t      *t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface);

void                t8dg_mortar_destroy (t8dg_mortar_t ** pmortar);

void                t8dg_mortar_get_idata_iface (t8dg_mortar_t * mortar, t8_locidx_t * pidata, int *piface, int side);

void                t8dg_mortar_invalidate (t8dg_mortar_t * mortar);

int                 t8dg_mortar_is_valid (t8dg_mortar_t * mortar);

int                 t8dg_mortar_get_side (t8dg_mortar_t * mortar, t8_locidx_t idata);

void                t8dg_mortar_fill (t8dg_mortar_t * mortar,
                                      sc_array_t * element_dof_values,
                                      t8dg_timestepping_data_t * time_data,
                                      t8dg_global_precomputed_values_t * global_values,
                                      t8dg_local_precomputed_values_t * local_values,
                                      t8dg_linear_flux_t * flux, t8dg_linear_numerical_flux_fn numerical_flux_fn);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_MORTAR_H_ */
