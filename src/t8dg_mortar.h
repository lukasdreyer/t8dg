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
#include "t8dg_geometry.h"

typedef struct t8dg_mortar t8dg_mortar_t;
typedef struct t8dg_mortar_array t8dg_mortar_array_t;

typedef struct t8dg_mortar_fill_data
{
  const t8dg_global_precomputed_values_t *global_values;
  const t8dg_local_precomputed_values_t *local_values;
  t8dg_geometry_transformation_data_t *geometry_data;
  const t8dg_flux_t  *flux;
  sc_array_t         *dof_values;
  const double        time;
} t8dg_mortar_fill_data_t;

//void                t8dg_mortar_fill (t8dg_mortar_t * mortar, t8dg_mortar_fill_data_t * mortar_fill_data);

void                t8dg_mortar_array_fill (t8dg_mortar_array_t * mortar_array, t8dg_mortar_fill_data_t * mortar_fill_data);

t8dg_mortar_array_t *t8dg_mortar_array_new_empty (t8_forest_t forest, int num_faces);

void                t8dg_mortar_array_invalidate_all (t8dg_mortar_array_t * mortar_array);

void                t8dg_mortar_array_destroy (t8dg_mortar_array_t ** pmortar_array);

/*TODO!!!!!*/
sc_array_t         *t8dg_mortar_array_get_oriented_flux (t8dg_mortar_array_t * mortar_array, t8_locidx_t idata, int iface);

#endif /* SRC_T8DG_MORTAR_H_ */
