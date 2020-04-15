/*
 * t8dg_geometry.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

/** @file t8dg_geometry.h */

#ifndef SRC_T8DG_COARSE_GEOMETRY_H_
#define SRC_T8DG_COARSE_GEOMETRY_H_

#include <sc_containers.h>
#include <t8.h>
#include <t8_forest.h>

#include "t8dg.h"
#include "t8dg_dmatrix.h"

/** returns the 3x3 jacobian of the geometry function*/
typedef void        (*t8dg_jacobian_fn_3D) (const double vertex[DIM3], t8dg_square_3D_matrix_t jacobian, t8_forest_t forest,
                                            t8_locidx_t itree);
/** returns the image vertex of the geometry function*/
typedef void        (*t8dg_geometry_fn_3D) (const double vertex[DIM3], double image_vertex[DIM3], t8_forest_t forest, t8_locidx_t itree);

typedef enum t8dg_coarse_geometry_data
{
  T8DG_TREE_VERTICES
} t8dg_coarse_geometry_data_t;

/** F_CE and DF_CE for the coarse geometry
 * TODO: How to chose unused dimensions so that DF_CE can be used to calculate the differential on the submanifold*/
typedef struct t8dg_coarse_geometry
{
  t8dg_jacobian_fn_3D jacobian; /**< Jacobian function*/
  t8dg_geometry_fn_3D geometry; /**< Geometry function*/
  t8dg_coarse_geometry_data_t data_type;
} t8dg_coarse_geometry_t;

/*TODO: Document In and Output*/

/**Create a new coarse geometry for the linear 1D case*/
t8dg_coarse_geometry_t *t8dg_coarse_geometry_new_1D_linear ();
/**Destroy a coarse geometry*/
void                t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_t ** pgeometry);

#endif /* SRC_T8DG_COARSE_GEOMETRY_H_ */
