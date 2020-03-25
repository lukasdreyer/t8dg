/*
 * t8dg_geometry.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

/** @file t8dg_geometry.h */

#ifndef SRC_T8DG_GEOMETRY_H_
#define SRC_T8DG_GEOMETRY_H_

#include <sc_containers.h>
#include "t8.h"
#include "t8_element.h"

#include "t8dg.h"
#include "t8dg_square_3D_matrix.h"

/** returns the 3x3 jacobian of the geometry function*/
typedef void        (*t8dg_jacobian_fn_3D) (t8dg_square_3D_matrix_t jacobian, const double vertex[DIM3], void *geometry_data);
/** returns the image vertex of the geometry function*/
typedef void        (*t8dg_geometry_fn_3D) (double image_vertex[DIM3], const double vertex[DIM3], void *geometry_data);

/** This information determines the transformation from the fine reference element to the subset of the coarse reference element
 * L(x) = hRx + x_0
 * */
typedef struct t8dg_element_fine_to_coarse_geometry_data
{
  double              scaling_factor;           /**< scaling factor h*/
  double              translation_vector[DIM3]; /**< translation vector x_0*/
  int                 idx_rotation_reflection;  /**< Index to determine which reflection/rotation matrix to use, TODO: change to enum? */
} t8dg_element_fine_to_coarse_geometry_data_t;

/** F_CE and DF_CE for the coarse geometry
 * TODO: How to chose unused dimensions so that DF_CE can be used to calculate the differential on the submanifold*/
typedef struct t8dg_coarse_geometry_3D
{
  t8dg_jacobian_fn_3D jacobian; /**< Jacobian function*/
  t8dg_geometry_fn_3D geometry; /**< Geometry function*/
} t8dg_coarse_geometry_3D_t;

/**TODO: If needed, export to t8_vec*/
void                t8dg_vec_print (const double vec_x[3]);

/*TODO: Document In and Output*/

/**Using the information of t8_code about location and rotation/reflection of the refined element in the coarse element,
 * save these in geometry_data*/
void                t8dg_element_set_geometry_data (t8dg_element_fine_to_coarse_geometry_data_t * geometry_data,
                                                    t8_element_t * element, t8_locidx_t idata, t8_eclass_scheme_c * scheme);

/**Create a new coarse geometry for the linear 1D case*/
t8dg_coarse_geometry_3D_t *t8dg_coarse_geometry_new_1D_linear ();
/**Destroy a coarse geometry*/
void                t8dg_coarse_geometry_destroy (t8dg_coarse_geometry_3D_t ** pgeometry);

/** Apply the transformation from fine reference element to the subset of the coarse reference element to a vertex*/
void                t8dg_fine_to_coarse_geometry (double coarse_element_vertex[DIM3],
                                                  double refined_element_vertex[DIM3],
                                                  t8dg_element_fine_to_coarse_geometry_data_t * element_data);

#endif /* SRC_T8DG_GEOMETRY_H_ */
