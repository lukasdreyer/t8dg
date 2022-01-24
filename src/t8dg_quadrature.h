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



#ifndef SRC_T8DG_QUADRATURE_H_
#define SRC_T8DG_QUADRATURE_H_

#include "t8dg.h"
#include "t8dg_vertexset.h"
#include <sc_containers.h>
#include "t8dg_common.h"

T8DG_EXTERN_C_BEGIN ();

typedef enum t8dg_quadrature_type
{
  T8DG_QUAD_LGL
} t8dg_quadrature_type_t;

/**Opaque handle typedef for quadrature*/
typedef struct t8dg_quadrature t8dg_quadrature_t;

/*TODO: document
 *
 * */
t8dg_quadrature_t  *t8dg_quadrature_new_vertexset (t8dg_vertexset_t * vertexset, int create_face_quadrature);

t8dg_quadrature_t  *t8dg_quadrature_new_tensor (t8dg_quadrature_t * tensor_first_quadrature, t8dg_quadrature_t * tensor_second_quadrature,
                                                int create_face_quadrature);

t8dg_quadrature_t  *t8dg_quadrature_new_hypercube (int dim, t8dg_vertexset_t * vertexset1D, int create_face_quadrature);

/*TODO: document
 *
 * */
void                t8dg_quadrature_destroy (t8dg_quadrature_t ** pquadrature);

void                t8dg_quadrature_reset (t8dg_quadrature_t ** pquadrature);

void                t8dg_quadrature_ref (t8dg_quadrature_t * quadrature);

void                t8dg_quadrature_unref (t8dg_quadrature_t ** pquadrature);

int                 t8dg_quadrature_get_num_face_quadrature (const t8dg_quadrature_t * quadrature);

/** Returns the number of element quadrature vertices of a quadrature
 * \param [in] quadrature            		quadrature
 * \return 					number of element quadrature vertices
 */
int                 t8dg_quadrature_get_num_element_vertices (const t8dg_quadrature_t * quadrature);

/** Returns the number of face quadrature vertices of a quadrature at face iface
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		index of the face
 * \return 					number of element quadrature vertices
 */
int                 t8dg_quadrature_get_num_face_vertices (const t8dg_quadrature_t * quadrature, int iface);

/** Returns the 3d - coordinates of the iquad^th element quadrature-vertex
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] quadrature            		quadrature
 * \param [in] iquad            		Index of the element quadrature point
*/
void                t8dg_quadrature_get_element_vertex (const t8dg_quadrature_t * quadrature, const int iquad, double vertex[3]);

/** Returns the quadrature weight of an element vertex. (TODO: sum equals 1 or VOL(element)?
 * \param [in] quadrature            		quadrature
 * \param [in] iquad            		index of the element quadrature vertex
 * \return 					quadrature weight of the vertex
 */
double              t8dg_quadrature_get_element_weight (const t8dg_quadrature_t * quadrature, const int iquad);

/** Returns the 3d - coordinates of the iquad^th face quadrature-vertex at face iface
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		Face index
 * \param [in] iquad            		Index of the element quadrature point
*/

 /*TODO*/
  void      t8dg_quadrature_get_face_vertex (const t8dg_quadrature_t * quadrature, const int iface, const int iquad, double vertex[3]);

/** Returns the quadrature weight of a face vertex. (TODO: sum equals 1 or VOL(element)?
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		index of the face
 * \param [in] iquad            		index of the face quadrature vertex
 * \return 					quadrature weight of the vertex
 */
double              t8dg_quadrature_get_face_weight (const t8dg_quadrature_t * quadrature, const int iface, const int iquad);

#if 0
t8_eclass_t         t8dg_quadrature_get_face_eclass (const t8dg_quadrature_t * quadrature, const int iface);
#endif

t8dg_quadrature_type_t t8dg_quadrature_get_type (const t8dg_quadrature_t * quadrature);

t8_eclass_t         t8dg_quadrature_get_eclass (const t8dg_quadrature_t * quadrature);

int                 t8dg_quadrature_get_dim (const t8dg_quadrature_t * quadrature);

double              t8dg_quadrature_integrate_reference_element (t8dg_quadrature_t * quadrature, t8dg_scalar_function_3d_fn integrand_fn,
                                                                 void *integrand_data);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_QUADRATURE_H_ */
