/*
 * t8dg_quadrature.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_QUADRATURE_H_
#define SRC_T8DG_QUADRATURE_H_

#include "t8dg.h"               /*for scalar_fn */
#include "t8dg_vertexset.h"
#include <sc_containers.h>

/** Index used to enumerate quadrature points*/
typedef int         t8dg_quad_idx_t;

typedef enum t8dg_quadrature_type
{
  T8DG_QUAD_LGL,
  T8DG_QUAD_GL,
  T8_QUAD_UNKNOWN
} t8dg_quadrature_type_t;

/**Opaque handle typedef for quadrature*/
typedef struct t8dg_quadrature t8dg_quadrature_t;

/** Returns a pointer to an array element indexed by a t8dg_quad_idx_t.
 * \param [in] index needs to be in [0]..[elem_count-1].
 */

void               *t8dg_sc_array_index_quadidx (const sc_array_t * array, t8dg_quad_idx_t iz);

/*TODO: document
 *
 * */
t8dg_quadrature_t  *t8dg_quadrature_new_vertexset (t8dg_vertexset_t * vertexset);

t8dg_quadrature_t  *t8dg_quadrature_new_tensor (int num_tensor, t8dg_quadrature_t * quad_tensors[3]);

t8dg_quadrature_t  *t8dg_quadrature_new_hypercube (int dim, t8dg_vertexset_t * vertexset1D);

/*TODO: document
 *
 * */
void                t8dg_quadrature_destroy (t8dg_quadrature_t ** pquadrature);

void                t8dg_quadrature_reset (t8dg_quadrature_t ** pquadrature);

void                t8dg_quadrature_ref (t8dg_quadrature_t * quadrature);

void                t8dg_quadrature_unref (t8dg_quadrature_t ** pquadrature);

int                 t8dg_quadrature_get_num_faces (const t8dg_quadrature_t * quadrature);

/** Returns the number of element quadrature vertices of a quadrature
 * \param [in] quadrature            		quadrature
 * \return 					number of element quadrature vertices
 */
t8dg_quad_idx_t     t8dg_quadrature_get_num_element_vertices (const t8dg_quadrature_t * quadrature);

/** Returns the number of face quadrature vertices of a quadrature at face iface
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		index of the face
 * \return 					number of element quadrature vertices
 */
t8dg_quad_idx_t     t8dg_quadrature_get_num_face_vertices (const t8dg_quadrature_t * quadrature, int iface);

/** Returns the 3d - coordinates of the iquad^th element quadrature-vertex
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] quadrature            		quadrature
 * \param [in] iquad            		Index of the element quadrature point
*/
void                t8dg_quadrature_get_element_vertex (double vertex[3],
                                                        const t8dg_quadrature_t * quadrature, const t8dg_quad_idx_t iquad);

/** Returns the quadrature weight of an element vertex. (TODO: sum equals 1 or VOL(element)?
 * \param [in] quadrature            		quadrature
 * \param [in] iquad            		index of the element quadrature vertex
 * \return 					quadrature weight of the vertex
 */
double              t8dg_quadrature_get_element_weight (const t8dg_quadrature_t * quadrature, const t8dg_quad_idx_t iquad);

/** Returns the 3d - coordinates of the iquad^th face quadrature-vertex at face iface
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		Face index
 * \param [in] iquad            		Index of the element quadrature point
*/
void                t8dg_quadrature_get_face_vertex (double vertex[3],
                                                     const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t iquad);

/** Returns the quadrature weight of a face vertex. (TODO: sum equals 1 or VOL(element)?
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		index of the face
 * \param [in] iquad            		index of the face quadrature vertex
 * \return 					quadrature weight of the vertex
 */
double              t8dg_quadrature_get_face_weight (const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t iquad);

t8dg_quadrature_type_t t8dg_quadrature_get_type (const t8dg_quadrature_t * quadrature);

t8_eclass_t         t8dg_quadrature_get_eclass (const t8dg_quadrature_t * quadrature);

int                 t8dg_quadrature_get_dim (const t8dg_quadrature_t * quadrature);

t8dg_quad_idx_t
    t8dg_quadrature_lgl_facequadix_lookup (const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t ifacequad);

double              t8dg_quadrature_integrate_reference_element (t8dg_quadrature_t * quadrature, t8dg_scalar_function_3d_fn integrand_fn,
                                                                 void *integrand_data);

#endif /* SRC_T8DG_QUADRATURE_H_ */
