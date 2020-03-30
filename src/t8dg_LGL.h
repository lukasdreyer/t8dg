#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

/** @file t8dg_LGL.h */

#include <sc_containers.h>
#include <sc_dmatrix.h>
#include <t8_forest.h>

#include "t8dg.h"

/** Index used to enumerate quadrature points*/
typedef int         t8dg_quad_idx_t;

/**Opaque handle typedef for quadrature*/
typedef struct t8dg_LGL_quadrature t8dg_LGL_quadrature_t;
/**Opaque handle typedef for functionbasis*/
typedef struct t8dg_LGL_functionbasis t8dg_LGL_functionbasis_t;

/** Returns a pointer to an array element indexed by a t8dg_quad_idx_t.
 * \param [in] index needs to be in [0]..[elem_count-1].
 */

static inline void *
t8dg_sc_array_index_quadidx (const sc_array_t * array, t8dg_quad_idx_t iz)
{
  T8DG_ASSERT (iz >= 0 && iz < (t8dg_quad_idx_t) array->elem_count);

  return (void *) (array->array + (array->elem_size * iz));
}

/** Creates a new 1D functionbasis and quadrature with the same LGL vertex set, with number_of_LGL many vertices
 *
 * \param [out] pquadrature			Pointer to the quadrature
 * \param [out] pfunctionbasis			Pointer to the functionbasis
 * \param [in] number_of_LGL 			number of 1D LGL vertices
 */

void                t8dg_LGL_quadrature_and_functionbasis_new_1D (t8dg_LGL_quadrature_t ** pquadrature,
                                                                  t8dg_LGL_functionbasis_t ** pfunctionbasis, int number_of_LGL);

/** Destroys a functionbasis and quadrature with the same LGL vertex set
 *
 * \param [in,out] pquadrature 			Pointer to the quadrature
 * \param [in,out] pfunctionbasis 		Pointer to the functionbasis
 */

void                t8dg_LGL_quadrature_and_functionbasis_destroy (t8dg_LGL_quadrature_t ** pquadrature,
                                                                   t8dg_LGL_functionbasis_t ** pfunctionbasis);

/** Applies Lookup-table for the transformation of an element-dof vector to a face-quadrature vector
 * TODO: implement
 * \param [in,out] face_quad_array 		sc_arry of values at face quadrature points
 * \param [in] element_dof_array     		sc_array of element dof values
 * \param [in] iface   				Index of the face
 * \param [in] quadrature 			quadrature
 * \param [in] functionbasis  			functionbasis
 */
void                t8dg_LGL_transform_element_dof_to_face_quad
  (sc_array_t * face_quad_array, const sc_array_t * element_dof_array, int iface, t8dg_LGL_quadrature_t * quadrature,
   t8dg_LGL_functionbasis_t * functionbasis);

/** Applies Lookup-table for the transformation of a face-quadrature vector to an element-dof vector
 * TODO: implement
 * \param [out] element_dof_array   		sc_array of element dof values
 * \param [in] face_quad_array 			sc_arry of values at face quadrature points
 * \param [in] iface 				Index of the face
 * \param [in] quadrature            		quadrature
 * \param [in] functionbasis            	functionbasis
 */
void                t8dg_LGL_transform_face_quad_to_element_dof (sc_array_t * element_dof_array,
                                                                 const sc_array_t * face_quad_array,
                                                                 int iface,
                                                                 t8dg_LGL_quadrature_t * quadrature,
                                                                 t8dg_LGL_functionbasis_t * functionbasis);
/** TODO:implement*/
void                t8dg_LGL_functionbasis_apply_derivative_matrix (sc_array_t * derivative_dof_values,
                                                                    const sc_array_t dof_values, t8dg_LGL_functionbasis_t * functionbasis);

/** TODO:implement*/
void                t8dg_LGL_functionbasis_apply_derivative_matrix_transpose (sc_array_t * dof_values,
                                                                              sc_array_t * derivative_dof_values,
                                                                              t8dg_LGL_functionbasis_t * functionbasis);

/** Returns the number of degrees of Freedom of a functionbasis
 * \param [in] functionbasis            	functionbasis
 * \return 					number of degrees of freedom
 */

t8dg_quad_idx_t     t8dg_LGL_functionbasis_get_num_dof (t8dg_LGL_functionbasis_t * functionbasis);

/** Returns the 3d - coordinates of the vertex corresponding to the idof^th nodal basis function
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] functionbasis            	functionbasis
 * \param [in] idof            			Index of the degree of freedom
*/
void                t8dg_LGL_functionbasis_get_vertex (double vertex[3], t8dg_LGL_functionbasis_t * functionbasis, int idof);

/** Returns the number of faces of the reference element
 * \param [in] quadrature            		quadrature
 * \return 					number of faces
 */
int                 t8dg_LGL_quadrature_get_num_faces (t8dg_LGL_quadrature_t * quadrature);

/** Returns the number of element quadrature vertices of a quadrature
 * \param [in] quadrature            		quadrature
 * \return 					number of element quadrature vertices
 */
t8dg_quad_idx_t     t8dg_LGL_quadrature_get_num_element_vertices (t8dg_LGL_quadrature_t * quadrature);

/** Returns the number of face quadrature vertices of a quadrature at face iface
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		index of the face
 * \return 					number of element quadrature vertices
 */
t8dg_quad_idx_t     t8dg_LGL_quadrature_get_num_face_vertices (t8dg_LGL_quadrature_t * quadrature, int iface);

/** Returns the 3d - coordinates of the iquad^th element quadrature-vertex
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] quadrature            		quadrature
 * \param [in] iquad            		Index of the element quadrature point
*/
void                t8dg_LGL_quadrature_get_element_vertex (double vertex[3], t8dg_LGL_quadrature_t * quadrature, t8dg_quad_idx_t iquad);

/** Returns the quadrature weight of an element vertex. (TODO: sum equals 1 or VOL(element)?
 * \param [in] quadrature            		quadrature
 * \param [in] iquad            		index of the element quadrature vertex
 * \return 					quadrature weight of the vertex
 */
double              t8dg_LGL_quadrature_get_element_weight (t8dg_LGL_quadrature_t * quadrature, t8dg_quad_idx_t iquad);

/** Returns the 3d - coordinates of the iquad^th face quadrature-vertex at face iface
 * \param [out] vertex            		3d vector to be filled with coordinates
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		Face index
 * \param [in] iquad            		Index of the element quadrature point
*/
void                t8dg_LGL_quadrature_get_face_vertex (double vertex[3], t8dg_LGL_quadrature_t * quadrature,
                                                         int iface, t8dg_quad_idx_t iquad);

/** Returns the quadrature weight of a face vertex. (TODO: sum equals 1 or VOL(element)?
 * \param [in] quadrature            		quadrature
 * \param [in] iface            		index of the face
 * \param [in] iquad            		index of the face quadrature vertex
 * \return 					quadrature weight of the vertex
 */
double              t8dg_LGL_quadrature_get_face_weight (t8dg_LGL_quadrature_t * quadrature, int iface, t8dg_quad_idx_t iquad);

#endif /* SRC_GLOBAL_H_ */
