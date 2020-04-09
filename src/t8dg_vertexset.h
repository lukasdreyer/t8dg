/*
 * t8dg_vertexset.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_VERTEXSET_H_
#define SRC_T8DG_VERTEXSET_H_

typedef enum t8dg_vertexset_type
{
  T8DG_LGL,
  T8DG_GL,
} t8dg_vertexset_type_t;

typedef struct t8dg_vertexset t8dg_vertexset_t;

t8dg_vertexset_type_t t8dg_vertexset_get_type (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_dim (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_num_element_vertices (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_num_faces (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_num_face_vertices (const t8dg_vertexset_t * vertexset, const int iface);

int                 t8dg_vertexset_get_LGL_facevertex_element_index (const t8dg_vertexset_t * vertexset, const int iface,
                                                                     const int ifacevertex);

void                t8dg_vertexset_get_vertex (double reference_vertex[3], const t8dg_vertexset_t * vertexset, const int ivertex);

t8dg_vertexset_t   *t8dg_vertexset_new_1D_LGL (const int number_of_LGL_vertices);

void                t8dg_vertexset_destroy (t8dg_vertexset_t ** pvertexset);

void                t8dg_vertexset_ref (t8dg_vertexset_t * vertexset);

void                t8dg_vertexset_unref (t8dg_vertexset_t ** vertexset);

#endif /* SRC_T8DG_VERTEXSET_H_ */
