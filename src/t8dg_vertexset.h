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
  T8DG_VERT_LGL
} t8dg_vertexset_type_t;

typedef struct t8dg_vertexset t8dg_vertexset_t;

t8dg_vertexset_type_t t8dg_vertexset_get_type (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_eclass_dim (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_num_vertices (const t8dg_vertexset_t * vertexset);

t8_eclass_t         t8dg_vertexset_get_eclass (const t8dg_vertexset_t * vertexset);

int                 t8dg_vertexset_get_embedded_dim (const t8dg_vertexset_t * vertexset);

double              t8dg_vertexset_get_first_coordinate (const t8dg_vertexset_t * vertexset, const int ivertex);

void                t8dg_vertexset_fill_vertex3D (const t8dg_vertexset_t * vertexset, const int ivertex, const int startdim,
                                                  double reference_vertex[3]);

t8dg_vertexset_t   *t8dg_vertexset_new_1D_LGL (const int number_of_LGL_vertices);

t8dg_vertexset_t   *t8dg_vertexset_new_childvertexset (const t8dg_vertexset_t * vertexset, int ichild);

t8dg_vertexset_t   *t8dg_vertexset_new_lgl_facevertexset (const t8dg_vertexset_t * vertexset, int iface);

void                t8dg_vertexset_destroy (t8dg_vertexset_t ** pvertexset);

void                t8dg_vertexset_ref (t8dg_vertexset_t * vertexset);

void                t8dg_vertexset_unref (t8dg_vertexset_t ** vertexset);

#endif /* SRC_T8DG_VERTEXSET_H_ */
