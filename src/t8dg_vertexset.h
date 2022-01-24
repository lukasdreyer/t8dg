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



#ifndef SRC_T8DG_VERTEXSET_H_
#define SRC_T8DG_VERTEXSET_H_

#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

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

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_VERTEXSET_H_ */
