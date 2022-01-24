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



#ifndef SRC_T8DG_TENSOR_H_
#define SRC_T8DG_TENSOR_H_

#include <sc_containers.h>
#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

void                t8dg_tensor_array_extract_vector (sc_array_t * tensor_array, const int ivector, const int stride, sc_array_t * vector);

void                t8dg_tensor_array_inject_vector (sc_array_t * vector, const int ivector, const int stride, sc_array_t * tensor_array);

/*deprecated?*/
int                 t8dg_tensor_mult_other_lengths (const int num_tensor, const int tensor_length[DIM3], const int itensor);

void                t8dg_tensor_transform_tensoridx (const int idx, const int tensordims[DIM3], int itensor[DIM3]);     /*MACRO? */

t8_eclass_t         t8dg_tensor_eclass (t8_eclass_t eclass_tensor1, t8_eclass_t eclass_tensor2);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_TENSOR_H_ */
