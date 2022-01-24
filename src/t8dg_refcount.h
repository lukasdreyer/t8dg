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



#ifndef SRC_T8DG_REFCOUNT_H_
#define SRC_T8DG_REFCOUNT_H_

#include "t8dg.h"
#include <sc_refcount.h>

T8DG_EXTERN_C_BEGIN ();

/** We can reuse the reference counter type from libsc. */
typedef sc_refcount_t t8dg_refcount_t;

/** Initialize a reference counter to 1.
 * It is legal if its status prior to this call is undefined.
 * \param [out] rc          The reference counter is set to one by this call.
 */
void                t8dg_refcount_init (t8dg_refcount_t * rc);

/** Create a new reference counter with count initialized to 1.
 * Equivalent to calling t8_refcount_init on a newly allocated refcount_t.
 * It is mandatory to free this with \ref t8_refcount_destroy.
 * \return An allocated reference counter whose count has been set to one.
 */
t8dg_refcount_t    *t8dg_refcount_new (void);

/** Destroy a reference counter that we allocated with \ref t8_refcount_new.
 * Its reference count must have decreased to zero.
 * \param [in,out] rc       Allocated, formerly valid reference counter.
 */
void                t8dg_refcount_destroy (t8dg_refcount_t * rc);

/** Increase the reference count by one.
 * It is not necessary to duplicate this functionality as a function. */
#define t8dg_refcount_ref(rc) sc_refcount_ref(rc)

/** Decrease the reference count by one.
 * The count must be greater zero on input.  If the reference count reaches
 * zero, which is indicated by the return value, the counter may NOT be used
 * furter with \ref t8_refcount_ref or \see t8_refcount_unref.  It IS legal to
 * query it with \ref t8_refcount_is_active and \ref t8_refcount_is_last and to
 * repurpose it later by calling \ref t8_refcount_init.
 * It is not necessary to duplicate this functionality as a function. */
#define t8dg_refcount_unref(rc) sc_refcount_unref(rc)

/** Query wether a reference counter is has a positive value. */
#define t8dg_refcount_is_active(rc) sc_refcount_is_active(rc)

/** Query wether a reference counter has value one. */
#define t8dg_refcount_is_last(rc) sc_refcount_is_last(rc)

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_REFCOUNT_H_ */
