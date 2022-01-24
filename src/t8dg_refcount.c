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

#include "t8dg_refcount.h"

void
t8dg_refcount_init (t8dg_refcount_t * rc)
{
  sc_refcount_init (rc, t8_get_package_id ());
}

t8dg_refcount_t    *
t8dg_refcount_new (void)
{
  t8dg_refcount_t    *rc;

  rc = T8_ALLOC (t8dg_refcount_t, 1);
  t8dg_refcount_init (rc);

  return rc;
}

void
t8dg_refcount_destroy (t8dg_refcount_t * rc)
{
  T8DG_ASSERT (!sc_refcount_is_active (rc));
  T8DG_FREE (rc);
}
