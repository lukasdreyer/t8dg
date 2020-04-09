#include "t8dg_refcount.h"

void
t8dg_refcount_init (t8dg_refcount_t * rc)
{
  sc_refcount_init (rc, t8_get_package_id ());
}

t8dg_refcount_t    *
t8_refcount_new (void)
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
