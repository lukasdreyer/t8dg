/*
 * t8dg_mptrac.h
 *
 *  Created on: October 15, 2021
 *      Author: Johannes Holke
 */

/** @file t8dg_mptrac.h */
#ifndef SRC_T8DG_MPTRAC_H_
#define SRC_T8DG_MPTRAC_H_

#include <t8/example/mptrac/t8_mptrac_interpolate.h>
#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

t8_mptrac_context_t* t8dg_mptrac_setup (const char *nc_filename);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_MPTRAC_H_ */
