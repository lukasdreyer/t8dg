/*
 * sc_array_access.h
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */

#ifndef SRC_SC_ARRAY_ACCESS_H_
#define SRC_SC_ARRAY_ACCESS_H_

#include "solver.hxx"

void t8dg_element_set_dofs_initial(t8dg_1D_advect_problem_t *problem,t8_locidx_t ielement);
void t8dg_element_set_jacobian(t8dg_1D_advect_problem_t *problem,t8_locidx_t ielement);
void t8dg_element_set_trafo_weights(t8dg_1D_advect_problem_t *problem,t8_locidx_t ielement);

#endif /* SRC_SC_ARRAY_ACCESS_H_ */
