/** \file t8dg_timestepping_dirk_coefficients.h
* This file keeps the coefficients from the Butcher Tableau of the DIRK(2,2) and DIRK(3,3) methods (=DIRK(order, stages)) 
* They descirbe an implicit RKV which is used in the time integration */

#ifndef SRC_T8DG_TIMESTEPPING_DIRK_COEFFICIENTS_H_
#define SRC_T8DG_TIMESTEPPING_DIRK_COEFFICIENTS_H_

#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

#if T8_WITH_PETSC
/* The resulting linear systems of the DIRK methods which has to be solved in each RKV step in order to obtain the intermediate/final solution is solved with the KSP implementation of PETSc. These coefficients are not needed otherwise, therefore they are surrunded by the T8_WITH_PETSC tags */

/* Number of the order/steps */
extern const int    dirk22_num_order_stages;
/* alpha coefficient of DIRK(2,2) */
extern const double dirk22_alpha;

/* Coefficient array of DIRK(2,2) */
extern double       dirk22_a_coeff[2][2];
extern double       dirk22_b_coeff[2];
extern double       dirk22_c_coeff[2];

/* Number of the order/steps */
extern const int    dirk33_num_order_stages;
/* Values of the Butcher Tableau of DIRK(3,3) */
extern const double dirk33_alpha;
extern const double dirk33_beta;
extern const double dirk33_b1;
extern const double dirk33_b2;
extern const double dirk33_tau1;
extern const double dirk33_tau2;

/* Coefficient arrays of DIRK(3,3) */
extern double       dirk33_a_coeff[3][3];
extern double       dirk33_b_coeff[3];
extern double       dirk33_c_coeff[3];

double              t8dg_timestepping_dirk_get_coeff_a (int order, int index_i, int index_j);

double              t8dg_timestepping_dirk_get_coeff_b (int order, int index);

double              t8dg_timestepping_dirk_get_coeff_c (int order, int index);

#endif
T8DG_EXTERN_C_END ();
#endif
