/** \file t8dg_timestepping_dirk_coefficients.c
* This file keeps the coefficients from the Butcher Tableau of the DIRK(2,2) and DIRK(3,3) methods (=DIRK(order, stages)) 
* They descirbe an implicit RKV which is used in the time integration */
#include "t8dg_timestepping_dirk_coefficients.h"

#if T8_WITH_PETSC
/* The resulting linear systems of the DIRK methods which has to be solved in each RKV step in order to obtain the intermediate/final solution is solved with the KSP implementation of PETSc. These coefficients are not needed otherwise, therefore they are surrunded by the T8_WITH_PETSC tags */

/* Number of the order/steps */
const int           dirk22_num_order_stages = 2;
/* alpha coefficient of DIRK(2,2) */
const double        dirk22_alpha = 0.2928932188;

/* Coefficient array of DIRK(2,2) */
double              dirk22_a_coeff[2][2] = { {dirk22_alpha, 0.0}, {(1.0 - dirk22_alpha), dirk22_alpha} };
double              dirk22_b_coeff[2] = { (1.0 - dirk22_alpha), dirk22_alpha };
double              dirk22_c_coeff[2] = { dirk22_alpha, 1.0 };

/* Number of the order/steps */
const int           dirk33_num_order_stages = 3;
/* Values of the Butcher Tableau of DIRK(3,3) */
const double        dirk33_alpha = 0.43586652;
const double        dirk33_beta = 0.28206673;
const double        dirk33_b1 = 1.2084966;
const double        dirk33_b2 = -0.64436317;
const double        dirk33_tau1 = 0.43586652;
const double        dirk33_tau2 = 0.71793326;

/* Coefficient arrays of DIRK(3,3) */
double              dirk33_a_coeff[3][3] =
  { {dirk33_alpha, 0.0, 0.0}, {dirk33_beta, dirk33_alpha, 0.0}, {dirk33_b1, dirk33_b2, dirk33_alpha} };
double              dirk33_b_coeff[3] = { dirk33_b1, dirk33_b2, dirk33_alpha };
double              dirk33_c_coeff[3] = { dirk33_tau1, dirk33_tau2, 1.0 };

double
t8dg_timestepping_dirk_get_coeff_a (int order, int index_i, int index_j)
{
  if (order == 2) {
    return dirk22_a_coeff[index_i][index_j];
  }
  else {
    return dirk33_a_coeff[index_i][index_j];
  }
}

double
t8dg_timestepping_dirk_get_coeff_b (int order, int index)
{
  if (order == 2) {
    return dirk22_b_coeff[index];
  }
  else {
    return dirk33_b_coeff[index];
  }
}

double
t8dg_timestepping_dirk_get_coeff_c (int order, int index)
{
  if (order == 2) {
    return dirk22_c_coeff[index];
  }
  else {
    return dirk33_c_coeff[index];
  }
}

#endif
