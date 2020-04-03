/*
 * t8dg_mortar.h
 *
 *  Created on: Apr 3, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_MORTAR_H_
#define SRC_T8DG_MORTAR_H_

#include <t8.h>
#include <sc_containers.h>
#include <t8_forest.h>

/** struct used to save the calculated numerical fluxes at quadrature points*/
typedef struct t8dg_mortar
{
  int                 number_face_quadrature_points;      /**< The number of face quadrature points*/
  t8_locidx_t         elem_idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_idata_plus;                    /**< Local index of the element corresponding to u_plus */

  int                 iface_plus, iface_minus;

  /*one value for each quadrature point */
  sc_array_t         *u_minus;                          /**< value of u on elem_minus at face quadrature points */
  sc_array_t         *u_plus;                           /**< value of u on elem_plus at face quadrature points */

  sc_array_t         *fluxes;                           /**< value of (cu)*.n at face quadrature points */
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/

} t8dg_mortar_t;                /*maybe change to opaque handle */

sc_array_t         *t8dg_mortar_get_flux (t8dg_mortar_t * mortar);

t8dg_mortar_t      *t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface);

void                t8dg_mortar_destroy (t8dg_mortar_t ** pmortar);

void                t8dg_mortar_get_idata_iface (t8dg_mortar_t * mortar, t8_locidx_t * pidata, int *piface, int side);

#endif /* SRC_T8DG_MORTAR_H_ */
