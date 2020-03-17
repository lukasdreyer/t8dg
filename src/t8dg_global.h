#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <sc_containers.h>
#include <sc_dmatrix.h>
#include <t8_forest.h>
#include <example/common/t8_example_common.h>



#define MAX_DIM 1
#define MAX_FACES (2*MAX_DIM)
#if 0
#define MAX_SUBFACES (1<<(MAX_DIM - 1))
#endif

typedef double      (*t8dg_scalar_function_MAX_DIMd_fn) (const double x[MAX_DIM]);
//#include "t8dg_geometry.h"


typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;


typedef void (*t8dg_matrix_application) (sc_array_t *dest, const sc_array_t *src, const void *application_data);
typedef void (*t8dg_time_matrix_application) (sc_array_t *dest, const sc_array_t *src, double t, const void *application_data);

typedef struct t8dg_quadrature t8dg_quadrature_t;
typedef struct t8dg_functionbasis t8dg_functionbasis_t;
typedef double (*t8dg_numerical_flux)(const double u_minus,const double u_plus, const void *application_data);

typedef struct t8dg_1D_advect_problem t8dg_1D_advect_problem_t;
typedef struct t8dg_1D_advect_element_precomputed_values t8dg_1D_advect_element_precomputed_values_t;
typedef struct t8dg_1D_advect_advance_element_data t8dg_1D_advect_advance_element_data_t;

struct t8dg_1D_advect_problem
{
  t8dg_scalar_function_MAX_DIMd_fn u_0;
  double flow_velocity;

  t8_forest_t			forest;

  sc_array_t			*element_values;		/*those get partitioned*/
  sc_array_t			*jacobian_invers_at_quad;	/* d*d*Q */
  sc_array_t			*element_trafo_quad_weight;	/* Q */
  sc_array_t			*face_trafo_quad_weight[MAX_FACES]; /* FQ */

  sc_array_t			*dof_values;			/*those get ghosted*/
  sc_array_t			*dof_new;			/*those are only needed locally to save the result of a rungekutta timestep*/
  sc_array_t			*advance_element_data;		/*those need to be recalculated for each time step, remain processor local*/

  int 				level;				/*uniform refinement level*/
  int 				dim;

  double			delta_t;/*time step*/
  double 			t;/*current time*/
  double 			T;/*end time*/
  double			cfl;

  t8dg_quadrature_t			*quadrature;
  t8dg_functionbasis_t			*functionbasis;
  t8dg_coarse_geometry_t		*coarse_geometry;
  t8dg_matrix_application		vandermonde, vandermonde_transpose;
  t8dg_matrix_application		face_vandermonde, face_vandermonde_transpose;
  t8dg_numerical_flux			numerical_flux;

  t8dg_time_matrix_application		evolution_matrix;

  sc_MPI_Comm				comm;
};


/*maybe obsolete, since t8code gives all these values*/
struct t8dg_1D_advect_element_precomputed_values
{
  int dim;
  int level;
  int num_faces;
  double diameter;
  double scaling_factor;
  double translation_vector[MAX_DIM];
  int idx_rotation_reflection;
  /* Jacobian, quadrature/trafo weights as seperate sc_array because of size dependency on number of quadrature points*/
};



struct t8dg_1D_advect_advance_element_data
{
  sc_array_t 			*fluxes[MAX_FACES];
  /*sc_array_t			*subfluxes[MAX_FACES][MAX_SUBFACES]; */
  t8_locidx_t			same_level_neighbour[MAX_FACES]; /**< Indices of the neighbor elements */
};




struct t8dg_quadrature
{
  int 			tensorflag;
  t8dg_quadrature_t 	*tensor1;
  t8dg_quadrature_t 	*tensor2;
  int			number_of_vertices;
  int			number_of_faces;
  int			number_of_facevertices[MAX_FACES];
  sc_array_t		*vertices;
  sc_array_t		*weights;
  /* how do i get to face-vertices/weights */
};

struct t8dg_functionbasis
{
  int number_of_dof;
  t8dg_matrix_application		directional_derivative_matrix;
};



void identity_matrix (sc_array_t *dest, const sc_array_t *src, const void *application_data);

void face_vandermonde_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data);
void face_vandermonde_transpose_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data);
void directional_derivative_1D_LGL2_matrix(sc_array_t *dest, const sc_array_t *src, const void *application_data);
double upwind_flux_1D(const double u_minus,const double u_plus, const void *application_data);

t8dg_quadrature_t * t8dg_1D_LGL_quadrature(int number_of_LGL);
void t8dg_quadrature_destroy(t8dg_quadrature_t **pquadrature);


t8dg_functionbasis_t * t8dg_1D_LGL_functionbasis(int number_of_LGL);
void t8dg_functionbasis_destroy(t8dg_functionbasis_t **pfunctionbasis);



#endif /* SRC_GLOBAL_H_ */
