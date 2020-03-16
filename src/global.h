#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <sc_containers.h>


#define MAX_FACES 2


typedef void (*t8dg_matrix_application) (sc_array_t *dest, const sc_array_t *src, const void *application_data);
typedef void (*t8dg_time_matrix_application) (sc_array_t *dest, const sc_array_t *src, double t, const void *application_data);

typedef struct t8dg_quadrature t8dg_quadrature_t;
typedef struct t8dg_functionbasis t8dg_functionbasis_t;
typedef double (*t8dg_numerical_flux)(const double u_minus,const double u_plus, const void *application_data);

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
