#ifndef SRC_T8DG_CMESH_H_
#define SRC_T8DG_CMESH_H_

#include <t8dg.h>

T8DG_EXTERN_C_BEGIN ();

int                 t8dg_cmesh_dim (int icmesh);

t8_cmesh_t          t8dg_cmesh_new_arg (int icmesh, sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_periodic_diagonal_line_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_half_moebius_more_trees (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_more_trees_different_size (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_moebius (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_square_tilted (sc_MPI_Comm comm);

t8_cmesh_t          t8dg_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm, int num_trees);

t8_cmesh_t          t8dg_cmesh_new_circle_ring (sc_MPI_Comm comm, double inner_radius, double outer_radius, int num_trees);

T8DG_EXTERN_C_END ();

#endif
