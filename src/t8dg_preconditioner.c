#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include "t8dg_preconditioner.h"

#if T8_WITH_PETSC
#if 0
typedef int         (*t8dg_corase_lvl_adapt_func_t) (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t itree, t8_locidx_t ielement,
                                                     t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

typedef struct
{
  t8dg_advect_diff_problem_t *problem;
  t8_forest_t forest  coarsened;
  t8dg_dof_values_t  *dof_values;
  t8dg_values_t      *dg_values;
  t8dg_adapt_data_t  *adapt_data;
  size_t              num_local_elements;
  PetscInt           *global_indexing;
  t8dg_mortar_array_t *tmp_mortar_coarse_lvl;
  t8_forest_t         problem_forest;
  t8dg_dof_values_t  *dof_values_adapt;
} t8dg_mg_coarse_lvl_t;

typedef struct
{
  t8dg_mg_coarse_lvl_t *coarse_lvl;
  void               *user_data;
  t8dg_dof_values_t  *problem_dofs;
  t8dg_dof_values_t  *problem_dofs_derivation;
  t8dg_time_matrix_application time_derivative_func;
  double              current_point_in_time;
} t8dg_coarse_matrix_ctx_t;

typedef struct
{
  t8dg_advect_diff_problem_t *problem;
  t8dg_mg_coarse_lvl_t *coarse_lvl;
  size_t              num_local_dofs_coarse_grid;
  size_t              num_local_dofs_fine_grid;
  PetscInt           *fine_lvl_global_indexing;
} t8dg_mg_interpolating_ctx_t;
#endif
extern PetscErrorCode MatMult_MF_Coarse_LVL (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_Prolongation (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_Restriction (Mat, Vec, Vec);

#if 0
extern PetscErrorCode CoarseLvlShellPCCreate (t8dg_timestepping_precon_jacobi_ctx_t **, double, double);
extern PetscErrorCode CoarseLvlShellPCSetUp (PC, Mat, Vec, t8dg_dof_values_t *, size_t);
extern PetscErrorCode CoarseLvlShellPCApply (PC, Vec x, Vec y);
extern PetscErrorCode CoarseLvlShellPCDestroy (PC);
#endif

void
t8dg_mg_set_up_two_lvl_precon (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                               t8dg_time_matrix_application time_derivative, t8dg_mg_coarse_lvl_t * coarse_lvl,
                               t8dg_corase_lvl_adapt_func_t coarsening_func, t8dg_mg_interpolating_ctx_t * res_prol_ctx,
                               t8dg_coarse_matrix_ctx_t * cmat_ctx, PetscInt * fine_forest_indexing)
{

  t8dg_debugf ("before call to coarse_lvl_construct\n");
  /* Build a coarse mesh which is used within the multigrid preconditioner */
  t8dg_mg_construct_coarse_lvl (problem, pdof_array, coarse_lvl, coarsening_func);
  t8dg_debugf ("after call to coarse_lvl_construct\n");
  /* Assign some data regarding the coarse matrix context */
  cmat_ctx->coarse_lvl = coarse_lvl;
  cmat_ctx->problem_dofs = *(coarse_lvl->dof_values);
  cmat_ctx->current_point_in_time = t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem));
  cmat_ctx->user_data = (void *) problem;;
  cmat_ctx->problem_dofs_derivation = t8dg_dof_values_duplicate (*(coarse_lvl->dof_values_adapt));
  cmat_ctx->time_derivative_func = time_derivative;

  /* Number of fine grid process-local degrees of freedom */
  res_prol_ctx->num_local_dofs_fine_grid =
    (size_t) (t8dg_dof_get_num_local_elements (*(coarse_lvl->dof_values)) * t8dg_dof_get_max_num_element_dof (*(coarse_lvl->dof_values)));
  t8dg_debugf ("before global indexing\n");
  /* Indexing vector for petsc entries on the initial/fine mesh */
  res_prol_ctx->fine_lvl_global_indexing = fine_forest_indexing;

  /* Number of coarse grid process-local degrees of freedom */
  res_prol_ctx->num_local_dofs_coarse_grid = coarse_lvl->num_local_dofs;

  /* Save the coarse level forest inside the restriction_prolongation context */
  res_prol_ctx->coarse_lvl = coarse_lvl;;

}

void
t8dg_mg_swap_coarse_fine_data (t8dg_mg_coarse_lvl_t * coarse_lvl)
{
  t8dg_dof_values_swap ((coarse_lvl->dof_values), (coarse_lvl->dof_values_adapt));
  t8dg_dof_values_swap (&(coarse_lvl->adapt_data->dof_values), &(coarse_lvl->adapt_data->dof_values_adapt));
  //t8dg_values_swap_to_adapt_data(coarse_lvl->dg_values, &(coarse_lvl->tmp_mortar_coarse_lvl));
}

PetscErrorCode
MatMult_MF_Restriction (Mat A, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_mg_interpolating_ctx_t *appctx;
  t8dg_mg_coarse_lvl_t *coarse_lvl;

  t8dg_debugf ("Restriction got called\n");

  /* Get the mtarix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Get the coarse level mesh */
  coarse_lvl = appctx->coarse_lvl;

  /* Write the entries of the 'in' vector in adapt_data->dof_values,  these are getting adapted/restricted */
  t8dg_precon_write_vec_to_dof (&in, coarse_lvl->adapt_data->dof_values);

  //t8dg_values_mg_allocate_adapt(coarse_lvl->dg_values, &coarse_lvl->forest_coarsened);
#if 0
  /* schon von Anfang an allokiert im ersten aufruf von forest_coarsened */
  /* Allocate space for the adaption/restriction step */
  t8dg_values_allocate_adapt (coarse_lvl->dg_values, coarse_lvl->forest_coarsened);

  /* Allocate inside the advect diff problem dofvalues meaning to hold the adapted/restricted dof_values */
  restricted_dofs = *(t8dg_advect_diff_problem_get_dof_values_adapt (coarse_lvl->problem));
  restricted_dofs = t8dg_dof_values_new (coarse_lvl->forest_coarsened, t8dg_values_get_global_values_array (coarse_lvl->dg_values));
#endif
  /* Therefore, the pointer in adapt_data has to be set to these dof_values_adapt */
  /* They are already allocated during construciton of coarse level */
  //coarse_lvl->adapt_data->dof_values_adapt;

  /* Calling iterate_replace executes the adaption/restriction */
  t8_forest_iterate_replace (coarse_lvl->forest_coarsened, coarse_lvl->problem_forest, t8dg_adapt_replace);

  /* Write the degrees of freedom to a petsc vec */
  t8dg_precon_write_dof_to_vec (coarse_lvl->adapt_data->dof_values_adapt, &out, coarse_lvl->global_indexing);

  /* Ununsed resources need to be cleaned up */
  /* The question is whether the problem should be set in the coarse level state or left in the fine level state */
  /* we choose to leave it in the fine level state, only the adapt_data has been filled with coarse grid values */
  /* no, it gets swapped */
  //t8dg_advect_diff_problem_swap_from_restriction_to_prolongation(coarse_lvl->problem, &coarse_lvl->tmp_mortar_coarse_lvl);
  //t8dg_mg_swap_coarse_fine_data(coarse_lvl);
  t8dg_dof_values_swap ((coarse_lvl->dof_values), (coarse_lvl->dof_values_adapt));

  t8dg_dof_values_swap (&(coarse_lvl->adapt_data->dof_values), &(coarse_lvl->adapt_data->dof_values_adapt));

  t8dg_values_mg_swap_instances_to_coarse_lvl (coarse_lvl->dg_values, &(coarse_lvl->tmp_mortar_coarse_lvl), &(coarse_lvl->forest_coarsened),
                                               &(coarse_lvl->problem_forest));

  //t8dg_values_mg_clean_restriction(coarse_lvl->dg_values);
  return 0;

#if 0
  /* Allocate space for dof_values getting restricted */
  restricting_dofs = t8dg_dof_values_duplicate ();

  coarse_lvl->adapt_data->dof_values =
    /* the 'in' vector should be restricted to the coarse level */
    pc_ctx->adapt_data->dof_values = pc_ctx->dofs_initial_problem;

  /* muss der forest vorher ge-ref-t werden */
  /* Call iterate_replace in order to interpolate the residual from the fine forest to the coarser forest mesh */
  t8_forest_iterate_replace (pc_ctx->forest_coarsened, pc_ctx->problem_forest, t8dg_adapt_replace);

  /* Change forest, values and degrees of freedom of the problem to the corresponding values of the coarse-level problem */
  t8dg_values_change_to_coarse_lvl_mg (pc, pc_ctx->dg_values);
#if 0
  /* muss alles ausgelagert in eine funktion in dg_values */
  pc_ctx->tmp_dg_values_forest = pc_ctx->problem->dg_values->forest;
  pc_ctx->problem->dg_values->forest = pc_ctx->problem->dg_values->forest_adapt;
  pc_ctx->problem->dg_values->forest_adapt = NULL;

  pc_ctx->tmp_local_values = pc_ctx->problem->dg_values->local_values;
  pc_ctx->problem->dg_values->local_values = pc_ctx->problem->dg_values->local_values_adapt;
  pc_ctx->problem->dg_values->local_values_adapt = NULL;

  t8dg_local_values_set_all_ghost_elements (pc_ctx->problem->dg_values->local_values);

  pc_ctx->tmp_mortar_array = pc_ctx->problem->dg_values->mortar_array;
  pc_ctx->problem->dg_values->mortar_array =
    t8dg_mortar_array_new_empty (pc_ctx->problem->dg_values->forest, pc_ctx->problem->dg_values->local_values);
  /* bis hier hin */
#endif
  pc_ctx->tmp_forest = pc_ctx->problem_forest;
  pc_ctx->problem_forest = pc_ctx->forest_coarsened;
  pc_ctx->forest_coarsened = NULL;
#endif

}

/* Prior to the prolongation, the restriction has been called, which swapped the sizes of 'normal' and 'adapt' data within the advect_diff_problem; therefore a pointer to dof_values is in fact now apointer to dof_values_adapt, etc. */
PetscErrorCode
MatMult_MF_Prolongation (Mat A, Vec in, Vec out)
{
  t8dg_debugf ("Prolongation got called\n");
  PetscErrorCode      ierr;
  t8dg_mg_interpolating_ctx_t *appctx;
  t8dg_mg_coarse_lvl_t *coarse_lvl;

  /* Get the mtarix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);
  t8dg_debugf ("Appctx got initialized\n");
  /* Get the coarse level mesh */
  coarse_lvl = appctx->coarse_lvl;
  t8dg_debugf ("before change vec to dof\n");
  //t8dg_values_swap_to_adapt_data(coarse_lvl->dg_values, &(coarse_lvl->tmp_mortar_coarse_lvl));
  //t8dg_mg_swap_coarse_fine_data(coarse_lvl);
  /* Write the entries of the 'in' vector in adapt_data->dof_values, these are getting adapted/restricted */
  t8dg_precon_write_vec_to_dof (&in, coarse_lvl->adapt_data->dof_values);
  t8_forest_set_user_data (coarse_lvl->problem_forest, coarse_lvl->adapt_data);
  t8dg_debugf ("before iterate replace\n");
  t8dg_values_mg_allocate_adapt (coarse_lvl->dg_values, coarse_lvl->problem_forest);
  t8dg_debugf ("forest coarsneed hat %d Elemente\n", t8_forest_get_global_num_elements (coarse_lvl->forest_coarsened));
  //t8dg_values_mg_allocate_prolongation(coarse_lvl->dg_values, coarse_lvl->problem_forest);
  t8dg_debugf ("problem forest hat %d Elemente\n", t8_forest_get_global_num_elements (coarse_lvl->problem_forest));
  /* Calling iterate_replace executes the adaption/prolongation */
  t8_forest_iterate_replace (coarse_lvl->problem_forest, coarse_lvl->forest_coarsened, t8dg_adapt_replace);
  t8dg_debugf ("nach iterate replace\n");
  /* Write the degrees of freedom to a petsc vec */
  t8dg_precon_write_dof_to_vec (coarse_lvl->adapt_data->dof_values_adapt, &out, appctx->fine_lvl_global_indexing);
  t8dg_debugf ("nach change dof to vec\n");
  /* set the problem into the intitial/fine-mesh state */
  //t8dg_advect_diff_problem_swap_from_prolongation_to_initial(coarse_lvl->);
  //t8dg_mg_swap_coarse_fine_data(coarse_lvl);
  t8dg_dof_values_swap ((coarse_lvl->dof_values), (coarse_lvl->dof_values_adapt));

  t8dg_dof_values_swap (&(coarse_lvl->adapt_data->dof_values), &(coarse_lvl->adapt_data->dof_values_adapt));

  t8dg_values_mg_swap_instances_to_fine_lvl (coarse_lvl->dg_values, &(coarse_lvl->tmp_mortar_coarse_lvl), &(coarse_lvl->forest_coarsened),
                                             &(coarse_lvl->problem_forest));
  t8dg_debugf ("nach back swap\n");

  /*Allocate for next restriction */
  t8dg_values_mg_allocate_adapt (coarse_lvl->dg_values, coarse_lvl->forest_coarsened);
  t8dg_debugf ("nach allocate adapt in prolongation\n");
  return 0;
}

PetscErrorCode
MatMult_MF_Coarse_LVL (Mat A, Vec in, Vec out)
{
  t8dg_debugf ("Coarse level iteration\n");
  PetscErrorCode      ierr;
  t8dg_coarse_matrix_ctx_t *c_appctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &c_appctx);
  CHKERRQ (ierr);

  /* Copy the petsc vec entries into the dof_values of the problem in order to compute their derivation */
  t8dg_precon_write_vec_to_dof (&in, *(c_appctx->coarse_lvl->dof_values));

  /* Calculate the time derivation of the degrees of freedom */
  (c_appctx->time_derivative_func) (*(c_appctx->coarse_lvl->dof_values), c_appctx->problem_dofs_derivation, c_appctx->current_point_in_time,
                                    c_appctx->user_data);

  /* Copy the dof_values back to a petsc_vector */
  t8dg_precon_write_dof_to_vec (c_appctx->problem_dofs_derivation, &out, c_appctx->coarse_lvl->global_indexing);

  return 0;

}

void
t8dg_mg_construct_coarse_lvl (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                              t8dg_mg_coarse_lvl_t * coarse_lvl_mesh, t8dg_corase_lvl_adapt_func_t coarsening_func)
{
  PetscErrorCode      ierr;
  t8_gloidx_t         global_offset_to_first_local_elem;
  size_t              iter;
  t8dg_debugf ("coarse_lvl_construct wurde aufgerufen\n");
  coarse_lvl_mesh->problem = problem;
  coarse_lvl_mesh->dof_values = pdof_array;     //  (t8dg_advect_diff_problem_get_dof_values(problem)); //*pdof_array;

  coarse_lvl_mesh->dg_values = *(t8dg_advect_diff_problem_get_dg_values (problem));
  coarse_lvl_mesh->adapt_data = *(t8dg_advect_diff_problem_get_adapt_data (problem));
  coarse_lvl_mesh->tmp_mortar_coarse_lvl = NULL;
  coarse_lvl_mesh->problem_forest = *(t8dg_advect_diff_problem_get_forest (problem));
  coarse_lvl_mesh->dof_values_adapt = (t8dg_advect_diff_problem_get_dof_values_adapt (problem));
  t8dg_debugf ("weitere zuweisungen\n");
  t8dg_adapt_data_set_time (coarse_lvl_mesh->adapt_data,
                            t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem)));
  t8dg_debugf ("nach set time\n");
  coarse_lvl_mesh->adapt_data->dof_values = t8dg_dof_values_clone (*(coarse_lvl_mesh->dof_values));
  t8dg_debugf ("vor source sink\n");
  /* If source and sink terms are considered */
  if (coarse_lvl_mesh->adapt_data->source_sink_fn != NULL) {
    t8dg_adapt_data_interpolate_source_fn (coarse_lvl_mesh->adapt_data);
  }
  t8dg_debugf ("nach source sink interpolation\n");
  /* Keep the original forest of the problem */
  t8_forest_ref (coarse_lvl_mesh->problem_forest);

  /* Initialize a coarsened forest conatining the elements of the coarser muligrid mesh */
  t8_forest_init (&(coarse_lvl_mesh->forest_coarsened));

  /* Set the user-data pointer of the new forest */
  t8_forest_set_user_data (coarse_lvl_mesh->forest_coarsened, coarse_lvl_mesh->adapt_data);

  /* Set the adapt function coarsening the forest */
  t8_forest_set_adapt (coarse_lvl_mesh->forest_coarsened, coarse_lvl_mesh->problem_forest, coarsening_func, 0);

  /* Ghost values are needed for solving the system on the coarser mesh */
  t8_forest_set_ghost (coarse_lvl_mesh->forest_coarsened, 1, T8_GHOST_FACES);

  /* Commit the pre-set forest, so that it will be adapted */
  t8_forest_commit (coarse_lvl_mesh->forest_coarsened);
  t8dg_debugf ("after construction forest_coarsened\n");
  /* If source and sink terms were considered in former calculations regarding the adapted degrees of freedom, they are going to be destroyed, because they are not needed anymore */
  if (coarse_lvl_mesh->adapt_data->source_sink_fn != NULL) {
    t8dg_dof_values_destroy (&(coarse_lvl_mesh->adapt_data->source_sink_dof));
  }

  /* Allocate space for new elements and ghost values needed by the coarsened forest */
  t8dg_values_allocate_adapt (coarse_lvl_mesh->dg_values, coarse_lvl_mesh->forest_coarsened);

  /* Forest gets references in allocate_adapt, therefore, it has be un-referenced one time */
  t8_forest_unref (&coarse_lvl_mesh->forest_coarsened);

  /* Create and allocate new degrees of freedom for the coarsened forest */
  *(coarse_lvl_mesh->dof_values_adapt) =
    t8dg_dof_values_new (coarse_lvl_mesh->forest_coarsened, t8dg_values_get_global_values_array (coarse_lvl_mesh->dg_values));

  /* Create and allocate new degrees of freedom of the same size for adapt_data->dof_values_adapt */
  coarse_lvl_mesh->adapt_data->dof_values_adapt =
    t8dg_dof_values_new (coarse_lvl_mesh->forest_coarsened, t8dg_values_get_global_values_array (coarse_lvl_mesh->dg_values));

  /* Number of local degrees of freedom */
  coarse_lvl_mesh->num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (*(coarse_lvl_mesh->dof_values_adapt)) *
              t8dg_dof_get_max_num_element_dof (*(coarse_lvl_mesh->dof_values_adapt)));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (coarse_lvl_mesh->num_local_dofs, &coarse_lvl_mesh->global_indexing);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (coarse_lvl_mesh->forest_coarsened);
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem =
      global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (*(coarse_lvl_mesh->dof_values_adapt));
  }

  /* Fill the array of global indices */
  for (iter = 0; iter < coarse_lvl_mesh->num_local_dofs; ++iter) {
    (coarse_lvl_mesh->global_indexing)[iter] = iter + global_offset_to_first_local_elem;
  }
  t8dg_debugf ("End of construction mehtod\n");
}

#if 0
void
newfunc ()
{

#if 0
  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  /* Store a pointer to the advect_diff_problem */
  pc_ctx->problem = problem;

  /* Store the function pointer to the advect_diff_time_derivative function */
  pc_ctx->time_derivative_func = time_derivative_func;

  /* Store the current time point in which this preconditioning happens */
  pc_ctx->current_time_point = current_time_point;

  /* Store a pointer to the problem related degrees of freedom */
  pc_ctx->dofs_initial_problem = *problem_dofs;
  //pc_ctx->adapt_data_dof_values = *adapt_data_dof_values;

  /* Store the initial problem's forest */
  pc_ctx->problem_forest = *initial_problem_forest;

  /* Store the adapt data which is needed to set up the coarse level mesh */
  pc_ctx->adapt_data = *adapt_data;

  /* Store a pointer to the dg_values */
  pc_ctx->dg_values = *initial_dg_values;
  /* Set the current time in adapt_data */
  t8dg_adapt_data_set_time (pc_ctx->adapt_data, pc_ctx->current_time_point);

  /* Set the degrees of freedom which are going to be adapted to the current approximation */
  pc_ctx->adapt_data->dof_values = pc_ctx->dofs_initial_problem;

  /* If source and sink terms are considered */
  if (pc_ctx->adapt_data->source_sink_fn != NULL) {
    t8dg_adapt_data_interpolate_source_fn (pc_ctx->adapt_data);
  }

  /* Keep the original forest of the problem */
  t8_forest_ref (pc_ctx->problem_forest);

  /* Initialize a coarsened forest conatining the elements of the coarser muligrid mesh */
  t8_forest_init (&pc_ctx->forest_coarsened);

  /* Set the user-data pointer of the new forest */
  t8_forest_set_user_data (pc_ctx->forest_coarsened, pc_ctx->adapt_data);

  /* Set the adapt function coarsening the forest */
  t8_forest_set_adapt (pc_ctx->forest_coarsened, pc_ctx->problem_forest, t8dg_adapt_multigrid_coarsen_finest_level, 0);

  /* Ghost values are needed for solving the system on the coarser mesh */
  t8_forest_set_ghost (pc_ctx->forest_coarsened, 1, T8_GHOST_FACES);

  /* Commit the pre-set forest, so that it will be adapted */
  t8_forest_commit (pc_ctx->forest_coarsened);

  /* If source and sink terms were considered in former calculations regarding the adapted degrees of freedom, they are going to be destroyed, because they are not needed anymore */
  if (pc_ctx->adapt_data->source_sink_fn != NULL) {
    t8dg_dof_values_destroy (&pc_ctx->adapt_data->source_sink_dof);
  }

  /* Allocate space for new elements and ghost values needed by the coarsened forest */
  t8dg_values_allocate_adapt (pc_ctx->dg_values, pc_ctx->forest_coarsened);

  /* Create and allocate new degrees of freedom for the coarsened forest */
  pc_ctx->problem_dofs_adapt = t8dg_dof_values_new (pc_ctx->forest_coarsened, t8dg_values_get_global_values_array (pc_ctx->dg_values));

  /* Set the degrees of freedom in adapt_data that will be interpolated/restricted to the coarse level forest */
  pc_ctx->adapt_data->dof_values = pc_ctx->problem_dofs_adapt;

  /* Save an indexing vector regarding this coarse forest mesh */
  /* Number of process-local degrees of freedom regarding the coarsened forest */
  pc_ctx->num_local_dofs_coarse_grid =
    (size_t) (t8dg_dof_get_num_local_elements (pc_ctx->problem_dofs_adapt) * t8dg_dof_get_max_num_element_dof (pc_ctx->problem_dofs_adapt));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (pc_ctx->num_local_dofs_coarse_grid, &pc_ctx->global_indexing_coarse_grid);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (pc_ctx->forest_coarsened);
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem = global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (pc_ctx->problem_dofs_adapt);
  }
  /* Fill the array of global indices */
  for (iter = 0; iter < pc_ctx->num_local_dofs_coarse_grid; ++iter) {
    (pc_ctx->global_indexing_coarse_grid)[iter] = iter + global_offset_to_first_local_elem;
  }

  /* Right now the coarsened forest functioning as the coarse mesh in the multigrid preconditioner with interpolated degrees of freedom (of the first residual) and further information needed should be set up */
#endif
  return 0;
}
#endif
void
t8dg_mg_create_coarse_lvl_matrix (Mat * A_coarse, t8dg_mg_coarse_lvl_t * coarse_lvl, t8dg_coarse_matrix_ctx_t * cmat_ctx)
{
  PetscErrorCode      ierr;
  /* Create the restricted/coarsened matrix which used to solve the problem on the coarse mesh */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, coarse_lvl->num_local_dofs, coarse_lvl->num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE, cmat_ctx,
                    A_coarse);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) * A_coarse, "Coarse-Matrix Application");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation on the coarse level mesh */
  ierr = MatShellSetOperation (*A_coarse, MATOP_MULT, (void (*)(void)) MatMult_MF_Coarse_LVL);
  CHKERRQ (ierr);
}

void
t8dg_mg_create_restriction_matrix (Mat * Restriction, t8dg_mg_interpolating_ctx_t * res_prol_ctx)
{
  PetscErrorCode      ierr;
  /* Define a Restriction matrix which interpolates a vector from the fine level onto the coarse level */
  /* During the restriction the initial problem gets changed in order to solve on the coarse level */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, res_prol_ctx->num_local_dofs_coarse_grid, res_prol_ctx->num_local_dofs_fine_grid, PETSC_DETERMINE,
                    PETSC_DETERMINE, res_prol_ctx, Restriction);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) * Restriction, "Restriction Matrix - interpolating a vector from the fine to the coarse level");
  CHKERRQ (ierr);
  /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
  ierr = MatShellSetOperation (*Restriction, MATOP_MULT, (void (*)(void)) MatMult_MF_Restriction);
  CHKERRQ (ierr);
}

void
t8dg_mg_create_prolongation_matrix (Mat * Prolongation, t8dg_mg_interpolating_ctx_t * res_prol_ctx)
{
  PetscErrorCode      ierr;
  /* Define a Prolongation matrix which interpolates a vector from the coarse level onto the fine level */
  /* During the prolongation the intital problem gets re-transformed from the coarse level properties */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, res_prol_ctx->num_local_dofs_fine_grid, res_prol_ctx->num_local_dofs_coarse_grid, PETSC_DETERMINE,
                    PETSC_DETERMINE, res_prol_ctx, Prolongation);
  CHKERRQ (ierr);
  ierr =
    PetscObjectSetName ((PetscObject) * Prolongation, "Prolongation Matrix - interpolating a vector from the coarse to the fine level");
  CHKERRQ (ierr);
  /* Define the (multiplicative) MatVec-Operation which resembles the application of the prolongation matrix */
  ierr = MatShellSetOperation (*Prolongation, MATOP_MULT, (void (*)(void)) MatMult_MF_Prolongation);
  CHKERRQ (ierr);
}

#if 0
void
ExampleRoutine ()
{
  PetscErrorCode      ierr;
  KSP                 ksp;
  KSP                 smoother;
  KSP                 coarse_solver;
  PC                  pc;
  PC                  coarse_pc;
  Mat                 Restriction, Prolongation;
  t8dg_coarse_matrix_ctx_t cmat_ctx;
  Mat                 A, A_coarse;
  t8dg_mg_interpolating_ctx_t res_prol_ctx;
  t8dg_mg_coarse_lvl_t c_lvl;

  /* Angenommen Mat A, Vec u und Vec f sind bereits auf dem fine level initialisiert */

  /* Number of fine grid process-local degrees of freedom */
  res_prol_ctx.num_local_dofs_fine_grid = t8dg_dof_get_num_local_elements (problem_dofs);

  /* Build a coarse mesh which is used within the multigrid preconditioner */
  t8dg_mg_construct_coarse_lvl (problem, &c_lvl, t8dg_adapt_multigrid_coarsen_finest_level);

  /* muessen noch zugewiesen werden */
  cmat_ctx.coarse_lvl = &c_lvl;
  cmat_ctx.problem_dofs =;
  cmat_ctx.current_point_in_time =;
  cmat_ctx.user_data =;
  cmat_ctx.problem_dofs_derivation =;
  cmat_ctx.time_derivative_func =;

  /* Number of coarse grid process-local degrees of freedom */
  res_prol_ctx.num_local_dofs_coarse_grid = c_lvl.num_local_dofs;

  /* Save the coarse level forest inside the restriction_prolongation context */
  res_prol_ctx.coarse_lvl = &c_lvl;
  /* Create the restricted matrix which used to solve the problem on the coarse mesh */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, c_lvl.num_local_dofs, c_lvl.num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE, &cmat_ctx, &A_coarse);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A_coarse, "Coarse-Matrix Application");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation on the coarse level mesh */
  ierr = MatShellSetOperation (A_coarse, MATOP_MULT, (void (*)(void)) MatMult_MF_Coarse_LVL);
  CHKERRQ (ierr);

  /* Define a Restriction matrix which interpolates a vector from the fine level onto the coarse level */
  /* During the restriction the initial problem gets changed in order to solve on the coarse level */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, res_prol_ctx.num_local_dofs_coarse_grid, res_prol_ctx.num_local_dofs_fine_grid, PETSC_DETERMINE,
                    PETSC_DETERMINE, &res_prol_ctx, &Restriction);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) Restriction, "Restriction Matrix - interpolating a vector from the fine to the coarse level");
  CHKERRQ (ierr);
  /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
  ierr = MatShellSetOperation (Restriction, MATOP_MULT, (void (*)(void)) MatMult_MF_Restriction);
  CHKERRQ (ierr);

  /* Define a Prolongation matrix which interpolates a vector from the coarse level onto the fine level */
  /* During the prolongation the intital problem gets re-transformed from the coarse level properties */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, res_prol_ctx.num_local_dofs_fine_grid, res_prol_ctx.num_local_dofs_coarse_grid, PETSC_DETERMINE,
                    PETSC_DETERMINE, &res_prol_ctx, &Prolongation);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) Prolongation, "Prolongation Matrix - interpolating a vector from the coarse to the fine level");
  CHKERRQ (ierr);
  /* Define the (multiplicative) MatVec-Operation which resembles the application of the prolongation matrix */
  ierr = MatShellSetOperation (Prolongation, MATOP_MULT, (void (*)(void)) MatMult_MF_Prolongation);
  CHKERRQ (ierr);

  /* KSP and PC for the initial problem */
  ierr = KSPCreate (PETSC_COMM_WORLD, &ksp);
  CHKERRQ (ierr);
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);

  ierr = PCSetType (pc, PCMG);
  CHKERRQ (ierr);
  ierr = PCMGSetLevels (pc, 2, NULL);
  CHKERRQ (ierr);
  ierr = PCMGSetType (pc, PC_MG_MULTIPLICATIVE);
  CHKERRQ (ierr);
  ierr = PCMGSetCycleType (pc, PC_MG_V_CYCLE);
  CHKERRQ (ierr);
  /* Define coarse solver */
  ierr = PCMGGetCoarseSolve (pc, &coarse_solver);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (coarse_solver, A_coarse, A_coarse);
  CHKERRQ (ierr);
  ierr = KSPSetTolerances (coarse_solver, 1.0e-8, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  ierr = KSPSetType (coarse_solver, KSPGMRES);
  CHKERRQ (ierr);
  /* Set the options in coarse_solver */
  /* --------- Define operators for coarse solvers  ----------- */
  /* Define smoother on level 1 (fine mesh) */
  /* Dont want to have a ksp smoother, therefore, a pc context needs to be extracted from the (ksp) smoother which then resembles for example a jacobi or richardson smoothing iteration */
  ierr = PCMGGetSmoother (pc, 1, &smoother);
  CHKERRQ (ierr);
  /* Set the options in smoother */
  /* --------Define operators for smoothers ------- */ */
    /* Smoother on level 1 (fine mesh) operates on the intial problem */
    ierr = KSPSetOperators (smoother, A, A);
  CHKERRQ (ierr);
  ierr = KSPSetType (smoother, KSPGMRES);
  CHKERRQ (ierr);
#if 0
  // das versuchen wir spaeter weiter, erstmal gmres weil weniger zu implmentieren 
  /* Get preconditioner on level 1 (fine mesh) */
  ierr = KSPGetPC (smoother, &coarse_pc);
  CHKERRQ (ierr);
  /* Set the coarse_pc to a smoother of choice - it has to be wrapped by a PCSHELL */
  ierr = PCSetType (coarse_pc, PCSHELL);
  CHKERRQ (ierr);
  //to be continued */
#endif

  /* Amount of pre- and post-smoothing steps */
  ierr = PCMGSetNumberSmooth (pc, 3);
  CHKERRQ (ierr);
  /* Set Restriction operator */
  ierr = PCMGSetRestriction (pc, 1, Restriction);
  CHKERRQ (ierr);
  /* Set Prolongation operator */
  ierr = PCMGSetInterpolation (pc, 1, Prolongation);
  CHKERRQ (ierr);
  /* Set residual calculation - my be default since only Mat's and Vec's are used */
  ierr = PCMGSetResidual (pc, 1, PCMGResidualDefault, Mat systemmatrix);
  CHKERRQ (ierr);
}
#endif

void
t8dg_precon_write_dof_to_vec (t8dg_dof_values_t * dofs, Vec * p_vec, PetscInt * indexing)
{
  PetscErrorCode      ierr;
  /* Sets the values of a petsc vector to the entries of a t8dg_dof_values_t */
  ierr =
    VecSetValues (*p_vec, (t8dg_dof_get_num_local_elements (dofs) * t8dg_dof_get_max_num_element_dof (dofs)), indexing,
                  t8dg_dof_values_get_double_pointer (dofs, 0), INSERT_VALUES);
  CHKERRQ (ierr);
  /* Assemble the petsc vector */
  ierr = VecAssemblyBegin (*p_vec);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (*p_vec);
  CHKERRQ (ierr);
}

void
t8dg_precon_write_vec_to_dof (Vec * p_vec, t8dg_dof_values_t * dofs)
{
  PetscErrorCode      ierr;
  PetscScalar        *vec_reader;
  double             *dof_pointer;
  int                 dof_iter;

  /* Retrieve the local part of a petsc vector */
  ierr = VecGetArrayRead (*p_vec, &vec_reader);
  CHKERRQ (ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  //dof_pointer = t8dg_dof_get_double_pointer_to_array (dofs);
  dof_pointer = t8dg_dof_values_get_double_pointer (dofs, 0);
  /* Overwrite the dof values with the entries of a petsc vector */
  /* it is assumed (not checked!) that the local petsc vector is of length dof_values->num_local_elements */
  for (dof_iter = 0; dof_iter < (t8dg_dof_get_num_local_elements (dofs) * t8dg_dof_get_max_num_element_dof (dofs)); ++dof_iter) {
    dof_pointer[dof_iter] = (double) vec_reader[dof_iter];
  }
  /* Restore the array entries in the petsc vector */
  ierr = VecRestoreArrayRead (*p_vec, &vec_reader);
  CHKERRQ (ierr);
}

#if 0
KSP                *
t8dg_mg_coarse_grid_solver_context (PC pc)
{

  KSP                 ksp;
  PC                  pc;
  double             *dof_values_ptr;
  PetscScalar        *approx_coarse_grid_correction;

  t8dg_debugf ("\n*********\n Multigrid Coarse Grid Solver has been called\n*********\n");

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  //pc_ctx.time_data = problem->time_data;
  //pc_ctx.user_data = user_data;
  //coarse_grid_dofs = *problem_dofs;
  pc_ctx->dof_current_coarse_grid_correction = t8dg_dof_values_duplicate (pc_ctx->dofs_initial_problem);

  /* Create the right-hand-side of the system resulting from coarsening the residual vector */
  ierr = VecCreateMPI (PETSC_COMM_WORLD, (PetscInt) pc_ctx->num_local_dofs_coarse_grid, PETSC_DETERMINE, &f);
  CHKERRQ (ierr);
  /* Fill the vector with the current degrees of freedom resulting from the former time step */
  ierr =
    VecSetValues (f, pc_ctx->num_local_dofs_coarse_grid, pc_ctx->global_indexing_coarse_grid,
                  t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem), INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (f);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (f);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side (coarsened fine-level residual)");
  CHKERRQ (ierr);

  /* Create an Approximation Vector storing the coarse grid correection of the residual */
  /* Use the size and allocation similar to vector f */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  ierr = VecCopy (f, u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation Coarse Grid Correction");
  CHKERRQ (ierr);

  /* Setting up matrix-free Matrix */
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context (in this case the preconditioning context) needed by the matrix-free MatVec multiplication */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, (PetscInt) pc_ctx->num_local_dofs_coarse_grid, (PetscInt) pc_ctx->num_local_dofs_coarse_grid,
                    PETSC_DETERMINE, PETSC_DETERMINE, (void *) &pc_pctx, &A);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix multigrid coarse level correction");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_MG_Coarse_Level);
  CHKERRQ (ierr);

  /* Setting up KSP Solver */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);
  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);
  /* Select GMRES as Solver; default: Restarted GMRES(30) */
  ierr = KSPSetType (ksp, KSPGMRES);
  CHKERRQ (ierr);
  /* Select a Preconditioner - without Preconditioning */
  ierr = PCSetType (pc, PCNONE);
  CHKERRQ (ierr);
  /* Set GMRES ieteration tolerances */
  ierr = KSPSetTolerances (ksp, 1.0e-7, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Set the Solver up */
  ierr = KSPSetUp (ksp);
  CHKERRQ (ierr);

  t8dg_debugf ("\nCoarse Grid GMRES got called\n");

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, f, u);
  CHKERRQ (ierr);

  t8dg_debugf ("\nCoarse Grid GMRES solve completed\n");

  /* View whether or not the GMRES did converge */
  ierr = KSPConvergedRateView (ksp, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

  /* Retrieve the local part of the 'solution' Vector u for reading purposes */
  ierr = VecGetArrayRead (u, &approx_coarse_grid_correction);
  CHKERRQ (ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem);

  /* Overwrite the coarse grid degrees of freedom with the newly calculated degrees of freedom containing the residual coarse grid correction */
  for (dof_iter = 0; dof_iter < pc_ctx->num_local_dofs_coarse_grid; ++dof_iter) {
    dof_values_ptr[dof_iter] = (double) approx_coarse_grid_correction[dof_iter];
  }
  /* Restore the array entries */
  ierr = VecRestoreArrayRead (u, &approx_coarse_grid_correction);
  CHKERRQ (ierr);

  /* Free used resources */
  t8dg_dof_values_destroy (&pc_ctx->dof_current_coarse_grid_correction);

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);

}

#endif

/*********************************

********* Lets try again *********

**********************************/

#if 0
#if 0
extern PetscErrorCode MGShellPCCreate (t8dg_preconditioner_mg_ctx_t **);
extern PetscErrorCode MGShellPCSetUp (PC, t8dg_linear_advection_diffusion_problem_t *, t8dg_dof_values_t **, t8dg_dof_values_t **,
                                      t8dg_adapt_data_t **, t8dg_dg_values_t **, t8_forest_t *, t8dg_time_matrix_application, double);
extern PetscErrorCode MGShellPCApply (PC, Vec, Vec);
extern PetscErrorCode MGShellPCDestroy (PC);
#endif
extern PetscErrorCode MatMult_MF_MG_Coarse_Level (Mat, Vec, Vec);
extern static void  t8dg_multigrid_coarse_level_solver (PC);
extern static void  t8dg_multigrid_smoother_richardson (PC pc, t8dg_dof_values_t * problem_dofs, t8dg_dof_values_t * rhs,
                                                        int smoothing_steps, double optimization_parameter);

/* Creates a two level multigrid preconditioner with corresponding 'preconditioning context' */
PetscErrorCode
MGShellPCCreate (t8dg_preconditioner_mg_ctx_t ** pcshell)
{
  PetscErrorCode      ierr;
  t8dg_preconditioner_mg_ctx_t *pc_ctx;

  /* Create a new two level multigrid preconditioning context */
  ierr = PetscNew (&pc_ctx);
  CHKERRQ (ierr);
  /* Initialize all members with NULL pointers */
  pc_ctx->problem = NULL;
  pc_ctx->forest_coarsened = NULL;
  pc_ctx->tmp_dg_values_forest = NULL;
  pc_ctx->tmp_forest = NULL;
  pc_ctx->tmp_local_values = NULL;
  pc_ctx->dof_values_fine_mesh = NULL;
  pc_ctx->dof_prolongated_coarse_grid_correction = NULL;
  pc_ctx->tmp_mortar_array = NULL;
  pc_ctx->mortar_array_coarsened = NULL;

  /* return the preconditioning context */
  *pcshell = pc_ctx;

  return 0;
}

/* Sets up the multigrid preconditioner - calculates a coarse mesh and other components */
PetscErrorCode
MGShellPCSetUp (PC pc, t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** problem_dofs,
                t8dg_dof_values_t ** problem_dofs_adapt, t8dg_adapt_data_t ** adapt_data, t8dg_dg_values_t ** initial_dg_values,
                t8_forest_t * initial_problem_forest, t8dg_time_matrix_application time_derivative_func, double current_time_point)
{
  PetscErrorCode      ierr;
  t8dg_preconditioner_mg_ctx_t *pc_ctx;
  t8_gloidx_t         global_offset_to_first_local_elem;
  size_t              iter;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  /* Store a pointer to the advect_diff_problem */
  pc_ctx->problem = problem;

  /* Store the function pointer to the advect_diff_time_derivative function */
  pc_ctx->time_derivative_func = time_derivative_func;

  /* Store the current time point in which this preconditioning happens */
  pc_ctx->current_time_point = current_time_point;

  /* Store a pointer to the problem related degrees of freedom */
  pc_ctx->dofs_initial_problem = *problem_dofs;
  //pc_ctx->adapt_data_dof_values = *adapt_data_dof_values;

  /* Store the initial problem's forest */
  pc_ctx->problem_forest = *initial_problem_forest;

  /* Store the adapt data which is needed to set up the coarse level mesh */
  pc_ctx->adapt_data = *adapt_data;

  /* Store a pointer to the dg_values */
  pc_ctx->dg_values = *initial_dg_values;
  /* Set the current time in adapt_data */
  t8dg_adapt_data_set_time (pc_ctx->adapt_data, pc_ctx->current_time_point);

  /* Set the degrees of freedom which are going to be adapted to the current approximation */
  pc_ctx->adapt_data->dof_values = pc_ctx->dofs_initial_problem;

  /* If source and sink terms are considered */
  if (pc_ctx->adapt_data->source_sink_fn != NULL) {
    t8dg_adapt_data_interpolate_source_fn (pc_ctx->adapt_data);
  }

  /* Keep the original forest of the problem */
  t8_forest_ref (pc_ctx->problem_forest);

  /* Initialize a coarsened forest conatining the elements of the coarser muligrid mesh */
  t8_forest_init (&pc_ctx->forest_coarsened);

  /* Set the user-data pointer of the new forest */
  t8_forest_set_user_data (pc_ctx->forest_coarsened, pc_ctx->adapt_data);

  /* Set the adapt function coarsening the forest */
  t8_forest_set_adapt (pc_ctx->forest_coarsened, pc_ctx->problem_forest, t8dg_adapt_multigrid_coarsen_finest_level, 0);

  /* Ghost values are needed for solving the system on the coarser mesh */
  t8_forest_set_ghost (pc_ctx->forest_coarsened, 1, T8_GHOST_FACES);

  /* Commit the pre-set forest, so that it will be adapted */
  t8_forest_commit (pc_ctx->forest_coarsened);

  /* If source and sink terms were considered in former calculations regarding the adapted degrees of freedom, they are going to be destroyed, because they are not needed anymore */
  if (pc_ctx->adapt_data->source_sink_fn != NULL) {
    t8dg_dof_values_destroy (&pc_ctx->adapt_data->source_sink_dof);
  }

  /* Allocate space for new elements and ghost values needed by the coarsened forest */
  t8dg_values_allocate_adapt (pc_ctx->dg_values, pc_ctx->forest_coarsened);

  /* Create and allocate new degrees of freedom for the coarsened forest */
  pc_ctx->problem_dofs_adapt = t8dg_dof_values_new (pc_ctx->forest_coarsened, t8dg_values_get_global_values_array (pc_ctx->dg_values));

  /* Set the degrees of freedom in adapt_data that will be interpolated/restricted to the coarse level forest */
  pc_ctx->adapt_data->dof_values = pc_ctx->problem_dofs_adapt;

  /* Save an indexing vector regarding this coarse forest mesh */
  /* Number of process-local degrees of freedom regarding the coarsened forest */
  pc_ctx->num_local_dofs_coarse_grid =
    (size_t) (t8dg_dof_get_num_local_elements (pc_ctx->problem_dofs_adapt) * t8dg_dof_get_max_num_element_dof (pc_ctx->problem_dofs_adapt));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (pc_ctx->num_local_dofs_coarse_grid, &pc_ctx->global_indexing_coarse_grid);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (pc_ctx->forest_coarsened);
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem = global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (pc_ctx->problem_dofs_adapt);
  }
  /* Fill the array of global indices */
  for (iter = 0; iter < pc_ctx->num_local_dofs_coarse_grid; ++iter) {
    (pc_ctx->global_indexing_coarse_grid)[iter] = iter + global_offset_to_first_local_elem;
  }

  /* Right now the coarsened forest functioning as the coarse mesh in the multigrid preconditioner with interpolated degrees of freedom (of the first residual) and further information needed should be set up */

  return 0;
}

/* muss in dg_values */
#if 0
#if T8_WITH_PETSC
void
t8dg_values_change_to_coarse_lvl_mg (PC pc, t8dg_values_t * values)
{
  t8dg_preconditioner_mg_ctx_t *pc_ctx;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  pc_ctx->tmp_dg_values_forest = values->forest;
  values->forest = values->forest_adapt;
  values->forest_adapt = NULL;

  pc_ctx->tmp_local_values = values->local_values;
  values->local_values = values->local_values_adapt;
  values->local_values_adapt = NULL;

  t8dg_local_values_set_all_ghost_elements (values->local_values);

  pc_ctx->tmp_mortar_array = values->mortar_array;
  values->mortar_array = t8dg_mortar_array_new_empty (values->forest, values->local_values);
}

void
t8dg_values_change_back_to_fine_lvl_mg (PC pc, t8dg_values_t * values)
{
  t8dg_preconditioner_mg_ctx_t *pc_ctx;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  values->forest = pc_ctx->tmp_dg_values_forest;
  pc_ctx->tmp_dg_values_forest = NULL;
  values->forest_adapt = NULL;

  t8dg_local_values_destroy (&values->local_values);
  values->local_values = pc_ctx->tmp_local_values;
  pc_ctx->tmp_local_values = NULL;
  t8dg_local_values_destroy (&values->local_values_adapt);
  values->local_values_adapt = NULL;

  t8dg_mortar_array_destroy (values->mortar_array);
  values->mortar_array = pc_ctx->tmp_mortar_array;

}
#endif

#endif

/* Applies the two level multigrid preconditioner to a vector */
PetscErrorCode
MGShellPCApply (PC pc, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_preconditioner_mg_ctx_t *pc_ctx;
  double             *dof_values_ptr;
  int                 smoothing_iterations = 3;
  PetscScalar        *current_local_approx;
  size_t              dof_iter;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  /* Retrieve the degrees of freedom of the 'in' vector */
  ierr = VecGetArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Get a pointer to the original degrees of freedom of the probblem */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem);

  /* Overwrite them with the 'in' Vector of the Matrix Multiplication */
  for (dof_iter = 0; dof_iter < t8dg_dof_get_num_local_elements (pc_ctx->dofs_initial_problem); ++dof_iter) {
    dof_values_ptr[dof_iter] = current_local_approx[dof_iter];
  }

  /* Restore the array entries of the 'in' Vector */
  ierr = VecRestoreArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* These are the dof_values and not the residual that is put into this function and i dont think i can calcualte residual newly because it may eventually differ from the one that has been put into this function */
  pc_ctx->dof_values_fine_mesh = t8dg_dof_values_clone (pc_ctx->dofs_initial_problem);

  /* Pre- Smoothing */
  /* Wird im Endeffekt nochmal das Residuum des Residuums berechenet, geglaettet und verwendet */

  /* Restriction */
  /* Set the degrees of freedom which are going to be adapted to the current approximation */
  /* adapt_data->fof_values should contain the degrees of freedom of the residual, which is going to be restricted onto the coarse level mesh */
  pc_ctx->adapt_data->dof_values = pc_ctx->dofs_initial_problem;

  /* muss der forest vorher ge-ref-t werden */
  /* Call iterate_replace in order to interpolate the residual from the fine forest to the coarser forest mesh */
  t8_forest_iterate_replace (pc_ctx->forest_coarsened, pc_ctx->problem_forest, t8dg_adapt_replace);

  /* Change forest, values and degrees of freedom of the problem to the corresponding values of the coarse-level problem */
  t8dg_values_change_to_coarse_lvl_mg (pc, pc_ctx->dg_values);
#if 0
  /* muss alles ausgelagert in eine funktion in dg_values */
  pc_ctx->tmp_dg_values_forest = pc_ctx->problem->dg_values->forest;
  pc_ctx->problem->dg_values->forest = pc_ctx->problem->dg_values->forest_adapt;
  pc_ctx->problem->dg_values->forest_adapt = NULL;

  pc_ctx->tmp_local_values = pc_ctx->problem->dg_values->local_values;
  pc_ctx->problem->dg_values->local_values = pc_ctx->problem->dg_values->local_values_adapt;
  pc_ctx->problem->dg_values->local_values_adapt = NULL;

  t8dg_local_values_set_all_ghost_elements (pc_ctx->problem->dg_values->local_values);

  pc_ctx->tmp_mortar_array = pc_ctx->problem->dg_values->mortar_array;
  pc_ctx->problem->dg_values->mortar_array =
    t8dg_mortar_array_new_empty (pc_ctx->problem->dg_values->forest, pc_ctx->problem->dg_values->local_values);
  /* bis hier hin */
#endif
  pc_ctx->tmp_forest = pc_ctx->problem_forest;
  pc_ctx->problem_forest = pc_ctx->forest_coarsened;
  pc_ctx->forest_coarsened = NULL;

  /* Coarse-Grid-Correction */
  /* After interpolating/restricting the problem to the coarse level problem, the exact solution to this coarse level problem needs to calculated/apprroximated */
  t8dg_multigrid_coarse_level_solver (pc_ctx->problem, pc);

  /* Now, when the coarse level solution is obtained, it has be interpolated/prolongated to the original (fine) forest mesh which is currently used in the initial problem */
  /* This can also be done by iterate_replace, because the finer and coarser forest are valid (and balanced) adaption of each other */

  /* Prolongation */
  /* Storage needs to be allocated, altough the old values are getting reassigned */
  t8dg_values_allocate_adapt (pc_ctx->dg_values, pc_ctx->tmp_forest);
  t8dg_dof_values_destroy (&pc_ctx->problem_dofs_adapt);        /* If they dont get destroyed a few lines earlier */
  pc_ctx->problem_dofs_adapt = t8dg_dof_values_new (pc_ctx->tmp_forest, t8dg_values_get_global_values_array (pc_ctx->dg_values));
  pc_ctx->adapt_data->dof_values_adapt = pc_ctx->problem_dofs_adapt;
  /* Call iterate_replace in order to interpolate/prolonagte the solution from the coarser forest to the initial (finer) forest mesh */
  t8_forest_iterate_replace (pc_ctx->tmp_forest, pc_ctx->problem_forest, t8dg_adapt_replace);

  /* Free some resources which are not needed */
  /* muss auch irgendwo anders gemacht werden */
  //t8dg_local_values_destroy(&pc_ctx->dg_values->local_values_adapt); /* wird in der methode weiter unten glaube ich schon gemacht */

  /* tmp_forest get referenced in allocate_adapt, therefore, it has to be un-referenced one time */
  t8_forest_unref (&pc_ctx->tmp_forest);

  /* Change values to finer forest */
  t8dg_values_change_back_to_fine_lvl_mg (pc, pc_ctx->dg_values);

#if 0
  //muss auch in values ausgelagert werden
  pc_ctx->problem->dg_values->forest = pc_ctx->tmp_dg_values_forest;
  pc_ctx->tmp_dg_values_forest = NULL;
  pc_ctx->problem->dg_values->forest_adapt = NULL;

  t8dg_local_values_destroy (&pc_ctx->problem->dg_values->local_values);
  pc_ctx->problem->dg_values->local_values = pc_ctx->tmp_local_values;
  pc_ctx->tmp_local_values = NULL;
  t8dg_local_values_destroy (&pc_ctx->problem->dg_values->local_values_adapt);

  t8dg_mortar_array_destroy (pc_ctx->problem->dg_values->mortar_array);
  pc_ctx->problem->dg_values->mortar_array = pc_ctx->tmp_mortar_array;
#endif
  /* Swap the forest pointers back to initial conditions */
  pc_ctx->forest_coarsened = pc_ctx->problem_forest;
  pc_ctx->problem_forest = pc_ctx->tmp_forest;
  pc_ctx->tmp_forest = NULL;

  /* have to work with problem->dof_values_adapt
     /* Post-Smoothing */

  /* Add the coarse level correction to the intial dofs of the fine mesh */
  /* Swap the values in order to leave the problem as it was before preconditioning, because the 'untouched' problem is needed by the timestepping methods */
  t8dg_dof_values_swap (&pc_ctx->dofs_initial_problem, &pc_ctx->dof_values_fine_mesh);
  t8dg_dof_values_axpy (1.0, pc_ctx->dof_values_fine_mesh, pc_ctx->dofs_initial_problem);

  /* Write out the 'in' vector with the applied coarse grid correction */
  ierr =
    VecSetValues (out, pc_ctx->num_local_dofs_coarse_grid, pc_ctx->global_indexing_coarse_grid,
                  t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem), INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (out);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (out);
  CHKERRQ (ierr);

  return 0;
}

/* Frees the resources used by the multigrid preconditioning method */
PetscErrorCode
MGShellPCDestroy (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_preconditioner_mg_ctx_t *pc_ctx;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  /* Free the coarse forest */
  t8_forest_unref (&pc_ctx->forest_coarsened);

  /* Free the degrees of freedom used for preconditioning purposes */
  t8dg_dof_values_destroy (&pc_ctx->problem->dof_values_adapt);

  /* Free the global indexing vector which enumerates the process-local coarse forest elements */
  ierr = PetscFree (pc_ctx->global_indexing_coarse_grid);
  CHKERRQ (ierr);

  /* Only if they are invoked by me */
  t8dg_dof_values_destroy (&pc_ctx->dof_values_fine_mesh);
  t8dg_dof_values_destroy (&pc_ctx->dof_prolongated_coarse_grid_correction);

  /* Free the prconditioning context */
  ierr = PetscFree (pc_ctx);
  CHKERRQ (ierr);

  return 0;
}

static void
t8dg_multigrid_coarse_level_solver (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_preconditioner_mg_ctx_t *pc_ctx;
  //t8dg_dof_values_t  *coarse_grid_dofs;
  t8_gloidx_t         global_offset_to_first_local_elem;
  Mat                 A;
  Vec                 u, f;
  KSP                 ksp;
  PC                  pc;
  double             *dof_values_ptr;
  PetscScalar        *approx_coarse_grid_correction;

  t8dg_debugf ("\n*********\n Multigrid Coarse Grid Solver has been called\n*********\n");

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  //pc_ctx.time_data = problem->time_data;
  //pc_ctx.user_data = user_data;
  //coarse_grid_dofs = *problem_dofs;
  pc_ctx->dof_current_coarse_grid_correction = t8dg_dof_values_duplicate (pc_ctx->dofs_initial_problem);

  /* Create the right-hand-side of the system resulting from coarsening the residual vector */
  ierr = VecCreateMPI (PETSC_COMM_WORLD, (PetscInt) pc_ctx->num_local_dofs_coarse_grid, PETSC_DETERMINE, &f);
  CHKERRQ (ierr);
  /* Fill the vector with the current degrees of freedom resulting from the former time step */
  ierr =
    VecSetValues (f, pc_ctx->num_local_dofs_coarse_grid, pc_ctx->global_indexing_coarse_grid,
                  t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem), INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (f);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (f);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side (coarsened fine-level residual)");
  CHKERRQ (ierr);

  /* Create an Approximation Vector storing the coarse grid correection of the residual */
  /* Use the size and allocation similar to vector f */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  ierr = VecCopy (f, u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation Coarse Grid Correction");
  CHKERRQ (ierr);

  /* Setting up matrix-free Matrix */
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context (in this case the preconditioning context) needed by the matrix-free MatVec multiplication */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, (PetscInt) pc_ctx->num_local_dofs_coarse_grid, (PetscInt) pc_ctx->num_local_dofs_coarse_grid,
                    PETSC_DETERMINE, PETSC_DETERMINE, (void *) &pc_pctx, &A);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix multigrid coarse level correction");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_MG_Coarse_Level);
  CHKERRQ (ierr);

  /* Setting up KSP Solver */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);
  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);
  /* Select GMRES as Solver; default: Restarted GMRES(30) */
  ierr = KSPSetType (ksp, KSPGMRES);
  CHKERRQ (ierr);
  /* Select a Preconditioner - without Preconditioning */
  ierr = PCSetType (pc, PCNONE);
  CHKERRQ (ierr);
  /* Set GMRES ieteration tolerances */
  ierr = KSPSetTolerances (ksp, 1.0e-7, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Set the Solver up */
  ierr = KSPSetUp (ksp);
  CHKERRQ (ierr);

  t8dg_debugf ("\nCoarse Grid GMRES got called\n");

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, f, u);
  CHKERRQ (ierr);

  t8dg_debugf ("\nCoarse Grid GMRES solve completed\n");

  /* View whether or not the GMRES did converge */
  ierr = KSPConvergedRateView (ksp, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

  /* Retrieve the local part of the 'solution' Vector u for reading purposes */
  ierr = VecGetArrayRead (u, &approx_coarse_grid_correction);
  CHKERRQ (ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem);

  /* Overwrite the coarse grid degrees of freedom with the newly calculated degrees of freedom containing the residual coarse grid correction */
  for (dof_iter = 0; dof_iter < pc_ctx->num_local_dofs_coarse_grid; ++dof_iter) {
    dof_values_ptr[dof_iter] = (double) approx_coarse_grid_correction[dof_iter];
  }
  /* Restore the array entries */
  ierr = VecRestoreArrayRead (u, &approx_coarse_grid_correction);
  CHKERRQ (ierr);

  /* Free used resources */
  t8dg_dof_values_destroy (&pc_ctx->dof_current_coarse_grid_correction);

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);

}

PetscErrorCode
MatMult_MF_MG_Coarse_Level (Mat A, Vec in, Vec out)
{
  t8dg_preconditioner_mg_ctx_t *pc_ctx;
  PetscErrorCode      ierr;
  PetscScalar        *current_local_approx;
  int                 dof_iter;
  double             *dof_values_ptr;

  /* Get the matrix-free application context equals the multigrid preconditioning context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Retrieve the local part of the 'in' Vector for reading purposes */
  ierr = VecGetArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Get a pointer to the original degrees of freedom of the probblem */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (pc_ctx->dofs_initial_problem);

  /* Overwrite them with the 'in' Vector of the Matrix Multiplication */
  /* Because this MatMult is only used within the GMRES, it is okay to overwrite the original problem dofs without resetting them */
  for (dof_iter = 0; dof_iter < pc_ctx->num_local_dofs_coarse_grid; ++dof_iter) {
    dof_values_ptr[dof_iter] = current_local_approx[dof_iter];
  }

  /* Restore the array entries of the 'in' Vector */
  ierr = VecRestoreArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Calculate the time derivative of the degrees of freedom passed to the MatMult (= 'in' Vector) and store their derivation inside 'dof_current_coarse_grid_correction' */
  (pc_pctx->time_derivative_func) (pc_ctx->dofs_initial_problem, pc_ctx->dof_current_coarse_grid_correction, pc_ctx->current_time_point,
                                   pc_ctx->problem);

  /* Write the result of the matrix-application to the 'in' Vector in the 'out' Vector */
  ierr =
    VecSetValues (out, pc_ctx->num_local_dofs_coarse_grid, pc_ctx->global_indexing_coarse_grid,
                  t8dg_dof_get_double_pointer_to_array (pc_ctx->dof_current_coarse_grid_correction), INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (out);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (out);
  CHKERRQ (ierr);

  return 0;
}

/* This functions performs a given number of richardson iterations */
static void
t8dg_multigrid_smoother_richardson (PC pc, t8dg_dof_values_t * problem_dofs, t8dg_dof_values_t * rhs, int smoothing_steps,
                                    double optimization_parameter)
{
  t8dg_preconditioner_mg_ctx_t *pc_ctx;
  t8dg_dof_values_t  *problem_dofs_derivative;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pc_ctx);
  CHKERRQ (ierr);

  /* Allocate space for derivation of dofs */
  problem_dofs_derivative = t8dg_dof_values_duplicate (problem_dofs);

  /* Iterate over the number of smoothing steps */
  for (iter = 0; iter < smoothing_steps; ++iter) {
    /* Calculate the derivation of the degrees of freedom and store them in problem_dofs_derivative */
    (pc_ctx->time_derivative_func) (problem_dofs, problem_dofs_derivative, pc_ctx->current_time_point, pc_ctx->problem);

    /* Apply one iteration of the richardson method */
    t8dg_dof_values_axpy (-optimization_parameter, problem_dofs_derivative, problem_dofs);
    t8dg_dof_values_axpy (optimization_parameter, rhs, problem_dofs);
  }

  /* Free used resources */
  t8dg_dof_values_destroy (problem_dofs_derivative);
}

#endif

#endif
