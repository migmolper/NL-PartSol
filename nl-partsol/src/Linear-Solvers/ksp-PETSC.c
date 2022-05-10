
#ifdef USE_PETSC
#include "Linear-Solvers/ksp-PETSC.h"

int krylov_PETSC(Mat *ptr_Tangent_Stiffness, Vec *ptr_Residual,
                 unsigned Nactivedofs) {

  Mat Tangent_Stiffness = *ptr_Tangent_Stiffness;
  Vec Residual = *ptr_Residual;

  // Solver start
  //    auto solver_begin = std::chrono::steady_clock::now();
  //    if (verbosity_ > 0 && mpi_rank == 0)
  //      console_->info("Type: \"{}\", Preconditioner: \"{}\", Begin!",
  //                     sub_solver_type_, preconditioner_type_);

  // Initialize PETSC matrix, vectors, and parameters
  KSP solver;
  KSPConvergedReason reason;
  PC pc;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  KSPCreate(PETSC_COMM_WORLD, &solver);

  KSPSetOperators(solver, Tangent_Stiffness, Tangent_Stiffness);

  /*     // Set tolerance
      KSPSetTolerances(solver, tolerance_, abs_tolerance_, div_tolerance_,
                       max_iter_); */

  KSPGetPC(solver, &pc);
  PCSetType(pc, PCJACOBI); // Default

  // Set solver_type
  KSPSetFromOptions(solver);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  KSPSolve(solver, Residual, Residual);
  KSPGetConvergedReason(solver, &reason);
  /*
      // Print residual in each iteration
      if (1) {
        PetscViewerAndFormat* vf;
        PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,
                                   PETSC_VIEWER_DEFAULT, &vf);
        PetscInt its;
        KSPGetIterationNumber(solver, &its);
        PetscPrintf(PETSC_COMM_WORLD, "\nConvergence in %d iterations.\n",
                    (int)its);
        PetscReal rnorm;
        if (1) {
          for (int i = 0; i < its; i++) {
            KSPMonitorTrueResidual(solver, i, rnorm, vf);
          }
        }
      }
  */

  /*     // Warn if solver does not converge
      if (reason < 0) {
        PetscPrintf(MPI_COMM_WORLD,
                    "\nKrylov PETSC solver \"%s\" with \"%s\" preconditioner "
                    "DIVERGED, try to modify the preconditioner, set tolerance "
                    "and maximum iteration.\n",
                    sub_solver_type_.c_str(), preconditioner_type_.c_str());
      } */

  /*     // Scatter and gather for cloning procedure
      if (mpi_size > 1) {
        // Initiate scatter arrays
        VecScatter ctx;
        Vec x_seq;
        PetscScalar* x_data;
        VecScatterCreateToAll(Residual, &ctx, &x_seq);
        VecScatterBegin(ctx, Residual, x_seq, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, Residual, x_seq, INSERT_VALUES, SCATTER_FORWARD);
        VecGetArray(x_seq, &x_data);

        // Copy petsc x to Eigen x
        for (unsigned i = 0; i < x.size(); i++) {
          const int global_index = rank_global_mapper_[i];
          x(i) = x_data[global_index];
        }

        // Destroy scatter arrays
        VecRestoreArray(x_seq, &x_data);
        VecScatterDestroy(&ctx);
        VecDestroy(&x_seq);
      } else {
        PetscScalar value;
        // Copy petsc x to Eigen x
        for (unsigned i = 0; i < x.size(); i++) {
          const int global_index = rank_global_mapper_[i];
          VecGetValues(Residual, 1, &global_index, &value);
          x(i) = value;
        }
      } */

  // Free work space
  KSPDestroy(&solver);
  /*
      // Solver End
      auto solver_end = std::chrono::steady_clock::now();
      if (verbosity_ > 0 && mpi_rank == 0)
        console_->info(
            "Type: \"{}\", Preconditioner: \"{}\", End! Duration: {} ms.",
            sub_solver_type_, preconditioner_type_,
            std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                                  solver_begin)
                .count()); */


  return 0;
}

#endif