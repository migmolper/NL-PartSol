/*!
  \mainpage NL-PartSol page

  \section Introduction
  This code is devoted to solve non-linear elastodynamic problems
  using the MPM approach. It is programed and mantained by
  Miguel Molinos Pérez (m.molinos@alumnos.upm.es)

  \section Installation

  \section Usage
 */

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_PETSC
#include <petscksp.h>
#endif

// clang-format off
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "InOutFun.h"
#include "Formulations/Displacements/U-Analisys.h"
#include "Formulations/Displacements/U-Discrete-Energy-Momentum.h"
#include "Formulations/Displacements/U-Forward-Euler.h"
#include "Formulations/Displacements/U-Generalized-Alpha.h"

#ifdef USE_PETSC
#include "Formulations/Displacements/U-Newmark-beta.h"
#endif

#include "Formulations/Displacements/U-Static.h"
#include "Formulations/Displacements/U-Verlet.h"
#include "Formulations/Displacements-Pressure/U-p-Analisys.h"
#include "Formulations/Displacements-Pressure/U-p-Newmark-beta.h"
#include "Formulations/Displacements-WaterPressure/U-pw-Analisys.h"
#include "Formulations/Displacements-WaterPressure/U-pw-Newmark-beta.h"
#include "Formulations/Displacements-WaterPressure/U-pw-Verlet.h"
// clang-format on

/*
  Call global variables
*/
char ShapeFunctionGP[MAXC];
char SimulationFile[MAXC];
char Static_conditons[MAXC];
char Formulation[MAXC];
char *TimeIntegrationScheme;
bool Flag_Print_Convergence;
Load gravity_field;
bool Driver_EigenErosion;
bool Driver_EigenSoftening;

//  Auxiliar functions for the main
static void nlpartsol_help_message();
static void free_nodes(Mesh);
static void free_particles(Particle);
static void standard_error(char *Error_message);

int main(int argc, char *argv[]) {

  char Error_message[10000];
  bool If_formulation = false;
  bool If_f_option = false;
  bool If_ff_option = false;
  int INFO_GramsSolid = 3;
  int STATUS = EXIT_SUCCESS;
  Mesh FEM_Mesh;
  Particle MPM_Mesh;
  Time_Int_Params Parameters_Solver;

  // OpenMP variables
#ifdef USE_OPENMP
  unsigned reqNumThreads = 1;
  const int maxNumThreads = omp_get_max_threads();
#endif

  // Default value for optional modulus
  Driver_EigenErosion = false;

  // Default values for the flags
  Flag_Print_Convergence = false;

  // Read simulation file and kind of simulation
  for (unsigned i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      nlpartsol_help_message();
      return EXIT_SUCCESS;
    }

    if (strcmp(argv[i], "--FORMULATION-U") == 0) {
      strcpy(Formulation, "-u");
      puts("" GREEN "U Formulation" RESET " ...");
      If_formulation = true;
    }

    if (strcmp(argv[i], "--FORMULATION-Up") == 0) {
      puts("" GREEN "U-p Formulation" RESET " ...");
      strcpy(Formulation, "-up");
      If_formulation = true;
    }

    if (strcmp(argv[i], "--FORMULATION-Upw") == 0) {
      puts("" GREEN "U-pw Formulation" RESET " ...");
      strcpy(Formulation, "-upw");
      If_formulation = true;
    }

    if (strcmp(argv[i], "--Fracture-Modulus") == 0) {
      i++;
      if (strcmp(argv[i], "Eigenerosion") == 0) {
        Driver_EigenErosion = true;
        puts("" GREEN "Activate Eigenerosion" RESET " ...");
      } else if (strcmp(argv[i], "Eigensoftening") == 0) {
        puts("" GREEN "Activate Eigensoftening" RESET " ...");
        Driver_EigenSoftening = true;
      } else {
        fprintf(stderr, "" RED "Wrong input for --Fracture-Modulus. Use "
                        "Eigenerosion or Eigensoftening" RESET " \n");
        return EXIT_FAILURE;
      }
    }

#ifdef USE_OPENMP
    if (strcmp(argv[i], "--OPENMP-CORES") == 0) {
      i++;
      reqNumThreads = atoi(argv[i]);
    }
#endif

    if (strcmp(argv[i], "--Print-Convergence") == 0) {
      puts("" GREEN "Print convergence activated" RESET " ...");
      Flag_Print_Convergence = true;
    }

    if (strcmp(argv[i], "-f") == 0) {
      i++;
      strcpy(SimulationFile, argv[i]);
      If_f_option = true;
      break;
    }

    if (strcmp(argv[i], "-ff") == 0) {
      i++;
      strcpy(Static_conditons, argv[i]);
      i++;
      strcpy(SimulationFile, argv[i]);
      If_ff_option = true;
      break;
    }
  }

  if ((If_f_option == false) && (If_ff_option == false)) {
    fprintf(stderr, "" RED "Wrong inputs : non input file" RESET " \n");
    return EXIT_FAILURE;
  }

  if (If_formulation == false) {
    fprintf(stderr, "" RED "Wrong inputs : select formulation " RESET " \n");
    return EXIT_FAILURE;
  }

// Initialize OpenMP
#ifdef USE_OPENMP
  fprintf(stderr, "" GREEN "Initialize OpenMP" RESET " ... \n");
  fprintf(stderr, "\t -> Threads requested : %i \n", reqNumThreads);
  fprintf(stderr, "\t -> Threads availables : %i \n", maxNumThreads);
  omp_set_num_threads(IMIN(maxNumThreads, reqNumThreads));
#endif

  // Initialize PETSc
#ifdef USE_PETSC
  puts("" GREEN "Initialize PETSc" RESET " ...");
  PetscInitialize(&argc, &argv, 0, 0);
#endif

  /* Select kinf of formulation  */
  if (strcmp(Formulation, "-u") == 0) {

    NumberDOF = NumberDimensions;

    if (If_ff_option) {
      puts("*************************************************");
      puts("" GREEN "Read solver" RESET " ...");
      Parameters_Solver = Solver_selector__InOutFun__(Static_conditons);

      puts("*************************************************");
      puts("" GREEN "Generating gravity field" RESET " ...");
      STATUS = Generate_Gravity_Field__InOutFun__(
          &gravity_field, Static_conditons, Parameters_Solver);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                        "Generate_Gravity_Field__InOutFun__() " RESET " \n");
        return EXIT_FAILURE;
      }

      puts("*************************************************");
      puts("" GREEN "Generating the background mesh" RESET " ...");
      FEM_Mesh = GramsBox(Static_conditons, Parameters_Solver);

      puts("*************************************************");
      puts("" GREEN "Generating new MPM simulation" RESET " ...");
      STATUS = Generate_One_Phase_Analysis__InOutFun__(
          &MPM_Mesh, Static_conditons, FEM_Mesh, Parameters_Solver);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED " Error in " RESET "" BOLDRED
                "Generate_One_Phase_Analysis__InOutFun__() " RESET " \n");
        printf("Computation " RED "abnormally" RESET " finished at : %s \n",
               __TIME__);
        puts("Exiting the program...");
        return EXIT_FAILURE;
      }

      puts("*************************************************");
      puts("" GREEN "Read outputs" RESET " ...");
      GramsOutputs(Static_conditons);
      NLPS_Out_nodal_path_csv__InOutFun__(Static_conditons,
                                          Parameters_Solver.NumTimeStep);
      NLPS_Out_particles_path_csv__InOutFun__(Static_conditons,
                                              Parameters_Solver.NumTimeStep);

      puts("*************************************************");
      printf("" GREEN "Start %s shape functions initialisation" RESET " ... \n",
             ShapeFunctionGP);
      initialise_shapefun__MeshTools__(MPM_Mesh, FEM_Mesh);
    }

    if (If_f_option) {
      puts("*************************************************");
      puts("" GREEN "Read solver" RESET " ...");
      puts("*************************************************");
      Parameters_Solver = Solver_selector__InOutFun__(SimulationFile);

      puts("*************************************************");
      puts("" GREEN "Generating gravity field" RESET " ...");
      STATUS = Generate_Gravity_Field__InOutFun__(
          &gravity_field, SimulationFile, Parameters_Solver);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in Generate_Gravity_Field__InOutFun__(,)" RESET
                " \n");
        return EXIT_FAILURE;
      }

      puts("*************************************************");
      puts("" GREEN "Generating the background mesh" RESET " ...");
      puts("*************************************************");
      FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

      puts("*************************************************");
      puts("" GREEN "Generating new MPM simulation" RESET " ...");
      puts("*************************************************");
      STATUS = Generate_One_Phase_Analysis__InOutFun__(
          &MPM_Mesh, SimulationFile, FEM_Mesh, Parameters_Solver);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED " Error in " RESET "" BOLDRED
                "Generate_One_Phase_Analysis__InOutFun__() " RESET " \n");
        puts("*************************************************");
        fprintf(stderr,
                "Computation " RED "abnormally" RESET " finished at : %s \n",
                __TIME__);
        fprintf(stderr, "Exiting the program...\n");
        return EXIT_FAILURE;
      }

      puts("*************************************************");
      puts("" GREEN "Read outputs" RESET " ...");
      GramsOutputs(SimulationFile);

      NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile,
                                          Parameters_Solver.NumTimeStep);
      NLPS_Out_particles_path_csv__InOutFun__(SimulationFile,
                                              Parameters_Solver.NumTimeStep);

      puts("*************************************************");
      printf("Start %s shape functions initialisation ... \n", ShapeFunctionGP);
      initialise_shapefun__MeshTools__(MPM_Mesh, FEM_Mesh);
    }

    puts("*************************************************");
    puts("" GREEN "Run dynamic simulation" RESET " ...");
    if (strcmp(Parameters_Solver.TimeIntegrationScheme, "FE") == 0) {
      U_Forward_Euler(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Generalized-alpha") == 0) {
      U_Generalized_alpha(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme, "NPC-FS") == 0) {
      STATUS = U_Verlet(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Discrete-Energy-Momentum") == 0) {
      U_Discrete_Energy_Momentum(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Newmark-beta-Finite-Strains") == 0) {
#ifdef USE_PETSC
      STATUS = U_Newmark_Beta(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in U_Newmark_Beta(,)" RESET " \n");
      }
#else
      fprintf(stderr,
              "" RED "To use Newmark-beta-Finite-Strains you need to use "
                     "USE_PETSC=true during compilation" RESET " \n");
      STATUS == EXIT_FAILURE;
#endif
    } else {
      sprintf(Error_message, "%s", "Wrong time integration scheme");
      standard_error(Error_message);
    }

  } else if (strcmp(Formulation, "-up") == 0) {

    NumberDOF = NumberDimensions;

    puts("*************************************************");
    puts("" GREEN "Read solver" RESET " ...");
    puts("*************************************************");
    Parameters_Solver = Solver_selector__InOutFun__(SimulationFile);

    puts("*************************************************");
    puts("" GREEN "Generating gravity field" RESET " ...");
    STATUS = Generate_Gravity_Field__InOutFun__(
        &gravity_field, Static_conditons, Parameters_Solver);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                      "Generate_Gravity_Field__InOutFun__() " RESET " \n");
      return EXIT_FAILURE;
    }

    puts("*************************************************");
    puts("" GREEN "Generating the background mesh" RESET " ...");
    puts("*************************************************");
    FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

    puts("*************************************************");
    puts("" GREEN "Generating new MPM simulation" RESET " ...");
    puts("*************************************************");
    STATUS = Generate_One_Phase_Analysis__InOutFun__(
        &MPM_Mesh, SimulationFile, FEM_Mesh, Parameters_Solver);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                      "Generate_One_Phase_Analysis__InOutFun__() " RESET " \n");
      return EXIT_FAILURE;
    }

    puts("*************************************************");
    puts("Read outputs ...");
    GramsOutputs(SimulationFile);
    NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile,
                                        Parameters_Solver.NumTimeStep);
    NLPS_Out_particles_path_csv__InOutFun__(SimulationFile,
                                            Parameters_Solver.NumTimeStep);

    puts("*************************************************");
    puts("Run simulation ...");
    if (strcmp(Parameters_Solver.TimeIntegrationScheme,
               "Newmark-beta-Finite-Strains") == 0) {
      Up_Newmark_beta_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else {
      sprintf(Error_message, "%s", "Wrong time integration scheme");
      standard_error(Error_message);
    }

  } else if (strcmp(Formulation, "-upw") == 0) {

    NumberDOF = NumberDimensions + 1;

    puts("*************************************************");
    puts("Read solver ...");
    Parameters_Solver = Solver_selector__InOutFun__(SimulationFile);

    puts("*************************************************");
    puts("" GREEN "Generating gravity field" RESET " ...");
    STATUS = Generate_Gravity_Field__InOutFun__(
        &gravity_field, Static_conditons, Parameters_Solver);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                      "Generate_Gravity_Field__InOutFun__() " RESET " \n");
      return EXIT_FAILURE;
    }

    puts("*************************************************");
    puts("Generating the background mesh ...");
    FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

    puts("*************************************************");
    puts("Generating new MPM simulation ...");
    MPM_Mesh = Generate_Soil_Water_Coupling_Analysis__InOutFun__(
        SimulationFile, FEM_Mesh, Parameters_Solver, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED " Error in " RESET "" BOLDRED
              "Generate_Soil_Water_Coupling_Analysis__InOutFun__() " RESET
              " \n");
      return EXIT_FAILURE;
    }

    puts("*************************************************");
    puts("Read VTK output directives ...");
    GramsOutputs(SimulationFile);

    puts("*************************************************");
    puts("Read nodal path output directives ...");
    NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile,
                                        Parameters_Solver.NumTimeStep);

    puts("*************************************************");
    puts("Read particle path output directives ...");
    NLPS_Out_particles_path_csv__InOutFun__(SimulationFile,
                                            Parameters_Solver.NumTimeStep);

    puts("*************************************************");
    puts("Run simulation ...");
    if (strcmp(Parameters_Solver.TimeIntegrationScheme, "NPC-FS") == 0) {
      upw_Verlet(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Newmark-beta-Finite-Strains") == 0) {
      upw_Newmark_beta_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else {
      sprintf(Error_message, "%s", "Wrong time integration scheme");
      standard_error(Error_message);
    }

  } else {
    sprintf(Error_message, "%s",
            "This formulation has not been yet implemented");
    standard_error(Error_message);
  }

#ifdef USE_PETSC
  // Finalize PETSc
  PetscFinalize();
#endif

  puts("*************************************************");
  puts("Free memory ...");
  free_nodes(FEM_Mesh);
  free_particles(MPM_Mesh);

  if (STATUS == EXIT_SUCCESS) {
    printf("Computation " GREEN "succesfully" RESET " finished at : %s \n",
           __TIME__);
    puts("Exiting of the program...");
    return EXIT_SUCCESS;
  } else {
    printf("Computation " RED "abnormally" RESET " finished at : %s \n",
           __TIME__);
    puts("Exiting the program...");
    return EXIT_FAILURE;
  }
}

/*********************************************************************/

static void nlpartsol_help_message() {

  puts("Non-Linear Particle Solver (NL-PartSol)");

  /* Read kind of system */
#ifdef __linux__
  puts("Linux version.");
#endif

#ifdef __APPLE__
  puts("Mac OSX version.");
#endif

#ifdef _WIN32
  puts("Windows version.");
#endif

#ifdef USE_PETSC
  puts("Usage : nl-partsol [NLPARTSOL Options] -f [commands.nlp] [PETSc "
       "Options]");
  puts("[NLPARTSOL Options]:");
  puts(" * --Print-Convergence: Display convergence stats.");
  puts(" * --FORMULATION-U : Displacement formulation");
  puts(" * --FORMULATION-Up : Velocity-Pressure formulation");
  puts(" * --FORMULATION-Upw : Soil-water mixture displacement-pressure "
       "formulation");
  puts(" * --Fracture-Modulus Eigenerosion/Eigensoftening : Phase Field "
       "fracture "
       "formulation");
  puts(" * --OPENMP-CORES i : Number of cores for OpenMP (i)");
#else
  puts("Usage : nl-partsol [NLPARTSOL Options] -f [commands.nlp]");
  puts("[NLPARTSOL Options]:");
  puts(" * --Print-Convergence: Display convergence stats.");
  puts(" * --FORMULATION-U : Displacement formulation");
  puts(" * --FORMULATION-Up : Velocity-Pressure formulation");
  puts(" * --FORMULATION-Upw : Soil-water mixture displacement-pressure "
       "formulation");
  puts(" * --Fracture-Modulus Eigenerosion/Eigensoftening : Phase Field "
       "fracture "
       "formulation");
  puts(" * --OPENMP-CORES i : Number of cores for OpenMP (i)");
#endif

  puts("The creator of NL-PartSol is Miguel Molinos");
  puts("mails to : m.molinos@alumnos.upm.es (Madrid-Spain)");
}

/*********************************************************************/

static void free_nodes(Mesh FEM_Mesh)
/*
  Function to free the reamaining memory
*/
{
  /* Free malloc in FEM_Mesh */
  free__MatrixLib__(FEM_Mesh.Coordinates);
  free(FEM_Mesh.NumNodesElem);
  free_table__SetLib__(FEM_Mesh.Connectivity, FEM_Mesh.NumElemMesh);
  free(FEM_Mesh.NumNeighbour);
  free_table__SetLib__(FEM_Mesh.NodeNeighbour, FEM_Mesh.NumNodesMesh);

  free(FEM_Mesh.SizeNodalLocality_0);
  free(FEM_Mesh.SizeNodalLocality);

  free_table__SetLib__(FEM_Mesh.NodalLocality_0, FEM_Mesh.NumNodesMesh);
  free_table__SetLib__(FEM_Mesh.NodalLocality, FEM_Mesh.NumNodesMesh);

  free(FEM_Mesh.ActiveNode);
  free(FEM_Mesh.BoundaryNode);

  if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
    free(FEM_Mesh.Num_Particles_Node);
    free_table__SetLib__(FEM_Mesh.List_Particles_Node, FEM_Mesh.NumNodesMesh);
  }

  if (FEM_Mesh.Locking_Control_Fbar) {
    free(FEM_Mesh.Idx_Patch);
    free(FEM_Mesh.Vol_Patch_n);
    free(FEM_Mesh.Vol_Patch_n1);
  }

  /* FEM_Mesh.Bounds */
}

/*********************************************************************/

static void free_particles(Particle MPM_Mesh)
/*
  Function to free the reamaining memory
*/
{

  /* Free malloc in MPM_Mesh */
  free(MPM_Mesh.I0);
  free(MPM_Mesh.Element_p);
  free(MPM_Mesh.NumberNodes);
  free_table__SetLib__(MPM_Mesh.ListNodes, MPM_Mesh.NumGP);

  if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
    free_table__SetLib__(MPM_Mesh.Beps, MPM_Mesh.NumGP);
  }

  if (strcmp(Formulation, "-u") == 0) {
    free_U_vars__Fields__(MPM_Mesh.Phi);
    free(MPM_Mesh.MatIdx);
  }

  if (strcmp(Formulation, "-up") == 0) {
    free_Up_vars__Fields__(MPM_Mesh.Phi);
    free(MPM_Mesh.MatIdx);
  }

  if (strcmp(Formulation, "-upw") == 0) {
    free_upw_vars__Fields__(MPM_Mesh.Phi);
    free(MPM_Mesh.MixtIdx);
  }

  /* Voxel size (Only uGIMP) */
  if (strcmp(ShapeFunctionGP, "uGIMP") == 0) {
    free__MatrixLib__(MPM_Mesh.lp);
  }

  /* Lagrange Multipliers / Beta (Only LME) */
  if (strcmp(ShapeFunctionGP, "LME") == 0) {
    free__MatrixLib__(MPM_Mesh.lambda);
    free__MatrixLib__(MPM_Mesh.Beta);
  }

  /* Lagrange Multipliers / Beta / Cut_off_Ellipsoid (Only aLME) */
  if (strcmp(ShapeFunctionGP, "aLME") == 0) {
    free__MatrixLib__(MPM_Mesh.lambda);
    free__MatrixLib__(MPM_Mesh.Beta);
    free__MatrixLib__(MPM_Mesh.Cut_off_Ellipsoid);
  }
}

/***************************************************************************/

static void standard_error(char *Error_message) {
  fprintf(stderr, "%s : %s !!! \n", "Error in nl-partsol", Error_message);
  exit(EXIT_FAILURE);
}

/*********************************************************************/
