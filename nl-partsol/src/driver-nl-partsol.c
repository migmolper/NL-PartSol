/*!
  \mainpage NL-PartSol page

  \section Introduction
  This code is devoted to solve non-linear elastodynamic problems
  using the MPM approach. It is programed and mantained by
  Miguel Molinos Pérez (m.molinos@alumnos.upm.es)

  \section Installation

  \section Usage
 */

#include <string.h>
#include "nl-partsol.h"
#include "petscsys.h"

/*
  Call global variables
*/
char *SimulationFile;
char *Static_conditons;
char *TimeIntegrationScheme;
char *Formulation;

static char help[] = "Appends to an ASCII file.\n\n";

/*
  Auxiliar functions for the main
*/
static void nlpartsol_help_message();
static void free_nodes(Mesh);
static void free_particles(Particle);
static void standard_error(char *Error_message);

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char*)0,help);

  char Error_message[MAXW];
  bool Is_Static_Initialization = false;
  bool Is_Restart_Simulation = false;
  int INFO_GramsSolid = 3;
  int STATUS = EXIT_SUCCESS;
  Mesh FEM_Mesh;
  Particle MPM_Mesh;
  Time_Int_Params Parameters_Solver;

  /*********************************************************************/
  /************ Read simulation file and kind of simulation ************/
  /*********************************************************************/
  if (argc == 2) {
    if ((strcmp(argv[1], "--help") == 0) || (strcmp(argv[1], "-h") == 0)) {
      nlpartsol_help_message();
    }
  }
  if (argc == 3) {
    Formulation = argv[1];
    SimulationFile = argv[2];
    Is_Static_Initialization = false;
  } else if (argc == 4) {
    Formulation = argv[1];
    Static_conditons = argv[2];
    SimulationFile = argv[3];
    Is_Static_Initialization = true;
  } else {
    sprintf(Error_message, "%s",
            "Wrong inputs, try to tip : nl-partsol --help");
    standard_error(Error_message);
  }

  /* Select kinf of formulation  */
  if (strcmp(Formulation, "-u") == 0) {

    NumberDOF = NumberDimensions;

    if (Is_Static_Initialization) {
      puts("*************************************************");
      puts("Read solver ...");
      Parameters_Solver = Solver_selector__InOutFun__(Static_conditons);

      puts("*************************************************");
      puts("Generating the background mesh ...");
      FEM_Mesh = GramsBox(Static_conditons, Parameters_Solver);

      puts("*************************************************");
      puts("Generating new MPM simulation ...");
      MPM_Mesh = Generate_One_Phase_Analysis__InOutFun__(
          Static_conditons, FEM_Mesh, Parameters_Solver);

      puts("*************************************************");
      puts("Read outputs ...");
      GramsOutputs(Static_conditons);
      NLPS_Out_nodal_path_csv__InOutFun__(Static_conditons);
      NLPS_Out_particles_path_csv__InOutFun__(Static_conditons);

      puts("*************************************************");
      printf("Start %s shape functions initialisation ... \n", ShapeFunctionGP);
      initialise_shapefun__MeshTools__(MPM_Mesh, FEM_Mesh);

    } else {
      puts("*************************************************");
      puts("Read solver ...");
      Parameters_Solver = Solver_selector__InOutFun__(SimulationFile);

      puts("*************************************************");
      puts("Generating the background mesh ...");
      FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

      puts("*************************************************");
      puts("Generating new MPM simulation ...");
      MPM_Mesh = Generate_One_Phase_Analysis__InOutFun__(
          SimulationFile, FEM_Mesh, Parameters_Solver);

      puts("*************************************************");
      puts("Read outputs ...");
      GramsOutputs(SimulationFile);
      NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile);
      NLPS_Out_particles_path_csv__InOutFun__(SimulationFile);

      puts("*************************************************");
      printf("Start %s shape functions initialisation ... \n", ShapeFunctionGP);
      initialise_shapefun__MeshTools__(MPM_Mesh, FEM_Mesh);
    }

    if (Is_Static_Initialization) {
      U_Static_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);

      puts("*************************************************");
      puts("Read solver ...");
      Parameters_Solver = Solver_selector__InOutFun__(SimulationFile);

      puts("*************************************************");
      puts("Generating the background mesh ...");
      free_nodes(FEM_Mesh);
      FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

      puts("*************************************************");
      printf("Start %s shape functions initialisation ... \n", ShapeFunctionGP);

      for (int p = 0; p < MPM_Mesh.NumGP; p++) {
        free__SetLib__(&MPM_Mesh.ListNodes[p]);
        MPM_Mesh.ListNodes[p] = NULL;
      }

      initialise_shapefun__MeshTools__(MPM_Mesh, FEM_Mesh);
    }

    puts("*************************************************");
    puts("Run dynamic simulation ...");
    if (strcmp(Parameters_Solver.TimeIntegrationScheme, "FE") == 0) {
      U_Forward_Euler(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Generalized-alpha") == 0) {
      U_Generalized_alpha(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme, "NPC") == 0) {
      U_Newmark_Predictor_Corrector(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme, "NPC-FS") == 0) {
      STATUS = U_Newmark_Predictor_Corrector_Finite_Strains(FEM_Mesh, MPM_Mesh,
                                                   Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Discrete-Energy-Momentum") == 0) {
      U_Discrete_Energy_Momentum(FEM_Mesh, MPM_Mesh, Parameters_Solver);
    } else if (strcmp(Parameters_Solver.TimeIntegrationScheme,
                      "Newmark-beta-Finite-Strains") == 0) {
      STATUS = U_Newmark_beta_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in U_Newmark_beta_Finite_Strains(,)"RESET" \n");
      }      
    } else {
      sprintf(Error_message, "%s", "Wrong time integration scheme");
      standard_error(Error_message);
    }

  } else if (strcmp(Formulation, "-up") == 0) {

    NumberDOF = NumberDimensions;

    puts("*************************************************");
    puts("Read solver ...");
    Parameters_Solver = Solver_selector__InOutFun__(SimulationFile);

    puts("*************************************************");
    puts("Generating the background mesh ...");
    FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

    puts("*************************************************");
    puts("Generating new MPM simulation ...");
    MPM_Mesh = Generate_One_Phase_Analysis__InOutFun__(SimulationFile, FEM_Mesh,
                                                       Parameters_Solver);

    puts("*************************************************");
    puts("Read outputs ...");
    GramsOutputs(SimulationFile);
    NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile);
    NLPS_Out_particles_path_csv__InOutFun__(SimulationFile);

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
    puts("Generating the background mesh ...");
    FEM_Mesh = GramsBox(SimulationFile, Parameters_Solver);

    puts("*************************************************");
    puts("Generating new MPM simulation ...");
    MPM_Mesh = Generate_Soil_Water_Coupling_Analysis__InOutFun__(
        SimulationFile, FEM_Mesh, Parameters_Solver);

    puts("*************************************************");
    puts("Read VTK output directives ...");
    GramsOutputs(SimulationFile);

    puts("*************************************************");
    puts("Read nodal path output directives ...");
    NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile);

    puts("*************************************************");
    puts("Read particle path output directives ...");
    NLPS_Out_particles_path_csv__InOutFun__(SimulationFile);

    puts("*************************************************");
    puts("Run simulation ...");
    if (strcmp(Parameters_Solver.TimeIntegrationScheme, "NPC-FS") == 0) {
      upw_Newmark_Predictor_Corrector_Finite_Strains(FEM_Mesh, MPM_Mesh,
                                                     Parameters_Solver);
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

  
  PetscFinalize();

  puts("*************************************************");
  puts("Free memory ...");
  free_nodes(FEM_Mesh);
  free_particles(MPM_Mesh);  

  if(STATUS == EXIT_SUCCESS)
  {
    printf("Computation "GREEN"succesfully"RESET" finished at : %s \n", __TIME__);
    puts("Exiting of the program...");
    return EXIT_SUCCESS;
  }
  else
  {
    printf("Computation "RED"abnormally"RESET" finished at : %s \n", __TIME__);
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

  puts("Usage : nl-partsol -Flag [commands.nlp]");
  puts("Flag values:");
  puts(" * -u   : Displacement formulation");
  puts(" * -up  : Velocity-Pressure formulation");
  puts(" * -upw : Soil-water mixture displacement-pressure formulation");
  puts(" * -uU  : Soil-water mixture velocity formulation (Not developed yet)");
  puts("The creator of NL-PartSol is Miguel Molinos");

  puts("mails to : m.molinos@alumnos.upm.es (Madrid-Spain)");

  exit(EXIT_SUCCESS);
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

  free(FEM_Mesh.Num_Particles_Node);
  free_table__SetLib__(FEM_Mesh.List_Particles_Node, FEM_Mesh.NumNodesMesh);
  //  free(FEM_Mesh.Num_Particles_Element);
  //  free_table__SetLib__(FEM_Mesh.List_Particles_Element,FEM_Mesh.NumElemMesh);

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
  free_table__SetLib__(MPM_Mesh.Beps, MPM_Mesh.NumGP);

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
