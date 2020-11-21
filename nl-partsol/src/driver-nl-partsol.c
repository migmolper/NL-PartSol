/*! 
  \mainpage NL-PartSol page
 
  \section Introduction
  This code is devoted to solve non-linear elastodynamic problems
  using the MPM approach. It is programed and mantained by
  Miguel Molinos Pérez (m.molinos@alumnos.upm.es)
  
  \section Installation
  
  \section Usage
 */

#include "nl-partsol.h"

/*
  Call global variables 
*/
char * SimulationFile;
char * RestartFile;
char * TimeIntegrationScheme;
char * Formulation;

/*
  Auxiliar functions for the main
*/
void nlpartsol_help_message();
void globalfree(Mesh,GaussPoint);


int main(int argc, char * argv[])
{  
  int InitialStep;
  bool Is_New_Simulation = false;
  bool Is_Restart_Simulation = false;
  Mesh FEM_Mesh;
  GaussPoint MPM_Mesh;

  /*********************************************************************/
  /************ Read simulation file and kind of simulation ************/
  /*********************************************************************/
  if(argc == 2)
    {      
      if((strcmp(argv[1],"--help") == 0) || (strcmp(argv[1],"-h") == 0))
      { 
        nlpartsol_help_message();
      }      
    }
  if(argc == 3)
    {      
      Formulation = argv[1];
      SimulationFile = argv[2];
      Is_New_Simulation = true;
      Is_Restart_Simulation = false;
      InitialStep = 0;
    }
  else if(argc == 4)
    {      
      Formulation = argv[1];
      SimulationFile = argv[2];
      RestartFile = argv[3];
      Is_New_Simulation = false;
      Is_Restart_Simulation = true;
      InitialStep = get_ResultStep(RestartFile);
    }
  else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in nl-partsol",
	      "Wrong inputs");
      exit(EXIT_FAILURE);
    }
    
  /* Assign degree of freedom */
  if((strcmp(Formulation,"-V") == 0) &&
     (NumberDimensions == 2))
    {
      NumberDOF = 2;
    }
  
  /*********************************************************************/
  /**************** DEFINE CALCULUS MESH and BCCs **********************/
  /*********************************************************************/
  puts("*************************************************");
  puts("Generating the background mesh ...");
  FEM_Mesh = GramsBox(SimulationFile);

  /*********************************************************************/
  /******************* TIME INTEGRATION SCHEME *************************/
  /*********************************************************************/
  puts("*************************************************");
  puts("Read time integration scheme ...");
  GramsTime(SimulationFile);

  /*********************************************************************/
  /******************* DEFINE GAUSS-POINT MESH *************************/
  /*********************************************************************/
  puts("*************************************************");
  if(Is_New_Simulation)
    {
      puts("Generating new MPM simulation ...");
      MPM_Mesh = GramsSolid2D(SimulationFile,FEM_Mesh);
    }
  if(Is_Restart_Simulation)
    {
      puts("Restarting old MPM simulation ...");
      MPM_Mesh = restart_Simulation(SimulationFile,RestartFile,FEM_Mesh);
    }

  /*********************************************************************/
  /*********************** OUTPUT VARIABLES ****************************/
  /*********************************************************************/
  puts("*************************************************");
  puts("Read outputs ...");
  GramsOutputs(SimulationFile);
  NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile);
    
  /*********************************************************************/
  /********************** RUN THE MPM CALCULUS *************************/
  /*********************************************************************/
  puts("*************************************************");
  puts("Run simulation ...");

  /* Forward Euler */
  if(strcmp(TimeIntegrationScheme,"FE") == 0 )
    { 
      U_Forward_Euler(FEM_Mesh, MPM_Mesh, InitialStep);
    }

  /* Generalized-alpha */
  if(strcmp(TimeIntegrationScheme,"GA") == 0 )
    { 
      U_GA(FEM_Mesh, MPM_Mesh, InitialStep);
    }
    
  /* Explicit Newmark predictor-corrector with infinitesimal strains */
  if(strcmp(TimeIntegrationScheme,"NPC") == 0 )
    { 
      U_Newmark_Predictor_Corrector(FEM_Mesh, MPM_Mesh, InitialStep);
    }

  /* Explicit Newmark predictor-corrector with finite strains */
  if(strcmp(TimeIntegrationScheme,"NPC-FS") == 0 )
    { 
      U_Newmark_Predictor_Corrector_Finite_Strains(FEM_Mesh, MPM_Mesh, InitialStep);
    }

  /* Discrete energy momentum method */
  if(strcmp(TimeIntegrationScheme,"Discrete-Energy-Momentum") == 0 )
    { 
      U_Discrete_Energy_Momentum(FEM_Mesh, MPM_Mesh, InitialStep);
    }

  /* Newmark-beta finite strains */
  if(strcmp(TimeIntegrationScheme,"Newmark-beta-Finite-Strains") == 0 )
    {
      U_Newmark_beta_Finite_Strains(FEM_Mesh, MPM_Mesh, InitialStep);
    }  
        
  /*********************************************************************/
  /************************* FREE ALL FIELDS ***************************/
  /*********************************************************************/
  puts("*************************************************");
  puts("Free memory ...");
  globalfree(FEM_Mesh, MPM_Mesh);
  
  /*********************************************************************/
  /******************** CRONOGRAPH CALCULUS : END **********************/
  /*********************************************************************/
  printf("Computation finished at : %s \n",__TIME__);   
         
  /*********************************************************************/
  /************************ CLOSE THE PROGRAM **************************/
  /*********************************************************************/
  puts("Exiting of the program...");
  exit(EXIT_SUCCESS);
   
   
}

/*********************************************************************/

void nlpartsol_help_message()
{

puts("Non-Linear Particle Solver (NL-PartSol)");

  /* Read kind of system */
#ifdef __linux__
  puts("This is the Linux version of the programa.");
#endif    

#ifdef __APPLE__
  puts("This is the Mac OSX version of the program."); 
#endif
  
#ifdef _WIN32 
  puts("This is the Windows version of the program."); 
#endif

puts("Usage : nl-partsol -V commands.nlp");

puts("The creator of NL-PartSol is Miguel Molinos");

puts("mails to m.molinos@alumnos.upm.es (Madrid-Spain)");

exit(EXIT_SUCCESS);
}


/*********************************************************************/

void globalfree(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*
  Function to free the reamaining memory
*/
{
  /* Free malloc in FEM_Mesh */
  free__MatrixLib__(FEM_Mesh.Coordinates);
  free(FEM_Mesh.NumNodesElem);
  free_table__SetLib__(FEM_Mesh.Connectivity,FEM_Mesh.NumElemMesh);
  free(FEM_Mesh.NumNeighbour);
  free_table__SetLib__(FEM_Mesh.NodeNeighbour,FEM_Mesh.NumNodesMesh);
  free(FEM_Mesh.SizeNodalLocality);
  free_table__SetLib__(FEM_Mesh.NodalLocality,FEM_Mesh.NumNodesMesh);
  free(FEM_Mesh.NumParticles);
  free_table__SetLib__(FEM_Mesh.I_particles,FEM_Mesh.NumNodesMesh);
  
  /* FEM_Mesh.Bounds */

  /* Free malloc in MPM_Mesh */
  free(MPM_Mesh.I0);
  free(MPM_Mesh.NumberNodes);
  free_table__SetLib__(MPM_Mesh.ListNodes,MPM_Mesh.NumGP);
  free_table__SetLib__(MPM_Mesh.Beps,MPM_Mesh.NumGP);
  free_Fields(MPM_Mesh.Phi);
  /* free_Fields(MPM_Mesh.Phi_n0);   */
  free(MPM_Mesh.MatIdx);

  /* Material * Mat; */

  /* Load * F; */
  /* Load * B; */

  if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    free__MatrixLib__(MPM_Mesh.lp);
  }
  /* Lagrange Multipliers / Beta (Only LME ) */
  if(strcmp(ShapeFunctionGP,"LME") == 0){
    free__MatrixLib__(MPM_Mesh.lambda);
    free__MatrixLib__(MPM_Mesh.Beta);
  }
      
}

/*********************************************************************/
