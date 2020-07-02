/*! 
  \mainpage NL-PartSol page
 
  \section Introduction
  This code is devoted to solve non-linear elastodynamic problems
  using the MPM approach. It is coded and mantained by
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

void globalfree(Mesh,GaussPoint);

int main(int argc, char * argv[])
{  
  /* Read kind of system */
  puts("*************************************************\n");
#ifdef linux
  puts("UNIX system");
#endif		

#ifdef _WIN32 
  puts("WIN32 system");	
#endif

  int InitialStep;
  bool Is_New_Simulation;
  bool Is_Restart_Simulation;
  Mesh FEM_Mesh;
  GaussPoint MPM_Mesh;

  /*********************************************************************/
  /************ Read simulation file and kind of simulation ************/
  /*********************************************************************/
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
    
  /*********************************************************************/
  /********************** RUN THE MPM CALCULUS *************************/
  /*********************************************************************/
  puts("*************************************************");
  puts("Run simulation ...");

  /* Forward Euler */
  if(strcmp(TimeIntegrationScheme,"FE") == 0 )
    { 
      U_FE(FEM_Mesh, MPM_Mesh, InitialStep);
    }

  /* Generalized-alpha */
  if(strcmp(TimeIntegrationScheme,"GA") == 0 )
    { 
      U_GA(FEM_Mesh, MPM_Mesh, InitialStep);
    }
    
  /* Explicit predictor-corrector */
  if(strcmp(TimeIntegrationScheme,"PCE") == 0 )
    { 
      U_PCE(FEM_Mesh, MPM_Mesh, InitialStep);
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

void globalfree(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*
  Function to free the reamaining memory
*/
{
  /* Free malloc in FEM_Mesh */
  FreeMat(FEM_Mesh.Coordinates);
  free(FEM_Mesh.NumNodesElem);
  free_SetTable(FEM_Mesh.Connectivity,FEM_Mesh.NumElemMesh);
  free(FEM_Mesh.NumNeighbour);
  free_SetTable(FEM_Mesh.NodeNeighbour,FEM_Mesh.NumNodesMesh);
  free(FEM_Mesh.SizeNodalLocality);
  free_SetTable(FEM_Mesh.NodalLocality,FEM_Mesh.NumNodesMesh);
  free(FEM_Mesh.NumParticles);
  free_SetTable(FEM_Mesh.I_particles,FEM_Mesh.NumNodesMesh);
  
  /* FEM_Mesh.Bounds */

  /* Free malloc in MPM_Mesh */
  free(MPM_Mesh.I0);
  free(MPM_Mesh.NumberNodes);
  free_SetTable(MPM_Mesh.ListNodes,MPM_Mesh.NumGP);
  free_SetTable(MPM_Mesh.Beps,MPM_Mesh.NumGP);
  free_Fields(MPM_Mesh.Phi);
  /* free_Fields(MPM_Mesh.Phi_n0);   */
  free(MPM_Mesh.MatIdx);

  /* Material * Mat; */

  /* Load * F; */
  /* Load * B; */

  if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    FreeMat(MPM_Mesh.lp);
  }
  /* Lagrange Multipliers / Beta (Only LME ) */
  if(strcmp(ShapeFunctionGP,"LME") == 0){
    FreeMat(MPM_Mesh.lambda);
    FreeMat(MPM_Mesh.Beta);
  }
      
}
