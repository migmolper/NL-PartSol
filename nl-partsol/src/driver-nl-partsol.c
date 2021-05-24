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
static void nlpartsol_help_message();
static void globalfree(Mesh,Particle);
static void standard_error(char * Error_message);


int main(int argc, char * argv[])
{  
  char Error_message[MAXW];
  bool Is_New_Simulation = false;
  bool Is_Restart_Simulation = false;
  Mesh FEM_Mesh;
  Particle MPM_Mesh;
  Time_Int_Params Parameters_Solver;

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
//      Is_Restart_Simulation = false;
//      InitialStep = 0;
    }
  else if(argc == 4)
    {      
      Formulation = argv[1];
      SimulationFile = argv[2];
      RestartFile = argv[3];
      Is_New_Simulation = false;
//      Is_Restart_Simulation = true;
//      InitialStep = get_ResultStep(RestartFile);
    }
  else
    {
      sprintf(Error_message,"%s","Wrong inputs, try to tip : nl-partsol --help");
      standard_error(Error_message); 
    }
    
  /* Select kinf of formulation  */
  if(strcmp(Formulation,"-u") == 0)
    {

      NumberDOF = NumberDimensions;

      puts("*************************************************");
      puts("Generating the background mesh ...");
      FEM_Mesh = GramsBox(SimulationFile);

      puts("*************************************************");
      puts("Read solver ...");
      Parameters_Solver = Solver_selector__InOutFun__(SimulationFile,FEM_Mesh.DeltaX);

      puts("*************************************************");
//      if(Is_New_Simulation)
//      {
      puts("Generating new MPM simulation ...");
      MPM_Mesh = GramsSolid(SimulationFile,FEM_Mesh);
//      }
//      if(Is_Restart_Simulation)
//      {
//        puts("Restarting old MPM simulation ...");
//        MPM_Mesh = restart_Simulation(SimulationFile,RestartFile,FEM_Mesh);
//      }

      puts("*************************************************");
      puts("Read outputs ...");
      GramsOutputs(SimulationFile);
      NLPS_Out_nodal_path_csv__InOutFun__(SimulationFile);
      NLPS_Out_particles_path_csv__InOutFun__(SimulationFile);

      puts("*************************************************");
      puts("Run simulation ...");
      if(strcmp(TimeIntegrationScheme,"FE") == 0 )
      { 
        U_Forward_Euler(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"GA") == 0 )
      { 
        U_Generalized_alpha(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"NPC") == 0 )
      { 
        U_Newmark_Predictor_Corrector(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"NPC-FS") == 0 )
      { 
        U_Newmark_Predictor_Corrector_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"Discrete-Energy-Momentum") == 0 )
      { 
        U_Discrete_Energy_Momentum(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"Newmark-beta-Finite-Strains") == 0 )
      {
        U_Newmark_beta_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"Newmark-beta-Finite-Strains-BDB") == 0 )
      {
        U_Newmark_beta_Finite_Strains_BDB(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      } 
      else
      {
        sprintf(Error_message,"%s","Wrong time integration scheme");
        standard_error(Error_message); 
      }

      puts("*************************************************");
      puts("Free memory ...");
      globalfree(FEM_Mesh, MPM_Mesh);       

      printf("Computation finished at : %s \n",__TIME__);  
      puts("Exiting of the program...");
      exit(EXIT_SUCCESS);

    }
    else if(strcmp(Formulation,"-upw") == 0)
    {

      NumberDOF = NumberDimensions + 1;

      puts("*************************************************");
      puts("Generating the background mesh ...");
      FEM_Mesh = GramsBox(SimulationFile);

      puts("*************************************************");
      puts("Read solver ...");
      Parameters_Solver = Solver_selector__InOutFun__(SimulationFile,FEM_Mesh.DeltaX);      

      puts("*************************************************");
      if(Is_New_Simulation)
      {
        puts("Generating new MPM simulation ...");
        MPM_Mesh = Generate_Soil_Water_Coupling_Analysis__InOutFun__(SimulationFile,FEM_Mesh);
      }
      if(Is_Restart_Simulation)
      {
        puts("Restarting old MPM simulation ...");
        MPM_Mesh = restart_Simulation(SimulationFile,RestartFile,FEM_Mesh);
      }

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
      if(strcmp(TimeIntegrationScheme,"NPC-FS") == 0)
      { 
        upw_Newmark_Predictor_Corrector_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else if(strcmp(TimeIntegrationScheme,"Newmark-beta-Finite-Strains") == 0)
      {
        upw_Newmark_beta_Finite_Strains(FEM_Mesh, MPM_Mesh, Parameters_Solver);
      }
      else
      {
        sprintf(Error_message,"%s","Wrong time integration scheme");
        standard_error(Error_message); 
      }

      puts("*************************************************");
      puts("Free memory ...");
      globalfree(FEM_Mesh, MPM_Mesh);       

      printf("Computation finished at : %s \n",__TIME__);  
      puts("Exiting of the program...");
      exit(EXIT_SUCCESS);

    }
    else if(strcmp(Formulation,"-GP") == 0)
    {
      /* Define variable */
      Particle PointAnalysis;

      puts("*************************************************");
      puts("Read solver ...");
      Parameters_Solver = Solver_selector__InOutFun__(SimulationFile, FEM_Mesh.DeltaX);

      puts("*************************************************");
      puts("Generating new Gauss Point simulation ...");
      PointAnalysis = Generate_Gauss_Point_Analysis__InOutFun__(SimulationFile);

      puts("*************************************************");
      puts("Run simulation ...");
      NonLinear_Gauss_Point_Analysis(PointAnalysis);

      puts("*************************************************");
      printf("Computation finished at : %s \n",__TIME__);  
      puts("Exiting of the program...");
      exit(EXIT_SUCCESS);
    }
    else
    {
      sprintf(Error_message,"%s","This formulation has not been yet implemented");
      standard_error(Error_message); 
    }

        
  
   
}

/*********************************************************************/

static void nlpartsol_help_message()
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

puts("Usage : nl-partsol -Flag [commands.nlp]");
puts("Flag values:");
puts(" * -GP  : Gauss Point Analysis");
puts(" * -u   : Displacement formulation");
puts(" * -up  : Velocity-Pressure formulation (Under development)");
puts(" * -upw : Soil-water mixture displacement-pressure formulation (Under development)");
puts(" * -uU  : Soil-water mixture velocity formulation (Under development)");
puts("The creator of NL-PartSol is Miguel Molinos");

puts("mails to : m.molinos@alumnos.upm.es (Madrid-Spain)");

exit(EXIT_SUCCESS);
}


/*********************************************************************/

static void globalfree(Mesh FEM_Mesh, Particle MPM_Mesh)
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

  if(strcmp(Formulation,"-u") == 0)
  {
    free_Fields(MPM_Mesh.Phi);
    free(MPM_Mesh.MatIdx);
  }

  if(strcmp(Formulation,"-upw") == 0)
  {
    free_upw_vars__Fields__(MPM_Mesh.Phi);
    free(MPM_Mesh.MixtIdx);
  }


  if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    free__MatrixLib__(MPM_Mesh.lp);
  }
  /* Lagrange Multipliers / Beta (Only LME ) */
  if(strcmp(ShapeFunctionGP,"LME") == 0){
    free__MatrixLib__(MPM_Mesh.lambda);
    free__MatrixLib__(MPM_Mesh.Beta);
  }
      
}

/***************************************************************************/

static void standard_error(char * Error_message)
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in nl-partsol",Error_message);
    exit(EXIT_FAILURE);
}

/*********************************************************************/
