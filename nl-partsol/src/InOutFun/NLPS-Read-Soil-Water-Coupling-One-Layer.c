#include "nl-partsol.h"


/*
  Call global variables 
*/
char * MPM_MeshFileName;

typedef struct
{
  bool Is_Soil_Water_Coupling;
  bool Is_ParticlesMesh;
  bool Is_GramsShapeFun;
  bool Is_GramsMaterials;
  bool Is_Particle_Initial;
  bool Is_Nodal_Initial;
  bool Is_GramsBodyForces;
  bool Is_GramsNeumannBC;

  int Counter_Materials;
  int Counter_BodyForces;
  int Counter_GramsNeumannBC;

} Simulation_Key_Words;

typedef struct
{

  int   GPxElement;
  char  Route_Mesh[MAXC];

} Mesh_Parameters;

/*
  Auxiliar functions and variables
*/
#ifdef _WIN32
static char * delimiters_1 = " (,)\r\n\t";
static char * delimiters_2 = " =\t\r\n"; 
#else
static char * delimiters_1 = " (,)\n\t";
static char * delimiters_2 = " =\t\n"; 
#endif
static char * delimiters_3 = "=";
static char * delimiters_4 = ";";

static char Error_message[MAXW];

static Simulation_Key_Words Check_Simulation_Key_Words(char *);
static Mesh_Parameters Read_Mesh_Parameters(char *);
static void initialise_2D_particles();
static bool Check_File(char *);
static void standard_error();
static FILE * Open_and_Check_simulation_file(char *);

/*********************************************************************/

GaussPoint Generate_Soil_Water_Coupling_Analysis__InOutFun__(char * Name_File, Mesh FEM_Mesh)
/*
 */
{
  int Ndim = NumberDimensions;
  int NumParticles;
  
  /* Parser num chars */
  int Num_words_parse;

  /* Special variables Soil_Water_Coupling_One_Layer */
  char * Parse_Soil_Water_Coupling_One_Layer[MAXW] = {NULL};
  char * Parse_Mesh_id[MAXW] = {NULL};
  char * Parse_Mesh_Properties[MAXW] = {NULL};

  Mesh MPM_GID_Mesh; /* GID mesh to define the material point mesh */
  GaussPoint MPM_Mesh; /* Material point mesh (Gauss-Points) */
  Simulation_Key_Words Sim_Params; /* Auxiliar varible to check key words */
  Mesh_Parameters Msh_Parms; /* Auxiliar variable to read mesh parameters */


  Sim_Params = Check_Simulation_Key_Words(Name_File);

  /*
    Define materials
  */
  puts("*************************************************");
  printf(" \t %s \n","* Read materials properties");
  if(Sim_Params.Is_GramsMaterials)
  {
    MPM_Mesh.NumberMaterials = Sim_Params.Counter_Materials;
    MPM_Mesh.Mat = Read_Materials__InOutFun__(Name_File, MPM_Mesh.NumberMaterials);
  }
  else
  {
    sprintf(Error_message,"%s","No materials were defined");
    standard_error(); 
  }
  printf("\t %s \n","DONE !!");

  
  /*
    Define particles mesh 
  */
  if(Sim_Params.Is_Soil_Water_Coupling)
  {
    
    /*
      Read particle mesh preliminar information
    */
    Msh_Parms = (Name_File);

    /*
      Read particles mesh 
    */
    MPM_GID_Mesh = ReadGidMesh(Msh_Parms.Route_Mesh);

    /*
      Define the number of particles
    */
    NumParticles = Msh_Parms.GPxElement*MPM_GID_Mesh.NumElemMesh;
    MPM_Mesh.NumGP = NumParticles;

    /*  
      Closest node to the particle 
    */
    MPM_Mesh.I0 = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* 
      Number of tributary nodes for each particle 
    */
    MPM_Mesh.NumberNodes = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* 
      Tributary nodes for each particle 
    */
    MPM_Mesh.ListNodes = alloc_table__SetLib__(NumParticles);

    /* 
      List of particles close to each particle 
    */
    MPM_Mesh.Beps = alloc_table__SetLib__(NumParticles);

    /*
      Define shape functions
    */
    if(Sim_Params.Is_GramsShapeFun)
    {

      /*
        Read Shape functions parameters 
      */
      GramsShapeFun(Name_File);
      /* 
        Lenght of the Voxel (Only GIMP) 
      */
      if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
      {
        MPM_Mesh.lp = allocZ__MatrixLib__(NumParticles,Ndim);
        strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
      }
      /*
         Lagrange Multipliers / Beta (Only LME ) 
      */
      if(strcmp(ShapeFunctionGP,"LME") == 0)
      {
        MPM_Mesh.lambda = allocZ__MatrixLib__(NumParticles,Ndim);
        strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
        MPM_Mesh.Beta = allocZ__MatrixLib__(NumParticles,Ndim);
        strcpy(MPM_Mesh.Beta.Info,"Beta parameter");
      }
      /*
         Anisotropic lagrange Multipliers / Beta (Only LME ) 
      */
      if(strcmp(ShapeFunctionGP,"aLME") == 0)
      {
        MPM_Mesh.lambda = allocZ__MatrixLib__(NumParticles,Ndim);
        strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
        MPM_Mesh.Beta = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
        strcpy(MPM_Mesh.Beta.Info,"Beta tensor");
        MPM_Mesh.Cut_Off_Ellipsoid = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
        strcpy(MPM_Mesh.Cut_Off_Ellipsoid.Info,"Cut-Off Ellipsoid");
      }
    }
    else
    {
      sprintf(Error_message,"%s","GramsShapeFun no defined");
      standard_error(); 
    }
    
    /*
      Allocate vectorial/tensorial fields 
    */
    MPM_Mesh.Phi = allocate_Fields(NumParticles);

    /*
      Assign material for each material point
    */
    MPM_Mesh.MatIdx = 

    /*
      Initialise particle 
    */    
    initial_position__Particles__(MPM_Mesh.Phi.x_GC,MPM_GID_Mesh,Msh_Parms.GPxElement); 
    if(Ndim == 2)
    {
      initialise_2D_particles(MPM_GID_Mesh,MPM_Mesh,Msh_Parms.GPxElement);
    }
     
    /*
      Initialise shape functions 
    */
    puts("*************************************************");
    if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){
      printf("\t * %s \n","Initialize MPMQ4 shape functions ...");
      initialize__Q4__(MPM_Mesh, FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
      printf("\t * %s \n","Initialize uGIMP shape functions ...");      
      initialize__GIMP__(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"LME") == 0){
      printf("\t * %s \n","Initialize LME shape functions ...");
      initialize__LME__(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"aLME") == 0){
      printf("\t * %s \n","Initialize aLME shape functions ...");
      initialize__aLME__(MPM_Mesh,FEM_Mesh);
    }
    printf("\t %s \n","DONE !!");


    /*
      Read initial values 
    */    
    if(Is_Particle_Initial)
    {
      Initial_condition_particles__InOutFun__(Name_File,MPM_Mesh,Msh_Parms.GPxElement);
    }
    else if(Is_Nodal_Initial)
    {
      Initial_condition_nodes__InOutFun__(Name_File,MPM_Mesh,FEM_Mesh);
    }
    else{
      puts("*************************************************");
      puts(" No initial conditions defined ");
    }

    /*
      Read external forces 
    */    
    if(Is_GramsNeumannBC)
    {
      MPM_Mesh.NumNeumannBC = Counter_GramsNeumannBC;
      MPM_Mesh.F = GramsNeumannBC(Name_File,Counter_GramsNeumannBC,Msh_Parms.GPxElement);
    }
    else
    {
      MPM_Mesh.NumNeumannBC = Counter_GramsNeumannBC;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
       "* No Neumann boundary conditions defined in",
       Name_File);
    }

    /*
      Read body forces 
    */    
    if(Is_GramsBodyForces)
    {
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      MPM_Mesh.B = GramsBodyForces(Name_File,Counter_BodyForces,Msh_Parms.GPxElement); 
    }
    else{
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
       "* No body forces defined in",Name_File);
    }
   
    /*
      Free the input data 
    */
    for(int i = 0 ; i<MPM_GID_Mesh.NumElemMesh ; i++)
      {
        free__SetLib__(&MPM_GID_Mesh.Connectivity[i]); 
      }   
    free(MPM_GID_Mesh.Connectivity);
    free__MatrixLib__(MPM_GID_Mesh.Coordinates);
    free(MPM_GID_Mesh.NumParticles);
    free(MPM_GID_Mesh.NumNeighbour);
    free(MPM_GID_Mesh.NodeNeighbour);

  } 
  else
  {
    fprintf(stderr,"%s : %s \n",
	    "Error in Soil-Water-Coupling-One-Layer()",
	    "Mesh file name and number of particles are required");
    exit(EXIT_FAILURE);
  }
  
  return MPM_Mesh;
}

/***************************************************************************/

static Simulation_Key_Words Check_Simulation_Key_Words(char * Name_File)
{

  Simulation_Key_Words Sim_Key_Wrds;
  char Line[MAXC] = {0}; 
  char * Words[MAXW] = {NULL};
  int Num_words = 0;
  int Num_line = 0;

  /*
    Default values
  */
  Sim_Key_Wrds.Is_Soil_Water_Coupling = false;
  Sim_Key_Wrds.Is_ParticlesMesh = false;
  Sim_Key_Wrds.Is_GramsShapeFun = false;
  Sim_Key_Wrds.Is_GramsMaterials = false;
  Sim_Key_Wrds.Is_Particle_Initial = false;
  Sim_Key_Wrds.Is_Nodal_Initial = false;
  Sim_Key_Wrds.Is_GramsBodyForces = false;
  Sim_Key_Wrds.Is_GramsNeumannBC = false;

  Sim_Key_Wrds.Counter_Materials = 0;
  Sim_Key_Wrds.Counter_BodyForces = 0;
  Sim_Key_Wrds.Counter_GramsNeumannBC = 0;
  

  /*
    Open and check file
  */
  FILE * Sim_dat = Open_and_Check_simulation_file(Name_File);
  

  while(fgets(Line,sizeof(Line),Sim_dat) != NULL)
  {

    /* Read the line with the space as separators */
    Num_words = parse(Words,Line,delimiters_1);

    /*
      Update line number
    */
    Num_line++;
   
    if (Num_words < 0)
    {
      sprintf(Error_message,"%s : %i","Parser failed in line",Num_line);
      standard_error(); 
    }
    else if (strcmp(Words[0],"Soil-Water-Coupling-One-Layer") == 0)
    {
      Sim_Key_Wrds.Is_Soil_Water_Coupling = true;
    }
    else if (strcmp(Words[0],"GramsShapeFun") == 0)
    {
      Sim_Key_Wrds.Is_GramsShapeFun = true;
    }
    else if (strcmp(Words[0],"GramsMaterials") == 0)
    {
      Sim_Key_Wrds.Is_GramsMaterials = true;
      Sim_Key_Wrds.Counter_Materials++;
    }
    else if (strcmp(Words[0],"GramsInitials") == 0)
    {
      Sim_Key_Wrds.Is_Particle_Initial = true;
    }
    else if (strcmp(Words[0],"Initial-nodal-values") == 0)
    {
      Sim_Key_Wrds.Is_Nodal_Initial = true;
    }
    else if (strcmp(Words[0],"GramsBodyForces") == 0)
    {
      Sim_Key_Wrds.Is_GramsBodyForces = true;
      Sim_Key_Wrds.Counter_BodyForces++;
    }
    else if (strcmp(Words[0],"GramsNeumannBC") == 0)
    {
      Sim_Key_Wrds.Is_GramsNeumannBC = true;
      Sim_Key_Wrds.Counter_GramsNeumannBC++;
    }    
  }

  /*
   Close  file 
  */
  fclose(Sim_dat);


  return Sim_Key_Wrds;
}

/***************************************************************************/

static Mesh_Parameters Read_Mesh_Parameters(char * Name_File)
{

  Mesh_Parameters Msh_Params;
  char Line[MAXC] = {0}; 
  char * Words[MAXW] = {NULL};
  char * File_Parameter[MAXW] = {NULL};
  char * GPxElement_Parameter[MAXW] = {NULL};
  int Num_words = 0;
  int Num_parameters = 0;
  int Num_line = 0;

  /*
    Open and check file
  */
  FILE * Sim_dat = Open_and_Check_simulation_file(Name_File);

  while(fgets(Line,sizeof(Line),Sim_dat) != NULL)
  {

    /*
     Read the line with the space as separators 
    */
    Num_words = parse(Words,Line,delimiters_1);

    /*
      Update line number
    */
    Num_line++;

    if ((strcmp(Words[0],"Soil-Water-Coupling-One-Layer") == 0) && (Num_words >= 3))
    {

     Num_parameters = parse(File_Parameter,Words[1],delimiters_3);
     if ((strcmp(File_Parameter[0],"File") == 0) && (Num_parameters == 2))
     {
        MPM_MeshFileName = File_Parameter[1];
        generate_route(GPxElement.Route_Mesh,Name_File);
        strcat(GPxElement.Route_Mesh,MPM_MeshFileName);
     }
     else
     {
        sprintf(Error_message,"%s (Line %i)","File parameter missed",Num_line);
        standard_error(); 
     }

      Num_parameters = parse(GPxElement_Parameter,Words[2],delimiters_3);
     if ((strcmp(GPxElement_Parameter[0],"GPxElement") == 0) && (Num_parameters == 2))
     {
        GPxElement.GPxElement = atoi(GPxElement_Parameter[1]);
     }
     else
     {
        sprintf(Error_message,"%s (Line %i)","GPxElement parameter missed",Num_line);
        standard_error(); 
     }

    }
    else
    {
      sprintf(Error_message,"%s (Line %i)","Insuficient number of parameters",Num_line);
      standard_error(); 
    }

  }

  /*
   Close  file 
  */
  fclose(Sim_dat);


  return Msh_Params;
}

/***************************************************************************/

static int * assign_material_to_particles(char * Name_File, int NumMaterials, int NumParticles, int GPxElement)
{
  int * MatIdx = (int *)calloc(NumParticles,sizeof(int));
  char Line[MAXC] = {0}; 
  char * Words[MAXW] = {NULL};
  char * File_Parameter[MAXW] = {NULL};
  char * MatIdx_Parameter[MAXW] = {NULL};
  int Num_words = 0;
  int Num_parameters = 0;
  int Num_line = 0;
 
  /*
    Open and check file
  */
  FILE * Sim_dat = Open_and_Check_simulation_file(Name_File);


  while(fgets(Line,sizeof(Line),Sim_dat) != NULL)
  {

    /*
     Read the line with the space as separators 
    */
    Num_words = parse(Words,Line,delimiters_1);

    /*
      Update line number
    */
    Num_line++;

    if ((strcmp(Words[0],"Assign-material-to-particles") == 0) && (Num_words >= 3))
    {

//      /* Read file with the nodes */
//      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Mat_id[1]);
//      printf("\t -> %s : %s \n","Material points",FileNodesRoute);

//      /* Get an array with the nodes */
//      Chain_Nodes = File2Chain(FileNodesRoute);
//      Num_Nodes = lenght__SetLib__(Chain_Nodes);
//      Array_Nodes = set_to_memory__SetLib__(Chain_Nodes,Num_Nodes);


      // Mat_GP.Id = atoi(Parse_Mat_Prop[1]);
      // if(Mat_GP.Id >= GP_Mesh.NumberMaterials)
      // {
      //   sprintf(Error_message,"%s %i !!! \n","Id should go from 0 to",GP_Mesh.NumberMaterials-1);
      // standard_error(Error_message);
      // }
      // for(int i = 0 ; i<Num_Nodes ; i++)
      // {
      //   for(int j = 0 ; j<GPxElement ; j++)
      //   {
      //   GP_Mesh.MatIdx[Array_Nodes[i]*GPxElement+j] = Mat_GP.Id;
      //   }
      // }


    }

  }

  /*
   Close  file 
  */
  fclose(Sim_dat);

  return MatIdx
}

/***************************************************************************/

static void initialise_2D_particles(Mesh MPM_GID_Mesh,GaussPoint MPM_Mesh, int GPxElement)
  /*
     Loop in the GID mesh to create particles from an element 
  */
{


  Matrix Poligon_Coordinates;
  double Area_Element, A_p, th_p, m_p, rho_p;
  int p;
  int MatIdx_p;

  for(int i = 0 ; i<MPM_GID_Mesh.NumElemMesh ; i++)
  {

    /* Get the coordinates of the element vertex */ 
    Poligon_Coordinates = get_nodes_coordinates__MeshTools__(MPM_GID_Mesh.Connectivity[i], MPM_GID_Mesh.Coordinates);

    /* Get the area (Poligon_Centroid.n) and the position of the centroid (Poligon_Centroid.nV) */
    Area_Element = area__MatrixLib__(Poligon_Coordinates);

    for(int j = 0 ; j<GPxElement ; j++)
    {

      /* 
        Get the index of the material point 
      */
      p = i*GPxElement+j;

      /*
        Get the index of the material 
      */
      MatIdx_p = MPM_Mesh.MatIdx[p];
    /* Get material properties */
      A_p = Area_Element/GPxElement;
      rho_p = MPM_Mesh.Mat[MatIdx_p].rho;
      th_p = MPM_Mesh.Mat[MatIdx_p].thickness;
      m_p = th_p*A_p*rho_p;
    /* Set the initial volume */
      MPM_Mesh.Phi.Vol_0.nV[p] = A_p;
    /* Set the initial density */
      MPM_Mesh.Phi.rho.nV[p] = rho_p; 
    /* Assign the mass parameter */
      MPM_Mesh.Phi.mass.nV[p] = m_p;
    /* Local coordinates of the element */
      MPM_Mesh.I0[p] = -999;
      MPM_Mesh.NumberNodes[p] = 4;

      if(strcmp(MPM_Mesh.Mat[MatIdx_p].Type,"Von-Mises") == 0)
      {
          MPM_Mesh.Phi.cohesion.nV[p] = MPM_Mesh.Mat[MatIdx_p].yield_stress_0;
      }
      if(strcmp(MPM_Mesh.Mat[MatIdx_p].Type,"Drucker-Prager-Plane-Strain") == 0)
      {
          MPM_Mesh.Phi.cohesion.nV[p] = MPM_Mesh.Mat[MatIdx_p].cohesion_reference;
      }
      if(strcmp(MPM_Mesh.Mat[MatIdx_p].Type,"Drucker-Prager-Outer-Cone") == 0)
      {
          MPM_Mesh.Phi.cohesion.nV[p] = MPM_Mesh.Mat[MatIdx_p].cohesion_reference;
      }

    }

      /* Free data */
    free__MatrixLib__(Poligon_Coordinates);

  }
}

/***************************************************************************/

static bool Check_File(char * Path_File)
{
  struct stat info;
  stat(Path_File,&info);
  bool status_check;

  if(S_ISREG(info.st_mode))
  {
    printf("\t \t -> %s : %s \n","File",Path_File);
    status_check = true;
  }
  else
  {
    sprintf(Error_message,"%s %s %s","File",Path_File,"does not exists");
    standard_error();
  } 

  return status_check;
}

/***************************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Generate-Soil-Water-Coupling-Analysis()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/

static FILE * Open_and_Check_simulation_file(char * Name_File)
{
  FILE * Simulation_file = fopen(Name_File,"r");  
  
  if (Simulation_file==NULL)
  {
    sprintf(Error_message,"%s %s","Incorrect lecture of",Name_File);
    standard_error(); 
  }  

  return Simulation_file;
}

/***************************************************************************/
