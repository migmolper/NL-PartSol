#include "nl-partsol.h"
#include <sys/stat.h>

/*
  Call global variables 
*/
char * MPM_MeshFileName;

int Number_Soil_Water_Mixtures; // Number of Soil-Water Mixtures in the sample
Mixture * Soil_Water_Mixtures; // Structure with the properties of the sample


typedef struct
{
  bool Is_Soil_Water_Coupling;
  bool Is_ParticlesMesh;
  bool Is_GramsShapeFun;
  bool Is_GramsMaterials;
  bool Is_Material_Mixtures;
  bool Is_Particle_Initial;
  bool Is_Nodal_Initial;
  bool Is_GramsBodyForces;
  bool Is_GramsNeumannBC;

  int Counter_Materials;
  int Counter_Mixtures;
  int Counter_BodyForces;
  int Counter_GramsNeumannBC;

} Simulation_Key_Words;

typedef struct
{

  int   GPxElement;
  char Route_Mesh[MAXC];

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
static int * assign_mixture_to_particles(char *, int, int, int);
static void initialise_2D_particles();
static void Check_File(char *);
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


  /*
    Loop in the file to find key words, some of them are
    mandatory, other no.
  */
  Sim_Params = Check_Simulation_Key_Words(Name_File);

  /*
    Define materials
  */
  puts("*************************************************");
  printf(" \t %s \n","* Read materials properties :");
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

  /*
    Define material mixtures
  */
  puts("*************************************************");
  printf(" \t %s \n","* Read mixtures :");
  if(Sim_Params.Is_Material_Mixtures && (MPM_Mesh.NumberMaterials > 1))
  {
    Number_Soil_Water_Mixtures = Sim_Params.Counter_Materials;
    Soil_Water_Mixtures = Read_Soil_Water_Mixtures__InOutFun__(Name_File, Number_Soil_Water_Mixtures);
  }
  else
  {
    sprintf(Error_message,"%s","No material mixtures were defined or insuficient number of material for the mixture");
    standard_error(); 
  }
  
  /*
    Define particles mesh 
  */
  if(Sim_Params.Is_Soil_Water_Coupling)
  {
    
    /*
      Read particle mesh preliminar information
    */
    puts("*************************************************");
    printf(" \t %s \n","* Read mesh properties for particles :");
    Msh_Parms = Read_Mesh_Parameters(Name_File);

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
    MPM_Mesh.Phi = allocate_upw_vars__Fields__(NumParticles);

    /*
      Assign mixture for each material point
    */
    puts("*************************************************");
    printf(" \t %s \n","* Start mixture assignement to particles ...");
    MPM_Mesh.MixtIdx = assign_mixture_to_particles(Name_File,Number_Soil_Water_Mixtures,NumParticles,Msh_Parms.GPxElement);
    printf(" \t %s \n","Success !!");

    /*
      Initialise particle 
    */    
    puts("*************************************************");
    printf(" \t %s \n","* Start particles initialisation ...");
    initial_position__Particles__(MPM_Mesh.Phi.x_GC,MPM_GID_Mesh,Msh_Parms.GPxElement); 
    if(Ndim == 2)
    {
      initialise_2D_particles(MPM_GID_Mesh,MPM_Mesh,Msh_Parms.GPxElement);
    }
    printf(" \t %s \n","Success !!");
 
    /*
      Initialise shape functions 
    */
    puts("*************************************************");
    if(strcmp(ShapeFunctionGP,"MPMQ4") == 0)
    {
      printf("\t * %s \n","Start MPMQ4 shape functions initialisation ...");
      initialize__Q4__(MPM_Mesh, FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
    {
      printf("\t * %s \n","Start uGIMP shape functions initialisation ...");      
      initialize__GIMP__(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"LME") == 0)
    {
      printf("\t * %s \n","Start LME shape functions initialisation ...");
      initialize__LME__(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"aLME") == 0)
    {
      printf("\t * %s \n","Start aLME shape functions initialisation ...");
      initialize__aLME__(MPM_Mesh,FEM_Mesh);
    }
    printf("\t %s \n","Success !!");

    /*
      Read initial values 
    */    
    if(Sim_Params.Is_Particle_Initial)
    {
      Initial_condition_particles__InOutFun__(Name_File,MPM_Mesh,Msh_Parms.GPxElement);
    }
    else if(Sim_Params.Is_Nodal_Initial)
    {
      Initial_condition_nodes__InOutFun__(Name_File,MPM_Mesh,FEM_Mesh);
    }
    else{
      puts("*************************************************");
      printf("\t * %s \n","No initial conditions defined");
    }

    /*
      Read external forces 
    */    
    if(Sim_Params.Is_GramsNeumannBC)
    {
      puts("*************************************************");
      printf("\t * %s \n","Read Newmann boundary conditions :");
      MPM_Mesh.Neumann_Contours = Read_upw_Neumann_Boundary_Conditions__InOutFun__(Name_File,Sim_Params.Counter_GramsNeumannBC,Msh_Parms.GPxElement);
    }
    else
    {
      MPM_Mesh.Neumann_Contours.NumBounds = 0;
      puts("*************************************************");
      printf(" \t %s \n","* No Neumann boundary conditions defined");
    }

    /*
      Read body forces 
    */    
    MPM_Mesh.b = alloc__TensorLib__(1); // Vector with the current value of the distance accelerations
    if(Sim_Params.Is_GramsBodyForces)
    {
      MPM_Mesh.NumberBodyForces = Sim_Params.Counter_BodyForces;
      MPM_Mesh.B = GramsBodyForces(Name_File,Sim_Params.Counter_BodyForces,Msh_Parms.GPxElement); 
    }
    else
    {
      MPM_Mesh.NumberBodyForces = Sim_Params.Counter_BodyForces;
      puts("*************************************************");
      printf(" \t %s \n","* No body forces defined");
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
    sprintf(Error_message,"Sintax error in file %s : Soil-Water-Coupling-One-Layer statement is required for a -u-pw analisis",Name_File);
    standard_error(); 
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
  Sim_Key_Wrds.Is_Material_Mixtures = false;
  Sim_Key_Wrds.Is_Particle_Initial = false;
  Sim_Key_Wrds.Is_Nodal_Initial = false;
  Sim_Key_Wrds.Is_GramsBodyForces = false;
  Sim_Key_Wrds.Is_GramsNeumannBC = false;

  Sim_Key_Wrds.Counter_Materials = 0;
  Sim_Key_Wrds.Counter_Mixtures = 0;
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
    else if ((Num_words > 0) && (strcmp(Words[0],"Soil-Water-Coupling-One-Layer") == 0))
    {
      Sim_Key_Wrds.Is_Soil_Water_Coupling = true;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"GramsShapeFun") == 0))
    {
      Sim_Key_Wrds.Is_GramsShapeFun = true;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"Define-Material") == 0))
    {
      Sim_Key_Wrds.Is_GramsMaterials = true;
      Sim_Key_Wrds.Counter_Materials++;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"Define-Mixture") == 0))
    {
      Sim_Key_Wrds.Is_Material_Mixtures = true;
      Sim_Key_Wrds.Counter_Mixtures++;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"GramsInitials") == 0))
    {
      Sim_Key_Wrds.Is_Particle_Initial = true;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"Initial-nodal-values") == 0))
    {
      Sim_Key_Wrds.Is_Nodal_Initial = true;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"GramsBodyForces") == 0))
    {
      Sim_Key_Wrds.Is_GramsBodyForces = true;
      Sim_Key_Wrds.Counter_BodyForces++;
    }
    else if ((Num_words > 0) && (strcmp(Words[0],"Define-Neumann-Boundary") == 0))
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
  char Route_Mesh[MAXC] = {0};
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

    if ((Num_words >= 3) && (strcmp(Words[0],"Soil-Water-Coupling-One-Layer") == 0))
    {

     Num_parameters = parse(File_Parameter,Words[1],delimiters_3);
     if ((Num_parameters == 2) && (strcmp(File_Parameter[0],"File") == 0))
     {
        MPM_MeshFileName = File_Parameter[1];
        generate_route(Route_Mesh,Name_File);
        strcat(Route_Mesh,MPM_MeshFileName);
        strcpy(Msh_Params.Route_Mesh,Route_Mesh);
        Check_File(Msh_Params.Route_Mesh);
        printf("\t -> %s : %s \n","File",Route_Mesh);
     }
     else
     {
        sprintf(Error_message,"Sintax error in line %i : %s",Num_line,"Soil-Water-Coupling-One-Layer (File=Mesh.msh, *)");
        standard_error(); 
     }

      Num_parameters = parse(GPxElement_Parameter,Words[2],delimiters_3);
     if ((Num_parameters == 2) && (strcmp(GPxElement_Parameter[0],"GPxElement") == 0))
     {
        Msh_Params.GPxElement = atoi(GPxElement_Parameter[1]);
        printf("\t -> %s : %i \n","Particles per element",Msh_Params.GPxElement);
     }
     else
     {
        sprintf(Error_message,"Sintax error in line %i : %s",Num_line,"Soil-Water-Coupling-One-Layer (*, GPxElement=int)");
        standard_error(); 
     }

    }

    

  }

  /*
   Close  file 
  */
  fclose(Sim_dat);


  return Msh_Params;
}

/***************************************************************************/

static int * assign_mixture_to_particles(char * Name_File, int NumMixtures, int NumParticles, int GPxElement)
{
  ChainPtr Chain_Nodes = NULL;
  int * Array_Nodes;
  int * MixtIdx = (int *)calloc(NumParticles,sizeof(int));
  char Line[MAXC] = {0}; 
  char FileNodesRoute[MAXC] = {0};
  char Route_Nodes[MAXC] = {0};
  char * Words[MAXW] = {NULL};
  char * File_Parameter[MAXW] = {NULL};
  char * MixtIdx_Parameter[MAXW] = {NULL};
  int Num_words = 0;
  int Num_parameters = 0;
  int Num_line = 0;
  int Num_Nodes_File = 0;
  int Mixture_Index;
 
  /*
    Open and check file
  */
  FILE * Sim_dat = Open_and_Check_simulation_file(Name_File);

  /*
    Generate route with the current possition of the command file
  */
  generate_route(Route_Nodes,Name_File);


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

    if ((Num_words >= 3) && (strcmp(Words[0],"Assign-mixture-to-particles") == 0))
    {

      /* 
        Read index of the mixture
      */
      Num_words = parse(MixtIdx_Parameter,Words[1],delimiters_3);
      if((Num_words == 2) && (strcmp(MixtIdx_Parameter[0],"MixtIdx") == 0))
      {

        Mixture_Index = atoi(MixtIdx_Parameter[1]);

        if(Mixture_Index >= NumMixtures)
        {
          sprintf(Error_message,"Sintax error in line %i : %s %i",Num_line,"MixtIdx should go from 0 to",NumMixtures-1);
          standard_error(); 
        }
      }
      else
      {
        sprintf(Error_message,"Sintax error in line %i : %s",Num_line,"Assign-mixture-to-particles (MixtIdx=Int, *)");
        standard_error(); 
      }

      /*
        Read file with the nodes 
      */
      Num_words = parse(File_Parameter,Words[2],delimiters_3);
      if((Num_words == 2) && (strcmp(File_Parameter[0],"Particles") == 0))
      {
        /*
          Generate array with the nodes
        */
        sprintf(FileNodesRoute,"%s%s",Route_Nodes,File_Parameter[1]);
        Check_File(FileNodesRoute);
        Chain_Nodes = File2Chain(FileNodesRoute);
        Num_Nodes_File = lenght__SetLib__(Chain_Nodes);
        Array_Nodes = set_to_memory__SetLib__(Chain_Nodes,Num_Nodes_File);
        free__SetLib__(&Chain_Nodes);

        /*
          Fill the MixtIdx array with the index
        */      
        for(int i = 0 ; i<Num_Nodes_File ; i++)
        {
          for(int j = 0 ; j<GPxElement ; j++)
          {
            MixtIdx[Array_Nodes[i]*GPxElement+j] = Mixture_Index;
          }
        }

        /*
          Information message
        */
        printf("\t -> Mixture %i has been assigned to %i particles \n",Mixture_Index,GPxElement*Num_Nodes_File);

        /*
          Free memory
        */
        free(Array_Nodes);

      }
      else
      {
        sprintf(Error_message,"Sintax error in line %i : %s",Num_line,"Assign-mixture-to-particles (*, Particles=List-Particles.txt)");
        standard_error(); 
      }

    }

  }

  /*
   Close  file 
  */
  fclose(Sim_dat);

  return MixtIdx;
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
  int Mixture_idx;
  int Material_Soil_idx;
  int Material_Water_idx;

  Material MatProp_Soil_p; /* Variable with the material properties of the solid phase */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
  double rho_0; /* Initial density of the mixture */
  double rho_s_0; /* Initial density of the fluid */
  double rho_f_0; /* Initial density of the fluid */
  double phi_s_0; /* Initial volume fraction (solid) */
  double phi_f_0; /* Initial volume fraction (fluid) */

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
        Get the index of the mixture 
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
      Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
      MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
      MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];

      /*
        Get the intrinsic density of each phase and
        assign it to the current intrinsic density of
        each phase
      */
      rho_s_0 = MatProp_Soil_p.rho;
      rho_f_0 = MatProp_Water_p.rho;
      MPM_Mesh.Phi.rho_s.nV[p] = rho_s_0; 
      MPM_Mesh.Phi.rho_f.nV[p] = phi_f_0;

      /*
        Get initial volume fraction of each phase and
        assign it to the current volume fraction of
        each phase
      */
      phi_s_0 = Soil_Water_Mixtures[Mixture_idx].phi_s_0;
      phi_f_0 = Soil_Water_Mixtures[Mixture_idx].phi_f_0;
      MPM_Mesh.Phi.phi_s.nV[p] = phi_s_0;
      MPM_Mesh.Phi.phi_f.nV[p] = phi_f_0;

      /*
        Compute material properties 
      */
      A_p   = Area_Element/GPxElement; // Material point area
      rho_p = phi_s_0*rho_s_0 + phi_f_0*rho_f_0; // Mixture density
      th_p  = MPM_Mesh.Mat[Mixture_idx].thickness; // Material point thickness
      m_p   = th_p*A_p*rho_p; // Material point mass
      
      /*
        Set the initial total volume volume fractions 
      */
      MPM_Mesh.Phi.Vol_0.nV[p] = A_p;

      /*
        Set the relative density of the mixture
      */
      MPM_Mesh.Phi.rho.nV[p] = rho_p;
    
      /*
        Assign the mass parameter
      */
      MPM_Mesh.Phi.mass.nV[p] = m_p;
      
      /*
        Local coordinates of the element
      */
      MPM_Mesh.I0[p] = -999;
      MPM_Mesh.NumberNodes[p] = 4;

      /*
        Initialise some plastic variables
      */
      if(strcmp(MPM_Mesh.Mat[Mixture_idx].Type,"Von-Mises") == 0)
      {
          MPM_Mesh.Phi.cohesion.nV[p] = MPM_Mesh.Mat[Mixture_idx].yield_stress_0;
      }
      if(strcmp(MPM_Mesh.Mat[Mixture_idx].Type,"Drucker-Prager-Plane-Strain") == 0)
      {
          MPM_Mesh.Phi.cohesion.nV[p] = MPM_Mesh.Mat[Mixture_idx].cohesion_reference;
      }
      if(strcmp(MPM_Mesh.Mat[Mixture_idx].Type,"Drucker-Prager-Outer-Cone") == 0)
      {
          MPM_Mesh.Phi.cohesion.nV[p] = MPM_Mesh.Mat[Mixture_idx].cohesion_reference;
      }

    }

      /* Free data */
    free__MatrixLib__(Poligon_Coordinates);

  }
}

/***************************************************************************/

static void Check_File(char * Path_File)
{
  struct stat info;
  stat(Path_File,&info);

  if(S_ISREG(info.st_mode) == 0)
  {
    sprintf(Error_message,"%s %s %s","File",Path_File,"does not exists");
    standard_error();
  } 
}

/***************************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s !!! \n",Error_message);
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
