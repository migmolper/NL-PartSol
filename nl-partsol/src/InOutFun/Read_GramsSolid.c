#include "nl-partsol.h"


/*
  Call global variables 
*/
double Thickness_Plain_Stress;
char * MPM_MeshFileName;
char   wrapper_LME[MAXC];


/*
  Auxiliar functions 
*/
static void initialise_particles(Mesh,Particle,int);

/*********************************************************************/

Particle GramsSolid(char * Name_File, Mesh FEM_Mesh, int * STATUS)
/*
 */
{

  int Ndim = NumberDimensions;
  int NumParticles;


  // Error signals
  int INFO_Hidrostatic = 0;

  /* Simulation file */
  FILE * Sim_dat;
  
  /* Parser num chars */
  int Num_words_parse;

  /* Special variables GramsSolid2D */
  char Line_GramsSolid2D[MAXC] = {0}; 
  char * Parse_GramsSolid2D[MAXW] = {NULL};
  char * Parse_Mesh_id[MAXW] = {NULL};
  char * Parse_Mesh_Properties[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char Route_Mesh[MAXC] = {0};
  
  /* Material point mesh (Gauss-Points) */
  Mesh MPM_GID_Mesh;
  Particle MPM_Mesh;

  /* Number of particles per GID element */
  int GPxElement = 1;

  /* Set to false check variables */
  bool Is_GramsSolid = false;
  bool Is_ParticlesMesh = false;
  bool Is_GPxElement = false;
  bool Is_GramsShapeFun = false;
  bool Is_GramsMaterials = false;
  bool Is_Particle_Initial = false;
  bool Is_Nodal_Initial = false;
  bool Is_Hidrostatic_Initial = false;
  bool Is_GramsBodyForces = false;
  bool Is_GramsNeumannBC = false;

  /* Initialize counters */
  int Counter_Materials = 0;
  int Counter_BodyForces = 0;
  int Counter_GramsNeumannBC = 0;
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL)
  {
    fprintf(stderr,"%s : \n\t %s %s","Error in GramsInitials()","Incorrect lecture of",Name_File);
    exit(EXIT_FAILURE);
  }
  
  /**************************************************/
  /********* Read and check GramsSolid2D ************/
  /**************************************************/
  while(fgets(Line_GramsSolid2D,sizeof(Line_GramsSolid2D),Sim_dat) != NULL)
  {

    /* Read the line with the space as separators */
    Num_words_parse = parse(Parse_GramsSolid2D,Line_GramsSolid2D," \n\t");
    if (Num_words_parse < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid()",
	      "Parser failed");
      exit(EXIT_FAILURE);
    }

    /*
      Read file mesh and number of material points per element
    */
    if((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"GramsSolid") == 0))
    {
      Is_GramsSolid = true;
      Num_words_parse = parse(Parse_Mesh_id, Parse_GramsSolid2D[1],"(,)");

      /*
	       Propertie 1
      */
      Num_words_parse = parse(Parse_Mesh_Properties, Parse_Mesh_id[0],"=");
      if((Num_words_parse == 2) && (strcmp(Parse_Mesh_Properties[0],"File") == 0))
      {
        MPM_MeshFileName = Parse_Mesh_Properties[1];
        generate_route(Route_Mesh,Name_File);
        strcat(Route_Mesh,MPM_MeshFileName);
        Is_ParticlesMesh = true;
      }
      
      /*
	      Propertie 2
      */
      Num_words_parse = parse(Parse_Mesh_Properties,Parse_Mesh_id[1],"=");
      if((Num_words_parse == 2) && (strcmp(Parse_Mesh_Properties[0],"GPxElement") == 0))
      {
        GPxElement = atoi(Parse_Mesh_Properties[1]);
        Is_GPxElement = true;
      }

      Thickness_Plain_Stress = 1.0;
      
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"GramsShapeFun") == 0))
    {
      Is_GramsShapeFun = true;
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"GramsMaterials") == 0))
    {
      Is_GramsMaterials = true;
      Counter_Materials++;
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"GramsInitials") == 0))
    {
      Is_Particle_Initial = true;
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"Initial-nodal-values") == 0))
    {
      Is_Nodal_Initial = true;
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"Hydrostatic-condition") == 0))
    {
      Is_Hidrostatic_Initial = true;
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"GramsBodyForces") == 0))
    {
      Is_GramsBodyForces = true;
      Counter_BodyForces++;
    }
    if ((Num_words_parse > 0) && (strcmp(Parse_GramsSolid2D[0],"Define-Neumann-Boundary") == 0))
    {
      Is_GramsNeumannBC = true;
      Counter_GramsNeumannBC++;
    }    
  }

  /**************************************************/
  /***************** Define MPM Mesh ****************/
  /**************************************************/
  /* Define GP Mesh */
  if(Is_GramsSolid && Is_GPxElement)
  {
    
    /* Read GP mesh */
    MPM_GID_Mesh = ReadGidMesh__MeshTools__(Route_Mesh);

    /* Define the number of particles */
    NumParticles = GPxElement*MPM_GID_Mesh.NumElemMesh;
    MPM_Mesh.NumGP = NumParticles;

    /* Closest node to the particle */
    MPM_Mesh.I0 = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* Element of the particle */
    MPM_Mesh.Element_p = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* Number of tributary nodes for each particle */
    MPM_Mesh.NumberNodes = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* Tributary nodes for each particle */
    MPM_Mesh.ListNodes = (ChainPtr *)malloc(NumParticles*sizeof(ChainPtr));
    if(MPM_Mesh.ListNodes == NULL)
    {
      printf("%s : %s \n","Define_GP_Mesh","Memory error for ListNodes");
      exit(EXIT_FAILURE);
    }
    for(int i = 0 ; i<NumParticles ; i++)
    {
      MPM_Mesh.ListNodes[i] = NULL;  
    }

    /* List of particles close to each particle */
    MPM_Mesh.Beps = (ChainPtr *)malloc(NumParticles*sizeof(ChainPtr));
    if(MPM_Mesh.Beps == NULL)
    {
      printf("%s : %s \n","Define_GP_Mesh","Memory error for Beps");
      exit(EXIT_FAILURE);
    }
    for(int i = 0 ; i<NumParticles ; i++)
    {
      MPM_Mesh.Beps[i] = NULL;  
    }

    /**************************************************/
    /********* Read Shape functions parameters ********/
    /**************************************************/
    if(Is_GramsShapeFun)
    {
      GramsShapeFun(Name_File);
      /* Lenght of the Voxel (Only GIMP) */
      if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
      {
        MPM_Mesh.lp = allocZ__MatrixLib__(NumParticles,Ndim);
        strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
      }
      /* Lagrange Multipliers / Beta (Only LME ) */
      if(strcmp(ShapeFunctionGP,"LME") == 0)
      {
        MPM_Mesh.lambda = allocZ__MatrixLib__(NumParticles,Ndim);
        strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
        MPM_Mesh.Beta = allocZ__MatrixLib__(NumParticles,1);
        strcpy(MPM_Mesh.Beta.Info,"Beta parameter");

        if(strcmp(wrapper_LME,"Newton-Raphson") == 0)
        {
          MPM_Mesh.update_lambda = update_lambda_Newton_Rapson__LME__;
        }
        else if(strcmp(wrapper_LME,"Nelder-Mead") == 0)
        {
          MPM_Mesh.update_lambda = update_lambda_Nelder_Mead__LME__;
        }
        else
        {
          fprintf(stderr,"%s : %s \n",
            "Error in GramsSolid()","Unrecognaised wrapper");
          exit(EXIT_FAILURE);      
        }
      }
    }
    else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid()","GramsShapeFun no defined");
      exit(EXIT_FAILURE);
    }


    /**************************************************/    
    /****** Allocate vectorial/tensorial fields *******/
    /**************************************************/
    MPM_Mesh.Phi = allocate_Fields(NumParticles);

    /**************************************************/
    /************ Read material properties ************/
    /**************************************************/
    puts("*************************************************");
    printf(" \t %s \n","* Read materials properties");
    if(Is_GramsMaterials)
    {
      MPM_Mesh.NumberMaterials = Counter_Materials;
      MPM_Mesh.MatIdx = (int *)malloc(NumParticles*sizeof(int));
      MPM_Mesh.Mat = GramsMaterials(Name_File,MPM_Mesh,GPxElement);
    }
    else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid()","GramsMaterials no defined");
      exit(EXIT_FAILURE);
    }

    /**************************************************/    
    /************** Initialise particle ***************/
    /**************************************************/    
    initial_position__Particles__(MPM_Mesh.Phi.x_GC,MPM_GID_Mesh,GPxElement);
    initialise_particles(MPM_GID_Mesh,MPM_Mesh,GPxElement);
     
    /**************************************************/
    /******** INITIALIZE SHAPE FUNCTIONS **************/
    /**************************************************/
    puts("*************************************************");
    if(strcmp(ShapeFunctionGP,"FEM") == 0)
    {
      if(strcmp(FEM_Mesh.TypeElem,"Triangle") == 0)
      {
        printf("\t * %s \n","Start FEM-T3 shape functions initialisation ...");
        initialize__T3__(MPM_Mesh, FEM_Mesh);
      }
      else if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0)
      {
        printf("\t * %s \n","Start FEM-Q4 shape functions initialisation ...");
        initialize__Q4__(MPM_Mesh, FEM_Mesh);
      }
      else if(strcmp(FEM_Mesh.TypeElem,"Tetrahedra") == 0)
      {
        printf("\t * %s \n","Start FEM-T4 shape functions initialisation ...");
        initialize__T4__(MPM_Mesh, FEM_Mesh);
      }
      else if(strcmp(FEM_Mesh.TypeElem,"Hexahedra") == 0)
      {
        printf("\t * %s \n","Start FEM-H8 shape functions initialisation ...");
        initialize__H8__(MPM_Mesh, FEM_Mesh);
      }
    }
    else if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
    {
      printf("\t * %s \n","Initialize uGIMP shape functions ...");      
      initialize__GIMP__(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"LME") == 0)
    {
      printf("\t * %s \n","Initialize LME shape functions ...");
      initialize__LME__(MPM_Mesh,FEM_Mesh);
    }
    else
    {
      fprintf(stderr,"%s : %s \n","Error in GramsSolid()","Undefined kind of shape function");
      exit(EXIT_FAILURE);
    }
    printf("\t %s \n","DONE !!");

    /**************************************************/    
    /************** Read initial values ***************/
    /**************************************************/
 
    if(Is_Particle_Initial)
    {
      Initial_condition_particles__InOutFun__(Name_File,MPM_Mesh,GPxElement);
    }
    else if(Is_Nodal_Initial)
    {
      Initial_condition_nodes__InOutFun__(Name_File,MPM_Mesh,FEM_Mesh);
    }
    else if(Is_Hidrostatic_Initial)
    {
      Hidrostatic_condition_particles__InOutFun__(Name_File,MPM_Mesh,GPxElement,&INFO_Hidrostatic);
      
      if(INFO_Hidrostatic)
      {
        fprintf(stderr,"Error in : %s\n","GramsSolid");
        (*STATUS) = 1;
        return MPM_Mesh;
      }

    }
    else
    {
      puts("*************************************************");
      puts(" No initial conditions defined ");
    }
    /**************************************************/    
    /************** Read external forces **************/
    /**************************************************/    
    if(Is_GramsNeumannBC)
    {
      MPM_Mesh.Neumann_Contours = Read_u_Neumann_Boundary_Conditions__InOutFun__(Name_File,Counter_GramsNeumannBC,GPxElement);
    }
    else
    {
      MPM_Mesh.Neumann_Contours.NumBounds = 0;
      puts("*************************************************");
      printf(" \t %s \n","* No Neumann boundary conditions defined");
    }
    /**************************************************/    
    /*************** Read body forces *****************/
    /**************************************************/    
    if(Is_GramsBodyForces)
    {
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      MPM_Mesh.B = GramsBodyForces(Name_File,Counter_BodyForces,GPxElement); 
    }
    else
    {
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
       "* No body forces defined in",Name_File);
    }
   
    /**************************************************/    
    /************* Free the input data ****************/
    /**************************************************/
    for(int i = 0 ; i<MPM_GID_Mesh.NumElemMesh ; i++)
    {
      free__SetLib__(&MPM_GID_Mesh.Connectivity[i]); 
    }   
    free(MPM_GID_Mesh.Connectivity);
    free__MatrixLib__(MPM_GID_Mesh.Coordinates);
    free(MPM_GID_Mesh.Num_Particles_Node);

  } 
  else
  {
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid()",
	    "Mesh file name and number of particles are required");
    exit(EXIT_FAILURE);
  }
  
  return MPM_Mesh;
}


/***************************************************************************/

static void initialise_particles(Mesh MPM_GID_Mesh, Particle MPM_Mesh, int GPxElement)
  /*
     Loop in the GID mesh to create particles from an element 
  */
{


  Matrix Element_Coordinates;
  double Vol_Element, V_p, m_p, rho_p;
  int p;
  int MatIdx_p;

  for(int i = 0 ; i<MPM_GID_Mesh.NumElemMesh ; i++)
  {

    /* Get the coordinates of the element vertex */ 
    Element_Coordinates = get_nodes_coordinates__MeshTools__(MPM_GID_Mesh.Connectivity[i], MPM_GID_Mesh.Coordinates);
    Vol_Element = MPM_GID_Mesh.volume_Element(Element_Coordinates);
    free__MatrixLib__(Element_Coordinates);

    if(Vol_Element <= 0.0)
    {
      fprintf(stderr,"%s : %s \n",
      "Error in GramsSolid()",
      "Element with negative volume");
      exit(EXIT_FAILURE);
    }

    for(int j = 0 ; j<GPxElement ; j++)
    {

    /* Get the index of the material point */
      p = i*GPxElement+j;

    /* Get the index of the material */
      MatIdx_p = MPM_Mesh.MatIdx[p];

    /* Get material properties */
      V_p = Vol_Element/GPxElement;
      rho_p = MPM_Mesh.Mat[MatIdx_p].rho;
      m_p = V_p*rho_p;

    /* Set the initial volume */
      MPM_Mesh.Phi.Vol_0.nV[p] = V_p;

    /* Set the initial density */
      MPM_Mesh.Phi.rho.nV[p] = rho_p; 

    /* Assign the mass parameter */
      MPM_Mesh.Phi.mass.nV[p] = m_p;

    }    

  }
}

/***************************************************************************/

