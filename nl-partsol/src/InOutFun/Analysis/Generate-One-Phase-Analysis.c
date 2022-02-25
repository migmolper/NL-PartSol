#include <sys/stat.h>
#include <string.h>
#include "nl-partsol.h"

/*
  Call global variables
*/
char *MPM_MeshFileName;

double Thickness_Plain_Stress; // For 2D cases

typedef struct {
  bool Is_One_Phase_Analysis;
  bool Is_ParticlesMesh;
  bool Is_GramsShapeFun;
  bool Is_Materials;
  bool Is_Particle_Initial;
  bool Is_Hydrostatic_conditions;
  bool Is_Nodal_Initial;
  bool Is_GramsBodyForces;
  bool Is_GramsNeumannBC;

  int Counter_Materials;
  int Counter_BodyForces;
  int Counter_GramsNeumannBC;

} Simulation_Key_Words;

typedef struct {

  int GPxElement;
  char Route_Mesh[MAXC];

} Mesh_Parameters;

/*
  Auxiliar functions and variables
*/
#ifdef _WIN32
static char *delimiters_1 = " (,)\r\n\t";
static char *delimiters_2 = " =\t\r\n";
#else
static char *delimiters_1 = " (,)\n\t";
static char *delimiters_2 = " =\t\n";
#endif
static char *delimiters_3 = "=";
static char *delimiters_4 = ";";

static char Error_message[MAXW];

static Simulation_Key_Words Check_Simulation_Key_Words(char *);
static Mesh_Parameters Read_Mesh_Parameters(char *);
static int *assign_material_to_particles(char *, int, int, int);
static void initialise_particles(Mesh, Particle, int);
static void Check_File(char *);
static void standard_error();
static FILE *Open_and_Check_simulation_file(char *);

/*********************************************************************/

Particle
Generate_One_Phase_Analysis__InOutFun__(char *Name_File, Mesh FEM_Mesh,
                                        Time_Int_Params Parameters_Solver) {
  int Ndim = NumberDimensions;
  int NumParticles;
  int status = 0;

  /* Parser num chars */
  int Num_words_parse;

  Mesh MPM_GID_Mesh; /* GID mesh to define the material point mesh */
  Particle MPM_Mesh; /* Material point mesh (Gauss-Points) */
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
  printf(" \t %s \n", "* Read materials properties :");
  if (Sim_Params.Is_Materials) {
    MPM_Mesh.NumberMaterials = Sim_Params.Counter_Materials;
    MPM_Mesh.Mat =
        Read_Materials__InOutFun__(Name_File, MPM_Mesh.NumberMaterials);
  } else {
    sprintf(Error_message, "%s", "No materials were defined");
    standard_error();
  }

  /*
    Define particles mesh
  */
  if (Sim_Params.Is_One_Phase_Analysis) {

    /*
      Read particle mesh preliminar information
    */
    printf(" \t %s \n", "* Read mesh properties for particles :");
    Msh_Parms = Read_Mesh_Parameters(Name_File);

    Thickness_Plain_Stress = 1.0;

    /*
      Read particles mesh
    */
    MPM_GID_Mesh = ReadGidMesh__MeshTools__(Msh_Parms.Route_Mesh);

    /*
      Define the number of particles
    */
    NumParticles = Msh_Parms.GPxElement * MPM_GID_Mesh.NumElemMesh;
    MPM_Mesh.NumGP = NumParticles;

    /*
      Closest node to the particle
    */
    MPM_Mesh.I0 = (int *)Allocate_ArrayZ(NumParticles, sizeof(int));

    /*
      Element of the particle
    */
    MPM_Mesh.Element_p = (int *)Allocate_ArrayZ(NumParticles, sizeof(int));

    /*
      Number of tributary nodes for each particle
    */
    MPM_Mesh.NumberNodes = (int *)Allocate_ArrayZ(NumParticles, sizeof(int));

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
    if (Sim_Params.Is_GramsShapeFun) {

      /*
        Read Shape functions parameters
      */
      GramsShapeFun(Name_File);

      if (strcmp(ShapeFunctionGP, "FEM") == 0) {

      } else if (strcmp(ShapeFunctionGP, "uGIMP") == 0) {
        MPM_Mesh.lp = allocZ__MatrixLib__(NumParticles, Ndim);
        strcpy(MPM_Mesh.lp.Info, "Voxel lenght GP");
      } else if (strcmp(ShapeFunctionGP, "LME") == 0) {
        MPM_Mesh.lambda = allocZ__MatrixLib__(NumParticles, Ndim);
        strcpy(MPM_Mesh.lambda.Info, "Lagrange Multiplier");
        MPM_Mesh.Beta = allocZ__MatrixLib__(NumParticles, 1);
        strcpy(MPM_Mesh.Beta.Info, "Beta");
      } else if (strcmp(ShapeFunctionGP, "aLME") == 0) {
        MPM_Mesh.lambda = allocZ__MatrixLib__(NumParticles, Ndim);
        strcpy(MPM_Mesh.lambda.Info, "Lagrange Multiplier");
        MPM_Mesh.Beta = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
        strcpy(MPM_Mesh.Beta.Info, "Beta parameter");
        MPM_Mesh.Cut_off_Ellipsoid =
            allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
        strcpy(MPM_Mesh.Cut_off_Ellipsoid.Info, "Cut-off Ellipsoid");
      } else {
        fprintf(stderr, "%s : %s \n", "Error in GramsShapeFun()",
                "Undefined kind of shape function");
        exit(EXIT_FAILURE);
      }
    } else {
      sprintf(Error_message, "%s", "GramsShapeFun no defined");
      standard_error();
    }

    /*
      Allocate variables for the shape function
    */

    /*
      Allocate vectorial/tensorial fields
    */
    if (strcmp(Formulation, "-u") == 0) {
      MPM_Mesh.Phi = allocate_U_vars__Fields__(NumParticles);
    } else if (strcmp(Formulation, "-up") == 0) {
      MPM_Mesh.Phi = allocate_Up_vars__Fields__(NumParticles);
    }

    /*
      Assign material for each material point
    */
    printf(" \t %s \n", "* Start material assignement to particles ...");
    MPM_Mesh.MatIdx =
        assign_material_to_particles(Name_File, MPM_Mesh.NumberMaterials,
                                     NumParticles, Msh_Parms.GPxElement);
    printf(" \t %s \n", "Success !!");

    /*
      Initialise particle
    */
    printf(" \t %s \n", "* Start particles initialisation ...");
    initial_position__Particles__(MPM_Mesh.Phi.x_GC, MPM_GID_Mesh,
                                  Msh_Parms.GPxElement);
    if (Ndim == 2) {
      initialise_particles(MPM_GID_Mesh, MPM_Mesh, Msh_Parms.GPxElement);
    }
    printf(" \t %s \n", "Success !!");

    /*
      Read initial values
    */
    if (Sim_Params.Is_Particle_Initial) {
      Initial_condition_particles__InOutFun__(Name_File, MPM_Mesh,
                                              Msh_Parms.GPxElement);
    } else if (Sim_Params.Is_Nodal_Initial) {
      Initial_condition_nodes__InOutFun__(Name_File, MPM_Mesh, FEM_Mesh);
    } else {
      printf("\t * %s \n", "No initial conditions defined");
    }

    /*
      Read hydrostatic conditions
    */
    if (Sim_Params.Is_Hydrostatic_conditions) {
      status = Hidrostatic_condition_particles__InOutFun__(
          Name_File, MPM_Mesh, Msh_Parms.GPxElement);
      if (status) {
        exit(0);
      }
    }

    /*
      Read external forces
    */
    if (Sim_Params.Is_GramsNeumannBC) {
      printf("\t * %s \n", "Read Newmann boundary conditions :");
      if (strcmp(Formulation, "-u") == 0) {
        MPM_Mesh.Neumann_Contours =
            Read_u_Neumann_Boundary_Conditions__InOutFun__(
                Name_File, Sim_Params.Counter_GramsNeumannBC,
                Msh_Parms.GPxElement, Parameters_Solver.NumTimeStep);
        Check_u_Neumann_Boundary_Conditions__InOutFun__(
            MPM_Mesh.Neumann_Contours, NumParticles);
      } else if (strcmp(Formulation, "-up") == 0) {
        MPM_Mesh.Neumann_Contours =
            Read_u_Neumann_Boundary_Conditions__InOutFun__(
                Name_File, Sim_Params.Counter_GramsNeumannBC,
                Msh_Parms.GPxElement, Parameters_Solver.NumTimeStep);
        Check_u_Neumann_Boundary_Conditions__InOutFun__(
            MPM_Mesh.Neumann_Contours, NumParticles);
      }
    } else {
      MPM_Mesh.Neumann_Contours.NumBounds = 0;
      printf(" \t %s \n", "* No Neumann boundary conditions defined");
    }

    /*
      Read body forces
    */
    MPM_Mesh.b = alloc__TensorLib__(
        1); // Vector with the current value of the distance accelerations
    if (Sim_Params.Is_GramsBodyForces) {
      MPM_Mesh.NumberBodyForces = Sim_Params.Counter_BodyForces;
      MPM_Mesh.B =
          GramsBodyForces(Name_File, Sim_Params.Counter_BodyForces,
                          Msh_Parms.GPxElement, Parameters_Solver.NumTimeStep);
    } else {
      MPM_Mesh.NumberBodyForces = Sim_Params.Counter_BodyForces;
      printf(" \t %s \n", "* No body forces defined");
    }

    /*
      Free the input data
    */
    for (int i = 0; i < MPM_GID_Mesh.NumElemMesh; i++) {
      free__SetLib__(&MPM_GID_Mesh.Connectivity[i]);
    }
    free(MPM_GID_Mesh.Connectivity);
    free__MatrixLib__(MPM_GID_Mesh.Coordinates);
    free(MPM_GID_Mesh.Num_Particles_Node);

  } else {
    sprintf(Error_message,
            "Sintax error in file %s : One-Phase-Analysis statement is "
            "required for a -u analisis",
            Name_File);
    standard_error();
  }

  return MPM_Mesh;
}

/***************************************************************************/

static Simulation_Key_Words Check_Simulation_Key_Words(char *Name_File) {

  Simulation_Key_Words Sim_Key_Wrds;
  char Line[MAXC] = {0};
  char *Words[MAXW] = {NULL};
  int Num_words = 0;
  int Num_line = 0;

  /*
    Default values
  */
  Sim_Key_Wrds.Is_One_Phase_Analysis = false;
  Sim_Key_Wrds.Is_ParticlesMesh = false;
  Sim_Key_Wrds.Is_GramsShapeFun = false;
  Sim_Key_Wrds.Is_Materials = false;
  Sim_Key_Wrds.Is_Particle_Initial = false;
  Sim_Key_Wrds.Is_Hydrostatic_conditions = false;
  Sim_Key_Wrds.Is_Nodal_Initial = false;
  Sim_Key_Wrds.Is_GramsBodyForces = false;
  Sim_Key_Wrds.Is_GramsNeumannBC = false;

  Sim_Key_Wrds.Counter_Materials = 0;
  Sim_Key_Wrds.Counter_BodyForces = 0;
  Sim_Key_Wrds.Counter_GramsNeumannBC = 0;

  /*
    Open and check file
  */
  FILE *Sim_dat = Open_and_Check_simulation_file(Name_File);

  while (fgets(Line, sizeof(Line), Sim_dat) != NULL) {

    /* Read the line with the space as separators */
    Num_words = parse(Words, Line, delimiters_1);

    /*
      Update line number
    */
    Num_line++;

    if (Num_words < 0) {
      sprintf(Error_message, "%s : %i", "Parser failed in line", Num_line);
      standard_error();
    } else if ((Num_words > 0) &&
               (strcmp(Words[0], "One-Phase-Analysis") == 0)) {
      Sim_Key_Wrds.Is_One_Phase_Analysis = true;
    } else if ((Num_words > 0) && (strcmp(Words[0], "GramsShapeFun") == 0)) {
      Sim_Key_Wrds.Is_GramsShapeFun = true;
    } else if ((Num_words > 0) && (strcmp(Words[0], "Define-Material") == 0)) {
      Sim_Key_Wrds.Is_Materials = true;
      Sim_Key_Wrds.Counter_Materials++;
    } else if ((Num_words > 0) && (strcmp(Words[0], "GramsInitials") == 0)) {
      Sim_Key_Wrds.Is_Particle_Initial = true;
    } else if ((Num_words > 0) &&
               (strcmp(Words[0], "Hydrostatic-condition") == 0)) {
      Sim_Key_Wrds.Is_Hydrostatic_conditions = true;
    } else if ((Num_words > 0) &&
               (strcmp(Words[0], "Initial-nodal-values") == 0)) {
      Sim_Key_Wrds.Is_Nodal_Initial = true;
    } else if ((Num_words > 0) && (strcmp(Words[0], "GramsBodyForces") == 0)) {
      Sim_Key_Wrds.Is_GramsBodyForces = true;
      Sim_Key_Wrds.Counter_BodyForces++;
    } else if ((Num_words > 0) &&
               (strcmp(Words[0], "Define-Neumann-Boundary") == 0)) {
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

static Mesh_Parameters Read_Mesh_Parameters(char *Name_File) {

  Mesh_Parameters Msh_Params;
  char Line[MAXC] = {0};
  char Route_Mesh[MAXC] = {0};
  char *Words[MAXW] = {NULL};
  char *File_Parameter[MAXW] = {NULL};
  char *GPxElement_Parameter[MAXW] = {NULL};
  int Num_words = 0;
  int Num_parameters = 0;
  int Num_line = 0;

  /*
    Open and check file
  */
  FILE *Sim_dat = Open_and_Check_simulation_file(Name_File);

  while (fgets(Line, sizeof(Line), Sim_dat) != NULL) {

    /*
     Read the line with the space as separators
    */
    Num_words = parse(Words, Line, delimiters_1);

    /*
      Update line number
    */
    Num_line++;

    if ((Num_words >= 3) && (strcmp(Words[0], "One-Phase-Analysis") == 0)) {

      Num_parameters = parse(File_Parameter, Words[1], delimiters_3);
      if ((Num_parameters == 2) && (strcmp(File_Parameter[0], "File") == 0)) {
        MPM_MeshFileName = File_Parameter[1];
        generate_route(Route_Mesh, Name_File);
        strcat(Route_Mesh, MPM_MeshFileName);
        strcpy(Msh_Params.Route_Mesh, Route_Mesh);
        Check_File(Msh_Params.Route_Mesh);
        printf("\t -> %s : %s \n", "File", Route_Mesh);
      } else {
        sprintf(Error_message, "Sintax error in line %i : %s", Num_line,
                "One-Phase-Analysis (File=Mesh.msh, *)");
        standard_error();
      }

      Num_parameters = parse(GPxElement_Parameter, Words[2], delimiters_3);
      if ((Num_parameters == 2) &&
          (strcmp(GPxElement_Parameter[0], "GPxElement") == 0)) {
        Msh_Params.GPxElement = atoi(GPxElement_Parameter[1]);
        printf("\t -> %s : %i \n", "Particles per element",
               Msh_Params.GPxElement);
      } else {
        sprintf(Error_message, "Sintax error in line %i : %s", Num_line,
                "One-Phase-Analysis (*, GPxElement=int)");
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

static int *assign_material_to_particles(char *Name_File, int NumMaterials,
                                         int NumParticles, int GPxElement) {
  ChainPtr Chain_Nodes = NULL;
  int *Array_Nodes;
  int *MatIdx = (int *)calloc(NumParticles, sizeof(int));
  char Line[MAXC] = {0};
  char FileNodesRoute[MAXC] = {0};
  char Route_Nodes[MAXC] = {0};
  char *Words[MAXW] = {NULL};
  char *File_Parameter[MAXW] = {NULL};
  char *MatIdx_Parameter[MAXW] = {NULL};
  int Num_words = 0;
  int Num_parameters = 0;
  int Num_line = 0;
  int Num_Nodes_File = 0;
  int Material_Index;

  /*
    Open and check file
  */
  FILE *Sim_dat = Open_and_Check_simulation_file(Name_File);

  /*
    Generate route with the current possition of the command file
  */
  generate_route(Route_Nodes, Name_File);

  while (fgets(Line, sizeof(Line), Sim_dat) != NULL) {

    /*
     Read the line with the space as separators
    */
    Num_words = parse(Words, Line, delimiters_1);

    /*
      Update line number
    */
    Num_line++;

    if ((Num_words >= 3) &&
        (strcmp(Words[0], "Assign-material-to-particles") == 0)) {

      /*
        Read index of the mixture
      */
      Num_words = parse(MatIdx_Parameter, Words[1], delimiters_3);
      if ((Num_words == 2) && (strcmp(MatIdx_Parameter[0], "MatIdx") == 0)) {

        Material_Index = atoi(MatIdx_Parameter[1]);

        if (Material_Index >= NumMaterials) {
          sprintf(Error_message, "Sintax error in line %i : %s %i", Num_line,
                  "MatIdx should go from 0 to", NumMaterials - 1);
          standard_error();
        }
      } else {
        sprintf(Error_message, "Sintax error in line %i : %s", Num_line,
                "Assign-material-to-particles (MatIdx=Int, *)");
        standard_error();
      }

      /*
        Read file with the nodes
      */
      Num_words = parse(File_Parameter, Words[2], delimiters_3);
      if ((Num_words == 2) && (strcmp(File_Parameter[0], "Particles") == 0)) {
        /*
          Generate array with the nodes
        */
        sprintf(FileNodesRoute, "%s%s", Route_Nodes, File_Parameter[1]);
        Check_File(FileNodesRoute);
        Chain_Nodes = File2Chain(FileNodesRoute);
        Num_Nodes_File = lenght__SetLib__(Chain_Nodes);
        Array_Nodes = set_to_memory__SetLib__(Chain_Nodes, Num_Nodes_File);
        free__SetLib__(&Chain_Nodes);

        /*
          Fill the MatIdx array with the index
        */
        for (int i = 0; i < Num_Nodes_File; i++) {
          for (int j = 0; j < GPxElement; j++) {
            MatIdx[Array_Nodes[i] * GPxElement + j] = Material_Index;
          }
        }

        /*
          Information message
        */
        printf("\t -> Material %i has been assigned to %i particles \n",
               Material_Index, GPxElement * Num_Nodes_File);

        /*
          Free memory
        */
        free(Array_Nodes);

      } else {
        sprintf(
            Error_message, "Sintax error in line %i : %s", Num_line,
            "Assign-material-to-particles (*, Particles=List-Particles.txt)");
        standard_error();
      }
    }
  }

  /*
   Close  file
  */
  fclose(Sim_dat);

  return MatIdx;
}

/***************************************************************************/

static void initialise_particles(Mesh MPM_GID_Mesh, Particle MPM_Mesh,
                                 int GPxElement)
/*
   Loop in the GID mesh to create particles from an element
*/
{

  Matrix Element_Coordinates;
  double Vol_Element, V_p, m_p, rho_p;
  int p;
  int MatIdx_p;

  for (int i = 0; i < MPM_GID_Mesh.NumElemMesh; i++) {

    /* Get the coordinates of the element vertex */
    Element_Coordinates = get_nodes_coordinates__MeshTools__(
        MPM_GID_Mesh.Connectivity[i], MPM_GID_Mesh.Coordinates);
    Vol_Element = MPM_GID_Mesh.volume_Element(Element_Coordinates);
    free__MatrixLib__(Element_Coordinates);

    if (Vol_Element <= 0.0) {
      fprintf(stderr, "%s : %s \n", "Error in One-Phase-Analysis()",
              "Element with negative volume");
      exit(EXIT_FAILURE);
    }

    for (int j = 0; j < GPxElement; j++) {

      /* Get the index of the material point */
      p = i * GPxElement + j;

      /* Get the index of the material */
      MatIdx_p = MPM_Mesh.MatIdx[p];

      /* Get material properties */
      V_p = Vol_Element / GPxElement;
      rho_p = MPM_Mesh.Mat[MatIdx_p].rho;
      m_p = V_p * rho_p;

      /* Set the initial volume */
      MPM_Mesh.Phi.Vol_0.nV[p] = V_p;

      /* Set the initial density */
      MPM_Mesh.Phi.rho.nV[p] = rho_p;

      /* Assign the mass parameter */
      MPM_Mesh.Phi.mass.nV[p] = m_p;

      /* Initialize frictional material */
      if (strcmp(MPM_Mesh.Mat[MatIdx_p].Type, "Granular") == 0) {
        Initialize_Frictional(&MPM_Mesh.Phi.Kappa_hardening.nV[p],
                              &MPM_Mesh.Phi.Equiv_Plast_Str.nV[p],
                              MPM_Mesh.Mat[MatIdx_p]);
      }
    }
  }
}

/***************************************************************************/

static void Check_File(char *Path_File) {
  struct stat info;
  stat(Path_File, &info);

  if (S_ISREG(info.st_mode) == 0) {
    sprintf(Error_message, "%s %s %s", "File", Path_File, "does not exists");
    standard_error();
  }
}

/***************************************************************************/

static void standard_error() {
  fprintf(stderr, "%s !!! \n", Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/

static FILE *Open_and_Check_simulation_file(char *Name_File) {
  FILE *Simulation_file = fopen(Name_File, "r");

  if (Simulation_file == NULL) {
    sprintf(Error_message, "%s %s", "Incorrect lecture of", Name_File);
    standard_error();
  }

  return Simulation_file;
}

/***************************************************************************/
