#include <sys/stat.h>
#include <string.h>
#include "nl-partsol.h"

/*
  Local structures and variables
*/
typedef struct {
  char Mesh_Generator[100];
  char Mesh_File[100];
  int Number_Boundaries;

} Nodes_Information;

#ifdef _WIN32
static char *delimiters = " =,()\r\n\t";
#else
static char *delimiters = " =,()\n\t";
#endif

static char Error_message[MAXW];

/*
  Auxiliar functions
*/
static Nodes_Information Read_Nodal_Set_Information(char *);
static void Check_Mesh_File(char *);
static void get_sourrounding_elements(Mesh);
static void fill_nodal_locality(Mesh, int);
static ChainPtr node_I_locality(int, Mesh);
static ChainPtr ring_search_nodal_locality(ChainPtr *, ChainPtr, Mesh);
static void compute_nodal_distance_local(Mesh);
static double mesh_size(Mesh);
static void standard_error(int, char *);
static void standard_output(char *);
static FILE *Open_and_Check_simulation_file(char *);
static void Check_Dirichlet_Boundary_Conditions(Boundaries, int);

/**********************************************************************/

Mesh GramsBox(char *Name_File, Time_Int_Params Parameters_Solver) {
  // Define mesh variable
  Mesh FEM_Mesh;

  int Num_nodal_rings = 4;

  // Read information in GramsBox and check sintax
  Nodes_Information Nodes_Info = Read_Nodal_Set_Information(Name_File);

  /*
    -> Read file with the mesh and allocate memory

    -> Get the sourrounding elements of each nodes

    -> Initialize nodal locality of each node :

      1º Define a pointer with the number of nodes close to each node

      2ª Define a table of sets with the list of nodes close to a node

      3ª Fill the table with the nodal locality and the pointer with the
        number of nodes close to each node

    -> Compute minimum nodal distance
  */
  puts("*************************************************");
  printf(" \t %s : %s \n", "* Read GID mesh in", Nodes_Info.Mesh_File);
  FEM_Mesh = ReadGidMesh__MeshTools__(Nodes_Info.Mesh_File);
  FEM_Mesh.NumNeighbour =
      (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh, sizeof(int));
  FEM_Mesh.NodeNeighbour = alloc_table__SetLib__(FEM_Mesh.NumNodesMesh);
  get_sourrounding_elements(FEM_Mesh);
  printf("\t \t %s : %s \n", "-> Compute sourrounding elements", "Done");

  FEM_Mesh.SizeNodalLocality_0 =
      (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh, sizeof(int));
  FEM_Mesh.NodalLocality_0 = alloc_table__SetLib__(FEM_Mesh.NumNodesMesh);
  fill_nodal_locality(FEM_Mesh, 1);
  printf("\t \t %s : %s \n", "-> Compute nodal neighborhood", "Done");

  if (Num_nodal_rings > 1) {
    FEM_Mesh.SizeNodalLocality =
        (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh, sizeof(int));
    FEM_Mesh.NodalLocality = alloc_table__SetLib__(FEM_Mesh.NumNodesMesh);
    fill_nodal_locality(FEM_Mesh, Num_nodal_rings);
    printf("\t \t %s : %s \n", "-> Compute extended nodal neighborhood",
           "Done");
  }

  FEM_Mesh.ActiveNode = (bool *)malloc(FEM_Mesh.NumNodesMesh * sizeof(bool));
  printf("\t \t %s : %s \n", "-> Allocate list of active nodes", "Done");

  FEM_Mesh.BoundaryNode = (bool *)malloc(FEM_Mesh.NumNodesMesh * sizeof(bool));
  printf("\t \t %s : %s \n", "-> Allocate list of boundary nodes", "Done");

  FEM_Mesh.Num_Particles_Node =
      (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh, sizeof(int));
  FEM_Mesh.List_Particles_Node = alloc_table__SetLib__(FEM_Mesh.NumNodesMesh);
  printf("\t \t %s : %s \n", "-> Allocate list of particles per node", "Done");

  //  FEM_Mesh.Num_Particles_Element = (int
  //  *)Allocate_ArrayZ(FEM_Mesh.NumElemMesh,sizeof(int));
  //  FEM_Mesh.List_Particles_Element = (ChainPtr
  //  *)malloc(FEM_Mesh.NumElemMesh*sizeof(ChainPtr)); printf("\t \t %s : %s
  //  \n","-> Allocate list of particle per element","Done");

  FEM_Mesh.h_avg =
      (double *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh, sizeof(double));
  compute_nodal_distance_local(FEM_Mesh);
  printf("\t \t %s : %s \n", "-> Compute local nodal distance", "Done");

  FEM_Mesh.DeltaX = mesh_size(FEM_Mesh);
  printf("\t \t %s : %f \n", "-> Compute mesh size", FEM_Mesh.DeltaX);

  /*
    -> Compute the boundary conditions
  */
  printf(" \t %s : %i \n", "* Number of boundaries",
         Nodes_Info.Number_Boundaries);
  FEM_Mesh.Bounds.NumBounds = Nodes_Info.Number_Boundaries;

  if (Nodes_Info.Number_Boundaries > 0) {
    if (strcmp(Formulation, "-u") == 0) {
      FEM_Mesh.Bounds = Read_u_Dirichlet_Boundary_Conditions__InOutFun__(
          Name_File, Nodes_Info.Number_Boundaries,
          Parameters_Solver.NumTimeStep);
    } else if (strcmp(Formulation, "-up") == 0) {
      FEM_Mesh.Bounds = Read_u_Dirichlet_Boundary_Conditions__InOutFun__(
          Name_File, Nodes_Info.Number_Boundaries,
          Parameters_Solver.NumTimeStep);
    } else if (strcmp(Formulation, "-upw") == 0) {
      FEM_Mesh.Bounds = Read_upw_Dirichlet_Boundary_Conditions__InOutFun__(
          Name_File, Nodes_Info.Number_Boundaries,
          Parameters_Solver.NumTimeStep);
    }

    Check_Dirichlet_Boundary_Conditions(FEM_Mesh.Bounds, FEM_Mesh.NumNodesMesh);
  }

  return FEM_Mesh;
}

/*********************************************************************/

static Nodes_Information Read_Nodal_Set_Information(char *Name_File) {

  Nodes_Information Nodes_Info;

  // Generate read route
  char Route_Mesh[MAXC] = {0};
  generate_route(Route_Mesh, Name_File);

  // Counter for the number of boundaries
  int Number_Boundaries = 0;

  // Variables to read file
  char line[MAXC] = {0};
  char *words[MAXW] = {NULL};
  int nwords;
  int Num_line = 0;
  int Node_i = 0;

  // Boolean variables to open and close braces
  bool Start_Reading_Nodes_Info = false;
  int Number_Left_Braces = 0;
  int Number_Right_Braces = 0;

  // Open the mesh file
  FILE *MeshFile = Open_and_Check_simulation_file(Name_File);

  // Read the file line by line
  while (fgets(line, sizeof(line), MeshFile) != NULL) {

    // Parse the conten of the line
    nwords = parse(words, line, delimiters);

    // Count the number of boundary conditions
    if (Start_Reading_Nodes_Info) {

      // Check number of open and closed braces
      for (int i = 0; i < nwords; i++) {
        if (strcmp(words[i], "{") == 0) {
          Number_Left_Braces++;
        }
        if (strcmp(words[i], "}") == 0) {
          Number_Right_Braces++;
        }
        if (strcmp(words[i], "GramsBoundary") == 0) {
          Number_Boundaries++;
        }
      }

    } else if (Start_Reading_Nodes_Info &&
               (Number_Left_Braces == Number_Right_Braces)) {
      break;
    }

    if ((nwords >= 5) && (strcmp(words[0], "GramsBox") == 0)) {

      if ((strcmp(words[1], "Type") == 0) && (strcmp(words[3], "File") == 0)) {
        Start_Reading_Nodes_Info = true;

        // Check kind of mesh generator
        if (strcmp(words[2], "GID") != 0) {
          sprintf(Error_message, "%s", "Unrecognised kind of mesh");
          standard_error(Num_line, Name_File);
        } else {
          strcpy(Nodes_Info.Mesh_Generator, words[2]);
        }

        // Check file with the mesh
        sprintf(Nodes_Info.Mesh_File, "%s%s", Route_Mesh, words[4]);
        Check_Mesh_File(Nodes_Info.Mesh_File);

        // Check number of open and closed braces
        for (int i = 0; i < nwords; i++) {
          if (strcmp(words[i], "{") == 0) {
            Number_Left_Braces++;
          }
          if (strcmp(words[i], "}") == 0) {
            Number_Right_Braces++;
          }
        }

      } else {
        sprintf(Error_message, "GramsBox has an unrecognised structure");
        standard_error(Num_line, Name_File);
      }
    }

    // Update the linea counter
    Num_line++;
  }

  if (!Start_Reading_Nodes_Info) {
    fprintf(stderr, "%s", "GramsBox keyword not founded");
    exit(EXIT_FAILURE);
  }

  if (Number_Left_Braces != Number_Right_Braces) {
    fprintf(stderr, "%s", "Umbalanced number of braces in GramsBox");
    exit(EXIT_FAILURE);
  }

  // Write the number of boundaries in the information structure
  Nodes_Info.Number_Boundaries = Number_Boundaries;

  //  At the end close the mesh file
  fclose(MeshFile);

  return Nodes_Info;
}

/*********************************************************************/

static void get_sourrounding_elements(Mesh FEM_Mesh) {

  /* Variable declaration */
  int *Element_Connectivity;
  int NumNodesElem;

  /* 1º Start the search of neighbour for each node */
  for (int i = 0; i < FEM_Mesh.NumNodesMesh; i++) {
    /* 2º Loop over all the elements in the mesh */
    for (int j = 0; j < FEM_Mesh.NumElemMesh; j++) {
      NumNodesElem = FEM_Mesh.NumNodesElem[j];
      Element_Connectivity =
          set_to_memory__SetLib__(FEM_Mesh.Connectivity[j], NumNodesElem);

      /* 3º Loop over the all the node in an element */
      for (int k = 0; k < NumNodesElem; k++) {
        /* 4º If my node belong to the element */
        if (Element_Connectivity[k] == i) {

          /* 5º Introduce the element in the chain */
          push__SetLib__(&FEM_Mesh.NodeNeighbour[i], j);

          /* 6º Update the counter */
          FEM_Mesh.NumNeighbour[i] += 1;
        }
      }

      /* Free memory */
      free(Element_Connectivity);
    }

    if (FEM_Mesh.NumNeighbour[i] == 0) {
      fprintf(stderr, "%s : %i \n",
              "Error computing the sourrounding elements of", i);
      exit(EXIT_FAILURE);
    }
  }
}

/*********************************************************************/

static void fill_nodal_locality(Mesh FEM_Mesh, int Num_nodal_rings) {
  /*
    Auxiliar variables for the nodal neighborhood reconstruction
  */
  int k_nodal_ring;    // Current search ring
  ChainPtr Search_Set; // Auxiliar set for recursive search

  for (int i = 0; i < FEM_Mesh.NumNodesMesh; i++) {

    if (Num_nodal_rings == 1) {
      FEM_Mesh.NodalLocality_0[i] = node_I_locality(i, FEM_Mesh);
      FEM_Mesh.SizeNodalLocality_0[i] =
          lenght__SetLib__(FEM_Mesh.NodalLocality_0[i]);
    } else if (Num_nodal_rings > 1) {

      k_nodal_ring = 0;
      Search_Set = NULL;
      push__SetLib__(&Search_Set, i);

      while (k_nodal_ring < Num_nodal_rings) {

        Search_Set = ring_search_nodal_locality(&FEM_Mesh.NodalLocality[i],
                                                Search_Set, FEM_Mesh);

        k_nodal_ring++;
      }

      FEM_Mesh.SizeNodalLocality[i] =
          lenght__SetLib__(FEM_Mesh.NodalLocality[i]);
    }
  }
}

/*********************************************************************/

static ChainPtr node_I_locality(int I, Mesh FEM_Mesh) {

  /* Define output */
  ChainPtr Nodes = NULL;

  /* Number of elements sourronding the node */
  int NumNeighbour = FEM_Mesh.NumNeighbour[I];

  /* Index of the elements sourronding the node */
  int *NodeNeighbour =
      set_to_memory__SetLib__(FEM_Mesh.NodeNeighbour[I], NumNeighbour);

  /* Table with the nodes of each element */
  ChainPtr *Table_ElemNodes = malloc(NumNeighbour * sizeof(ChainPtr));

  /* Fill each position of the table with a list of nodes in the element */
  for (int i = 0; i < NumNeighbour; i++) {
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[NodeNeighbour[i]];
  }

  /* Free table with elements */
  free(NodeNeighbour);

  /* Get the union of this nodes */
  Nodes = union__SetLib__(Table_ElemNodes, NumNeighbour);

  /* Return nodes close to the node I */
  return Nodes;
}

/*********************************************************************/

static ChainPtr ring_search_nodal_locality(ChainPtr *Set_k, ChainPtr Search_Set,
                                           Mesh FEM_Mesh) {

  /*
    Variables
  */
  ChainPtr aux_Set = NULL;
  ChainPtr new_Search_Set = NULL;
  ChainPtr i_Search_Set = Search_Set;

  /*
    Loop in the search set
  */
  while (i_Search_Set != NULL) {

    /*
      For each node in the set, get the closest node
    */
    aux_Set = node_I_locality(i_Search_Set->Idx, FEM_Mesh);

    /*
      Loop in the closest set of nodes of the original set
    */
    while (aux_Set != NULL) {
      /*
        Update the set with the new node and the search set
      */
      if (!inout__SetLib__(*Set_k, aux_Set->Idx)) {
        push__SetLib__(Set_k, aux_Set->Idx);
        push__SetLib__(&new_Search_Set, aux_Set->Idx);
      }

      /*
        Update second iterator
      */
      aux_Set = aux_Set->next;
    }

    /*
      Free auxiliar set
    */
    free__SetLib__(&aux_Set);

    /*
      Update first iterator
    */
    i_Search_Set = i_Search_Set->next;
  }

  /*
    Free old search list
  */
  free__SetLib__(&Search_Set);

  return new_Search_Set;
}

/*********************************************************************/

static void compute_nodal_distance_local(Mesh FEM_Mesh) {
  int Ndim = NumberDimensions;
  int A, B;
  int NumCloseNodes_A;

  ChainPtr NodalLocality_A = NULL;

  Matrix h_AB = alloc__MatrixLib__(1, Ndim);
  double avg_h_A;

  for (A = 0; A < FEM_Mesh.NumNodesMesh; A++) {

    NodalLocality_A = node_I_locality(A, FEM_Mesh);

    avg_h_A = 0.0;
    NumCloseNodes_A = 0;

    /*
      Loop in the closest set of nodes
    */
    while (NodalLocality_A != NULL) {
      /*
        Get the index of the node
      */
      B = NodalLocality_A->Idx;

      if (A != B) {

        for (int i = 0; i < Ndim; i++) {
          h_AB.nV[i] =
              FEM_Mesh.Coordinates.nM[B][i] - FEM_Mesh.Coordinates.nM[A][i];
        }

        avg_h_A += norm__MatrixLib__(h_AB, 2);

        NumCloseNodes_A++;
      }

      NodalLocality_A = NodalLocality_A->next;
    }

    FEM_Mesh.h_avg[A] = avg_h_A / (double)NumCloseNodes_A;

    free__SetLib__(&NodalLocality_A);
  }

  free__MatrixLib__(h_AB);
}

/*********************************************************************/

static double mesh_size(Mesh FEM_Mesh)
/*
  Function to get the minimum mesh size.
*/
{

  /* Auxiliar variables of the function */
  int Ndim = NumberDimensions;
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem;              /* Number of nodes of each element */
  ChainPtr Element_Connectivity; /* Connectivity of the element */
  ChainPtr Element_Connectivity_Circular;
  double MinElementSize = 10e16;

  /* Loop over the elements in the mesh */
  for (int i = 0; i < NumElemMesh; i++) {
    /* Connectivity of the element */
    NumNodesElem = FEM_Mesh.NumNodesElem[i];
    Element_Connectivity = copy__SetLib__(FEM_Mesh.Connectivity[i]);
    Element_Connectivity_Circular =
        create_circular_set__SetLib__(Element_Connectivity);

    /* Get the gradient of the element for each node */
    if ((Ndim == 2) && (strcmp(FEM_Mesh.TypeElem, "Triangle") == 0) &&
        (NumNodesElem == 3)) {
      MinElementSize =
          DMIN(MinElementSize, min_DeltaX__T3__(Element_Connectivity_Circular,
                                                FEM_Mesh.Coordinates));
    } else if ((Ndim == 2) && (strcmp(FEM_Mesh.TypeElem, "Triangle") == 0) &&
               (NumNodesElem == 6)) {
      MinElementSize =
          DMIN(MinElementSize, min_DeltaX__T3__(Element_Connectivity_Circular,
                                                FEM_Mesh.Coordinates));
    } else if ((Ndim == 2) &&
               (strcmp(FEM_Mesh.TypeElem, "Quadrilateral") == 0) &&
               (NumNodesElem == 4)) {
      MinElementSize =
          DMIN(MinElementSize, min_DeltaX__Q4__(Element_Connectivity_Circular,
                                                FEM_Mesh.Coordinates));
    } else if ((Ndim == 3) && (strcmp(FEM_Mesh.TypeElem, "Tetrahedra") == 0) &&
               (NumNodesElem == 4)) {
      MinElementSize =
          DMIN(MinElementSize, min_DeltaX__T4__(Element_Connectivity_Circular,
                                                FEM_Mesh.Coordinates));
    } else if ((Ndim == 3) && (strcmp(FEM_Mesh.TypeElem, "Hexahedra") == 0) &&
               (NumNodesElem == 8)) {
      MinElementSize =
          DMIN(MinElementSize, min_DeltaX__H8__(Element_Connectivity_Circular,
                                                FEM_Mesh.Coordinates));
    } else {
      printf("%s : %s %i %s \n", "Error in mesh_size", "Element with ",
             NumNodesElem, "nodes is not implemented !!!");
      exit(EXIT_FAILURE);
    }

    /* Free memory */
    free__SetLib__(&Element_Connectivity);
  }

  return MinElementSize;
}

/*********************************************************************/

static void standard_error(int Num_line, char *Name_File) {
  fprintf(stderr, "GramsBox : Error in line %i while reading %s \n", Num_line,
          Name_File);
  fprintf(stderr, "%s \n", Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/

static void standard_output(char *Status_message) {
  fprintf(stdout, "%s \n", Status_message);
}

/***************************************************************************/

static void Check_Mesh_File(char *PATH_Name) {
  struct stat info;
  stat(PATH_Name, &info);
  char Error_message[MAXW];

  if (!S_ISREG(info.st_mode)) {
    fprintf(stderr, "GramsBox : Mesh file %s does not exists \n", PATH_Name);
    exit(EXIT_FAILURE);
  }
}

/**********************************************************************/

static FILE *Open_and_Check_simulation_file(char *Name_File) {
  FILE *Simulation_file = fopen(Name_File, "r");

  if (Simulation_file == NULL) {
    fprintf(stderr, "GramsBox : Incorrect lecture of %s \n", Name_File);
    exit(EXIT_FAILURE);
  }

  return Simulation_file;
}

/***************************************************************************/

static void Check_Dirichlet_Boundary_Conditions(Boundaries Dirichlet_BCC,
                                                int NumNodesMesh) {
  /*
  Define auxilar variables
*/
  int Number_of_BCC = Dirichlet_BCC.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */

  for (int i = 0; i < Number_of_BCC; i++) {

    NumNodesBound = Dirichlet_BCC.BCC_i[i].NumNodes;

    for (int j = 0; j < NumNodesBound; j++) {
      if ((Dirichlet_BCC.BCC_i[i].Nodes[j] > NumNodesMesh) ||
          (Dirichlet_BCC.BCC_i[i].Nodes[j] < 0)) {
        fprintf(stderr, "%s : %s 0 - %i \n", "Error in GramsBox()",
                "The index of the nodes with Dirichlet BCC should be between",
                NumNodesMesh);
        exit(EXIT_FAILURE);
      }
    }
  }
}

/***************************************************************************/