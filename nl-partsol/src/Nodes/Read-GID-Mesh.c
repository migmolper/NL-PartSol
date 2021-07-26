#include "nl-partsol.h"


/*
  Local structures
*/
typedef struct
{
  int   Dimension;
  char  ElemType [100];
  int   NumNodesElem;
  int   NumNodesMesh;
  int   NumElemMesh;

} Mesh_Information;

/*
  Auxiliar functions and variables
*/

#ifdef __linux__
  static char * delimiters = " \t\n"; 
#endif    

#ifdef __APPLE__
  static char * delimiters = " \t\n"; 
#endif
  
#ifdef _WIN32 
  static char * delimiters = " \r\n\t";
#endif


static char Error_message[MAXW];

static Mesh_Information Read_Mesh_Information(char *);
static void Fill_Coordinates(char * , Mesh_Information, Matrix);
static void Fill_Linear_Conectivity(char * , Mesh_Information, ChainPtr *, int *);
static void Fill_Quadratic_Conectivity(char *,Mesh_Information,ChainPtr *, int *, int *);
static void standard_error(int, char *);
static void standard_output(char *);
static FILE * Open_and_Check_simulation_file(char *);
/***************************************************************************/

Mesh ReadGidMesh__MeshTools__(char * MeshName)
{
  /* Output element mesh */
  Mesh GID_Mesh;

  /* Pointer to the file */
  FILE * MeshFile;
     
  /*
    Read information of the mesh such as number of nodes, number of elements, etc
  */
  Mesh_Information Mesh_Info = Read_Mesh_Information(MeshName);
  


  /*
    Select the kind of element of the mesh and assign functions
  */
  if((strcmp(Mesh_Info.ElemType,"Triangle") == 0) && (Mesh_Info.NumNodesElem == 3))
  {
    GID_Mesh.N_ref = N__T3__;
    GID_Mesh.dNdX_ref = dN_Ref__T3__;
    GID_Mesh.dNdX = dN__T3__;
    GID_Mesh.volume_Element = volume__T3__;
    GID_Mesh.In_Out_Element = in_out__T3__;
    GID_Mesh.NumNodesMesh = Mesh_Info.NumNodesMesh;
    GID_Mesh.NumElemMesh = Mesh_Info.NumElemMesh;
    GID_Mesh.Dimension = Mesh_Info.Dimension;
    strcpy(GID_Mesh.TypeElem,Mesh_Info.ElemType);
    GID_Mesh.Coordinates  = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
    GID_Mesh.Connectivity = alloc_table__SetLib__(GID_Mesh.NumElemMesh);  
    GID_Mesh.NumNodesElem = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
    GID_Mesh.Num_Particles_Node = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
    GID_Mesh.Locking_Control_Fbar = false;
    Fill_Coordinates(MeshName,Mesh_Info,GID_Mesh.Coordinates);
    Fill_Linear_Conectivity(MeshName,Mesh_Info,GID_Mesh.Connectivity,GID_Mesh.NumNodesElem);
  }
  else if((strcmp(Mesh_Info.ElemType,"Triangle") == 0) && (Mesh_Info.NumNodesElem == 6))
  {
    GID_Mesh.N_ref = N__T3__;
    GID_Mesh.dNdX_ref = dN_Ref__T3__;
    GID_Mesh.dNdX = dN__T3__;
    GID_Mesh.volume_Element = volume__T3__;
    GID_Mesh.In_Out_Element = in_out__T3__;
    GID_Mesh.NumNodesMesh = Mesh_Info.NumNodesMesh;
    GID_Mesh.NumElemMesh = Mesh_Info.NumElemMesh*4;
    GID_Mesh.Num_Patch_Mesh = Mesh_Info.NumElemMesh;
    GID_Mesh.Dimension = Mesh_Info.Dimension;
    strcpy(GID_Mesh.TypeElem,Mesh_Info.ElemType);
    GID_Mesh.Coordinates  = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
    GID_Mesh.Connectivity = alloc_table__SetLib__(GID_Mesh.NumElemMesh);  
    GID_Mesh.NumNodesElem = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
    GID_Mesh.Num_Particles_Node = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
    GID_Mesh.Locking_Control_Fbar = true;
    GID_Mesh.Idx_Patch = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
    GID_Mesh.Vol_Patch_n = (double *)Allocate_ArrayZ(Mesh_Info.NumElemMesh,sizeof(double));
    GID_Mesh.Vol_Patch_n1 = (double *)Allocate_ArrayZ(Mesh_Info.NumElemMesh,sizeof(double));
    Fill_Coordinates(MeshName,Mesh_Info,GID_Mesh.Coordinates);
    Fill_Quadratic_Conectivity(MeshName,Mesh_Info,GID_Mesh.Connectivity,GID_Mesh.Idx_Patch,GID_Mesh.NumNodesElem);
  }
  else if((strcmp(Mesh_Info.ElemType,"Quadrilateral") == 0) && (Mesh_Info.NumNodesElem == 4))
  {
    GID_Mesh.N_ref = N__Q4__;
    GID_Mesh.dNdX_ref = dN_Ref__Q4__;
    GID_Mesh.dNdX = dN__Q4__;
    GID_Mesh.volume_Element = volume__Q4__;
    GID_Mesh.In_Out_Element = in_out__Q4__;
    GID_Mesh.NumNodesMesh = Mesh_Info.NumNodesMesh;
    GID_Mesh.NumElemMesh = Mesh_Info.NumElemMesh;
    GID_Mesh.Dimension = Mesh_Info.Dimension;
    strcpy(GID_Mesh.TypeElem,Mesh_Info.ElemType);
    GID_Mesh.Coordinates  = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
    GID_Mesh.Connectivity = alloc_table__SetLib__(GID_Mesh.NumElemMesh);  
    GID_Mesh.NumNodesElem = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
    GID_Mesh.Num_Particles_Node = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
    GID_Mesh.Locking_Control_Fbar = false;
    Fill_Coordinates(MeshName,Mesh_Info,GID_Mesh.Coordinates);
    Fill_Linear_Conectivity(MeshName,Mesh_Info,GID_Mesh.Connectivity,GID_Mesh.NumNodesElem);
  }
  else if((strcmp(Mesh_Info.ElemType,"Tetrahedra") == 0) && (Mesh_Info.NumNodesElem == 4))
  {
    GID_Mesh.N_ref = N__T4__;
    GID_Mesh.dNdX_ref = dN_Ref__T4__;
    GID_Mesh.dNdX = dN__T4__;
    GID_Mesh.volume_Element = volume__T4__;
    GID_Mesh.In_Out_Element = in_out__T4__;
    GID_Mesh.NumNodesMesh = Mesh_Info.NumNodesMesh;
    GID_Mesh.NumElemMesh = Mesh_Info.NumElemMesh;
    GID_Mesh.Dimension = Mesh_Info.Dimension;
    strcpy(GID_Mesh.TypeElem,Mesh_Info.ElemType);
    GID_Mesh.Coordinates  = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
    GID_Mesh.Connectivity = alloc_table__SetLib__(GID_Mesh.NumElemMesh);  
    GID_Mesh.NumNodesElem = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
    GID_Mesh.Num_Particles_Node = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
    GID_Mesh.Locking_Control_Fbar = false;
    Fill_Coordinates(MeshName,Mesh_Info,GID_Mesh.Coordinates);
    Fill_Linear_Conectivity(MeshName,Mesh_Info,GID_Mesh.Connectivity,GID_Mesh.NumNodesElem);
  }
  else if((strcmp(Mesh_Info.ElemType,"Hexahedra") == 0) && (Mesh_Info.NumNodesElem == 8))
  {
    GID_Mesh.N_ref = N__H8__;
    GID_Mesh.dNdX_ref = dN_Ref__H8__;
    GID_Mesh.dNdX = dN__H8__;
    GID_Mesh.volume_Element = volume__H8__;
    GID_Mesh.In_Out_Element = in_out__H8__;
    GID_Mesh.NumNodesMesh = Mesh_Info.NumNodesMesh;
    GID_Mesh.NumElemMesh = Mesh_Info.NumElemMesh;
    GID_Mesh.Dimension = Mesh_Info.Dimension;
    strcpy(GID_Mesh.TypeElem,Mesh_Info.ElemType);
    GID_Mesh.Coordinates  = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
    GID_Mesh.Connectivity = alloc_table__SetLib__(GID_Mesh.NumElemMesh);  
    GID_Mesh.NumNodesElem = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
    GID_Mesh.Num_Particles_Node = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
    GID_Mesh.Locking_Control_Fbar = false;
    Fill_Coordinates(MeshName,Mesh_Info,GID_Mesh.Coordinates);
    Fill_Linear_Conectivity(MeshName,Mesh_Info,GID_Mesh.Connectivity,GID_Mesh.NumNodesElem);
  }
  else
  {
    sprintf(Error_message,"Non suported element -> TypeElem : %s ; NumNodesElem : %i\n",
     Mesh_Info.ElemType,Mesh_Info.NumNodesElem);
    standard_error(0, MeshName);
  }


  return GID_Mesh;
}

/***************************************************************************/

static Mesh_Information Read_Mesh_Information(char * MeshName)
{

  Mesh_Information Mesh_Info;

  /* Variables to read file */
  char line[MAXC] = {0};
  char * words[MAXW] = {NULL};
  int nwords;
  int Num_line = 0;

  /*
    Auxiliar variable to count the number of 
    nodes and elements of the mesh
  */
  bool Count_Nodes = false;
  bool Count_Elements = false;
  
  /* Open the mesh file */
  FILE * MeshFile = Open_and_Check_simulation_file(MeshName);

  /*
    Read the file line by line
  */
  while(fgets(line, sizeof(line), MeshFile) != NULL)
  {

    // Parse the conten of the line
    nwords = parse (words, line, delimiters);
 
    //  Read the first line with the header
    if(Num_line == 0)
    {
      if((nwords == 7) && 
        (strcmp(words[0],"MESH") == 0) &&
        (strcmp(words[1],"dimension") == 0) &&
        (strcmp(words[3],"ElemType") == 0) &&
        (strcmp(words[5],"Nnode") == 0))
      {
        Mesh_Info.Dimension = atoi(words[2]);
        strcpy(Mesh_Info.ElemType,words[4]);
        Mesh_Info.NumNodesElem = atoi(words[6]);
      }
      else
      {
        sprintf(Error_message,"The header of a GID has a non-suported structure");
        standard_error(Num_line, MeshName);
      }
    }

    // Start counting the number of nodes
    if((nwords == 1) && (strcmp(words[0],"Coordinates") == 0))
    {
      Mesh_Info.NumNodesMesh = -2;
      Count_Nodes = true;
    }

    if(Count_Nodes)
    {
      Mesh_Info.NumNodesMesh++;
    }

    if((nwords == 2) && (strcmp(words[0],"End") == 0) && (strcmp(words[1],"Coordinates") == 0))
    {
      Count_Nodes = false;
    }

    // Start counting the number of elements
    if((nwords == 1) && (strcmp(words[0],"Elements") == 0))
    {
      Mesh_Info.NumElemMesh = -2;
      Count_Elements = true;
    }

    if(Count_Elements)
    {
      Mesh_Info.NumElemMesh++;
    }

    if((nwords == 2) && (strcmp(words[0],"End") == 0) && (strcmp(words[1],"Elements") == 0))
    {
      Count_Elements = false;
    }

    // Update the linea counter
    Num_line++;
  }

  /* At the end close the mesh file */
  fclose(MeshFile);


  return Mesh_Info;
}

/***************************************************************************/

static void Fill_Coordinates(
  char * MeshName,
  Mesh_Information Mesh_Info,
  Matrix Coordinates)
{

  /* Variables to read file */
  char line[MAXC] = {0};
  char * words[MAXW] = {NULL};
  int nwords;
  int Num_line = 0;
  int Node_i = 0;

  /*
    Auxiliar variable to count the number of 
    nodes and elements of the mesh
  */
  bool Read_Coordinates_Nodes = false;

  /* Open the mesh file */
  FILE * MeshFile = Open_and_Check_simulation_file(MeshName);

  /*
    Read the file line by line
  */
  while(fgets(line, sizeof(line), MeshFile) != NULL)
  {

    // Parse the conten of the line
    nwords = parse (words, line, delimiters);

    // Start reading the coordinates of the nodes
    if((nwords == 1) && (strcmp(words[0],"Coordinates") == 0))
    {
      Read_Coordinates_Nodes = true;
    }

    if((Read_Coordinates_Nodes == true) && (nwords == 4))
    {

      for(int j = 0 ; j<Mesh_Info.Dimension ; j++)
      {
        Coordinates.nM[Node_i][j] = atof(words[j+1]);
      } 

      Node_i++;
    }

    if((nwords == 2) && (strcmp(words[0],"End") == 0) && (strcmp(words[1],"Coordinates") == 0))
    {
      Read_Coordinates_Nodes = false;
    }

    /* 
      Update the linea counter
    */
    Num_line++;
  }

  /* 
    At the end close the mesh file
  */
  fclose(MeshFile);
}

/***************************************************************************/

static void Fill_Linear_Conectivity(
  char * MeshName,
  Mesh_Information Mesh_Info,
  ChainPtr * Connectivity, 
  int * NumNodesElem)
{

  /* Variables to read file */
  char line[MAXC] = {0};
  char * words[MAXW] = {NULL};
  int nwords;
  int Num_line = 0;
  int Element_i = 0;

  /*
    Auxiliar variable to count the number of 
    nodes and elements of the mesh
  */
  bool Read_Element_Connectivity = false;

  /* Open the mesh file */
  FILE * MeshFile = Open_and_Check_simulation_file(MeshName);

  /*
    Read the file line by line
  */
  while(fgets(line, sizeof(line), MeshFile) != NULL)
  {

    // Parse the conten of the line
    nwords = parse (words, line, delimiters);

    // Start reading the coordinates of the nodes
    if((nwords == 1) && (strcmp(words[0],"Elements") == 0))
    {
      Read_Element_Connectivity = true;
    }

    if((Read_Element_Connectivity == true) && (nwords == (Mesh_Info.NumNodesElem + 1)))
    {

      for(int j = 0 ; j<Mesh_Info.NumNodesElem ; j++)
      {
        push__SetLib__(&Connectivity[Element_i], atoi(words[j+1]) - 1);
      } 

      NumNodesElem[Element_i] = Mesh_Info.NumNodesElem;

      Element_i++;
    }

    if((nwords == 2) && (strcmp(words[0],"End") == 0) && (strcmp(words[1],"Elements") == 0))
    {
      Read_Element_Connectivity = false;
    }

    /* 
      Update the linea counter
    */
    Num_line++;
  }

  /* 
    At the end close the mesh file
  */
  fclose(MeshFile);

}

/***************************************************************************/

static void Fill_Quadratic_Conectivity(
  char * MeshName,
  Mesh_Information Quadratic_Mesh_Info,
  ChainPtr * Connectivity, 
  int * Idx_Patch,
  int * NumNodesElem)
{

  ChainPtr * Quadratic_Connectivity = alloc_table__SetLib__(Quadratic_Mesh_Info.NumElemMesh);

  int NumInternalElements = 4;
  int Num_Linear_Element_Mesh = NumInternalElements*Quadratic_Mesh_Info.NumElemMesh;
  int * Quadratic_Connectivity_i;

  Fill_Linear_Conectivity(MeshName,Quadratic_Mesh_Info,Quadratic_Connectivity,NumNodesElem);

  for(int i = 0 ; i<Quadratic_Mesh_Info.NumElemMesh ; i++)
  {

    Quadratic_Connectivity_i = set_to_memory__SetLib__(Quadratic_Connectivity[i], 6);

    // Internal element 1
    push__SetLib__(&Connectivity[i*NumInternalElements + 0],Quadratic_Connectivity_i[5]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 0],Quadratic_Connectivity_i[2]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 0],Quadratic_Connectivity_i[0]);
    Idx_Patch[i*NumInternalElements + 0] = i;

    // Internal element 2
    push__SetLib__(&Connectivity[i*NumInternalElements + 1],Quadratic_Connectivity_i[4]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 1],Quadratic_Connectivity_i[2]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 1],Quadratic_Connectivity_i[1]);
    Idx_Patch[i*NumInternalElements + 1] = i;

    // Internal element 3
    push__SetLib__(&Connectivity[i*NumInternalElements + 2],Quadratic_Connectivity_i[3]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 2],Quadratic_Connectivity_i[1]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 2],Quadratic_Connectivity_i[0]);
    Idx_Patch[i*NumInternalElements + 2] = i;

    // Internal element 4
    push__SetLib__(&Connectivity[i*NumInternalElements + 3],Quadratic_Connectivity_i[2]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 3],Quadratic_Connectivity_i[1]);
    push__SetLib__(&Connectivity[i*NumInternalElements + 3],Quadratic_Connectivity_i[0]);
    Idx_Patch[i*NumInternalElements + 3] = i;

    free(Quadratic_Connectivity_i);

  }

  /*
    Fill number of nodes per element
  */
  for(int i = 0 ; i<Num_Linear_Element_Mesh ; i++)
  {
    NumNodesElem[i] = 3;
  }
  
  free_table__SetLib__(Quadratic_Connectivity,Quadratic_Mesh_Info.NumElemMesh); 

}

/***************************************************************************/

static void standard_error(
  int Num_line, 
  char * MeshName)
{
  fprintf(stderr,"ReadGidMesh__MeshTools__ : Error in line %i while reading %s \n",Num_line,MeshName);
  fprintf(stderr,"%s \n",Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/

static void standard_output(
  char * Status_message)
{
  fprintf(stdout,"%s \n",Status_message);
}

/**********************************************************************/

static FILE * Open_and_Check_simulation_file(
  char * Name_File)
{
  FILE * Simulation_file = fopen(Name_File,"r");  
  
  if (Simulation_file==NULL)
  {
    fprintf(stderr,"ReadGidMesh__MeshTools__  : Incorrect lecture of %s \n",Name_File);
    exit(EXIT_FAILURE);
  }  

  return Simulation_file;
}

/***************************************************************************/
