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

#ifdef _WIN32
static char * delimiters = " \r\n\t";
#else
static char * delimiters = " \t\n"; 
#endif

static char Error_message[MAXW];

static Mesh_Information Read_Mesh_Information(char * );
static void Fill_Coordinates(char * , Mesh_Information, Matrix);
static void Fill_Conectivity(char * , Mesh_Information, ChainPtr *, int *);
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
    Fill information of the mesh
  */
  GID_Mesh.NumNodesMesh = Mesh_Info.NumNodesMesh;
  GID_Mesh.NumElemMesh = Mesh_Info.NumElemMesh;
  GID_Mesh.Dimension = Mesh_Info.Dimension;
  strcpy(GID_Mesh.TypeElem,Mesh_Info.ElemType);
  
  /*
    Select the kind of element of the mesh
  */
  if((strcmp(GID_Mesh.TypeElem,"Triangle") == 0) && (Mesh_Info.NumNodesElem == 3))
  {
    GID_Mesh.N_ref = N__T3__;
    GID_Mesh.dNdX_ref = dN_Ref__T3__;
  }
  else if((strcmp(GID_Mesh.TypeElem,"Quadrilateral") == 0) && (Mesh_Info.NumNodesElem == 4))
  {
    GID_Mesh.N_ref = N__Q4__;
    GID_Mesh.dNdX_ref = dN_Ref__Q4__;
  }
  else if((strcmp(GID_Mesh.TypeElem,"Tetrahedra") == 0) && (Mesh_Info.NumNodesElem == 4))
  {
    GID_Mesh.N_ref = N__T4__;
    GID_Mesh.dNdX_ref = dN_Ref__T4__;
  }
  else if((strcmp(GID_Mesh.TypeElem,"Hexahedra") == 0) && (Mesh_Info.NumNodesElem == 8))
  {
    GID_Mesh.N_ref = N__H8__;
    GID_Mesh.dNdX_ref = dN_Ref__H8__; 
  }
  else
  {
    sprintf(Error_message,"This kind of element is not suported");
    standard_error(0, MeshName);
  }

  /* 
    Allocate data 
  */
  GID_Mesh.Coordinates  = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
  GID_Mesh.Connectivity = alloc_table__SetLib__(GID_Mesh.NumElemMesh);  
  GID_Mesh.NumNodesElem = (int *)Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
  GID_Mesh.NumParticles = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));


  /*  
    Fill coordinates and connectivity
  */
  Fill_Coordinates(MeshName,Mesh_Info,GID_Mesh.Coordinates);
  Fill_Conectivity(MeshName,Mesh_Info,GID_Mesh.Connectivity,GID_Mesh.NumNodesElem);

  for(int i = 0; i<GID_Mesh.NumElemMesh ; i++)
  {
    print__SetLib__(GID_Mesh.Connectivity[i]); 
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
    if(strcmp(words[0],"Coordinates") == 0)
    {
      Mesh_Info.NumNodesMesh = -2;
      Count_Nodes = true;
    }

    if(Count_Nodes)
    {
      Mesh_Info.NumNodesMesh++;
    }

    if((strcmp(words[0],"End") == 0) && (strcmp(words[1],"Coordinates") == 0))
    {
      Count_Nodes = false;
    }

    // Start counting the number of elements
    if(strcmp(words[0],"Elements") == 0)
    {
      Mesh_Info.NumElemMesh = -2;
      Count_Elements = true;
    }

    if(Count_Elements)
    {
      Mesh_Info.NumElemMesh++;
    }

    if((strcmp(words[0],"End") == 0) && (strcmp(words[1],"Elements") == 0))
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
    if(strcmp(words[0],"Coordinates") == 0)
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

    if((strcmp(words[0],"End") == 0) && (strcmp(words[1],"Coordinates") == 0))
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

static void Fill_Conectivity(
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
    if(strcmp(words[0],"Elements") == 0)
    {
      Read_Element_Connectivity = true;
    }

    if((Read_Element_Connectivity == true) && (nwords == (Mesh_Info.NumNodesElem + 1)))
    {

      for(int j = 0 ; j<Mesh_Info.NumNodesElem ; j++)
      {
        push__SetLib__(&Connectivity[Element_i], atoi(words[j+1]));
      } 

      NumNodesElem[Element_i] = Mesh_Info.NumNodesElem;

      Element_i++;
    }

    if((strcmp(words[0],"End") == 0) && (strcmp(words[1],"Elements") == 0))
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
