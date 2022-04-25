#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "InOutFun.h"

/*
  Local structures
*/
typedef struct {
  int *Dir;
  Curve *Value;

} BCC_Properties;

/*
  Auxiliar functions and variables
*/

#ifdef _WIN32
static char *delimiters_1 = " ,()\r\n\t";
static char *delimiters_2 = " =\t\r\n";
#else
static char *delimiters_1 = " ,()\n\t";
static char *delimiters_2 = " =\t\n";
#endif
static char *delimiters_3 = "=";

static char Error_message[MAXW];

static BCC_Properties Read_Boundary_Conditions_Properties(FILE *, char *, int,
                                                          int);
static void Check_Curve_File(char *);
static void standard_error();
static void standard_output(char *);
static FILE *Open_and_Check_simulation_file(char *);
static void active_direction(int *, int);

/**********************************************************************/

Boundaries
Read_u_Neumann_Boundary_Conditions__InOutFun__(char *Name_File, int NumBounds,
                                               int GPxElement, int NumTimeStep)
/*

        Define-Neumann-Boundary(File=Right_contour.txt)
        {
                T.x CurveConstant.txt
                T.y NULL
        }

*/
{

  /*
          Boundaries
  */
  Boundaries Bounds;

  /* Simulation file */
  FILE *Sim_dat;

  /*
          Number of words
  */
  int Num_words_line;

  /*
          Parse lines of GramsBoundary
  */
  char Line_GramsBoundary[MAXC] = {0};
  char *Parse_GramsBoundary[MAXW] = {NULL};
  char *Parse_Nodes[MAXW] = {NULL};

  /*
          Boundaries iterator
  */
  int IndexBoundary = 0;

  /*
          Auxiliar boundary to generate the list of loaded material points
  */
  int NumElements;
  int Element_idx;
  int *ElementList;

  /*
          Read current line
  */
  int Num_line = 0;

  /*
          Parse file name with the list of nodes
  */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];
  ChainPtr Chain_Nodes = NULL;

  /*
          Parse GramsBoundary properties
  */
  char Line_Properties[MAXC] = {0};
  char *Parse_Properties[MAXW] = {NULL};
  char FileLoadRoute[MAXC];

  /*
          Assign the number of bounds
  */
  Bounds.NumBounds = NumBounds;

  /*
          Allocate boundaries
  */
  Bounds.BCC_i = (Load *)Allocate_Array(NumBounds, sizeof(Load));

  /*
          Open and check file
  */
  Sim_dat = Open_and_Check_simulation_file(Name_File);

  /*
          Generate route
  */
  generate_route(Route_Nodes, Name_File);

  /*
          Read Define-Neumann-Boundary line
  */
  while (fgets(Line_GramsBoundary, sizeof(Line_GramsBoundary), Sim_dat) !=
         NULL) {

    /*
            Update line number
    */
    Num_line++;

    /*
    Read the line with the space as separators
    */
    Num_words_line =
        parse(Parse_GramsBoundary, Line_GramsBoundary, delimiters_1);

    /*
            Find Define-Neumann-Boundary line
    */
    if ((Num_words_line > 0) &&
        (strcmp(Parse_GramsBoundary[0], "Define-Neumann-Boundary") == 0)) {

      /*
              File with the nodes
      */
      Num_words_line = parse(Parse_Nodes, Parse_GramsBoundary[1], delimiters_3);

      /*
              Check the sintax of the parameters
      */
      if ((Num_words_line != 2) || (strcmp(Parse_Nodes[0], "File") != 0)) {
        sprintf(
            Error_message,
            "Sintax error in line %i : Define-Neumann-Boundary(File=Nodes.txt)",
            Num_line);
        standard_error();
      }

      /*
      Read file with the nodes and get an array with the nodes
*/
      sprintf(FileNodesRoute, "%s%s", Route_Nodes, Parse_Nodes[1]);
      Chain_Nodes = File2Chain(FileNodesRoute);
      NumElements = lenght__SetLib__(Chain_Nodes);
      ElementList = set_to_memory__SetLib__(Chain_Nodes, NumElements);
      free__SetLib__(&Chain_Nodes);

      /*
              Use the auxilar list to fill the list of material points loaded
      */
      Bounds.BCC_i[IndexBoundary].NumNodes = GPxElement * NumElements;
      Bounds.BCC_i[IndexBoundary].Nodes =
          (int *)Allocate_ArrayZ(GPxElement * NumElements, sizeof(int));
      for (int i = 0; i < NumElements; i++) {
        Element_idx = ElementList[i];

        for (int j = 0; j < GPxElement; j++) {
          Bounds.BCC_i[IndexBoundary].Nodes[i * GPxElement + j] =
              Element_idx * GPxElement + j;
        }
      }

      free(ElementList);

      /*
      Number of dimensions of the BCC
*/
      Bounds.BCC_i[IndexBoundary].Dim = NumberDOF;

      /*
              Read parameters
      */
      BCC_Properties Parameters = Read_Boundary_Conditions_Properties(
          Sim_dat, Name_File, Bounds.BCC_i[IndexBoundary].NumNodes,
          NumTimeStep);

      /*
              Direction of the BCC
      */
      Bounds.BCC_i[IndexBoundary].Dir = Parameters.Dir;

      /*
              Curve for each dimension
      */
      Bounds.BCC_i[IndexBoundary].Value = Parameters.Value;

      /*
              Information of the BCC
      */
      strcpy(Bounds.BCC_i[IndexBoundary].Info, "Newmann boundary conditions");

      /*
          Increment the index
      */
      IndexBoundary++;
    }
  }

  return Bounds;
}

/**********************************************************************/

void Check_u_Neumann_Boundary_Conditions__InOutFun__(
    Boundaries Neumann_Contours, int NumParticles) {

  int NumContactForces = Neumann_Contours.NumBounds;
  int NumNodesLoad;

  for (int i = 0; i < NumContactForces; i++) {

    NumNodesLoad = Neumann_Contours.BCC_i[i].NumNodes;

    for (int j = 0; j < NumNodesLoad; j++) {
      if ((Neumann_Contours.BCC_i[i].Nodes[j] > NumParticles) ||
          (Neumann_Contours.BCC_i[i].Nodes[j] < 0)) {
        fprintf(stderr, "%s : %s 0 - %i \n",
                "Error in Define-Neumann-Boundary()",
                "The index of the particle with Neumann BCC should be between",
                NumParticles);
        exit(EXIT_FAILURE);
      }
    }
  }
}

/**********************************************************************/

static BCC_Properties Read_Boundary_Conditions_Properties(FILE *Simulation_file,
                                                          char *Name_File,
                                                          int NumNodes,
                                                          int NumTimeStep) {
  int Ndim = NumberDimensions;

  /*
          Output structure with the properties
  */
  BCC_Properties Properties;

  /*
          Variables for reading purposes
  */
  char Parameter_line[MAXC] = {0};
  char *Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /*
          Check variables for sintax
  */
  bool Is_Open = false;
  bool Is_Close = false;

  /*
          Count lines
  */
  int Num_line = 0;

  /*
          Generate route
  */
  char FileLoadRoute[MAXC] = {0};
  char Route_Nodes[MAXC] = {0};
  generate_route(Route_Nodes, Name_File);

  /*
          Initialise and allocate Dir vector for each DOF
  */
  Properties.Dir = (int *)calloc(NumTimeStep * NumberDOF, sizeof(int));

  /*
  Initialise and allocate curve for each DOF
*/
  Properties.Value = (Curve *)Allocate_Array(NumberDOF, sizeof(Curve));

  while (fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) !=
         NULL) {
    /*
            Update line number
    */
    Num_line++;

    /*
            Parse each line
    */
    Parser_status = parse(Parameter_pars, Parameter_line, delimiters_2);

    if ((Parser_status == 1) && (strcmp(Parameter_pars[0], "{") == 0)) {
      Is_Open = true;
    } else if ((Parser_status == 2) &&
               (strcmp(Parameter_pars[0], "T.x") == 0)) {
      if (strcmp(Parameter_pars[1], "NULL") != 0) {
        sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parameter_pars[1]);
        Properties.Value[0] = ReadCurve(FileLoadRoute);
        active_direction(&Properties.Dir[NumTimeStep * 0],
                         IMIN(NumTimeStep, Properties.Value[0].Num));
        printf(" \t %s (%s) : \n \t \t Number of particles = %i \n \t \t File "
               "curve %s \n",
               "-> BcNeumann ", Parameter_pars[0], NumNodes, FileLoadRoute);
      }
    } else if ((Parser_status == 2) &&
               (strcmp(Parameter_pars[0], "T.y") == 0)) {
      if (strcmp(Parameter_pars[1], "NULL") != 0) {
        sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parameter_pars[1]);
        Properties.Value[1] = ReadCurve(FileLoadRoute);
        active_direction(&Properties.Dir[NumTimeStep * 1],
                         IMIN(NumTimeStep, Properties.Value[1].Num));
        printf(" \t %s (%s) : \n \t \t Number of particles = %i \n \t \t File "
               "curve %s \n",
               "-> BcNeumann ", Parameter_pars[0], NumNodes, FileLoadRoute);
      }
    } else if ((Parser_status == 2) &&
               (strcmp(Parameter_pars[0], "T.z") == 0)) {
      if (strcmp(Parameter_pars[1], "NULL") != 0) {
        sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parameter_pars[1]);
        Properties.Value[Ndim - 1] = ReadCurve(FileLoadRoute);
        active_direction(&Properties.Dir[NumTimeStep * 2],
                         IMIN(NumTimeStep, Properties.Value[2].Num));
        printf("\t \t %s (%s) : \n \t \t Number of particles = %i \n \t \t "
               "File curve %s \n",
               "-> BcNeumann ", Parameter_pars[0], NumNodes, FileLoadRoute);
      }
    } else if ((Parser_status == 1) && (strcmp(Parameter_pars[0], "}") == 0)) {
      Is_Close = true;
      break;
    }
  }

  return Properties;
}

/***************************************************************************/

static void Check_Curve_File(char *PATH_Name) {
  struct stat info;
  stat(PATH_Name, &info);
  char Error_message[MAXW];

  if (S_ISREG(info.st_mode)) {
    printf("\t -> %s : %s \n", "Curve", PATH_Name);
  } else {
    sprintf(Error_message, "\t -> %s : %s %s \n", "List of nodes", PATH_Name,
            "does not exists");
    standard_error(Error_message);
  }
}

/**********************************************************************/

static void standard_error() {
  fprintf(stderr, "%s : %s !!! \n",
          "Error in Read_u_Neumann_Boundary_Conditions__InOutFun__()",
          Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/

static void standard_output(char *Status_message) {
  fprintf(stdout, "%s \n", Status_message);
}

/**********************************************************************/

static FILE *Open_and_Check_simulation_file(char *Name_File) {
  FILE *Simulation_file = fopen(Name_File, "r");

  if (Simulation_file == NULL) {
    sprintf(Error_message, "%s %s", "Incorrect lecture of", Name_File);
    standard_error();
  }

  return Simulation_file;
}

/***************************************************************************/

static void active_direction(int *Dir, int Num) {
  for (unsigned i = 0; i < Num; i++) {
    Dir[i] = 1;
  }
}

/**********************************************************************/
