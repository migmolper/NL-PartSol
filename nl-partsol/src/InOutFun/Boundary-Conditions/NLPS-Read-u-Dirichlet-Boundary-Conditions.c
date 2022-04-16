#include <sys/stat.h>
#include <string.h>
#include "nl-partsol.h"

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

Boundaries Read_u_Dirichlet_Boundary_Conditions__InOutFun__(char *Name_File,
                                                            int NumBounds,
                                                            int NumTimeStep)
/*

GramsBoundary (File=Right_contour.txt)
{
        BcDirichlet V.x CurveConstant.txt
        BcDirichlet V.y NULL
        BcDirichlet V.z NULL
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
          Read GramsBoundary line
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
            Find GramsBoundary line
    */
    if ((Num_words_line > 0) &&
        (strcmp(Parse_GramsBoundary[0], "GramsBoundary") == 0)) {

      /*
              File with the nodes
      */
      Num_words_line = parse(Parse_Nodes, Parse_GramsBoundary[1], delimiters_3);

      /*
              Check the sintax of the parameters
      */
      if ((Num_words_line != 2) || (strcmp(Parse_Nodes[0], "File") != 0)) {
        sprintf(Error_message,
                "Sintax error in line %i : GramsBoundary(File=Nodes.txt)",
                Num_line);
        standard_error();
      }

      /*
      Read file with the nodes and get an array with the nodes
*/
      sprintf(FileNodesRoute, "%s%s", Route_Nodes, Parse_Nodes[1]);
      Chain_Nodes = File2Chain(FileNodesRoute);
      Bounds.BCC_i[IndexBoundary].NumNodes = lenght__SetLib__(Chain_Nodes);
      Bounds.BCC_i[IndexBoundary].Nodes = set_to_memory__SetLib__(
          Chain_Nodes, Bounds.BCC_i[IndexBoundary].NumNodes);
      free__SetLib__(&Chain_Nodes);

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
      strcpy(Bounds.BCC_i[IndexBoundary].Info, "Displacements");

      /*
          Increment the index
      */
      IndexBoundary++;
    }
  }

  return Bounds;
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

    if ((strcmp(Parameter_pars[0], "{") == 0) && (Parser_status == 1)) {
      Is_Open = true;
    } else if ((strcmp(Parameter_pars[0], "BcDirichlet") == 0) &&
               (Parser_status == 3)) {
      if (strcmp(Parameter_pars[1], "V.x") == 0) {
        if (strcmp(Parameter_pars[2], "NULL") != 0) {
          sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parameter_pars[2]);
          Properties.Value[0] = ReadCurve(FileLoadRoute);
          active_direction(&Properties.Dir[NumTimeStep * 0],
                           IMIN(NumTimeStep, Properties.Value[0].Num));
          printf(" \t %s (%s) : \n \t \t Number of nodes = %i \n \t \t File "
                 "curve %s \n",
                 "-> BcDirichlet ", Parameter_pars[1], NumNodes, FileLoadRoute);
        }
      } else if (strcmp(Parameter_pars[1], "V.y") == 0) {
        if (strcmp(Parameter_pars[2], "NULL") != 0) {
          sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parameter_pars[2]);
          Properties.Value[1] = ReadCurve(FileLoadRoute);
          active_direction(&Properties.Dir[NumTimeStep * 1],
                           IMIN(NumTimeStep, Properties.Value[1].Num));
          printf(" \t %s (%s) : \n \t \t Number of nodes = %i \n \t \t File "
                 "curve %s \n",
                 "-> BcDirichlet ", Parameter_pars[1], NumNodes, FileLoadRoute);
        }
      } else if (strcmp(Parameter_pars[1], "V.z") == 0) {
        if (strcmp(Parameter_pars[2], "NULL") != 0) {
          sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parameter_pars[2]);
          Properties.Value[2] = ReadCurve(FileLoadRoute);
          active_direction(&Properties.Dir[NumTimeStep * 2],
                           IMIN(NumTimeStep, Properties.Value[2].Num));
          printf("\t \t %s (%s) : \n \t \t Number of nodes = %i \n \t \t File "
                 "curve %s \n",
                 "-> BcDirichlet ", Parameter_pars[1], NumNodes, FileLoadRoute);
        }
      } else {
        sprintf(Error_message,
                "Sintax error in line %i : BcDirichlet ->%s<- %s", Num_line,
                Parameter_pars[2], Parameter_pars[3]);
        standard_error();
      }
    } else if ((strcmp(Parameter_pars[0], "}") == 0) && (Parser_status == 1)) {
      Is_Close = true;
      break;
    } else if (Parser_status > 0) {
      sprintf(Error_message, "Sintax error in line %i : Undefined parameter %s",
              Num_line, Parameter_pars[0]);
      standard_error();
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
  fprintf(stderr, "%s : %s !!! \n", "Error in Define-Material()",
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
