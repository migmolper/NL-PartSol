#include "nl-partsol.h"
#include <ctype.h>

/*
  Local structures
*/
typedef struct {
  int N_Rows;
  int N_Cols;
  char *Parser_string;

} Param_Data;

/*
  Auxiliar functions and variables
*/
#ifdef _WIN32
static char *delimiters_1 = " \r\n";
#else
static char *delimiters_1 = " \n";
#endif
static char *delimiters_2 = "=";

static char Error_message[MAXW];

static Param_Data Read_header(FILE *);
static Matrix Read_data(Param_Data, FILE *);
static int Read_NROWS_param(char *);
static int Read_NCOLS_param(char *);
static char *Read_PARSER_param(char *, int);
static void Check_header_parameters(bool, bool, bool);
static void standard_error();
static FILE *Open_and_Check_simulation_file(char *);

/**************************************************************/

Matrix Read_Delimited_File__InOutLib__(char *Name_File)
/*
        # NROWS=4000 NCOLS=2 PARSER=%d,%d,%d,%d
*/
{
  Matrix Data;

  FILE *Data_File = Open_and_Check_simulation_file(Name_File);

  Param_Data Header_Parameters = Read_header(Data_File);

  Data = Read_data(Header_Parameters, Data_File);

  fclose(Data_File);

  return Data;
}

/**************************************************************/

static Param_Data Read_header(FILE *Data_File) {
  Param_Data Parameters;

  /* Special variables for line-reading */
  char line[MAXC] = {0};       /* Variable for reading the lines in the files */
  char *kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */
  char *parameters[MAXW] = {NULL};
  int nparameters;

  bool Is_NROWS = false;
  bool Is_NCOLS = false;
  bool Is_PARSER = false;

  /* Read the file line by line */
  while (fgets(line, sizeof line, Data_File) != NULL) {
    /* Read the line with the space as separators */
    nkwords = parse(kwords, line, delimiters_1);

    /* Read header */
    if (strcmp(kwords[0], "#") == 0) {
      for (int i = 1; i < nkwords; i++) {
        nparameters = parse(parameters, kwords[i], delimiters_2);

        if (strcmp(parameters[0], "NROWS") == 0) {
          Parameters.N_Rows = Read_NROWS_param(parameters[1]);
          Is_NROWS = true;
        }

        if (strcmp(parameters[0], "NCOLS") == 0) {
          Parameters.N_Cols = Read_NCOLS_param(parameters[1]);
          Is_NCOLS = true;
        }

        if (strcmp(parameters[0], "PARSER") == 0) {
          Parameters.Parser_string =
              Read_PARSER_param(parameters[1], Parameters.N_Cols);
          Is_PARSER = true;
        }
      }

      Check_header_parameters(Is_NROWS, Is_NCOLS, Is_PARSER);

      /* Stop reading	*/
      break;
    }
  }

  return Parameters;
}

/**************************************************************/

static Matrix Read_data(Param_Data Header_Parameters, FILE *Data_File) {

  int N_Rows = Header_Parameters.N_Rows;
  int N_Cols = Header_Parameters.N_Cols;

  Matrix Data = allocZ__MatrixLib__(N_Rows, N_Cols);

  /* Special variables for line-reading */
  char *STATUS_line = {NULL};
  char line[MAXC] = {0};       /* Variable for reading the lines in the files */
  char *kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  for (int i = 0; i < N_Rows; i++) {
    STATUS_line = fgets(line, sizeof line, Data_File);

    if (STATUS_line == NULL) {
      sprintf(Error_message, "%s", "Incorrect number of rows");
      standard_error();
    }

    nkwords = parse(kwords, line, " ,\n");

    if (nkwords != N_Cols) {
      sprintf(Error_message, "%s %i %s", "The line", i + 1,
              "does not match with the parser");
      standard_error();
    }

    for (int j = 0; j < N_Cols; j++) {
      Data.nM[i][j] = atof(kwords[j]);
    }
  }

  return Data;
}

/**************************************************************/

static int Read_NROWS_param(char *NROWS_string) {
  int NROWS = atoi(NROWS_string);

  if (NROWS < 1) {
    sprintf(Error_message, "%s", "Number of rows is < 1");
    standard_error();
  }

  return NROWS;
}

/**************************************************************/

static int Read_NCOLS_param(char *NCOLS_string) {
  int NCOLS = atoi(NCOLS_string);

  if (NCOLS < 1) {
    sprintf(Error_message, "%s", "Number of columns is < 1");
    standard_error();
  }

  return NCOLS;
}

/**************************************************************/

static char *Read_PARSER_param(char *PARSER_string, int N_Cols) {
  char *PARSER = (char *)malloc((strlen(PARSER_string)) * sizeof(char));

  strcpy(PARSER, PARSER_string);

  //  	nkwords = parse (kwords, line, Header_Parameters.Parser_string);

  return PARSER;
}

/**************************************************************/

static void Check_header_parameters(bool Is_NROWS, bool Is_NCOLS,
                                    bool Is_PARSER) {
  if (!Is_NROWS || !Is_NCOLS || !Is_PARSER) {
    fprintf(stderr, "%s\n",
            "Missing parameters in the header for Read_Delimited_File()");
    fputs(Is_NROWS ? "NROWS : true \n" : "NROWS : false \n", stdout);
    fputs(Is_NCOLS ? "NCOLS : true \n" : "NCOLS : false \n", stdout);
    fputs(Is_PARSER ? "PARSER : true \n" : "PARSER : false \n", stdout);
    exit(EXIT_FAILURE);
  }
}

/**************************************************************/

static void standard_error() {
  fprintf(stderr, "%s : %s !!! \n", "Error in Read_Delimited_File()",
          Error_message);
  exit(EXIT_FAILURE);
}

/**************************************************************/

static FILE *Open_and_Check_simulation_file(char *Name_File) {
  FILE *Simulation_file = fopen(Name_File, "r");
  char Error_message[MAXW];

  if (Simulation_file == NULL) {
    sprintf(Error_message, "%s %s", "Incorrect lecture of", Name_File);
    standard_error();
  }

  return Simulation_file;
}

/**************************************************************/
