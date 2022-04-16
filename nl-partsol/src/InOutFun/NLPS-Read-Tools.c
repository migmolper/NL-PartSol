#include <sys/stat.h>
#include <string.h>
#include "nl-partsol.h"

static char *vector_delimiter = " {;}";

/***************************************************************************/

int Parse_string__InOutFun__(char **words, char *str, char *delims,
                             int *STATUS) {
  char *delims_cpy = (char *)malloc((strlen(delims)) * sizeof(char));
  strcpy(delims_cpy, delims);

  // Set to zero the number of words
  int n = 0;

  // string pointer
  char *p;

  for (p = strtok(str, delims_cpy); p; p = strtok(NULL, delims_cpy)) {
    // allocate/copy
    words[n++] = strdup(p);

    // limit reached - realloc/break
    if (n == MAXW) {
      fprintf(stderr, "Error in : %s. %s \n", "Parse_string__InOutFun__",
              "MAXW reached");
      (*STATUS) = 1;
      return n;
    }
  }

  free(delims_cpy);

  if (n < 0) {
    fprintf(stderr, "Error in : %s. %s \n", "Parse_string__InOutFun__",
            "Parser failed");
    (*STATUS) = 1;
  }

  return n;
}

/***************************************************************************/

FILE *Open_and_Check_File__InOutFun__(char *Name_File, int *STATUS) {
  FILE *Simulation_file = NULL;
  struct stat info;
  stat(Name_File, &info);

  Simulation_file = fopen(Name_File, "r");

  if ((S_ISREG(info.st_mode) == false) || (Simulation_file == NULL)) {
    fprintf(stderr, "%s %s : %s \n", "Open_and_Check_File__InOutFun__",
            "Incorrect lecture of file", Name_File);
    (*STATUS) = 1;
  }

  return Simulation_file;
}

/***************************************************************************/

ChainPtr File_to_Chain__InOutFun__(char *Name_File, int *STATUS) {
  // Output
  ChainPtr File_Chain = NULL;

  // Open and check file
  int INFO_OpenFile = 0;

  FILE *Sim_dat = Open_and_Check_File__InOutFun__(Name_File, &INFO_OpenFile);

  if (INFO_OpenFile) {
    fprintf(stderr, "%s %s : %s \n", "File_to_Chain__InOutFun__",
            "Incorrect lecture of file", Name_File);
    (*STATUS) = 1;

    return File_Chain;
  }

  // Read file linea by line
  int Numers_Line;
  char Line_File[MAXC] = {0};
  char *Parse_File[MAXW] = {NULL};

  while (fgets(Line_File, sizeof(Line_File), Sim_dat) != NULL) {
    Numers_Line = parse(Parse_File, Line_File, " \n\t");
    push__SetLib__(&File_Chain, atoi(Parse_File[0]));
  }

  // Close file
  fclose(Sim_dat);

  return File_Chain;
}

/***************************************************************************/

Tensor Read_Vector__InOutFun__(char *String_vector, int *STATUS)
/*!
 * Vector = {x, y}
 * */
{
  int Ndim = NumberDimensions;

  Tensor Output = alloc__TensorLib__(1);

  /* Variables for reading purposes */
  char *Parser_vector[MAXW] = {NULL};
  int Dimensions_vector;

  Dimensions_vector = parse(Parser_vector, String_vector, vector_delimiter);

  if (Dimensions_vector != Ndim) {
    fprintf(stderr, "%s : %s \n", "Error in", "Read_Vector__InOutFun__");
    (*STATUS) = 1;
    return Output;
  }

  for (int i = 0; i < Ndim; i++) {
    Output.n[i] = atof(Parser_vector[i]);
  }

  return Output;
}

/***************************************************************************/
