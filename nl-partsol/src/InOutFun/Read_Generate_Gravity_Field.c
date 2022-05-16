/**
 * @file Read_Generate_Gravity_Field.c
 * @author Miguel Molinos (@migmolper)
 * @brief Read the gravity field of the simulation
 * use the following sintax for constant gravity:
 *
 * generate-gravity-field-constant
 * {
 *  g.x DOUBLE
 *  g.y DOUBLE
 * }
 *
 * use the following sintax for curve-value gravity:
 *
 * generate-gravity-field-curve
 * {
 *  g ROUTE/FILE.csv
 * }
 *
 * @version 0.1
 * @date 2022-05-16
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "nl-partsol.h"
#include <string.h>
#include <sys/stat.h>

typedef struct {
  double *x;
  double *y;
#if NumberDimensions == 3
  double *z;
#endif
  int NumTimeStep;

} Gravity_curve;

/*
  Call global variables
*/
int NumberDOF;

/*
  Auxiliar functions and variables
*/
static bool __gravity_field_definiton(FILE *Sim_dat, bool *g_is_constant,
                                      bool *g_is_curve, int *STATUS);

static int __fill_gravity_constant(FILE *Sim_dat, Gravity_curve G);

static int __fill_gravity_curve(FILE *Sim_dat, Gravity_curve G);

/**********************************************************************/

int Generate_Gravity_Field__InOutFun__(Load *gravity_field,
                                       const char *Name_File,
                                       Time_Int_Params Parameters_Solver) {

  int STATUS = EXIT_SUCCESS;
  unsigned NumTimeStep = Parameters_Solver.NumTimeStep;

  /* Simulation file */
  FILE *Sim_dat;
  Gravity_curve G;
  bool g_is_constant = false;
  bool g_is_curve = false;

  /* Open and check file */
  Sim_dat = fopen(Name_File, "r");
  if (Sim_dat == NULL) {
    fprintf(stderr, "%s : \n\t %s %s", "Error in fopen()",
            "Incorrect lecture of", Name_File);
    return EXIT_FAILURE;
  }

  // Number of dimensions of the BCC
  (*gravity_field).Dim = NumberDOF;

  // Direction of the BCC
  (*gravity_field).Dir = (int *)calloc(NumTimeStep * NumberDOF, sizeof(int));
  if ((*gravity_field).Dir == NULL) {
    fprintf(stderr, "" RED "Error in calloc()" RESET " \n");
    return EXIT_FAILURE;
  }

  // Curve for each dimension
  (*gravity_field).Value = (Curve *)malloc(NumberDOF * sizeof(Curve));
  if ((*gravity_field).Value == NULL) {
    fprintf(stderr, "" RED "Error in malloc" RESET " \n");
    return EXIT_FAILURE;
  }

  // Information of the BCC
  strcpy((*gravity_field).Info, "Gravity field");

  // Read Generate_Gravity_Field__InOutFun__ line
  if (__gravity_field_definiton(Sim_dat, &g_is_constant, &g_is_curve,
                                &STATUS)) {

    // Activate gravity field
    (*gravity_field).STATUS = true;

    // Allocate memory for the gravity curve
    G.x = (double *)calloc(NumTimeStep, __SIZEOF_DOUBLE__);
    G.y = (double *)calloc(NumTimeStep, __SIZEOF_DOUBLE__);
#if NumberDimensions == 3
    G.z = (double *)calloc(NumTimeStep, __SIZEOF_DOUBLE__);
#endif
    G.NumTimeStep = NumTimeStep;

    //
    if (g_is_constant == true) {
      STATUS = __fill_gravity_constant(Sim_dat, G);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __fill_gravity_constant(,)" RESET " \n");
        return EXIT_FAILURE;
      }
    }

    if (g_is_curve == true) {
      STATUS = __fill_gravity_curve(Sim_dat, G);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __fill_gravity_curve(,)" RESET " \n");
        return EXIT_FAILURE;
      }
    }

    (*gravity_field).Value[0].Fx = G.x;
    (*gravity_field).Value[0].Num = NumTimeStep;
    strcpy((*gravity_field).Value[0].Info, "x");
    (*gravity_field).Dir[0] = 1;

    (*gravity_field).Value[1].Fx = G.y;
    (*gravity_field).Value[1].Num = NumTimeStep;
    strcpy((*gravity_field).Value[1].Info, "y");
    (*gravity_field).Dir[1] = 1;

#if NumberDimensions == 3
    (*gravity_field).Value[2].Fx = G.z;
    (*gravity_field).Value[2].Num = NumTimeStep;
    strcpy((*gravity_field).Value[2].Info, "z");
    (*gravity_field).Dir[2] = 1;
#endif
  }

  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __gravity_field_definiton(,)" RESET " \n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**********************************************************************/

static bool __gravity_field_definiton(FILE *Sim_dat, bool *g_is_constant,
                                      bool *g_is_curve, int *STATUS) {
  unsigned Ndim = NumberDimensions;
  int NWL;

  char Line[MAXC] = {0};
  char *Parse_Line[MAXW] = {NULL};
  char *Parse_Nodes[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];

  char File[MAXC];

  while (fgets(Line, sizeof(Line), Sim_dat) != NULL) {

    /* Read the line with the space as separators */
    NWL = parse(Parse_Line, Line, " \n\t");
    if (NWL < 0) {
      fprintf(stderr, "" RED "Error in parse()" RESET " \n");
      *STATUS = EXIT_FAILURE;
      return false;
    }

    if ((NWL > 0) &&
        (strcmp(Parse_Line[0], "generate-gravity-field-constant") == 0)) {
      *g_is_constant = true;
      fprintf(stderr, "\t * generate-gravity-field-constant : " GREEN
                      "True" RESET " \n");
      return true;
    }

    if ((NWL > 0) &&
        (strcmp(Parse_Line[0], "generate-gravity-field-curve") == 0)) {
      *g_is_curve = true;
      fprintf(stderr,
              "\t * generate-gravity-field-curve : " GREEN "True" RESET " \n");
      return true;
    }
  }

  return false;
}

/**********************************************************************/

static int __fill_gravity_constant(FILE *Sim_dat, Gravity_curve G) {

  int NWL;
  char Line_Properties[MAXC] = {0};
  char *Parse_Properties[MAXW] = {NULL};

  bool Is_Open_brace = false;
  bool Is_Close_brace = false;

  double gx_value, gy_value, gz_value;

  /* Read the curve to impose the boundary condition */
  while (fgets(Line_Properties, sizeof(Line_Properties), Sim_dat) != NULL) {

    NWL = parse(Parse_Properties, Line_Properties, "= \n\t");
    if (NWL < 0) {
      fprintf(stderr, "" RED "Error in parse()" RESET " \n");
      return EXIT_FAILURE;
    }

    if ((NWL > 0) && (strcmp(Parse_Properties[0], "{") == 0)) {
      Is_Open_brace = true;
    }

    if ((NWL > 0) && (Is_Open_brace == true) &&
        (strcmp(Parse_Properties[0], "}") == 0)) {
      Is_Close_brace = true;
      break;
    }

    if ((NWL > 0) && (Is_Open_brace == true) &&
        (strcmp(Parse_Properties[0], "g.x") == 0)) {

      gx_value = atof(Parse_Properties[1]);

      for (unsigned i = 0; i < G.NumTimeStep; i++)
        G.x[i] = gx_value;
    }

    if ((NWL > 0) && (Is_Open_brace == true) &&
        (strcmp(Parse_Properties[0], "g.y") == 0)) {

      gy_value = atof(Parse_Properties[1]);

      for (unsigned i = 0; i < G.NumTimeStep; i++)
        G.y[i] = gy_value;
    }

#if NumberDimensions == 3
    if ((NWL > 0) && (Is_Open_brace == true) &&
        (strcmp(Parse_Properties[0], "g.z") == 0)) {

      gz_value = atof(Parse_Properties[1]);

      for (unsigned i = 0; i < G.NumTimeStep; i++)
        G.z[i] = gz_value;
    }
#endif
  }

  if ((Is_Open_brace == false) || (Is_Close_brace == false)) {
    fprintf(stderr, "" RED "Unbalanced braces {...}" RESET " \n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**********************************************************************/

static int __fill_gravity_curve(FILE *Sim_dat, Gravity_curve G) {

  unsigned Ndim = NumberDimensions;
  int NWL;
  char Line_Properties[MAXC] = {0};
  char *Parse_Properties[MAXW] = {NULL};

  FILE *CSV_file;
  char *Aux_Line_CSV;
  char Line_CSV[MAXC] = {0};
  char *Columns_CSV[MAXW] = {NULL};
  int NCOLS;

  bool Is_Open_brace = false;
  bool Is_Close_brace = false;

  char g_file[MAXC];

  struct stat info;

  /* Read the curve to impose the boundary condition */
  while (fgets(Line_Properties, sizeof(Line_Properties), Sim_dat) != NULL) {

    NWL = parse(Parse_Properties, Line_Properties, "= \n\t");
    if (NWL < 0) {
      fprintf(stderr, "" RED "Error in parse()" RESET " \n");
      return EXIT_FAILURE;
    }

    if ((NWL > 0) && (strcmp(Parse_Properties[0], "{") == 0)) {
      Is_Open_brace = true;
    }

    if ((NWL > 0) && (Is_Open_brace == true) &&
        (strcmp(Parse_Properties[0], "}") == 0)) {
      Is_Close_brace = true;
      break;
    }

    if ((NWL > 0) && (Is_Open_brace == true) &&
        (strcmp(Parse_Properties[0], "g") == 0)) {

      strcpy(g_file, Parse_Properties[1]);

      stat(g_file, &info);

      if (S_ISREG(info.st_mode) == 0) {
        fprintf(stderr, "" RED "The file %s does not exists" RESET " \n",
                g_file);
        return EXIT_FAILURE;
      }

      CSV_file = fopen(g_file, "r");

      if (CSV_file == NULL) {
        perror(g_file);
        return EXIT_FAILURE;
      }

      for (unsigned i = 0; i < G.NumTimeStep; i++) {
        Aux_Line_CSV = fgets(Line_CSV, sizeof Line_CSV, CSV_file);

        if (Aux_Line_CSV == NULL) {
          fprintf(stderr, "" RED "Error in file %s during line %i" RESET " \n",
                  g_file, i+1);
          return EXIT_FAILURE;
        }

        NCOLS = parse(Columns_CSV, Line_CSV, ",\n");

        if (NCOLS != Ndim) {
          fprintf(stderr,
                  "" RED "Error in file %s during line %i, NCOLS != Ndim" RESET
                  " \n",
                  g_file, i+1);
          return EXIT_FAILURE;
        }

        G.x[i] = atof(Parse_Properties[0]);
        G.y[i] = atof(Parse_Properties[1]);
#if NumberDimensions == 3
        G.z[i] = atof(Parse_Properties[2]);
#endif
      }
    }
  }

  if ((Is_Open_brace == false) || (Is_Close_brace == false)) {
    fprintf(stderr, "" RED "Unbalanced braces {...}" RESET " \n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**********************************************************************/