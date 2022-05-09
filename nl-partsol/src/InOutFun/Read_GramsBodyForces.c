#include "nl-partsol.h"
#include <string.h>

typedef struct {
  double *x;
  double *y;
#if NumberDimensions == 3
  double *z;
#endif

} Gravity_curve;

/*
  Call global variables
*/
int NumberDOF;

/*
  Auxiliar functions and variables
*/
static int __search_for_gravity_field_definiton(Gravity_curve *G, FILE *Sim_dat,
                                                unsigned NumTimeStep);

/**********************************************************************/

int GramsBodyForces(Load *gravity_field, const char *Name_File,
                    Time_Int_Params Parameters_Solver) {

  int STATUS = EXIT_SUCCESS;
  unsigned NumTimeStep = Parameters_Solver.NumTimeStep;

  /* Simulation file */
  FILE *Sim_dat;
  Gravity_curve G;

  /* Open and check file */
  Sim_dat = fopen(Name_File, "r");
  if (Sim_dat == NULL) {
    fprintf(stderr, "%s : \n\t %s %s", "Error in GramsBodyForces()",
            "Incorrect lecture of", Name_File);
    return EXIT_FAILURE;
  }

  /* Generate route */
  generate_route(Route_Nodes, Name_File);

  // Number of dimensions of the BCC
  (*gravity_field).Dim = NumberDOF;

  // Direction of the BCC
  (*gravity_field).Dir = (int *)calloc(NumTimeStep * NumberDOF, sizeof(int));

  // Curve for each dimension
  (*gravity_field).Value = (Curve *)Allocate_Array(NumberDOF, sizeof(Curve));

  // Information of the BCC
  strcpy((*gravity_field).Info, "Gravity field");

  // Read GramsBodyForces line
  STATUS = __search_for_gravity_field_definiton(&G, Sim_dat);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in __search_for_gravity_field_definiton(,)" RESET
            " \n");
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

  return STATUS;
}

/**********************************************************************/

static int __search_for_gravity_field_definiton(Gravity_curve *G, FILE *Sim_dat,
                                                unsigned NumTimeStep) {
  unsigned Ndim = NumberDimensions;
  int NWL;
  bool g_is_constant = false;
  bool g_is_curve = false;
  bool x_component = false;
  bool y_component = false;
#if NumberDimensions == 3
  bool z_component = false;
#endif

  bool Is_Key_word = false;
  bool Is_Open_braces = false;
  bool Is_Close_braces = false;

  char Line[MAXC] = {0};
  char *Parse_Line[MAXW] = {NULL};
  char *Parse_Nodes[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];

  char File[MAXC];

  /* Parse generate-gravity-field properties */
  char Line_Properties[MAXC] = {0};
  char *Parse_Properties[MAXW] = {NULL};

  while (fgets(Line, sizeof(Line), Sim_dat) != NULL) {

    /* Read the line with the space as separators */
    NWL = parse(Parse_Line, Line, " \n\t");
    if (NWL < 0) {
      fprintf(stderr, "%s : %s \n", "Error in generate-gravity-field ()",
              "Parser failed");
      return EXIT_FAILURE;
    }

    if (strcmp(Parse_Line[0], "generate-gravity-field-constant")) {

      g_is_constant = true;

      NWL = parse(Parse_Nodes, Parse_Line[1], "(,)");

    } else if (strcmp(Parse_Line[0], "generate-gravity-field-curve")) {

      g_is_curve = true;

    }

    /* Find generate-gravity-field line */
    if ((NWL > 0) && (strcmp(Parse_Line[0], "generate-gravity-field") == 0) &&
        ((strcmp(Parse_Line[2], "{") == 0))) {

      /* File with the nodes */

      // Allocate memory for the gravity curve
      G->x = (double *)calloc(NumTimeStep, __SIZEOF_DOUBLE__);
      G->y = (double *)calloc(NumTimeStep, __SIZEOF_DOUBLE__);
#if NumberDimensions == 3
      G->z = (double *)calloc(NumTimeStep, __SIZEOF_DOUBLE__);
#endif

      /* Read file with the nodes */
      sprintf(FileNodesRoute, "%s%s", Route_Nodes, Parse_Nodes[1]);

      /* Read the curve to impose the boundary condition */
      while (fgets(Line_Properties, sizeof(Line_Properties), Sim_dat) != NULL) {

        NWL = parse(Parse_Properties, Line_Properties, " \n\t");
        if (NWL < 0) {
          fprintf(stderr, "%s : %s \n", "Error in generate-gravity-field ()",
                  "Parser failed");
          return EXIT_FAILURE;
        }
        if ((NWL > 0) && (strcmp(Parse_Properties[0], "}") != 0)) {

          if (strcmp(Parse_Properties[0], "b.x") == 0) {
            x_component = true;

            if (g_is_constant) {

            } else if (g_is_curve) {
              if (strcmp(Parse_Properties[1], "NULL") != 0) {
                sprintf(FileLoadRoute, "%s%s", Route_Nodes,
                        Parse_Properties[1]);
              }
            }

          }

          else if (strcmp(Parse_Properties[0], "b.y") == 0) {
            y_component = true;

            if (g_is_constant) {

            } else if (g_is_curve) {
              if (strcmp(Parse_Properties[1], "NULL") != 0) {
                sprintf(FileLoadRoute, "%s%s", Route_Nodes,
                        Parse_Properties[1]);
              }
            }

          }
#if NumberDimensions == 3
          else if (strcmp(Parse_Properties[0], "b.z") == 0) {
            z_component = true;

            if (g_is_constant) {

            } else if (g_is_curve) {
            }

            if (strcmp(Parse_Properties[1], "NULL") != 0) {
              sprintf(FileLoadRoute, "%s%s", Route_Nodes, Parse_Properties[1]);
            } else {
              (*gravity_field).Dir[2] = 0;
            }
          }
#endif
          else {
            fprintf(stderr, "%s : %s \n", "Error in generate-gravity-field()",
                    "Use -> b.x, b.y, b.z");
            return EXIT_FAILURE;
          }
        } else {
          break;
        }
      }
    }
  }
}

/**********************************************************************/
