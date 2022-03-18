#include <sys/stat.h>
#include <string.h>
#include "nl-partsol.h"

/*
  Call global variables
*/

int ResultsTimeStep;
// int NumTimeStep //

char Output_Backup_File[MAXC];

bool Backup_damage = false;
bool Backup_plastic_deformation_gradient = false;
bool Backup_EPS = false;

/*
  Auxiliar functions
*/
static FILE *Open_and_Check_simulation_file(char *);
static bool Is_Output_Activate(char *, char *);
static void standard_error(char *);
static void read_Output_intervals(char *);
static bool Check_Output_directory(char *);

/**********************************************************************/

void Read_Backup(char *Name_File)
/*
  Example :
  Backup (i=100) {
        DIR=test/Sulsky_MPM
  }
*/
{
  /* Simulation file */
  FILE *Sim_dat;

  char *Error_message;

  /* Temporal integator */
  int Aux_Out_id;
  char *Parse_Out_id[MAXW] = {NULL};

  /* Temporal integrator properties */
  char Line_Out_Prop[MAXC] = {0};
  char *Parse_Out_Prop[MAXW] = {NULL};

  /* Parse file for the route */
  bool Is_OutputDir = false;
  char Route_Outs[MAXC] = {0};

  /* Special variables for line-reading */
  char line[MAXC] = {0};       /* Variable for reading the lines in the files */
  char *kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char *STATUS_LINE;

  /* Open and check file */
  Sim_dat = Open_and_Check_simulation_file(Name_File);

  /* Set output file names to default */
  strcpy(Output_Backup_File, "Backup");

  /* Read the file line by line */
  while (fgets(line, sizeof line, Sim_dat) != NULL) {

    /* Read the line with the space as separators */
    nkwords = parse(kwords, line, " \n\t");
    if (nkwords < 0) {
      standard_error("Parser failed");
    }

    if ((nkwords > 0) && (strcmp(kwords[0], "GramsOutputs") == 0)) {

      /* Read output period */
      read_Output_intervals(kwords[1]);

      /* Look for the curly brace { */
      if (strcmp(kwords[2], "{") == 0) {
        /* Initial line */
        STATUS_LINE = fgets(Line_Out_Prop, sizeof(Line_Out_Prop), Sim_dat);
        if (STATUS_LINE == NULL) {
          standard_error("Unspected EOF");
        }
        Aux_Out_id = parse(Parse_Out_Prop, Line_Out_Prop, " =\t\n");
        if (strcmp(Parse_Out_Prop[0], "}") == 0) {
          /* Check output dir */
          if (!Is_OutputDir) {
            standard_error("Non output dir defined");
          }
          break;
        }
        while (STATUS_LINE != NULL) {

          if (Aux_Out_id != 2) {
            standard_error("Use this format -> Propertie = value");
          }

          if (strcmp(Parse_Out_Prop[0], "DIR") == 0) {
            generate_route(Route_Outs, Name_File);
            sprintf(OutputDir, "%s%s", Route_Outs, Parse_Out_Prop[1]);
            Is_OutputDir = Check_Output_directory(OutputDir);
          } else if (strcmp(Parse_Out_Prop[0], "Particles-file") == 0) {
            strcpy(Output_Backup_File, Parse_Out_Prop[1]);
            printf("\t -> %s : %s \n", "Particles files name",
                   OutputParticlesFile);
          } else if (strcmp(Parse_Out_Prop[0], "Damage") == 0) {
            Backup_damage =
                Is_Output_Activate(Parse_Out_Prop[0], Parse_Out_Prop[1]);
          } else if (strcmp(Parse_Out_Prop[0],
                            "Plastic-Deformation-Gradient") == 0) {
            Backup_plastic_deformation_gradient =
                Is_Output_Activate(Parse_Out_Prop[0], Parse_Out_Prop[1]);
          } else if (strcmp(Parse_Out_Prop[0], "Equivalent-Plastic-Strain") ==
                     0) {
            Backup_EPS =
                Is_Output_Activate(Parse_Out_Prop[0], Parse_Out_Prop[1]);
          } else {
            sprintf(Error_message, "%s %s", "Undefined", Parse_Out_Prop[0]);
            standard_error(Error_message);
          }
          /* Read next line and check */
          STATUS_LINE = fgets(Line_Out_Prop, sizeof(Line_Out_Prop), Sim_dat);
          Aux_Out_id = parse(Parse_Out_Prop, Line_Out_Prop, " =\t\n");
          if (strcmp(Parse_Out_Prop[0], "}") == 0) {
            break;
          }
        }
        if (STATUS_LINE == NULL) {
          standard_error("you forget to put a }");
        }

        if (strcmp(Parse_Out_Prop[0], "}") == 0) {
          /* Check output dir */
          if (!Is_OutputDir) {
            standard_error("No output dir was defined");
          }
          break;
        }
      } else {
        standard_error("Use this format -> GramsOutputs (Type=string) { ");
      }
    }
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);
}

/***************************************************************************/

static FILE *Open_and_Check_simulation_file(char *Name_File) {
  FILE *Simulation_file = fopen(Name_File, "r");
  char *Error_message;

  if (Simulation_file == NULL) {
    sprintf(Error_message, "%s %s", "Incorrect lecture of", Name_File);
    standard_error(Error_message);
  }

  return Simulation_file;
}

/***************************************************************************/

static bool Is_Output_Activate(char *output_field, char *status_text) {
  bool status;
  char *Error_message;

  if (strcmp(status_text, "true") == 0) {
    printf("\t -> %s : true \n", output_field);
    return true;
  } else if (strcmp(status_text, "false") == 0) {
    printf("\t -> %s : False \n", output_field);
    return false;
  } else {
    sprintf(Error_message, "The input was %s. Please, use : true/false",
            status_text);
    standard_error(Error_message);
  }

  return status;
}

/***************************************************************************/

static void standard_error(char *Error_message) {
  fprintf(stderr, "%s : %s !!! \n", "Error in GramsOutputs()", Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/

static void read_Output_intervals(char *Interval_message) {
  int Interval_status;
  char *Parse_message[MAXW] = {NULL};

  Interval_status = parse(Parse_message, Interval_message, "(=)");

  /* Check format */
  if ((Interval_status != 2) || (strcmp(Parse_message[0], "i") != 0)) {
    standard_error("Use this format -> (i=int)");
  }

  /* Get interval output and store in a global variable */
  ResultsTimeStep = atoi(Parse_message[1]);
  //    if(ResultsTimeStep > NumTimeStep)
  //      {
  //      	standard_error("The result time step should be less than the
  //      total time steps");
  //      }

  /* Print in screen some information */
  printf("\t -> %s : %i \n", "Output values each", ResultsTimeStep);
}

/***************************************************************************/

static bool Check_Output_directory(char *Output_directory) {
  struct stat info;
  stat(Output_directory, &info);
  if (S_ISDIR(info.st_mode)) {
    printf("\t -> %s : %s \n", "Output directory", OutputDir);
    return true;
  } else {
    printf("\t -> %s : %s %s \n", "Output directory", OutputDir,
           "does not exists");
    return false;
  }
}

/***************************************************************************/