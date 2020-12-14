#include "nl-partsol.h"
#include <ctype.h>

/*
  Local structures
*/
typedef struct
{
  int   N_Rows;
  int   N_Cols;
  char  Parser_string[MAXC];

} Param_Data;

/*
  Auxiliar functions and variables
*/
static char Error_message[MAXW];

static void   standard_error();
static FILE * Open_and_Check_simulation_file(char *);

/**************************************************************/

Matrix Read_Delimited_File(char * Name_File)
/*
	NROWS=4000 NCOLS=2 PARSER=%d,%d,%d,%d
*/
{
	Matrix Data;

	// FILE * Data_File = Open_and_Check_simulation_file(Name_File);

	return Data;
}

/**************************************************************/

// Param_Data Read_header(FILE * Data_File)
// {
// 	Param_Data Parameters;

// 	/* Special variables for line-reading */
//   	char line[MAXC] = {0}; /* Variable for reading the lines in the files */
//   	char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
//   	int nkwords; /* Number of element in the line , just for check */

// 	 Read the file line by line 
//   	while(fgets(line, sizeof line, Sim_dat) != NULL )
//   	{
//   		/* Read the line with the space as separators */
//     	nkwords = parse (kwords, line," \n\t");

//   		/* Read header */
//     	if (strcmp(kwords[0],"#") == 0 )
//     	{
//     		for(int i = 0 ; i<nkwords ; i++)
//     		{
//     			if (strcmp(kwords[0],"#") == 0 )
//     		}

//     	}
//   	}

// 	Parameters;
// }


/**************************************************************/

// static void standard_error()
// {
//   fprintf(stderr,"%s : %s !!! \n",
//      "Error in Read_Delimited_File()",Error_message);
//     exit(EXIT_FAILURE);
// }

/**************************************************************/

// static FILE * Open_and_Check_simulation_file(char * Name_File)
// {
//   FILE * Simulation_file = fopen(Name_File,"r");  
//   char Error_message[MAXW];
  
//   if (Simulation_file==NULL)
//   {
//   	sprintf(Error_message,"%s %s","Incorrect lecture of",Name_File);
// 	  standard_error(); 
//   }  

//   return Simulation_file;
// }

/**************************************************************/
