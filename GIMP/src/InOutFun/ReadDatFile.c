
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/


void ReadDatFile(char * Name_File)
/*
  Read data from the .DAT file and initialize the variables
  
  Inputs
  - Name_file : Name of the file
  - 
  Outputs
  - Conectivity matrix
  - Coordenates of the nodes
  - phi_n : Values of the variables in the n step, 
        initializerd with the Initial conditions
  - DeltaT : Time-step
  - A_el : Area of the element
  - type_elem : Type of the element (1)
  - 
  - N_nodes : Number of nodes
  - N_elem : Number of elements
  - N_steps : Number of time steps
*/
{

  /* Simulation file */
  FILE * Sim_dat;

  /* Char */
  char line[80];

  /* Number of element in the line , just for check */
  int nwords;
  int aux;
  char * words[20];
  int Element_i,Nodes_i;

  

  
  /* Array with the possible index fields */
  int AuxiliarArrayFields[14];

  /* Variable for the index in nodes and conectivitie */
  int ix;

  /* Auxiliar variable for reading the lines in the files */
  char line[200];
  
  /* Allocate a string with enough space for the extensions. */
  char *Name_Simulation = malloc(strlen(Name_File+5));
  
  /* Copy the name with extensions into fn. */
  sprintf(Name_Simulation, "%s.dat", Name_File); 
  
  /* Open and check .dat file */
  Sim_dat = fopen(Name_Simulation,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .dat file");
    exit(0);
  }

  
    
  // Close .dat file
  free(Name_Simulation); // Free the memory gained from malloc. 
  fclose(Sim_dat);
  

} /* void read_dat(char * Name_File) */



/***************************************************************************/
