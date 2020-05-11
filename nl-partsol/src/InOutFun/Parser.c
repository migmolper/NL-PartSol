#include "nl-partsol.h"

/***************************************************************************/

int parse(char ** words, char * str, char * delims)
/*
  Parser function :
  Inputs :
  - str -> Input string
  - delims -> parser symbol o list of symbols
  Outputs :
  - words -> Output string with the words parsed
  - n -> Number of words in str
*/
{
  /* Set to zero the number of words */
  int n = 0;
  /* string pointer */
  char * p;
  
  for (p = strtok (str, delims); p; p = strtok (NULL, delims))  {

    /* allocate/copy */
    words[n++] = strdup (p);
    
    /* limit reached - realloc/break */    
    if (n == MAXW) { 
      fprintf (stderr, "warning: MAXW reached.\n");
      break;
    }
    
  }
  
  return n;
}

/***************************************************************************/


void generate_route(char * Route_Nodes, char * Name_File)
/*
  Use the route of the .gdf file to generate a route to other files
*/
{

  char * Name_File_cpy = malloc(strlen(Name_File));
  int Num_words_route;
  char * Name_Parse[MAXW] = {NULL};

  /* Copy the file name to avoid lost data */
  strcpy(Name_File_cpy, Name_File);  

  #ifdef linux
  Num_words_route = parse(Name_Parse,Name_File_cpy,"(/)") - 1;
  strcat(Route_Nodes,"./");
  for(int i = 0 ; i<Num_words_route ; i++){
    strcat(Route_Nodes, Name_Parse[i]);
    strcat(Route_Nodes,"/");
  }
  #endif

  #ifdef _WIN32 
  strcpy(Name_File_cpy, Name_File);
  Num_words_route = parse(Name_Parse,Name_File_cpy,"\\");
  for(int i = 0 ; i<Num_words_route-1 ; i++){
     strcat(Route_Nodes, Name_Parse[i]);
	 strcat(Route_Nodes, "\\");
   }
  #endif 

  /* Free auxiliar table */
  free(Name_File_cpy);

}
