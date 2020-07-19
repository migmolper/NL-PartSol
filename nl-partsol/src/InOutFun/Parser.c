#include "nl-partsol.h"

/***************************************************************************/

int parse(char ** words, char * str, char * delims)
{
  /*!
    Set to zero the number of words 
  */
  int n = 0;
  
  /*!
    string pointer 
  */
  char * p;
  
  for (p = strtok (str, delims); p; p = strtok (NULL, delims))
    {

      /*!
	allocate/copy 
      */
      words[n++] = strdup (p);
    
      /*!
	limit reached - realloc/break
      */    
      if (n == MAXW)
	{ 
	  fprintf (stderr, "warning: MAXW reached.\n");
	  break;
	}
    
    }
  
  return n;
}

/***************************************************************************/


void generate_route(char * Route_Nodes, char * Name_File)
{

  char * Name_File_cpy = (char *)malloc((strlen(Name_File) + 1)*sizeof(char));
  int Num_words_route;
  char * Name_Parse[MAXW] = {NULL};

  /*!
    Copy the file name to avoid lost data 
  */
  strcpy(Name_File_cpy, Name_File);  

#ifdef linux
  Num_words_route = parse(Name_Parse,Name_File_cpy,"(/)") - 1;
  strcat(Route_Nodes,"./");
  for(int i = 0 ; i<Num_words_route ; i++)
    {
      strcat(Route_Nodes, Name_Parse[i]);
      strcat(Route_Nodes,"/");
    }
#endif

#ifdef _WIN32 
  strcpy(Name_File_cpy, Name_File);
  Num_words_route = parse(Name_Parse,Name_File_cpy,"\\");
  for(int i = 0 ; i<Num_words_route-1 ; i++)
    {
      strcat(Route_Nodes, Name_Parse[i]);
      strcat(Route_Nodes, "\\");
    }
#endif 

  /*!
    Free auxiliar table 
  */
  free(Name_File_cpy);

}

/***************************************************************************/

int get_ResultStep(char * File_Result)
{
  char File_Result_Copy[80];
  char * Parse_1[MAXW] = {NULL};
  char * Parse_2[MAXW] = {NULL};
  FILE * Ptr_File_Result;
  char Line[MAXC] = {0};
  char * kwords[MAXW] = {NULL};
  int n_words;
  int step;
  
  strcpy(File_Result_Copy,File_Result);

#ifdef linux
  n_words = parse(Parse_1, File_Result_Copy, "_\n");
  n_words = parse(Parse_2, Parse_1[n_words-1], ".\n");
#endif

#ifdef _WIN32
  n_words = parse(Parse_1, File_Result_Copy, "_\r\n");
  n_words = parse(Parse_2, Parse_1[n_words-1], ".\r\n");
#endif
 
  step = atoi(Parse_2[0]);

  /* Open and check file */
  Ptr_File_Result = fopen(File_Result,"r");  
  if (Ptr_File_Result == NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in get_ResultStep()",
	    "Incorrect lecture of",File_Result);
    exit(EXIT_FAILURE);
  }

  /* Read the first line */
  fgets(Line, sizeof Line, Ptr_File_Result);
  /* Read the second line */
  fgets(Line, sizeof Line, Ptr_File_Result);
  
  /* Read the line with the space as separators */
  n_words = parse (kwords, Line," \n\t");

  if ((n_words == 4) &&
      (strcmp(kwords[0],"Results") == 0) &&
      (strcmp(kwords[1],"time") == 0) &&
      (strcmp(kwords[2],"step") == 0))
    {
      step += atoi(kwords[3]);
    }

  /* Close file */
  fclose(Ptr_File_Result);

  return step;
}

/***************************************************************************/
