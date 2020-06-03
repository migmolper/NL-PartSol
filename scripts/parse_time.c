#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>

#define MAXW 100
#define MAXC 1000

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


void main()
{
  char Name_file_t[80];
  char Name_file_copy[80];
  char * Parse_Line_1[MAXW] = {NULL};
  char * Parse_Line_2[MAXW] = {NULL};
  int n_words;
  
  sprintf(Name_file_t,"%s/MPM_%s_%i.vtk","Resultados","MPM_VALUES",50);
  strcpy(Name_file_copy,Name_file_t);

  printf("antes %s \n",Name_file_t);

  n_words = parse(Parse_Line_1, Name_file_copy, "_\n");

  n_words = parse(Parse_Line_2, Parse_Line_1[n_words-1], ".\n");

  printf("%i , %s, %i \n",n_words,Parse_Line_2[0],atoi(Parse_Line_2[0]));

  printf("despues %s \n",Name_file_t);
}
