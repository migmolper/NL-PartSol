#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include "../GRAMS/TypeDefinitions.h"

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
  
  for (p = strtok (str, delims); p; p = strtok (NULL, delims)) 
    {
      words[n++] = strdup (p);    /* allocate/copy */
      
      if (n == MAXW) { /* limit reached - realloc/break */
	fprintf (stderr, "warning: MAXW reached.\n");
            break;
      }
    }

  return n;
}
