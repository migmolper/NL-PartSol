#include "grams.h"

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
