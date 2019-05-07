#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/

ParserDictionary InitParserDictionary(void)
/*

*/
{

  char * sep [8] = {"#","=","@","&",";",","," \n","%"};

  char * KeyWords [28] = {"NUM_NODES","NUM_GAUSSPOINTS",
			  "ELEM_TYPE","DOF","RESTART",
			  "TRUE","FALSE","KIND_ANALYSIS",
			  "U","U_P","SIGMA_V","G","RHO",
			  "MATERIAL","AIR","WATER","SOIL",
			  "X_GP","U_X","V_X","SIGMA_X",
			  "TIME_STEP","NUM_STEP","MESH_FILE",
			  "COND_INIT","BOUND_COND","2STG","MASS"};

  ParserDictionary Dict;
  
  Dict.sep = sep;
  Dict.NumberSeparators = 7;
  Dict.KeyWords = KeyWords;
  Dict.NumberKeyWords = 28;

  return Dict;
  
}

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
