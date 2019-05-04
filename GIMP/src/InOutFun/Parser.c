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

  char * KeyWords [27] = {"NUM_NODES","NUM_GAUSSPOINTS",
     "ELEM_TYPE","DOF","RESTART","TRUE","FALSE","KIND_ANALYSIS",
     "U","U_P","SIGMA_V","G","RHO","MATERIAL","AIR","WATER",
     "SOIL","X_GP","U_X","V_X","SIGMA_X","TIME_STEP","NUM_STEP",
     "MESH_FILE","COND_INIT","BOUND_COND","2STG"};

  ParserDictionary Dict;
  
  Dict.sep = sep;
  Dict.NumberSeparators = 7;
  Dict.KeyWords = KeyWords;
  Dict.NumberKeyWords = 27;

  return Dict;
  
}

/***************************************************************************/


/* int GetWords(char * line, char * words []  , int separator, int maxwords) */
/* /\* */
/*  *\/   */
/* {   */
/*   char *p = line; */
/*   int nwords = 0; */

/*   while(1){ */
    
/*     while( (int)*p == separator ){ */
/*       p++; */
/*     } */
    
/*     if( *p == '\0' ) return nwords; */
    
/*     words[nwords++] = p; */
      
/*     while( ( (int)*p != separator ) && */
/* 	   ( *p != '\0')  ){ */
/*       p++; */
/*     } */
      
/*     if(*p == '\0') return nwords; */
      
/*     *p++ = '\0'; */
      
/*     if(nwords >= maxwords) */
/*       return nwords; */
/*   } */
  
  
/* } */

/***************************************************************************/


int parse (char **words, char *str, char * delims)
{
  int n = 0;
  char *p;
  
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
