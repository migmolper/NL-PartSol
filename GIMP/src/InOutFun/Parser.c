#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/

void InitParserDictionary(void)
/*
*/
{
 
  char *hashtag = "#";
  int ascii_hashtag = (int)*hashtag;
  char *equal = "=";
  int ascii_equal = (int)*equal;  
  char *at = "@";
  int ascii_at = (int)*at;
  char *ampersand = "&";
  int ascii_ampersand = (int)*ampersand;  
  char *semicolon = ";";
  int ascii_semicolon = (int)*semicolon;  
  char *comma = ",";
  int ascii_comma = (int)*comma; 
  char *white = " ";
  int ascii_white = (int)*white;

  int ascii_sep [7] =
    {ascii_hashtag,
     ascii_equal,
     ascii_at,
     ascii_ampersand,
     ascii_semicolon,
     ascii_comma,
     ascii_white};

  char * KeyWords [21] =
    {"NUM_NODES",
     "NUM_GAUSSPOINTS",
     "ELEM_TYPE",
     "DOF",
     "RESTART",
     "TRUE",
     "FALSE",  
     "KIND_ANALYSIS",
     "U",
     "U_P",
     "SIGMA_V",
     "G",
     "RHO",
     "MATERIAL",
     "AIR",
     "WATER",
     "SOIL",
     "X_GP",
     "U_X",
     "V_X",
     "SIGMA_X",
     "TIME_STEP",
     "NUM_STEP",
     "MESH_FILE",
     "COND_INIT",
     "BOUND_COND",
     "2STG"};

  Dict.ascii_sep = ascii_sep;
  Dict.NumberSeparators = 7;
  Dict.KeyWords = KeyWords;
  Dict.NumberKeyWords = 21;
  
}

/***************************************************************************/


int GetWords(char *line, char *words[], int separator ,int maxwords)
/*
*/
{  
 
  char *p = line;
  int nwords = 0;    
 

  while(1)
    {
      while( isspace(*p) ){
	p++;
      }
	  
      if(*p == '\0') return nwords;
	  
      words[nwords++] = p;
	  
      while( ((int)*p != separator ) &&
	     (*p != '\0') ){
	p++;
      }
	  
      if(*p == '\0') return nwords;
	  
      *p++ = '\0';
	  
      if(nwords >= maxwords)
	return nwords;
    }
}

