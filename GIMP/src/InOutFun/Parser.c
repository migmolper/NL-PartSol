#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include "../ToolsLib/TypeDefinitions.h"

/***************************************************************************/

ParserDictionary InitParserDictionary(void)
/*
*/
{

  ParserDictionary Dict;
  
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

  char * KeyWords [14] =
    {"NumberNodes",
     "NumberGaussPoints",
     "ElemType",
     "DegreeFreedom",
     "KindOfAnalisys",
     "U",
     "U_P",
     "Sigma_V",
     "Gravity",
     "Density",
     "Materials",
     "Air",
     "Water",
     "Soil"};

  Dict.ascii_sep = ascii_sep;
  Dict.NumberSeparators = 7;
  Dict.KeyWords = KeyWords;
  Dict.NumberKeyWords = 14;
  
  return Dict;
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

/***************************************************************************/

