#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>



int getwords(char *line, char *words[], int separator ,int maxwords)
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

int main(void){

  /* Send to global variables */

  
  
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

  


  FILE * DatFile = fopen("ejemplo.dat","r");
  if( DatFile != NULL ){

    char line[80];
    int nwords;
    char * words[20];

    while ( fgets(line, sizeof line, DatFile) != NULL ){
      nwords = getwords(line, words,ascii_white,20);
      for(int i  = 0 ; i<nwords ; i++){

	
	printf("%s ", words[i]);

	
      }
      printf("\n");      
    } 
    fclose(DatFile);
    
  }
  else{
    perror("ejemplo.dat");
  }
  
  return 0;
}
