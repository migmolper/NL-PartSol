#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif



/* Dictionay definition : */
ParserDictionary Dict;

/* Functions definitions : */
void InitParserDictionary(void);
int GetWords(char *, char * [], int,int);
void ReadGidMesh(char *);
Matrix Read_CSV(char *,int);
/* void ReadData(char *); */

