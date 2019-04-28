
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* In ReadGidMesh.c : */
void ReadGidMesh(char *,
		 ElementProperties *,
		 MeshProperties *,
		 ParserDictionary);

/* /\* In ReadDatFile.c : *\/ */
/* void ReadData(char *); */

/* In Parser.c : */
ParserDictionary InitParserDictionary(void);

int GetWords(char *, char * [], int,int);
