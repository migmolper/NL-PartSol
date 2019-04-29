
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* In Parser.c : */
ParserDictionary InitParserDictionary(void);

int GetWords(char *, char * [], int,int);


/* In ReadGidMesh.c : */
void ReadGidMesh(char *,
		 ElementProperties *,
		 MeshProperties *,
		 ParserDictionary *);

/* /\* In ReadDatFile.c : *\/ */
/* void ReadData(char *); */

