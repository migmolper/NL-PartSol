#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Word Parser */
enum { MAXW = 100, MAXC = 1000 };

/* Functions definitions : */
ParserDictionary InitParserDictionary(void);
/* int GetWords(char *, char * [], int, int); */
int parse (char **words, char *str, char *delims);
Mesh ReadGidMesh(char *);
Matrix Read_CSV(char *,int);
void ReadDatFile(char *);
void ReadBCC(char *);
void GnuplotOutput1D(Matrix, Matrix, double, double, int, int, char[20]);

/* Outputs */
void WriteVtk_MPM(char *,GaussPoint,Matrix,int);
