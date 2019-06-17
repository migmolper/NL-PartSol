#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Word Parser */
enum { MAXW = 100, MAXC = 1000 };

/* Auxiliar functions */
ParserDictionary InitParserDictionary(void);
Curve ReadCurve(char *);
int parse(char **words, char *str, char *delims);

/* Inputs */
Mesh ReadGidMesh(char *);
Matrix Read_CSV(char *, int);
void ReadDatFile(char *);
void ReadBCC(char * , Mesh);

/* Outputs */
void WriteGnuplot(Matrix, Matrix, double, double, int, int, char[20]);
void WriteVtk_MPM(char * , GaussPoint, Matrix, int);
void WriteVtk_FEM(char * , Mesh, Matrix, int);
