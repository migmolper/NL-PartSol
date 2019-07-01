/***************************************/
/******** Name of inputs files *********/
/***************************************/

char * FEM_MeshFileName;
char * MPM_MeshFileName;
char * InitCondFileName;
char * BCC_FEM_FileName;
char * OutputDir;

/***************************************/
/*********** Find of analysis **********/
/***************************************/

char * KindAnalysis;
char * FieldsAnalysis;
int NumberDimensions;
int NumberDOF;
char * TimeIntegration;

/***************************************/
/********* Time integration ************/
/***************************************/

double DeltaTimeStep;
int NumTimeStep;

/***************************************/
/******* Constitutive parameters *******/
/***************************************/
double ElasticModulus;
double PoissonModulus;
double Density;

/***************************************/
/******** Boundary conditions **********/
/***************************************/

char Field[10];
int * Node;
double * Value;


