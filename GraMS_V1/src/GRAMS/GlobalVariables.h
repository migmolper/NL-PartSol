/***************************************/
/******** Name of inputs files *********/
/***************************************/

char * FEM_MeshFileName;
char * MPM_MeshFileName;
char * OutputDir;

/***************************************/
/*********** Find of analysis **********/
/***************************************/

char * ShapeFunctionGP;
char * Formulation;
int NumberDimensions;
int NumberDOF;
char * TimeIntegration;

/***************************************/
/********* Time integration ************/
/***************************************/

double DeltaTimeStep;
int NumTimeStep;
int ResultsTimeStep;

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


