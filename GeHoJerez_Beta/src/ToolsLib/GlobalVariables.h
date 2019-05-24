
/* Name of inputs files */
char * FEM_MeshFileName;
char * MPM_MeshFileName;
char * InitCondFileName;
char * BounCondFileName;

/* Find of analysis */
char * KindAnalysis;
char * FieldsAnalysis;
int NumberDimensions;
char * TimeIntegration;

/* Time integration */
double DeltaTimeStep;
int NumTimeStep;

/* Constitutive parameters */
double g;
double ElasticModulus;
double PoissonModulus;
double Density;

/* Boundary conditions */
char Field[10];
int * Node;
double * Value;


