
/* Name of inputs files */
char * MeshFileName;
char * InitCondFileName;
char * BounCondFileName;

/* Numerical parameters */
int NumberDimensions;
char * KIND_ANALYSIS[100];

/* Time integration */
double DeltaTimeStep;
int NumTimeStep;

/* Constitutive parameters */
double g;
double ElasticModulus;
double Density;

/* Boundary conditions */
char Field[10];
int * Node;
double * Value;

