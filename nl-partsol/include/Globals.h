
#define MAXW 100
#define MAXC 1000
#define NumberDimensions 2
#define TOL_InOut 10E-23
#define TOL_NR 10E-6
#define TOL_zero 10E-10

extern char * SimulationFile;
extern char * RestartFile;
extern char * FEM_MeshFileName;
extern char * MPM_MeshFileName;
extern char * TimeIntegrationScheme;
extern char * ShapeFunctionGP;
extern char * Formulation;
extern char   OutputDir[MAXC];
extern char   Field[10];

extern double gamma_LME;
extern double TOL_lambda;
extern double CFL; /* Courant number (0-1) */
extern double DeltaTimeStep;
extern double SpectralRadius;


/* Convergence parameters for radia returning algorithm */
extern double TOL_Radial_Returning;
extern int Max_Iterations_Radial_Returning;


extern int NumTimeStep;
extern int ResultsTimeStep;
extern int NumberDOF;



/* int * Node; */
/* double * Value; */
