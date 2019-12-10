
#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Loads definition :
 An important aclaration about this code, a force is a load,
 or even a boundary condition is defined as a load from the
 point of view of the code (as a structure). The idea is
 a force and a boundary condition both of them has a direction, 
 a value, a list of nodes/GPs where it is applied,... So it's
 a good idea to use this structure for both porpouses.
*/
typedef struct {

  /* Number of nodes/GP with this load */
  int NumNodes;
  /* Number of dimensions of the load */
  int Dim;
  /* Direction of the load {0,0} {1,0} {0,1} {1,1} */
  int * Dir;
  /* List of nodes with this load */
  int * Nodes;
  /* Curve for each dimension with the evolution
     value with the time */
  Curve * Value;
  /* Some information about this load */
  char Info [100];

} Load;

/*******************************************************/

/* Boundary conditions definition */
typedef struct {

  /* Number of boundaries of the domain */
  int NumBounds;
  /* Table with all the boundaries and its values */
  Load * BCC_i;
  /* Some information about the boundaries */
  char Info [100];
  
} Boundaries;


/*******************************************************/

/* Load case : multiple loads */
typedef struct {

  /* Number of nodes/GP with this load */
  int NumLoads;
  /* List of loads */
  Load * Load_i;
  /* Some information about this load */
  char Info [100];

} LoadCase;


/*******************************************************/

/* Physical fields */
typedef struct {

  /* Position in global coordinates */
  Matrix x_GC;  
  /* Position in element coordiantes */
  Matrix x_EC;
  /* Voxel shape */
  Matrix lp;
  /* Density field */
  Matrix rho;
  /* Mass field */
  Matrix mass;
  /* Displacement field */
  Matrix dis;
  /* Velocity field */
  Matrix vel;
  /* Acceleration field */
  Matrix acc;
  /* Stress field */
  Matrix Stress;
  /* Strain field */
  Matrix Strain;
  
} Fields;

/*******************************************************/

/* Gauss points definitions */
typedef struct {

  /* Total number of the GP */
  int NumGP;
  /* Identification of the element where it is */
  int * Element_id;
  /* Kind of shape function for the GP mesh */
  char ShapeFunctionGP [20];
  /* Tributary nodes variables */
  int * NumberNodes;
  ChainPtr * ListNodes;
  /* List of Fields */
  Fields Phi; /* Values from the actual step */
  Fields Phi_n0; /* Values from the previous step */
  /* Constitutive response */
  int * MatIdx;
  ConstLib D;
  /* External forces */
  LoadCase F;
  /* Body forces */
  LoadCase B;

} GaussPoint;

/*******************************************************/

/* Mesh definition */
typedef struct {

  /*** GENERAL MESH PROPERTIES ***/
  /* Number of nodes in the mesh */
  int NumNodesMesh;
  /* Number of elements in the mesh */
  int NumElemMesh;
  /* Table with the coordinates of the nodes of the mesh */
  /* double ** Coordinates; */
  Matrix Coordinates;
  /* List of nodes for each element (Connectivity) */
  int * NumNodesElem;
  ChainPtr * Connectivity;
  /* Active node : 
     Boolean variable that set the node ative (1>) or not (0) */
  int * ActiveNode;
  /* Number of elements that shares a node and list 
     of elements that share a node */
  int * NumNeighbour;
  ChainPtr * NodeNeighbour;

  /*** BOUNDARIES ***/
  Boundaries Bounds;
  
  /*** ELEMENT PROPERTIES ***/
  /* Number of dimensions of the element */
  int Dimension;
  /* Size of the element (COURANT) */
  double DeltaX;
  /* Name of the element */
  char TypeElem [20];
  /* Shape function for the mesh, for interpolation porpouses */
  Matrix (* N_ref)(Matrix );
  /* Derivative shape function, for interpolation porpouses */
  Matrix (* dNdX_ref)(Matrix );
    
} Mesh;

/*******************************************************/

/* Define and initialize both mesh */
GaussPoint Define_GP_Mesh(char *, double);
GaussPoint InitializeGP(char *, Mesh);
Mesh InitializeMesh(char *);

/* int SearchGaussPoint(int, Matrix, Matrix, Mesh); */
Matrix GetNodalMassMomentum(GaussPoint, Mesh);
Matrix GetNodalVelocity(Mesh, Matrix, Matrix);

void BCC_Nod_VALUE(Mesh, Matrix, int);
void UpdateGaussPointStrain(GaussPoint, Mesh, Matrix);
double UpdateGaussPointDensity(double, double);
void UpdateGaussPointStress(GaussPoint);
Matrix GetNodalForces(GaussPoint, Mesh, int);
void UpdateGridNodalMomentum(Mesh, Matrix, Matrix);
void UpdateVelocityAndPositionGP(GaussPoint, Mesh,
				 Matrix, Matrix, Matrix);

void GetNodalConnectivity(Mesh);
double GetMinElementSize(Mesh);
void GlobalSearchGaussPoints(GaussPoint, Mesh);
void LocalSearchGaussPoints(GaussPoint, Mesh);
Matrix Get_B_GP(Matrix);
