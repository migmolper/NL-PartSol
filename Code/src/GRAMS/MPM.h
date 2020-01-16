
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

/* Curve definition */
typedef struct{

  /* Number of items in the curve */
  int Num; 
  /* Values for each time */
  double * Fx; 
  /* Aditional information */
  char Info [100];
  
} Curve;

/*******************************************************/

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

  /* Density field */
  Matrix rho;
  /* Mass field */
  Matrix mass;
  /* Position in global coordinates */
  Matrix x_GC;  
  /* Position in element coordiantes */
  Matrix x_EC;
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
  /* Deformation Energy */
  Matrix W;
  /* Damage parameter (Fracture) */
  Matrix ji;
  
} Fields;

/*******************************************************/

/* Gauss points definitions */
typedef struct {

  /* Total number of the GP */
  int NumGP;
  /* Identification of the element where it is */
  int * Element_id;

  /* Tributary nodes variables */
  int * NumberNodes;
  ChainPtr * ListNodes;

  /* GPs near to each GP */
  ChainPtr * Beps;

  /* List of Fields */
  Fields Phi; /* Values from the actual step */
  Fields Phi_n0; /* Values from the previous step */
  
  /* Constitutive response */
  int * MatIdx; /* Index of the material for each GP */
  Material * Mat; /* Array wit the number of materials */

  /* Forces over the GP */
  LoadCase F; /* External forces */
  LoadCase B; /* Body forces */

  /* Shape functions variables */
  /* GIMP */
  Matrix lp; /* Voxel shape  */
  /* LME */
  Matrix lambda; /* Lagrange multiplier */
  Matrix Beta; /* Norm parameter */

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

  /* List of GPs in a element */
  ChainPtr * GPsElements;

  /*** BOUNDARIES ***/
  Boundaries Bounds;
  
  /*** ELEMENT PROPERTIES ***/
  /* Number of dimensions of the element */
  int Dimension;
  /* Size of the element */
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
GaussPoint Define_GP_Mesh(char *, char *);
GaussPoint InitializeGP(char *, Mesh);
Mesh InitializeMesh(char *);

/* int SearchGaussPoint(int, Matrix, Matrix, Mesh); */
Matrix GetNodalMassMomentum(GaussPoint, Mesh);
Matrix GetNodalVelocity(Mesh, Matrix, Matrix);

/* Boundary conditions */
Curve BcDirichlet(char *);
void BCC_Nod_VALUE(Mesh, Matrix, int);
Matrix Eval_Body_Forces(LoadCase, int, int);
Matrix Eval_Contact_Forces(LoadCase, int, int);

/* MPM functions  */
void UpdateGaussPointStrain(GaussPoint, Mesh, Matrix);
double UpdateGaussPointDensity(double, double);
void UpdateGaussPointStress(GaussPoint);
Matrix GetNodalForces(GaussPoint, Mesh, int);
void UpdateGridNodalMomentum(Mesh, Matrix, Matrix);
void UpdateVelocityAndPositionGP(GaussPoint, Mesh,
				 Matrix, Matrix, Matrix);

Boundaries GetBoundaryBox(Mesh);
double GetMinElementSize(Mesh);
void GetNodalConnectivity(Mesh);
Matrix ElemCoordinates(Mesh, int *, int);
/* void GlobalSearchGaussPoints(GaussPoint, Mesh); */
ChainPtr DiscardElements(ChainPtr, Matrix, Matrix, Mesh);
void LocalSearchGaussPoints(GaussPoint, Mesh);
void UpdateBeps(GaussPoint, Mesh);
void GPinCell(ChainPtr *, ChainPtr *, Matrix, int, double);
