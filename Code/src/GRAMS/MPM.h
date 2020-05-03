
/*! \file MPM.h
 * Definitions of structures
 */


#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/*! \struct Curve
 *  Curve definition
 */
typedef struct{

  /*! Number of items in the curve */
  int Num; 
  /*! Values for each time */
  double * Fx; 
  /*! Aditional information */
  char Info [100];
  
} Curve;

/*******************************************************/

/*! \struct Load
 *   An important aclaration about this code, a force is a load,
 * or even a boundary condition is defined as a load from the
 * point of view of the code (as a structure). The idea is
 * a force and a boundary condition both of them has a direction, 
 * a value, a list of nodes/GPs where it is applied,... So it's
 * a good idea to use this structure for both porpouses.
 */
typedef struct {

  /*! Number of nodes/GP with this load */
  int NumNodes;
  /*! Number of dimensions of the load */
  int Dim;
  /*! Direction of the load {0,0} {1,0} {0,1} {1,1} */
  int * Dir;
  /*! List of nodes with this load */
  int * Nodes;
  /*! Curve for each dimension with the evolution
     value with the time */
  Curve * Value;
  /*! Some information about this load */
  char Info [100];

} Load;

/*******************************************************/

/*! \struct Boundaries
 * Boundary conditions definition 
 */
typedef struct {

  /*! Number of boundaries of the domain */
  int NumBounds;
  /*! Table with all the boundaries and its values */
  Load * BCC_i;
  
} Boundaries;

/*******************************************************/

/*! \struct Fields
 * Physical fields 
*/
typedef struct {

  /*! Density field */
  Matrix rho;
  /*! Mass field */
  Matrix mass;
  /*! Position in global coordinates */
  Matrix x_GC;  
  /*! Position in element coordiantes */
  Matrix x_EC;
  /*! Displacement field */
  Matrix dis;
  /*! Velocity field */
  Matrix vel;
  /*! Acceleration field */
  Matrix acc;
  /*! Stress field */
  Matrix Stress;
  /*! Strain field */
  Matrix Strain;
  /*! Strain during crack */
  Matrix StrainF;
  /*! Deformation Energy */
  Matrix W;
  /*! Damage parameter (Fracture) */
  Matrix ji;
  
} Fields;

/*******************************************************/

/*! \struct GaussPoint
 * Gauss points definitions
 */
typedef struct {

  /*! Total number of the GP */
  int NumGP;
  /*! Identification of the element where it is */
  int * I0;

  /*! Tributary nodes variables */
  int * NumberNodes;
  ChainPtr * ListNodes;

  /*! GPs near to each GP */
  ChainPtr * Beps;

  /*! Values from the actual step */
  Fields Phi;
  /*! Values from the previous step */
  Fields Phi_n0; 
  
  /*! Number of materials */
  int NumberMaterials;
  /*! Index of the material for each GP */
  int * MatIdx;
  /*! Array wit the number of materials */
  Material * Mat;

  /*! Neumann boundary conditions */
  int NumNeumannBC;
  Load * F;

  /*! Body forces */
  int NumberBodyForces;
  Load * B; 

  /*! Shape functions variables */
  Matrix lp; /*! Voxel shape  */
  Matrix lambda; /*! Lagrange multiplier */
  Matrix Beta; /*! Norm parameter */

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
  /* Number of elements that shares a node and list 
     of elements that share a node */
  int * NumNeighbour;
  ChainPtr * NodeNeighbour;

  /* Number of particles close to this node */
  int * NumParticles;
  /* List of particles in a node */
  ChainPtr * I_particles;

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

/*! \struct Element
 *  Individual mesh properties of a GPs
 */
typedef struct{

  /*! Index of the GP */ 
  int i_GP;
  /*! Number of nodes near to the GP */
  int NumberNodes;
  /*! Nodal connectivity of the GP */ 
  int * Connectivity;
  
} Element;

/*******************************************************/

/*! \struct Time_Int_Params
 */
typedef struct {

  /*! Generalized alpha parameters */
  double GA_alpha;
  double GA_beta;
  double GA_gamma;
  
} Time_Int_Params;

/*******************************************************/

/* Boundary conditions */
Curve BcDirichlet(char *);
void imposse_NodalMomentum(Mesh, Matrix, int);
void imposse_NodalVelocity(Mesh, Matrix, int);
Matrix Eval_Body_Forces(Load *, int, int, int);
Matrix Eval_Contact_Forces(Load *, int, int, int);
Matrix compute_Reactions(Mesh, Matrix);

/* Forces-stress-strain functions */
void ComputeDamage(int, GaussPoint, Mesh);
Tensor compute_RateOfStrain(Matrix, Matrix);
Tensor update_Strain(Tensor, Tensor, double);
double update_Density(double, double, Tensor);
Tensor compute_Stress(Tensor, Tensor, Material);
void update_LocalState(Matrix, GaussPoint, Mesh, double);
double compute_InternalEnergy(Tensor, Tensor);
Matrix compute_InternalForces(Matrix, GaussPoint,Mesh);
Matrix compute_BodyForces(Matrix, GaussPoint, Mesh, int);
Matrix compute_ContacForces(Matrix, GaussPoint, Mesh, int);

/* Forward-Euler */
Matrix compute_NodalMomentumMass(GaussPoint, Mesh);
Matrix compute_NodalVelocity(Mesh, Matrix);
void update_NodalMomentum(Mesh, Matrix, Matrix);
void update_Particles_FE(GaussPoint, Mesh, Matrix, Matrix, double);

/* Explicit Predictor-Corrector */
Matrix compute_NodalMass(GaussPoint, Mesh);
Matrix compute_VelocityPredictor(GaussPoint, Mesh, Matrix,
				 Matrix, Time_Int_Params,double);
Matrix compute_VelocityCorrector(Mesh, Matrix, Matrix,
				 Matrix, Time_Int_Params,double);
void update_Particles_PCE(GaussPoint, Mesh, Matrix,
			  Matrix, Matrix, double);

/* Generalized-alpha */
void GA_UpdateNodalKinetics(Mesh, Matrix, Matrix, Time_Int_Params);
Matrix GetNodalKinetics(GaussPoint, Mesh);
Matrix GetNodalVelocityDisplacement(GaussPoint, Mesh);
void update_Particles_GA(GaussPoint, Mesh, Matrix, Time_Int_Params);


/* Mesh utilities */
Matrix GetInitialGaussPointPosition(Mesh, int);
double GetMinElementSize(Mesh);
void GetNodalConnectivity(Mesh);

ChainPtr DiscardElements(ChainPtr, Matrix, Matrix, Mesh);

void LocalSearchGaussPoints(GaussPoint, Mesh);

void ComputeBeps(int, GaussPoint, Mesh);

void GPinCell(ChainPtr *, ChainPtr *, Matrix, int, double);

Element get_Element(int, ChainPtr, int);
Matrix get_set_Coordinates(ChainPtr, Matrix, Matrix);
Matrix get_Element_Field(Matrix, Element);
ChainPtr get_locality_of_node(int, Mesh);
int get_closest_node_to(Matrix, ChainPtr, Matrix);
bool InOut_Element(Matrix, ChainPtr, Matrix);
int search_particle_in(int, Matrix, ChainPtr, Mesh);
Matrix ElemCoordinates(ChainPtr, Matrix);
void asign_particle_to_nodes(int, ChainPtr, Mesh);
