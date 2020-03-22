
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
  int * Element_id;

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

/* Variational recovery process */
Matrix GetNodalMassMomentum(GaussPoint, Mesh);
Matrix GetNodalVelocity(Mesh, Matrix, Matrix);
Matrix GetNodalKinetics(GaussPoint, Mesh);
Matrix GetNodalVelocityDisplacement(GaussPoint, Mesh);

/* Boundary conditions */
Curve BcDirichlet(char *);
void BCC_Nod_VALUE(Mesh, Matrix, int);
Matrix Eval_Body_Forces(Load *, int, int, int);
Matrix Eval_Contact_Forces(Load *, int, int, int);

/* Eulerian functions  */
void UpdateGaussPointStrain(GaussPoint, Mesh, Matrix);
double UpdateGaussPointDensity(double, double);
void UpdateGaussPointStress(GaussPoint);
Matrix GetNodalForces(GaussPoint, Mesh, int);
void UpdateGridNodalMomentum(Mesh, Matrix, Matrix);

/* Forward Euler */
void UpdateVelocityAndPositionGP(GaussPoint, Mesh, Matrix, Matrix, Matrix);

/* Generalized-alpha */
void GA_UpdateNodalKinetics(Mesh, Matrix, Matrix, Time_Int_Params);
void GA_AdvectionKinetics(GaussPoint, Mesh, Matrix, Time_Int_Params);


Matrix initialize_NodalVelocity(GaussPoint, Mesh);
Matrix GetNodalMass(GaussPoint, Mesh);
Matrix PredictorNodalVelocity(GaussPoint, Mesh, Matrix,
			      Matrix, Time_Int_Params,double);
Matrix CorrectorNodalVelocity(Mesh, Matrix, Matrix,
			      Matrix, Time_Int_Params,double);
void Update_Lagrangian_PCE(GaussPoint,
			   Mesh, Matrix,
			   Matrix, Matrix,
			   double);

/* Predicto-corrector explicit */
void PCE_Predictor(GaussPoint, Mesh, Matrix, Matrix,
		   Matrix, Time_Int_Params);
void PCE_Corrector(GaussPoint, Mesh, Matrix,
		   Matrix, Time_Int_Params);



Matrix GetInitialGaussPointPosition(Matrix, Mesh, int);
double GetMinElementSize(Mesh);
void GetNodalConnectivity(Mesh);
Matrix ElemCoordinates(Mesh, int *, int);
ChainPtr DiscardElements(ChainPtr, Matrix, Matrix, Mesh);
void LocalSearchGaussPoints(GaussPoint, Mesh);
void UpdateBeps(GaussPoint, Mesh);
void GPinCell(ChainPtr *, ChainPtr *, Matrix, int, double);
Element GetElementGP(int, ChainPtr, int);
Matrix GetElementField(Matrix, Element);
