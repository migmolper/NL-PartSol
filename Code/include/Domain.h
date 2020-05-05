#ifndef _DOMAIN_H_
#define _DOMAIN_H_

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

#endif
