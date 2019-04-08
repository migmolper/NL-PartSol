#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************** Mesh parameters ******************/


/* Number of spatial-dimensions
   1D->1 ; 2D->2 ; 3D->3
*/
int SpatialDimensions;

/* Number of nodes of the mesh*/
int NumberNodes;

/* Array with the nodes coordinates, it is allocated
   with a size of NumberNodes*SpatialDimensions  */
double * NodesCoordinates;

/* Number of elements in the mesh */
int NumberElements;

/* Type of element adopted for the simulation 
   2-nodes-linear->1

   Gives you the number of nodes in a element 
*/
int KindOfElement;

/* Array with the conectivity of the mesh, it is allocated
   with a size of NumberElements*NumberNodesPerElement */
int * MeshConectivity;

/*****************************************************/



/**************** Physical parameters ****************/

/* Number of fields in the simulation */
int NumberFields;

/* Table of pointer with the fields,
   the size depends of the number of fields and the number of nodes 
   
   Fields = [ Field 1 ; Field 2 ; ... ; Field N | ... | Field 1 ; Field 2 ; ... ; Field N ]
             <----------- Element 1 ----------->       <----------- Element N ----------->
   
   Possible fields, each one are allocated with the 
   number of elements, all of them are called from Fields :
   U_x, U_y, U_z,
   V_x, V_y, V_z, 
   Sigma_x, Sigma_y, Sigma_z, 
   Taub_xy, Taub_xz, Taub_yz,
   P_w, P_a
*/
double * Fields;

/* Gravity constant */
double g;


/*****************************************************/



/************* Simulation parameters *****************/


/* Number of time-steps of the simulation */
int NumberTimeSteps;


/* Present time-step in the simulation
   minval : 0
   maxval : NumberTimeSteps - 1 */
int ActualTimeStep;

/* Time increment between the steps*/
double DeltaTime;

/* Temporal integration scheme 
   0->Forward-Euler
   1->Predictor-corrector (Two-steps)
   2->Runge-Kutta 4
   3->Runge-Kuta 45
 */
int TemporalIntegrationSheme;

/* Name of the .dat file (Name of the simulation) */
char NameFileDat[80];

/* Represent the .dat file */
FILE * SimulationDat;

/*****************************************************/




/******* General Structure of the element ************/

struct Element {

  /* Global index of the element in the mesh, its depends in 
     general of the mesh generator that we are using  */
  int Index;
  
  /* Nodes of the element, it is a pointer to the MeshConectivity
   table */
  int * Node;

  /* As all the elements are joined between them using a chain,
   it is necesary to define where it is placed in memory the next 
   and the previous element, it is done with a pointer */
  struct Element * NextElement;
  struct Element * PreviousElement;

  /* Table of pointer where each component points to a function 
   where the basis function is defined, the half of the table is
  filled with the basis function for each node, and the other 
  half is filled with its derivatives */
  double(** BasisFunction)(double);


}


/*****************************************************/
