#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

void UpdateVelocityAndPositionGP(GaussPoint MPM_Mesh,
				 Mesh FEM_Mesh,
				 Matrix Nodal_MASS,
				 Matrix Nodal_MOMENTUM,
				 Matrix Nodal_TOT_FORCES){

  Matrix N_GP; /* Value of the shape-function in the GP */
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I; /* Index of each tributary node for the GP */
  double mass_I; /* Value of the nodal mass */
  double N_I_GP; /* Nodal value for the GP */

  /* 1º iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Define element of the GP */
    GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i],
			      MPM_Mesh.NumberNodes[i]);

    /* 3º Evaluate shape function in the GP i */
    N_GP = Get_Operator("N",GP_Element,
			MPM_Mesh,FEM_Mesh);
    
    /* 4º Iterate over the nodes of the element */
    for(int j = 0; j<GP_Element.NumberNodes; j++){
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if(fabs(N_I_GP) <= TOL_zero){
	continue;
      }
      /* Get the nodal mass */
      mass_I = Nodal_MASS.nV[GP_I];
      /* Update GP cuantities with nodal values */
      for(int k = 0 ; k<NumberDimensions ; k++){
	/* Update the GP velocities */
	MPM_Mesh.Phi.vel.nM[i][k] +=
	  DeltaTimeStep*N_I_GP*
	  Nodal_TOT_FORCES.nM[k][GP_I]/mass_I;	
	/* Update the GP position */
	MPM_Mesh.Phi.x_GC.nM[i][k] +=
	  DeltaTimeStep*N_I_GP*
	  Nodal_MOMENTUM.nM[k][GP_I]/mass_I;
      } 
    }
    
    /* 5º Free memory */
    free(GP_Element.Connectivity), FreeMat(N_GP);
  }  
}


/*******************************************************/
