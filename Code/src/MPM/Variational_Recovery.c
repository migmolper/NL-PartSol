#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*********************************************************************/

Matrix GetNodalMassMomentum(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  /* Output */
  Matrix Nodal_FIELDS;
  Matrix N_GP;  /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double GP_mass; /* Mass of the GP */ 
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I;

  /* 1º Allocate the output list of fields */
  Nodal_FIELDS = MatAllocZ(NumberDimensions + 1,FEM_Mesh.NumNodesMesh);
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Define element of the GP */
    GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i],
			      MPM_Mesh.NumberNodes[i]);
    
    /* 4º Evaluate the shape function in the coordinates of the GP */
    N_GP = Get_Operator("N",GP_Element,
			MPM_Mesh,FEM_Mesh);
   
    /* 5º Get the mass of the GP */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];

    /* 6º Get the nodal mass and mommentum */
    for(int k = 0 ; k<GP_Element.NumberNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Element.Connectivity[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(N_GP_I == 0){
	continue;
      }
      /* Nodal mass */
      Nodal_FIELDS.nM[0][GP_I] += GP_mass*N_GP_I;
      /* Nodal momentum */
      for(int l = 0 ; l<NumberDimensions ; l++){
	Nodal_FIELDS.nM[l+1][GP_I] +=
	  GP_mass*MPM_Mesh.Phi.vel.nM[i][l]*N_GP_I;
      }   
    }

    /* 7º Free the value of the shape functions */
    FreeMat(N_GP), free(GP_Element.Connectivity);
  }
 
  return Nodal_FIELDS;
  
}

/*******************************************************/

Matrix GetNodalVelocity(Mesh FEM_Mesh,
			Matrix Nodal_MOMENTUM,
			Matrix Nodal_MASS){
  /*
    Get the nodal velocity using : 
    v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
    Initialize nodal velocities 
  */
  Matrix Vel_Mesh
    = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Vel_Mesh.Info,"VELOCITY");
  
  /* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      if(Nodal_MASS.nV[i] > 0){
	Vel_Mesh.nM[j][i] = (double)Nodal_MOMENTUM.nM[j][i]/Nodal_MASS.nV[i];
      }
    }    
  }
  
  return Vel_Mesh;
}

/*******************************************************/
