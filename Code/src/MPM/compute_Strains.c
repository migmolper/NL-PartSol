#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*************************************************************/

Tensor compute_RateOfStrain(Matrix Velocity, Matrix Gradient)
{
  Tensor Rate_Strain = alloc_Tensor(2);
  Tensor Velocity_I;
  Tensor Gradient_I;
  Tensor VoG_I;

  int NodesElem = Gradient.N_rows;

  /* Compute rate of strain */
  for(int I = 0 ; I<NodesElem ; I++){
    /* Assign from matrix to tensor */
    Velocity_I = memory_to_Tensor(Velocity.nM[I], 1);
    Gradient_I = memory_to_Tensor(Gradient.nM[I], 1);
    /* Compute the dyadic product of the nodal velocity and the
       gradient of the shape functions */
    VoG_I = get_dyadicProduct_Of(Velocity_I, Gradient_I);
    /* Ad the nodal contribution to the train tensor */
    for(int i = 0 ; i<3 ; i++){
      for(int j = 0 ; j<3 ; j++){
	Rate_Strain.N[i][j] +=
	  0.5*(VoG_I.N[i][j] + VoG_I.N[j][i]);
      }
    }
    /* Free memory */
    free_Tensor(VoG_I);
  }
  
  return Rate_Strain;
}

/*******************************************************/

Tensor update_Strain(Tensor Strain, Tensor Rate_Strain, double TimeStep)
{
  /* Check in the input its is ok */
  if ((Strain.Order == 2) && (Rate_Strain.Order == 2)){
    /* Update strain tensor with the rate of strain tensor */
    for(int i = 0 ; i<3 ; i++){
      for(int j = 0 ; j<3 ; j++){
	Strain.N[i][j] += TimeStep*Rate_Strain.N[i][j];
      }
    }
  }
  else{
    fprintf(stderr,"%s : %s %s !!! \n",
	    "Error in update_Strain()",
	    "The input should be",
	    "two tensors of 2nd order and a scalar");
    exit(EXIT_FAILURE);    
  }

  return Strain;
}

/*******************************************************/

/* void UpdateGaussPointStrain(GaussPoint MPM_Mesh, */
/* 			    Mesh FEM_Mesh, */
/* 			    Matrix Mesh_Vel) */
/* /\*! */
/*  * \brief Brief description of UpdateGaussPointStrain. */
/*  *        Update the strain state of the body with the */
/*  *        infinetesimal strain tensor.  */
/*  * */
/*  * \f[  */
/*  \Delta\Epsilon_{ij,p}^{k-1/2} =  */
/*  \frac{\Delta t}{2} \cdot \sum_{I=0}^{Nn}(N^{k}_{Ip,j} */
/*  \cdot v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot v_{jI}^{k-1/2}) */
/*   \f] */
/*  * */
/*  *  The parameters for this functions are  : */
/*  *  @param MPM_Mesh : Mesh with the material points. */
/*  *  @param FEM_Mesh : Background mesh for calculations. */
/*  *  @param Mesh_Vel : Nodal values of the velocity field. */
/*  * */
/*  *\/ */
/* { */
  
/*   /\* Variable definition *\/ */
/*   Element GP_Element; /\* Element for each Gauss-Point *\/ */
/*   Matrix Element_Velocity; /\* Array with the nodal velocities *\/ */
/*   Matrix dNdx_GP; /\*  Matrix with the nodal derivatives *\/ */
/*   Matrix B; /\* B matrix to get the deformation *\/ */
/*   Matrix Delta_Strain_GP; /\* Vectoriced Strain tensor *\/ */
/*   double Delta_TraceStrain_GP; /\* Increment of the trace of the Stress tensor *\/ */
 
/*   for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){ */
   
/*     /\* 1º Define element of the GP *\/ */
/*     GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i], */
/* 			      MPM_Mesh.NumberNodes[i]); */
    
/*     /\* 2º Get the element gradient *\/ */
/*     dNdx_GP = Get_Operator("dNdx", GP_Element, */
/* 			   MPM_Mesh, FEM_Mesh); */
	    
/*     /\* 3º Calcule the B matrix and free the gradient *\/ */
/*     B = Get_B_GP(dNdx_GP); */
/*     FreeMat(dNdx_GP); */

/*     /\* 4º Get the nodal velocities in the element and free connectivity *\/ */
/*     Element_Velocity = GetElementField(Mesh_Vel, GP_Element); */
/*     free(GP_Element.Connectivity); */
   
/*     /\* 5º Multiply B by the velocity array and by the time step to get */
/*        the increment stress tensor *\/ */
/*     Delta_Strain_GP = Scalar_prod(B,Element_Velocity); */
/*     Delta_Strain_GP = Matrix_x_Scalar(Delta_Strain_GP, */
/* 					  DeltaTimeStep); */

/*     /\* 6º Free the array with the nodal velocity of the element and the B matrix *\/ */
/*     FreeMat(Element_Velocity); */
/*     FreeMat(B); */

/*     /\* 7º Udate the Gauss-Point strain tensor *\/ */
/*     Delta_TraceStrain_GP = 0; /\* Set to zero the trace increment *\/ */
/*     for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){ */
/*       MPM_Mesh.Phi.Strain.nM[i][j] += Delta_Strain_GP.nV[j]; */
/*       /\* Get the trace of the stress tensor *\/ */
/*       if(j<NumberDimensions){ */
/* 	Delta_TraceStrain_GP += Delta_Strain_GP.nV[j]; */
/*       } */
/*     } */
/*     FreeMat(Delta_Strain_GP); */

/*     /\* 8º Update the density of the GP *\/ */
/*     MPM_Mesh.Phi.rho.nV[i] = UpdateGaussPointDensity(MPM_Mesh.Phi.rho.nV[i], */
/* 						     Delta_TraceStrain_GP);     */
/*   } */
/* } */

/* /\*******************************************************\/ */
