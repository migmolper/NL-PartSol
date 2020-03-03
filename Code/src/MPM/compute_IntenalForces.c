#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*************************************************************/

Matrix compute_InternalForces(Matrix F_I, Matrix V_I,
			      GaussPoint MPM_Mesh,
			      Mesh FEM_Mesh, double TimeStep)
{  
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix Gradient_p; /* Shape functions gradients */
  Matrix Nodal_Velocity_p; /* Velocity of the element nodes */
  Material Material_p; /* Properties of the Gauss-Point material */
  Tensor Rate_Strain_p; /* Increment of strain tensor */
  Tensor Strain_p; /*  Strain tensor */
  Tensor Stress_p; /* Stress tensor */
  Tensor Gradient_pI;
  Tensor InternalForcesDensity_Ip;
  double W_p; /* Internal energy of the Gauss-Point */
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V_p; /* Volumen of the Gauss-Point */
  int Np = MPM_Mesh.NumGP;
  int Ip;
  int Nn;
  int Idx_Mat_p;

  /* Loop in the GPs */
  for(int p = 0 ; p<Np ; p++){

    /* Get the value of the density */
    rho_p = MPM_Mesh.Phi.rho.nV[p];

    /* Get the value of the mass */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Asign memory to tensors */
    Strain_p = memory_to_Tensor(MPM_Mesh.Phi.Strain.nM[p], 2);
    Stress_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nodes_p =
      get_Element(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* Get the velocity of the nodes of the element */
    Nodal_Velocity_p = get_Element_velocity(Nodes_p, V_I);

    /* Compute gradient of the shape function in each node */
    Gradient_p =
      compute_ShapeFunction_Gradient(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the material properties */
    Idx_Mat_p = MPM_Mesh.MatIdx[p];
    Material_p = MPM_Mesh.Mat[Idx_Mat_p];

    /* Update Strain tensor */
    Rate_Strain_p = compute_RateOfStrain(Nodal_Velocity_p,Gradient_p);
    Strain_p = update_Strain(Strain_p, Rate_Strain_p, TimeStep);

    /* Update density field */
    rho_p = update_Density(rho_p, TimeStep, Rate_Strain_p);

    /* Compute stress tensor */
    Stress_p = compute_Stress(Strain_p,Stress_p,Material_p);

    /* Compute deformation energy */
    W_p = 0.5*get_innerProduct_Of(Strain_p, Stress_p);

    /* Compute the volume of the Gauss-Point */
    V_p = m_p/rho_p;

    /* Compute nodal forces */
    for(int I = 0 ; I<Nn ; I++){
      /* Pass by reference the nodal gradient to the tensor */
      Gradient_pI = memory_to_Tensor(Gradient_p.nM[I], 1);
      
      /* Compute the nodal forces of the Gauss-Point */
      InternalForcesDensity_Ip =
	get_firstOrderContraction_Of(Stress_p, Gradient_pI);
      
      /* Get the node of the mesh for the contribution */
      Ip = Nodes_p.Connectivity[I];
      
      /* Asign the nodal forces contribution to the node */
      for(int i = 0 ; i<3 ; i++){
	F_I.nM[Ip][i] -= InternalForcesDensity_Ip.n[i]*V_p;
      }

      /* Free the internal forces density */
      free_Tensor(InternalForcesDensity_Ip);
    }

    /* Update memory */
    MPM_Mesh.Phi.rho.nV[p] = rho_p;
    MPM_Mesh.Phi.W.nV[p] = W_p;
    
    
    /* Free the matrix with the nodal velocity of the element */
    FreeMat(Nodal_Velocity_p);
    
    /* Free the matrix with the nodal gradient of the element */
    FreeMat(Gradient_p);
    
  }

  return F_I;
  
}

/*******************************************************/

/* Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep) */
/* /\*! */
/*  * \brief Brief description of GetNodalForces. */
/*  *        Compute the nodal contribution of each GP to the total forces. */
/*  * */
/*  *  The parameters for this functions are  : */
/*  *  @param MPM_Mesh : Mesh with the material points. */
/*  *  @param FEM_Mesh : Background mesh for calculations. */
/*  *  @param TimeStep : Parameter to evaluate the loads curve. */
/*  * */
/*  *\/ */
/* { */
/*   double GP_mass; /\* Mass of the GP *\/ */
/*   double Vol_GP; /\* Gauss-Point volumen *\/ */
/*   double ji_GP; /\* Damage parameter *\/ */
/*   Matrix N_GP; /\* Matrix with the nodal shape functions *\/ */
/*   double N_GP_I; /\* Evaluation of the GP in the node *\/ */
/*   Matrix dNdx_GP; /\* Matrix with the nodal derivatives *\/ */
/*   Matrix N_dNdx_GP; /\* Operator Matrix *\/ */
/*   Matrix B, B_T; /\* B matrix *\/ */
/*   Element GP_Element; /\* Element for each Gauss-Point *\/ */
/*   int GP_I; /\* Node of the GP *\/ */
/*   Matrix Stress_GP = /\* Stress tensor of a Gauss-Point *\/ */
/*     MatAssign(MPM_Mesh.Phi.Stress.N_cols,1,NAN,NULL,NULL); */
/*   Matrix D_Stress_GP; /\* Divergence of the stress tensor *\/ */
/*   Matrix Nodal_TOT_FORCES = /\* Total forces *\/ */
/*     MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh); */
/*   strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES"); */
/*   Matrix Body_Forces_t = /\* Matrix with the body forces for TimeStep *\/ */
/*     Eval_Body_Forces(MPM_Mesh.B,MPM_Mesh.NumberBodyForces,MPM_Mesh.NumGP,TimeStep); */
/*   Matrix Contact_Forces_t = /\* Matrix with the contact forces for TimeStep *\/ */
/*     Eval_Contact_Forces(MPM_Mesh.F,MPM_Mesh.NumNeumannBC,MPM_Mesh.NumGP,TimeStep); */

/*   for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){ */

/*     /\* 1º Define element of the GP *\/ */
/*     GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i], MPM_Mesh.NumberNodes[i]); */

/*     /\* 2º Evaluate the shape function and its gradient in the GP *\/ */
/*     N_dNdx_GP = Get_Operator("N_dNdx", GP_Element, */
/* 			     MPM_Mesh, FEM_Mesh); */
/*     /\* 3º Asign values to the pointer structures *\/ */
/*     N_GP = MatAssign(1,N_dNdx_GP.N_cols,NAN,N_dNdx_GP.nM[0],NULL); */
/*     dNdx_GP = MatAssign(2,N_dNdx_GP.N_cols,NAN,NULL, */
/* 			(double **)malloc(2*sizeof(double *))); */
/*     dNdx_GP.nM[0] = N_dNdx_GP.nM[1]; */
/*     dNdx_GP.nM[1] = N_dNdx_GP.nM[2]; */
/*     free(N_dNdx_GP.nM); */
           
/*     /\* 4º Get the B_T matrix for the derivates *\/ */
/*     B = Get_B_GP(dNdx_GP); */
/*     FreeMat(dNdx_GP); */
/*     B_T = Transpose_Mat(B); */
/*     FreeMat(B); */
    
/*     /\* 5º Asign to an auxiliar variable the value of the stress tensor *\/ */
/*     Stress_GP.nV = MPM_Mesh.Phi.Stress.nM[i]; */

/*     /\* 6º Get the divergence stress tensor evaluates in the Gauss-Point */
/*        and free the B_T matrix *\/ */
/*     D_Stress_GP = Scalar_prod(B_T,Stress_GP), FreeMat(B_T); */
    
/*     /\* 7º Calcule the volumen of the Gauss-Point *\/ */
/*     GP_mass = MPM_Mesh.Phi.mass.nV[i]; */
/*     Vol_GP = GP_mass/MPM_Mesh.Phi.rho.nV[i]; */

/*     /\* 8º Damage parameter for the Gauss-point (fracture) *\/ */
/*     ji_GP = MPM_Mesh.Phi.ji.nV[i]; */

/*     /\* 9º Acumulate this forces to the total array with the internal forces *\/ */
/*     for(int k = 0; k<GP_Element.NumberNodes; k++){ */
/*       /\* Get the node for the GP *\/ */
/*       GP_I = GP_Element.Connectivity[k]; */
/*       /\* Evaluate the GP function in the node *\/ */
/*       N_GP_I = N_GP.nV[k]; */
/*       /\* If this node has a null Value of the SHF continue *\/ */
/*       if(fabs(N_GP_I) <= TOL_zero) continue; */
/*       /\* Loop in the dimensions *\/ */
/*       for(int l = 0; l<NumberDimensions; l++){ */
/* 	/\* Add the internal forces with */
/* 	   damage variable option *\/ */
/* 	Nodal_TOT_FORCES.nM[l][GP_I] -= (1-ji_GP)* */
/* 	  D_Stress_GP.nV[k*NumberDimensions+l]*Vol_GP; */
/* 	/\* Add the body forces *\/ */
/* 	Nodal_TOT_FORCES.nM[l][GP_I] += */
/* 	  N_GP_I*Body_Forces_t.nM[l][i]*GP_mass; */
/* 	/\* Add the contact forces *\/ */
/* 	Nodal_TOT_FORCES.nM[l][GP_I] += */
/* 	  N_GP_I*Contact_Forces_t.nM[l][i]*Vol_GP; */
/*       } */
/*     } */
    
/*     /\* 10 º Free memory *\/ */
/*     free(GP_Element.Connectivity); */
/*     FreeMat(D_Stress_GP); */
/*     FreeMat(N_GP); */

/*   } */

/*   /\* 11º Free memory *\/ */
/*   FreeMat(Contact_Forces_t); */
/*   FreeMat(Body_Forces_t); */
  
/*   return Nodal_TOT_FORCES; */
  
/* } */

/*******************************************************/

