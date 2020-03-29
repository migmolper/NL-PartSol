#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*******************************************************/

void UpdateGaussPointStrain(GaussPoint MPM_Mesh,
			    Mesh FEM_Mesh,
			    Matrix Mesh_Vel)
/*!
 * \brief Brief description of UpdateGaussPointStrain.
 *        Update the strain state of the body with the
 *        infinetesimal strain tensor. 
 *
 * \f[ 
 \Delta\Epsilon_{ij,p}^{k-1/2} = 
 \frac{\Delta t}{2} \cdot \sum_{I=0}^{Nn}(N^{k}_{Ip,j}
 \cdot v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot v_{jI}^{k-1/2})
  \f]
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param FEM_Mesh : Background mesh for calculations.
 *  @param Mesh_Vel : Nodal values of the velocity field.
 *
 */
{
  
  /* Variable definition */
  Element GP_Element; /* Element for each Gauss-Point */
  Matrix Element_Velocity; /* Array with the nodal velocities */
  Matrix dNdx_GP; /*  Matrix with the nodal derivatives */
  Matrix B; /* B matrix to get the deformation */
  Matrix Delta_Strain_GP; /* Vectoriced Strain tensor */
  double Delta_TraceStrain_GP; /* Increment of the trace of the Stress tensor */
 
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
   
    /* 1º Define element of the GP */
    GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i],
			      MPM_Mesh.NumberNodes[i]);
    
    /* 2º Get the element gradient */
    dNdx_GP = Get_Operator("dNdx", GP_Element,
			   MPM_Mesh, FEM_Mesh);
	    
    /* 3º Calcule the B matrix and free the gradient */
    B = Get_B_GP(dNdx_GP);
    FreeMat(dNdx_GP);

    /* 4º Get the nodal velocities in the element and free connectivity */
    Element_Velocity = GetElementField(Mesh_Vel, GP_Element);
    free(GP_Element.Connectivity);
   
    /* 5º Multiply B by the velocity array and by the time step to get
       the increment stress tensor */
    Delta_Strain_GP = Scalar_prod(B,Element_Velocity);
    Delta_Strain_GP = Matrix_x_Scalar(Delta_Strain_GP,
					  DeltaTimeStep);

    /* 6º Free the array with the nodal velocity of the element and the B matrix */
    FreeMat(Element_Velocity);
    FreeMat(B);

    /* 7º Udate the Gauss-Point strain tensor */
    Delta_TraceStrain_GP = 0; /* Set to zero the trace increment */
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      MPM_Mesh.Phi.Strain.nM[i][j] += Delta_Strain_GP.nV[j];
      /* Get the trace of the stress tensor */
      if(j<NumberDimensions){
	Delta_TraceStrain_GP += Delta_Strain_GP.nV[j];
      }
    }
    FreeMat(Delta_Strain_GP);

    /* 8º Update the density of the GP */
    MPM_Mesh.Phi.rho.nV[i] = UpdateGaussPointDensity(MPM_Mesh.Phi.rho.nV[i],
						     Delta_TraceStrain_GP);    
  }
}

/*******************************************************/


double UpdateGaussPointDensity(double rho_n_GP,
			       double Delta_TraceStrain_GP)
/*!
 * \brief Brief description of UpdateGaussPointDensity.
 *        Update the density field of the Gauss Point . 
 *
 *  The parameters for this functions are  :
 *  @param rho_n_GP : Density of the previous step.
 *  @param Delta_TraceStrain_GP : Increment of the trace of the strain tensor.
 *
 */
{
  double rho_n1_GP; /* Density for the next step */

  /* Update the density */
  rho_n1_GP = (double)rho_n_GP/(1 + Delta_TraceStrain_GP);

  return rho_n1_GP;  
}

/*******************************************************/


void UpdateGaussPointStress(GaussPoint MPM_Mesh)
/*!
 * \brief Brief description of UpdateGaussPointStress.
 *        Update the stress state of the body evaluating the
 *        strains in each Gauss-Point to get the stress state.
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *
 */
{
  /* Variable definition  */
  Matrix Strain_n1; /* Strain tensor of the GP in the next step */
  Matrix Stress_n; /* Stress tensor of the GP in the previous step */
  Matrix Stress_n1; /* Stress tensor of the GP in the next step */
  Material Material_GP; /* Material properties of the GP */
  int N_GP = MPM_Mesh.NumGP;
  int Mat_GP;

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<N_GP ; i++){
    if(MPM_Mesh.Phi.ji.nV[i] != 1.0){
      
      /* 3º Asign materials of the GP */
      Mat_GP = MPM_Mesh.MatIdx[i];
      Material_GP = MPM_Mesh.Mat[Mat_GP]; 
      
      /* 4º Use pointers to memory manage */
      Strain_n1.nV = MPM_Mesh.Phi.Strain.nM[i];
      Stress_n.nV = MPM_Mesh.Phi.Stress.nM[i];
      Stress_n1.nV = MPM_Mesh.Phi.Stress.nM[i];
    
      /* 5º Get the new stress tensor */
      if(strcmp(Material_GP.Type,"LE") == 0){
	Stress_n1 = LinearElastic(Strain_n1,Stress_n,Material_GP);
      }
      else{
	exit(0);
      }
      
      /* 6º Get the deformation energy */
      MPM_Mesh.Phi.W.nV[i] =
	W_LinearElastic(Strain_n1,Stress_n1,MPM_Mesh.Phi.ji.nV[i]);
    }
    else{
      for(int j = 0 ; j<MPM_Mesh.Phi.Stress.N_cols ; j++){
	MPM_Mesh.Phi.Stress.nM[i][j] = 0.0;
      }
    }
  }
}

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep)
/*!
 * \brief Brief description of GetNodalForces.
 *        Compute the nodal contribution of each GP to the total forces. 
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param FEM_Mesh : Background mesh for calculations.
 *  @param TimeStep : Parameter to evaluate the loads curve.
 *
 */
{
  double GP_mass; /* Mass of the GP */
  double Vol_GP; /* Gauss-Point volumen */
  double ji_GP; /* Damage parameter */
  Matrix N_GP; /* Matrix with the nodal shape functions */
  double N_GP_I; /* Evaluation of the GP in the node */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix N_dNdx_GP; /* Operator Matrix */
  Matrix B, B_T; /* B matrix */
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I; /* Node of the GP */
  Matrix Stress_GP = /* Stress tensor of a Gauss-Point */
    MatAssign(MPM_Mesh.Phi.Stress.N_cols,1,NAN,NULL,NULL);
  Matrix D_Stress_GP; /* Divergence of the stress tensor */
  Matrix Nodal_TOT_FORCES = /* Total forces */
    MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES");
  Matrix Body_Forces_t = /* Matrix with the body forces for TimeStep */
    Eval_Body_Forces(MPM_Mesh.B,MPM_Mesh.NumberBodyForces,MPM_Mesh.NumGP,TimeStep);
  Matrix Contact_Forces_t = /* Matrix with the contact forces for TimeStep */
    Eval_Contact_Forces(MPM_Mesh.F,MPM_Mesh.NumNeumannBC,MPM_Mesh.NumGP,TimeStep);

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 1º Define element of the GP */
    GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i], MPM_Mesh.NumberNodes[i]);

    /* 2º Evaluate the shape function and its gradient in the GP */
    N_dNdx_GP = Get_Operator("N_dNdx", GP_Element,
			     MPM_Mesh, FEM_Mesh);
    /* 3º Asign values to the pointer structures */
    N_GP = MatAssign(1,N_dNdx_GP.N_cols,NAN,N_dNdx_GP.nM[0],NULL);
    dNdx_GP = MatAssign(2,N_dNdx_GP.N_cols,NAN,NULL,
			(double **)malloc(2*sizeof(double *)));
    dNdx_GP.nM[0] = N_dNdx_GP.nM[1];
    dNdx_GP.nM[1] = N_dNdx_GP.nM[2];
    free(N_dNdx_GP.nM);
           
    /* 4º Get the B_T matrix for the derivates */
    B = Get_B_GP(dNdx_GP);
    FreeMat(dNdx_GP);
    B_T = Transpose_Mat(B);
    FreeMat(B);
    
    /* 5º Asign to an auxiliar variable the value of the stress tensor */
    Stress_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

    /* 6º Get the divergence stress tensor evaluates in the Gauss-Point 
     and free the B_T matrix */
    D_Stress_GP = Scalar_prod(B_T,Stress_GP), FreeMat(B_T);
    
    /* 7º Calcule the volumen of the Gauss-Point */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];
    Vol_GP = GP_mass/MPM_Mesh.Phi.rho.nV[i];

    /* 8º Damage parameter for the Gauss-point (fracture) */
    ji_GP = MPM_Mesh.Phi.ji.nV[i];

    /* 9º Acumulate this forces to the total array with the internal forces */  
    for(int k = 0; k<GP_Element.NumberNodes; k++){
      /* Get the node for the GP */
      GP_I = GP_Element.Connectivity[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(fabs(N_GP_I) <= TOL_zero) continue;
      /* Loop in the dimensions */
      for(int l = 0; l<NumberDimensions; l++){
	/* Add the internal forces with 
	 damage variable option */
	Nodal_TOT_FORCES.nM[l][GP_I] -= (1-ji_GP)* 
	  D_Stress_GP.nV[k*NumberDimensions+l]*Vol_GP;
	/* Add the body forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Body_Forces_t.nM[l][i]*GP_mass;
	/* Add the contact forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Contact_Forces_t.nM[l][i]*Vol_GP;
      }      
    }
    
    /* 10 º Free memory */
    free(GP_Element.Connectivity);
    FreeMat(D_Stress_GP);
    FreeMat(N_GP);

  }

  /* 11º Free memory */
  FreeMat(Contact_Forces_t);
  FreeMat(Body_Forces_t);
  
  return Nodal_TOT_FORCES;
  
}

/*******************************************************/

void FE_Update_Momentum(Mesh FEM_Mesh,
			Matrix Nodal_MOMENTUM,
			Matrix Nodal_TOT_FORCES)
/*!
 * \brief Brief description of UpdateGridNodalMomentum.
 *        Compute the nodal contribution of each GP to the total forces.
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param Nodal_MOMENTUM : Nodal value of the momentum.
 *  @param Nodal_TOT_FORCES : Nodal value of the total forces.
 *
 */
{
  /* Update the grid nodal momentum */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      if(FEM_Mesh.ActiveNode[i] > 0){
	Nodal_MOMENTUM.nM[j][i] +=
	  DeltaTimeStep*Nodal_TOT_FORCES.nM[j][i];
      }
    }
  }  
}

/*******************************************************/

void GA_UpdateNodalKinetics(Mesh FEM_Mesh,
			    Matrix Nodal_Kinetics,
			    Matrix Nodal_Forces,
			    Time_Int_Params Params)
/*!
 * \brief Brief description of Update_Nodal_Acceleration_Velocity.
 *
 *  The parameters for this functions are  :
 *  @param FEM_Mesh
 *  @param Nodal_Kinetics = {m, a0, a1, v}
 *  @param Nodal_Forces
 *  @param Params
 */
{
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = NumberDimensions;
  double SizeTable = N_dim*sizeof(double *);
  double Mass_I;

  /* Time integration parameters */
  double alpha = Params.GA_alpha;
  double gamma = Params.GA_gamma;

  /* Asign forces to an auxiliar variable */
  Matrix F = Nodal_Forces;
  
  /* Nodal values the fields */
  Matrix Nodal_Acceleration_t0 =
    MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));
  Matrix Nodal_Acceleration_t1 =
    MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));
  Matrix Nodal_Velocity =
    MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));

  /* 1º Asign Kinetics values to matricial tables */
  for(int i = 0 ; i<N_dim ; i++){
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1+i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1+N_dim+i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+2*N_dim+i];
  }
  /* 2º Update the grid nodal variales */
  for(int i = 0 ; i<N_Nodes ; i++){
    Mass_I = Nodal_Kinetics.nM[0][i];
    if(Mass_I > 0){
      for(int j = 0 ; j<N_dim ; j++){
	/* Get the nodal acceleration t1 */
	Nodal_Acceleration_t1.nM[j][i] =
	  (F.nM[j][i]/Mass_I - alpha*Nodal_Acceleration_t0.nM[j][i])/(1-alpha);
	/* Update nodal velocity */
	Nodal_Velocity.nM[j][i] +=
	  ((1-gamma)*Nodal_Acceleration_t0.nM[j][i] +
	   gamma*Nodal_Acceleration_t1.nM[j][i])*DeltaTimeStep;
      }
    }
  }

  /* 3º Free tables */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);
  
}
