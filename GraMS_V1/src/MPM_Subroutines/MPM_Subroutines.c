#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MeshTools/MeshTools.h"
#include "../InOutFun/InOutFun.h"
#include "MPM_Subroutines.h"

/*********************************************************************/

Matrix GetNodalMassMomentum(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  /* 0º Variable declaration */

  /* Output */
  Matrix Nodal_FIELDS;
  
  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix N_GP; /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double GP_mass; /* Mass of the GP */

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  int GP_I;

  /* 1º Allocate the output list of fields */
  Nodal_FIELDS = MatAllocZ(NumberDimensions + 1,FEM_Mesh.NumNodesMesh);
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 4º Evaluate the shape function in the coordinates of the GP */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Fill the poligon */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	/* Get the node for the GP */
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the lenght voxel */
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      /* Evaluate the shape function */
      N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
   
    /* 5º Get the mass of the GP */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];

    /* 6º Get the nodal mass and mommentum */
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
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
    FreeMat(N_GP);
    free(GP_Connect);
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
  Matrix Vel_Mesh;
  Vel_Mesh = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
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

void UpdateGaussPointStrain(GaussPoint MPM_Mesh,
			    Mesh FEM_Mesh,
			    Matrix Mesh_Vel)
/*
  Calcule the particle stress increment :

  \Delta\Epsilon_{ij,p}^{k-1/2} = 
  \frac{\Delta t}{2} \cdot
  \sum_{I=0}^{Nn}(N^{k}_{Ip,j} \cdot
  v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot 
  v_{jI}^{k-1/2})
*/
{
  /* Gauss-Point properties */
  Matrix X_GP; /* Element coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  
  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */
  
  /* Mesh variables */
  Matrix Elem_Vel; /* Array with the nodal velocities */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix B; /* B marix to get the deformation */
  Matrix Increment_Strain_GP; /* Vectoriced Strain tensor */
  double Incr_TraceStrain; /* Increment of the trace of the Stress tensor */
 
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
   
    /* 1º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);
    
    /* 2º Get the element gradient */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      /* Fill the poligon */
      GP_ElemCoord = MatAllocZ(GP_NumNodes,NumberDimensions);
      for(int k = 0; k<GP_NumNodes ; k++){
    	/* Get the node for the GP */
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  GP_ElemCoord.nM[k][l] =
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the element coordinates of the GP */
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      /* Evaluate the shape function gradient */
      dNdx_GP = Get_dNdX_Q4(X_GP,GP_ElemCoord);
      /* Free memory */
      FreeMat(GP_ElemCoord);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Generate a matrix with the distances to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	/* Get the node for the GP */
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the GP voxel */
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      /* Evaluate the shape function gradient */
      dNdx_GP = dGIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
	    
    /* Calcule the B matrix */
    B = Get_B_GP(dNdx_GP);
    /* Free shape-function derivatives */
    FreeMat(dNdx_GP);

    /* 3º Allocate and array with the velocities of the element */
    Elem_Vel = MatAllocZ(GP_NumNodes*NumberDimensions,1);

    /* 4º Get the nodal velocities in the element */
    for(int k = 0; k<GP_NumNodes; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<NumberDimensions ; l++){
	Elem_Vel.nV[k*NumberDimensions + l] =
	  Mesh_Vel.nM[l][GP_I];
      }
    }
    /* Free data */
    free(GP_Connect);
   
    /* 5º Multiply B by the velocity array and by the time step to get
       the increment stress tensor */
    Increment_Strain_GP = Scalar_prod(B,Elem_Vel);    
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      Increment_Strain_GP.nV[j] *= DeltaTimeStep;
    }

    /* Free memory */
    FreeMat(Elem_Vel);
    FreeMat(B);

    /* 6º Set to zero the trace of the stress tensor */
    Incr_TraceStrain = 0;

    /* 7º Udate the Gauss-Point strain tensor */
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      MPM_Mesh.Phi.Strain.nM[i][j] += Increment_Strain_GP.nV[j];
      /* 7º Get the trace of the stress tensor */
      if(j<NumberDimensions){
	Incr_TraceStrain += Increment_Strain_GP.nV[j];
      }
    }
    FreeMat(Increment_Strain_GP);

    /* 8º Update the density of the GP */
    MPM_Mesh.Phi.rho.nV[i] =
      UpdateGaussPointDensity(MPM_Mesh.Phi.rho.nV[i],
			      Incr_TraceStrain);    
  }
}


/*******************************************************/


double UpdateGaussPointDensity(double rho_n,
			       double Incr_TraceStrain){

  /* 1º Density for the next step */
  double rho_n1;

  /* 2º Update the density */
  rho_n1 = (double)rho_n/(1 + Incr_TraceStrain);

  /* 3º Return density updated */
  return rho_n1;  
}


/*******************************************************/

void UpdateGaussPointStress(GaussPoint MPM_Mesh){

  /* 0º Variable declaration  */
  Matrix StrainTensor_GP;
  Matrix StressTensor_GP;

  /* 1º Switch the dimensions of the aulixiar strain tensor */
  switch(NumberDimensions){
  case 1:
    StrainTensor_GP.N_rows = 1;
     StrainTensor_GP.N_cols = 1;
    StrainTensor_GP.nM = NULL;
    break;
  case 2:
    StrainTensor_GP.N_rows = 3;
    StrainTensor_GP.N_cols = 1;
    StrainTensor_GP.nM = NULL;
    break;
  case 3:
    StrainTensor_GP.N_rows = 6;
    StrainTensor_GP.N_cols = 1;
    StrainTensor_GP.nM = NULL;
    break;
  default :
    printf("%s : %s \n","Error in UpdateGaussPointStress()",
	   "Wrong number of dimensions !!! ");
    exit(0);
  }

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    /* 3º Store in an auxiliar variable the strain tensor in the GP */
    StrainTensor_GP.nV = MPM_Mesh.Phi.Strain.nM[i];
    /* 4º Get the new stress tensor (2D Linear elastic) */
    StressTensor_GP =
      MPM_Mesh.D.LE2D(StrainTensor_GP,PoissonModulus,ElasticModulus);
    /* 5º Update the stress tensor with the new-one */
    for(int j = 0 ; j<StrainTensor_GP.N_rows ; j++){
      MPM_Mesh.Phi.Stress.nM[i][j] = StressTensor_GP.nV[j];
    }
    /* 6º Free memory */
    FreeMat(StressTensor_GP);
  }
  
}

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep)
{
  
  /* 0º Auxiliar variable declaration */

  /* Properties of each Gauss-Point */
  Matrix X_GP; /* Coordinate for each Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  double GP_mass; /* Mass of the GP */
  double Vol_GP; /* Gauss-Point volumen */

  /* Mesh properties evaluated in Gauss-Point coords */
  Matrix N_GP; /* Matrix with the nodal shape functions */
  double N_GP_I; /* Evaluation of the GP in the node */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix B, B_T;

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Node of the GP */

  /* Stress tensor of a Gauss-Point */
  Matrix StressTensor_GP; 
  StressTensor_GP.N_rows = MPM_Mesh.Phi.Stress.N_cols;
  StressTensor_GP.N_cols = 1;
  StressTensor_GP.nM = NULL;

  /* Internal forces for each node in a element (by a GP)*/  
  Matrix Div_Stress_Tensor;

  /* Total forces */
  Matrix Nodal_TOT_FORCES; 
  Nodal_TOT_FORCES = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES");

  /* Contact forces */
  Matrix Contact_Forces_t;
  Contact_Forces_t = MatAllocZ(NumberDimensions,MPM_Mesh.NumGP);

  /* Body Forces */
  Matrix Body_Forces_t;
  Body_Forces_t = MatAllocZ(NumberDimensions,MPM_Mesh.NumGP);
  int GP_Force;
  
  /* 1º Fill matrix with the body forces for TimeStep */ 
  if(MPM_Mesh.B.NumLoads>0){
    for(int i = 0 ; i<MPM_Mesh.B.NumLoads; i++){
      for(int j = 0 ; j<MPM_Mesh.B.Load_i[i].NumNodes ; j++){
	GP_Force = MPM_Mesh.B.Load_i[i].Nodes[j];
	for(int k = 0 ; k<MPM_Mesh.B.Load_i[i].Dim ; k++){
	  if( (MPM_Mesh.B.Load_i[i].Dir[k] == 1) ||
	      (MPM_Mesh.B.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > MPM_Mesh.B.Load_i[i].Value[k].Num)){
	      printf("%s : %s\n","Error in GetNodalForces()",
		      "The time step is out of the curve !!");
	      exit(0);
	    }
	    Body_Forces_t.nM[k][GP_Force] +=
	      MPM_Mesh.B.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)MPM_Mesh.B.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  
  /* 2º Fill matrix with the contact forces for TimeStep */
  if(MPM_Mesh.F.NumLoads>0){
    for(int i = 0 ; i<MPM_Mesh.F.NumLoads; i++){
      for(int j = 0 ; j<MPM_Mesh.F.Load_i[i].NumNodes ; j++){
	GP_Force = MPM_Mesh.F.Load_i[i].Nodes[j];
	for(int k = 0 ; k<MPM_Mesh.F.Load_i[i].Dim ; k++){
	  if( (MPM_Mesh.F.Load_i[i].Dir[k] == 1) ||
	      (MPM_Mesh.F.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > MPM_Mesh.F.Load_i[i].Value[k].Num)){
	      printf("%s : %s \n","Error in GetNodalForces()",
		     "The time step is out of the curve !!");
	      exit(0);
	    }
	    Contact_Forces_t.nM[k][GP_Force] +=
	      MPM_Mesh.F.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)MPM_Mesh.F.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  /* 3º Iterate over all the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 4º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 5º Evaluate the shape function in the GP */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      /* Get the coordinates of the element */
      GP_ElemCoord = MatAllocZ(GP_NumNodes,NumberDimensions);
      for(int k = 0; k<GP_NumNodes; k++){
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  GP_ElemCoord.nM[k][l] =
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the natural coordinates of the GP */
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      /* Evaluate the shape function and it gradient */
      N_GP = Q4(X_GP);
      dNdx_GP = Get_dNdX_Q4(X_GP,GP_ElemCoord);
      /* Free memory */
      FreeMat(GP_ElemCoord);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Get the distance of the GP to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the voxel of the GP */
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      /* Evaluate the shape function and it gradient */
      N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      dNdx_GP = dGIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
       
    /* 6º Get the B_T matrix for the derivates */
    B = Get_B_GP(dNdx_GP);    
    FreeMat(dNdx_GP);
    B_T = Transpose_Mat(B);
    FreeMat(B);
    
    /* 7º Asign to an auxiliar variable the value of the stress tensor */
    StressTensor_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

    /* 8º Get the divergence stress tensor evaluates in the Gauss-Point 
     and free the B_T matrix */
    Div_Stress_Tensor = Scalar_prod(B_T,StressTensor_GP);    
    FreeMat(B_T);
    
    /* 9º Calcule the volumen of the Gauss-Point */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];
    Vol_GP = GP_mass/MPM_Mesh.Phi.rho.nV[i];

    /* 10º Acumulate this forces to the total array with the internal forces */  
    for(int k = 0; k<GP_NumNodes; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(N_GP_I == 0){
	continue;
      }
      for(int l = 0; l<NumberDimensions; l++){
	/* 10aº Add the internal forces */
	Nodal_TOT_FORCES.nM[l][GP_I] -=
	  Div_Stress_Tensor.nV[k*NumberDimensions + l]*Vol_GP;
	/* 10bº Add the body forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  GP_mass*N_GP_I*Body_Forces_t.nM[l][i];	
	/* 10cº Add the contact forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Contact_Forces_t.nM[l][i]*Vol_GP;	
      }      
    }
	
    /* 10dº Free memory */
    free(GP_Connect);
    FreeMat(Div_Stress_Tensor);
    FreeMat(N_GP);

  }

  /* 11º Free memory */
  FreeMat(Contact_Forces_t);
  FreeMat(Body_Forces_t);  
  
  return Nodal_TOT_FORCES;
  
}

/*******************************************************/

void UpdateGridNodalMomentum(Mesh FEM_Mesh,
			     Matrix Nodal_MOMENTUM,
			     Matrix Nodal_TOT_FORCES)
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

void UpdateVelocityAndPositionGP(GaussPoint MPM_Mesh,
				 Mesh FEM_Mesh,
				 Matrix Nodal_MASS,
				 Matrix Nodal_MOMENTUM,
				 Matrix Nodal_TOT_FORCES){

  /* 0º Variable declaration */

  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix N_GP; /* Value of the shape-function in the GP */

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  int GP_I;

  /* Mesh properties */
  double mass_I; /* Value of the nodal mass */
  double N_I_GP; /* Nodal value for the GP */

  /* 1º iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 3º Evaluate MPM-Q4 shape function  */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
    }
    /* 4º Evaluate GIMP shape function */
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Get the distance of the GP to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_Connect[k]][l];
    	}
      }
      /* Get the GP voxel */
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      /* Evaluate shape function */
      N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
    /* 5º Iterate over the nodes of the element */
    for(int j = 0; j<GP_NumNodes; j++){
      /* Node of the GP */
      GP_I = GP_Connect[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if(N_I_GP == 0){
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
    
    /* 6º Free memory */
    free(GP_Connect);
    FreeMat(N_GP);
    
  }  
}


/*******************************************************/


  