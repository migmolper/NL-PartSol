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

Matrix GetNodalKinetics(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*!
 *  Nodal_Kinetics = {mass, a0, a1, v}
 */
{
  /* */
  int N_Nodes_Mesh = FEM_Mesh.NumNodesMesh;
  int N_GPs = MPM_Mesh.NumGP;
  
  /* Sizes */
  int N_dim = NumberDimensions;
  int N_Acceleration_dim = N_dim;
  int N_Velocity_dim = N_dim;
  int N_Kinetics_dim = 1 + 2*N_Acceleration_dim + N_Velocity_dim;
  double SizeKinetics = N_Kinetics_dim*sizeof(double *);

  /* Nodal values the fields */
  Matrix Nodal_Mass =
    MatAllocZ(1,N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t0 =
    MatAllocZ(N_Acceleration_dim,N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t1 =
    MatAllocZ(N_Acceleration_dim,N_Nodes_Mesh);
  Matrix Nodal_Velocity =
    MatAllocZ(N_Velocity_dim,N_Nodes_Mesh);
  Matrix Nodal_Kinetics =
    MatAssign(N_Kinetics_dim,N_Nodes_Mesh,NAN,NULL,
	      (double**)malloc(SizeKinetics));
  strcpy(Nodal_Kinetics.Info,
	 "MASS;ACCELERATION_t0;ACCELERATION_t1;VELOCITY");

  /* GPs values of the fields */
  double Mass_GP; /* Mass of the GP */
  Matrix Vel_GP; /* Velocity of the GP */
  Matrix Acc_GP; /* Stress of the GP */
  
  /* Shape function auxiliar variables */
  Element GP_Element; /* Element for each Gauss-Point */
  Matrix N_GP; /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double Mass_GP_I; /* Nodal contribution of each GP */  
  int GP_I;
  
  /* 1º Asign Kinetics values to matricial tables */
  Nodal_Kinetics.nM[0] = Nodal_Mass.nV;
  for(int i = 0 ; i<N_dim ; i++){
    Nodal_Kinetics.nM[1+i] = Nodal_Acceleration_t0.nM[i];
    Nodal_Kinetics.nM[1+N_dim+i] = Nodal_Acceleration_t1.nM[i];
    Nodal_Kinetics.nM[1+2*N_dim+i] = Nodal_Velocity.nM[i];
  }
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<N_GPs ; i++){

    /* 3º Define element of the GP */
    GP_Element = GetElementGP(i, MPM_Mesh.ListNodes[i],
			      MPM_Mesh.NumberNodes[i]);
    
    /* 4º Evaluate the shape function in the coordinates of the GP */
    N_GP = Get_Operator("N",GP_Element,
			MPM_Mesh,FEM_Mesh);
   
    /* 5º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Acc_GP.nV = MPM_Mesh.Phi.acc.nM[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];

    /* 6º Get the nodal kinetics (I) */
    for(int k = 0 ; k<GP_Element.NumberNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Element.Connectivity[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(N_GP_I == 0){
	continue;
      }
      /* Nodal constribution */
      Mass_GP_I = Mass_GP*N_GP_I;
      /* Nodal mass */
      Nodal_Mass.nV[GP_I] += Mass_GP_I;
      /* Nodal acceleration t0 */
      for(int l = 0 ; l<N_Acceleration_dim ; l++){
	Nodal_Acceleration_t0.nM[l][GP_I] += Acc_GP.nV[l]*Mass_GP_I;
      }
      /* Nodal velocity */
      for(int l = 0 ; l<N_Velocity_dim ; l++){
	Nodal_Velocity.nM[l][GP_I] += Vel_GP.nV[l]*Mass_GP_I;
      }
    }
    /* 7º Free the value of the shape functions */
    FreeMat(N_GP);
    free(GP_Element.Connectivity);
  }

  /* 8º Get the nodal kinetics (II) */
  for(int i = 0 ; i<N_Nodes_Mesh ; i++){
    if(Nodal_Mass.nV[i] > 0){
      /* Get the nodal acceleration */
      for(int j = 0 ; j<N_Acceleration_dim ; j++){
	Nodal_Acceleration_t0.nM[j][i] =
	  Nodal_Acceleration_t0.nM[j][i]*(double)1/Nodal_Mass.nV[i];
      }
      /* Get the nodal velocity */    
      for(int j = 0 ; j<N_Velocity_dim ; j++){
	Nodal_Velocity.nM[j][i] =
	  Nodal_Velocity.nM[j][i]*(double)1/Nodal_Mass.nV[i];
      }
    }
  }

  /* Free table of pointers */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);

  /* Return the kinetics variables in the nodes */
  return Nodal_Kinetics;
}

/*******************************************************/

void UpdateNodalMass(GaussPoint MPM_Mesh, Mesh FEM_Mesh, Matrix Nodal_Kinetics)
{

  int N_Nodes = FEM_Mesh.NumNodesMesh;
  Matrix N_GP;  /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double GP_mass; /* Mass of the GP */ 
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I;

  /* 1º Allocate the output list of fields */
  for(int i = 0 ; i<N_Nodes ; i++){
    Nodal_Kinetics.nM[0][i] = 0.0;
  }
  
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
      Nodal_Kinetics.nM[0][GP_I] += GP_mass*N_GP_I;
    }

    /* 7º Free the value of the shape functions */
    FreeMat(N_GP), free(GP_Element.Connectivity);
  } 
}
