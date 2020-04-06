#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*********************************************************************/

Matrix compute_NodalFields(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{
  /* Define some dimensions */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Np = MPM_Mesh.NumGP;
  int Ip;

  /* Auxiliar variables */
  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pI; /* Evaluation of the GP in the node */
  double m_p; /* Mass of the GP */
  double M_Ip;
  Element Nodes_p; /* Element for each Gauss-Point */

  /* Output */
  Matrix Phi_I;

  /* Allocate the output list of fields */
  Phi_I = MatAllocZ(Nnodes, Ndim + 1);
  
  /* Iterate over the GP to get the nodal values */
  for(int p = 0 ; p<Np ; p++){

    /* Define element of the GP */
    Nodes_p = get_Element(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);
   
    /* Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Get the nodal mass and mommentum */
    for(int I = 0 ; I<Nodes_p.NumberNodes ; I++){
      
      /* Get the node for the GP */
      Ip = Nodes_p.Connectivity[I];
      
      /* Evaluate the GP function in the node */
      ShapeFunction_pI = ShapeFunction_p.nV[I];
      
      /* If this node has a null Value of the SHF continue */
      if(ShapeFunction_pI == 0){
	continue;
      }

      /* Gauss-point contribution to node I */
      M_Ip = m_p*ShapeFunction_pI;

      /* Nodal mass */
      Phi_I.nM[Ip][Ndim] += M_Ip;
      
      /* Nodal momentum */
      for(int i = 0 ; i<Ndim ; i++){
	Phi_I.nM[Ip][i] += M_Ip*MPM_Mesh.Phi.vel.nM[p][i];
      }   
    }

    /* Free the value of the shape functions */
    FreeMat(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }
 
  return Phi_I;
  
}

/*******************************************************/

Matrix compute_NodalVelocity(Mesh FEM_Mesh, Matrix Phi_I)
{
  /*
    Get the nodal velocity using : 
    v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
    Initialize nodal velocities 
  */

  /* Define some dimensions */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /* Define output */
  Matrix V_I = MatAllocZ(Nnodes,Ndim);
  strcpy(V_I.Info,"VELOCITY");

  /* Value of the nodal mass */
  double M_I;
  
  /* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<Nnodes ; i++){
    M_I = Phi_I.nM[i][Ndim];
    if(M_I > 0){
      for(int j = 0 ; j<Ndim ; j++){      
	V_I.nM[i][j] = Phi_I.nM[i][j]/M_I;
      }
    }    
  }
  
  return V_I;
}

/*******************************************************/

Matrix PCE_Predictor(GaussPoint MPM_Mesh,
		     Mesh FEM_Mesh,
		     Matrix Nodal_VEL,
		     Matrix Nodal_MASS,
		     Time_Int_Params Params,
		     double DeltaTimeStep){
  /*
    Get the nodal velocity using : 
    v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
    Initialize nodal velocities 
  */

  /* Time integration parameters */
  double gamma = Params.GA_gamma;

  Matrix N_GP;  /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double GP_mass; /* Mass of the GP */ 
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I;
  
  /* Matrix Vel_Mesh */
  Nodal_VEL = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_VEL.Info,"VELOCITY");
 
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
      /* Nodal momentum */
      if(Nodal_MASS.nV[GP_I] > 0){
	for(int l = 0 ; l<NumberDimensions ; l++){
	  Nodal_VEL.nM[l][GP_I] +=
	    GP_mass*N_GP_I*
	    (MPM_Mesh.Phi.vel.nM[i][l] +
	     MPM_Mesh.Phi.acc.nM[i][l]*(1-gamma)*DeltaTimeStep);
	}
      }
    }

    /* 7º Free the value of the shape functions */
    FreeMat(N_GP), free(GP_Element.Connectivity);
  }

  /* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    if(Nodal_MASS.nV[i] > 0){
      for(int j = 0 ; j<NumberDimensions ; j++){
	Nodal_VEL.nM[j][i] = (double)Nodal_VEL.nM[j][i]/Nodal_MASS.nV[i];
      }
    }    
  }
  
  return Nodal_VEL;
}

/*******************************************************/

Matrix PCE_Corrector(Mesh FEM_Mesh,
		     Matrix Nodal_VEL,
		     Matrix Nodal_TOT_FORCES,
		     Matrix Nodal_MASS,
		     Time_Int_Params Params,
		     double DeltaTimeStep){
  /*
    Get the nodal velocity using : 
    v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
    Initialize nodal velocities 
  */

  /* Time integration parameters */
  double gamma = Params.GA_gamma;
  
  /* Matrix Vel_Mesh */
  /*   = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh); */
  /* strcpy(Vel_Mesh.Info,"VELOCITY"); */
  
  /* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    if((FEM_Mesh.ActiveNode[i] > 0) && (Nodal_MASS.nV[i] > 0)){
      for(int j = 0 ; j<NumberDimensions ; j++){
	Nodal_VEL.nM[j][i] +=
	  gamma*DeltaTimeStep*
	  Nodal_TOT_FORCES.nM[j][i]/Nodal_MASS.nV[i];
      }
    }
  }
  
  return Nodal_VEL;
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
  int Ndim = NumberDimensions;
  int N_Acceleration_dim = Ndim;
  int N_Velocity_dim = Ndim;

  /* Nodal values the fields */
  Matrix Nodal_Mass = MatAllocZ(1,N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t0 = MatAllocZ(N_Acceleration_dim,N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t1 = MatAllocZ(N_Acceleration_dim,N_Nodes_Mesh);
  Matrix Nodal_Velocity = MatAllocZ(N_Velocity_dim,N_Nodes_Mesh);
 
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
  
  /* 1º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<N_GPs ; i++){

    /* 2º Define element of the GP */
    GP_Element = get_Element(i, MPM_Mesh.ListNodes[i], MPM_Mesh.NumberNodes[i]);
    
    /* 3º Evaluate the shape function in the coordinates of the GP */
    N_GP = compute_ShapeFunction(GP_Element, MPM_Mesh,FEM_Mesh);
   
    /* 4º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Acc_GP.nV = MPM_Mesh.Phi.acc.nM[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];

    /* 5º Get the nodal kinetics (I) */
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
    /* 6º Free the value of the shape functions */
    FreeMat(N_GP);
    free(GP_Element.Connectivity);
  }

  /* 7º Get the nodal kinetics (II) */
  for(int i = 0 ; i<N_Nodes_Mesh ; i++){
    if(Nodal_Mass.nV[i] > 0){
      /* Get the nodal acceleration */
      for(int j = 0 ; j<N_Acceleration_dim ; j++){
	Nodal_Acceleration_t0.nM[j][i] =
	  Nodal_Acceleration_t0.nM[j][i]/Nodal_Mass.nV[i];
      }
      /* Get the nodal velocity */    
      for(int j = 0 ; j<N_Velocity_dim ; j++){
	Nodal_Velocity.nM[j][i] =
	  Nodal_Velocity.nM[j][i]/Nodal_Mass.nV[i];
      }
    }
  }

  /* 8º Asign Kinetics values to matricial tables */

  int N_Kinetics_dim = 1 + 2*N_Acceleration_dim + N_Velocity_dim;
  double SizeKinetics = N_Kinetics_dim*sizeof(double *);

  Matrix Nodal_Kinetics =
    MatAssign(N_Kinetics_dim,N_Nodes_Mesh,NAN,NULL,
	      (double**)malloc(SizeKinetics));
  strcpy(Nodal_Kinetics.Info,
	 "MASS;ACCELERATION_t0;ACCELERATION_t1;VELOCITY");
  
  Nodal_Kinetics.nM[0] = Nodal_Mass.nV;
  for(int i = 0 ; i<Ndim ; i++){
    Nodal_Kinetics.nM[1+i] = Nodal_Acceleration_t0.nM[i];
    Nodal_Kinetics.nM[1+Ndim+i] = Nodal_Acceleration_t1.nM[i];
    Nodal_Kinetics.nM[1+2*Ndim+i] = Nodal_Velocity.nM[i];
  }

  /* Free table of pointers */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);

  /* Return the kinetics variables in the nodes */
  return Nodal_Kinetics;
}

/*******************************************************/


Matrix GetNodalVelocityDisplacement(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*!
 *  Nodal_Kinetics = {m, d, v}
 */
{
  /* */
  int N_Nodes_Mesh = FEM_Mesh.NumNodesMesh;
  int N_GPs = MPM_Mesh.NumGP;
  
  /* Sizes */
  int Ndim = NumberDimensions;

  /* Nodal values the fields */
  Matrix Nodal_Mass = MatAllocZ(1,N_Nodes_Mesh);
  Matrix Nodal_Displacement = MatAllocZ(Ndim,N_Nodes_Mesh);
  Matrix Nodal_Velocity = MatAllocZ(Ndim,N_Nodes_Mesh);
  Matrix Nodal_Acceleration = MatAllocZ(Ndim,N_Nodes_Mesh);
 
  /* GPs values of the fields */
  double Mass_GP; /* Mass of the GP */
  Matrix Disp_GP; /* Stress of the GP */
  Matrix Vel_GP; /* Velocity of the GP */
  Matrix Acel_GP; /* Velocity of the GP */
  
  /* Shape function auxiliar variables */
  Element GP_Element; /* Element for each Gauss-Point */
  Matrix N_GP; /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double Mass_GP_I; /* Nodal contribution of each GP */  
  int GP_I;
  
  /* 1º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<N_GPs ; i++){

    /* 2º Define element of the GP */
    GP_Element = get_Element(i, MPM_Mesh.ListNodes[i], MPM_Mesh.NumberNodes[i]);
    
    /* 3º Evaluate the shape function in the coordinates of the GP */
    N_GP = compute_ShapeFunction(GP_Element, MPM_Mesh,FEM_Mesh);
   
    /* 4º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];
    Disp_GP.nV = MPM_Mesh.Phi.dis.nM[i];
    Acel_GP.nV = MPM_Mesh.Phi.acc.nM[i];

    /* 5º Get the nodal kinetics (I) */
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

      for(int l = 0 ; l<Ndim ; l++){
	/* Nodal displacement */
	Nodal_Displacement.nM[l][GP_I] += Disp_GP.nV[l]*Mass_GP_I;
	/* Nodal velocity */
	Nodal_Velocity.nM[l][GP_I] += Vel_GP.nV[l]*Mass_GP_I;
	/* Nodal acceleration */
	Nodal_Acceleration.nM[l][GP_I] += Acel_GP.nV[l]*Mass_GP_I;
      }
      
    }
    /* 6º Free the value of the shape functions */
    FreeMat(N_GP);
    free(GP_Element.Connectivity);
  }

  /* 7º Get the nodal kinetics (II) */
  for(int i = 0 ; i<N_Nodes_Mesh ; i++){
    if(Nodal_Mass.nV[i] > 0){
      /* Get the nodal displacement */
      for(int j = 0 ; j<Ndim ; j++){
	Nodal_Displacement.nM[j][i] = Nodal_Displacement.nM[j][i]/Nodal_Mass.nV[i];
      }
      /* Get the nodal velocity */    
      for(int j = 0 ; j<Ndim ; j++){
	Nodal_Velocity.nM[j][i] = Nodal_Velocity.nM[j][i]/Nodal_Mass.nV[i];
      }
      /* Get the nodal acceleration */    
      for(int j = 0 ; j<Ndim ; j++){
	Nodal_Acceleration.nM[j][i] = Nodal_Acceleration.nM[j][i]/Nodal_Mass.nV[i];
      }
    }
  }

  /* 8º Asign Kinetics values to matricial tables */

  int N_Kinetics_dim = 1 + 3*Ndim;
  double SizeKinetics = N_Kinetics_dim*sizeof(double *);

  Matrix Nodal_Kinetics =
    MatAssign(N_Kinetics_dim,N_Nodes_Mesh,NAN,NULL,
	      (double**)malloc(SizeKinetics));
  strcpy(Nodal_Kinetics.Info,
	 "MASS;DISPLACEMENT;VELOCITY;ACCELERATION");
  
  Nodal_Kinetics.nM[0] = Nodal_Mass.nV;
  for(int i = 0 ; i<Ndim ; i++){
    Nodal_Kinetics.nM[1+i] = Nodal_Displacement.nM[i];
    Nodal_Kinetics.nM[1+Ndim+i] = Nodal_Velocity.nM[i];
    Nodal_Kinetics.nM[1+2*Ndim+i] = Nodal_Acceleration.nM[i];
  }

  /* Free table of pointers */
  free(Nodal_Displacement.nM);
  free(Nodal_Velocity.nM);
  free(Nodal_Acceleration.nM);

  /* Return the kinetics variables in the nodes */
  return Nodal_Kinetics;
}

/*******************************************************/
