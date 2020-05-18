#include "nl-partsol.h"

/*******************************************************/

void update_Particles_FE(GaussPoint MPM_Mesh, Mesh FEM_Mesh,
			 Matrix Phi_I, Matrix F_I, double Dt){

  int Ndim = NumberDimensions;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p; /* Value of the shape-function in the GP */
  double N_pI; /* Nodal value for the GP */
  double M_I; /* Value of the nodal mass */
  int Np = MPM_Mesh.NumGP;
  int Nnodes;
  int Ip; /* Index of each tributary node for the GP */

  /* 1º iterate over the Gauss-Points */
  for(int p = 0 ; p<Np ; p++){

    /* 2º Define element of the GP */
    Nnodes = MPM_Mesh.NumberNodes[p];
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nnodes);

    /* 3º Evaluate shape function in the GP i */
    N_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);
    
    /* 4º Iterate over the nodes of the element */
    for(int I = 0; I<Nnodes; I++){
      /* Node of the GP */
      Ip = Nodes_p.Connectivity[I];
      /* Evaluate the GP function in the node */
      N_pI = N_p.nV[I];
      /* If this node has a null Value of the SHF continue */
      if(fabs(N_pI) <= TOL_zero){
	continue;
      }
      /* Get the nodal mass */
      M_I = Phi_I.nM[Ip][Ndim];
      /* Update GP cuantities with nodal values */
      for(int i = 0 ; i<Ndim ; i++){
	/* Update the GP velocities */
	MPM_Mesh.Phi.vel.nM[p][i] += Dt*N_pI*F_I.nM[Ip][i]/M_I;
	/* Update the GP position */
	MPM_Mesh.Phi.x_GC.nM[p][i] += Dt*N_pI*Phi_I.nM[Ip][i]/M_I;	  
      } 
    }
    
    /* 5º Free memory */
    free(Nodes_p.Connectivity);
    FreeMat(N_p);
  }  
}


/*******************************************************/

void update_Particles_PCE(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
			   Matrix M_I,Matrix V_I,
			   Matrix F_I,double DeltaTimeStep){
  
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  Matrix ShapeFunction_p; /* Value of the shape-function in the GP */
  Element Nodes_p; /* Element for each Gauss-Point */
  int Nnodes;
  int Ip; /* Index of each tributary node for the GP */
  double mass_I; /* Value of the nodal mass */
  double ShapeFunction_pI; /* Nodal value for the GP */

  /* 1º iterate over the Gauss-Points */
  for(int p = 0 ; p<Np ; p++){

    /* 2º Define element of the GP */
    Nnodes = MPM_Mesh.NumberNodes[p];
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nnodes);

    /* 3º Evaluate shape function in the GP i */
    ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);

    for(int i = 0 ; i<Ndim ; i++){
      MPM_Mesh.Phi.acc.nM[p][i] = 0.0;
    }
    
    /* 4º Iterate over the nodes of the element */
    for(int I = 0; I<Nodes_p.NumberNodes; I++){
      /* Node of the GP */
      Ip = Nodes_p.Connectivity[I];
      /* Evaluate the GP function in the node */
      ShapeFunction_pI = ShapeFunction_p.nV[I];
      /* If this node has a null Value of the SHF continue */
      if(fabs(ShapeFunction_pI) <= TOL_zero){
	continue;
      }
      /* Get the nodal mass */
      mass_I = M_I.nV[Ip];
      /* Update GP cuantities with nodal values */
      if(mass_I>0){
	for(int i = 0 ; i<Ndim ; i++){
	  /* Update the GP velocities */
	  MPM_Mesh.Phi.acc.nM[p][i] +=
	    ShapeFunction_pI*F_I.nM[Ip][i]/mass_I;
	  /* Update the GP velocities */
	  MPM_Mesh.Phi.vel.nM[p][i] +=
	    ShapeFunction_pI*DeltaTimeStep*F_I.nM[Ip][i]/mass_I;
	  /* Update the GP displacement */
	  MPM_Mesh.Phi.x_GC.nM[p][i] +=
	    ShapeFunction_pI*DeltaTimeStep*V_I.nM[Ip][i] +
	    ShapeFunction_pI*0.5*pow(DeltaTimeStep,2)*F_I.nM[Ip][i]/mass_I;
	}
      } 
    }
    
    /* 5º Free memory */
    free(Nodes_p.Connectivity);
    FreeMat(ShapeFunction_p);
  }  
}


/*******************************************************/

void update_Particles_GA(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
			  Matrix Nodal_Kinetics,Time_Int_Params GA_Params)
{

  
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = NumberDimensions;
  double SizeTable = N_dim*sizeof(double *);
  
  /* Shape function nodal parameters */
  Matrix N_GP; /* Value of the shape-function in the GP */
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I; /* Index of each tributary node for the GP */
  double N_I_GP; /* Nodal value for the GP */

  /* Time integration parameters */
  double beta = GA_Params.GA_beta;
  double gamma = GA_Params.GA_gamma;

  /* Asign Kinetics values to matricial tables */
  Matrix  Nodal_Acceleration_t0;
  /* = MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Acceleration_t1;
  /* = MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Velocity;
   /* = MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  for(int i = 0 ; i<N_dim ; i++){
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1+i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1+N_dim+i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+2*N_dim+i];
  }
  
  /* 1º iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Define element of the GP */
    GP_Element = get_particle_Set(i, MPM_Mesh.ListNodes[i],MPM_Mesh.NumberNodes[i]);

    /* 3º Evaluate shape function in the GP i */
    N_GP = compute_ShapeFunction(GP_Element, MPM_Mesh, FEM_Mesh);


    /* 4º Set to zero the GPs acceleration and velocity of the previous step */
    for(int k = 0 ; k<N_dim ; k++){
      MPM_Mesh.Phi.acc.nM[i][k] = 0.0;
    }


    /* 5º Iterate over the nodes of the element */
    for(int j = 0; j<GP_Element.NumberNodes; j++){
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if(fabs(N_I_GP) <= TOL_zero){
	continue;
      }
      /* Update GP cuantities with nodal values */
      for(int k = 0 ; k<N_dim ; k++){
	/* Get nodal values
	   Nodal_Kinetics = {m, a0, a1, v}
	*/
	/* Get the GP accelerations */
	MPM_Mesh.Phi.acc.nM[i][k] +=
	  N_I_GP*Nodal_Acceleration_t1.nM[k][GP_I];
      } 
    }
        
    /* 5º Iterate over the nodes of the element */
    for(int j = 0; j<GP_Element.NumberNodes; j++){
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if(fabs(N_I_GP) <= TOL_zero){
	continue;
      }
      /* Update GP cuantities with nodal values */
      for(int k = 0 ; k<N_dim ; k++){
	/* Get nodal values
	   Nodal_Kinetics = {m, a0, a1, v}
	 */
	/* Update the GP velocities */
	MPM_Mesh.Phi.vel.nM[i][k] +=
	  N_I_GP*((1-gamma)*Nodal_Acceleration_t0.nM[k][GP_I] +
		  gamma*Nodal_Acceleration_t1.nM[k][GP_I])*DeltaTimeStep;
      } 
    }

    /* 5º Iterate over the nodes of the element */
    for(int j = 0; j<GP_Element.NumberNodes; j++){
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if(fabs(N_I_GP) <= TOL_zero){
	continue;
      }
      /* Update GP cuantities with nodal values */
      for(int k = 0 ; k<N_dim ; k++){
	/* Get nodal values
	   Nodal_Kinetics = {m, a0, a1, v}
	 */
	/* Update the GP position I */
	MPM_Mesh.Phi.x_GC.nM[i][k] +=
	  N_I_GP*Nodal_Velocity.nM[k][GP_I]*DeltaTimeStep +
	  N_I_GP*((0.5 - beta)*Nodal_Acceleration_t0.nM[k][GP_I] +
		  beta*Nodal_Acceleration_t1.nM[k][GP_I])*
	  DeltaTimeStep*DeltaTimeStep;
      } 
    }
    
    /* 6º Free memory */
    free(GP_Element.Connectivity);
    FreeMat(N_GP);
  }

  /* Update the grid nodal variales */
  for(int i = 0 ; i<N_Nodes ; i++){
      for(int j = 0 ; j<N_dim ; j++){
	Nodal_Acceleration_t0.nM[j][i] = Nodal_Acceleration_t1.nM[j][i];
	Nodal_Acceleration_t1.nM[j][i] = 0.0;
      }
  }

  /* Free tables */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);
  
}

/*******************************************************/
