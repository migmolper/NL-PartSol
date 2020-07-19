#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Matrix compute_Nodal_Kinetics(GaussPoint, Mesh);
static void update_Kinetics(Mesh, Matrix, Matrix, Time_Int_Params);
static Matrix GetNodalVelocityDisplacement(GaussPoint, Mesh);
static void update_Particles(GaussPoint, Mesh, Matrix, Time_Int_Params);

/**************************************************************/

void U_GA(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /*!
    Time step 
  */
  int TimeStep;

  /*!
    Control parameters of the generalized-alpha algorithm 
    all the parameters are controled by a simple parameter :
    SpectralRadius 
  */
  Time_Int_Params Params;
  Params.GA_alpha = (2*SpectralRadius-1)/(1+SpectralRadius);
  Params.GA_beta = (5-3*SpectralRadius)/(pow((1+SpectralRadius),2)*(2-SpectralRadius));
  Params.GA_gamma = 3/2 - Params.GA_alpha;

  /*!
    Auxiliar variable for the nodal kinetics
    Nodal_Kinetics = {mass, a0, a1, v}
  */
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  Matrix V_I;
  Matrix Nodal_Kinetics;

  /*!
    Nodal forces for the balance 
  */
  Matrix F_I = memory_to_matrix__MatrixLib__(Ndim,Nnodes,NULL);
  Matrix R_I = memory_to_matrix__MatrixLib__(Ndim,Nnodes,NULL);

  
  puts("*************************************************");
  puts(" First step : Get the nodal kinetics");
  puts(" \t WORKING ...");
  Nodal_Kinetics = compute_Nodal_Kinetics(MPM_Mesh,FEM_Mesh);
  /* V_I =  MatAssign(Nnodes,Ndim,NAN,NULL, */
  /* 		   (double**)malloc(Nnodes*sizeof(double *))); */
  for(int i = 0 ; i<Ndim ; i++)
    {
      V_I.nM[i] = Nodal_Kinetics.nM[1+2*Ndim+i];
    }
  puts(" \t DONE !!! \n");


  for(TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("First step : Compute nodal kinetics ... WORKING",TimeStep);
      /* BCC_Nod_VALUE(FEM_Mesh, V_I, TimeStep); */
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
	{
	  /*! 
	    Print Nodal values after appling the BCCs 
	  */
	  WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_Kinetics,
		       (int)TimeStep/ResultsTimeStep);
	  /*!
	    Print GPs results
	  */
	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,"ALL",
		       (int)TimeStep/ResultsTimeStep,ResultsTimeStep);
	}
      
      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute equilibrium ... WORKING",TimeStep);
      F_I = compute_equilibrium_U(V_I,MPM_Mesh,FEM_Mesh,TimeStep); 
      R_I = compute_Reactions(FEM_Mesh, F_I);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Third step : Update kinetics ... WORKING",TimeStep);
      update_Kinetics(FEM_Mesh, Nodal_Kinetics, F_I, Params);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update particles lagrangian ... WORKING",TimeStep);
      update_Particles(MPM_Mesh, FEM_Mesh, Nodal_Kinetics, Params);
      LocalSearchGaussPoints(MPM_Mesh, FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Reset nodal values ... WORKING",TimeStep);
      free__MatrixLib__(F_I);
      print_Status("DONE !!!",TimeStep);

    }
  
}


/*******************************************************/

static void update_Kinetics(Mesh FEM_Mesh,
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
  int Ndim = NumberDimensions;
  double SizeTable = Ndim*sizeof(double *);
  double Mass_I;

  /* Time integration parameters */
  double alpha = Params.GA_alpha;
  double gamma = Params.GA_gamma;

  /* Asign forces to an auxiliar variable */
  Matrix F = Nodal_Forces;
  
  /* Nodal values the fields */
  Matrix Nodal_Acceleration_t0;
    /* MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Acceleration_t1;
    /* MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Velocity;
    /* MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */

  /* 1º Asign Kinetics values to matricial tables */
  for(int i = 0 ; i<Ndim ; i++){
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1+i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1+Ndim+i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+2*Ndim+i];
  }
  /* 2º Update the grid nodal variales */
  for(int i = 0 ; i<N_Nodes ; i++){
    Mass_I = Nodal_Kinetics.nM[0][i];
    if(Mass_I > 0){
      for(int j = 0 ; j<Ndim ; j++){
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

/*********************************************************************/


static Matrix compute_Nodal_Kinetics(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
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
  Matrix Nodal_Mass = allocZ__MatrixLib__(1,N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t0 = allocZ__MatrixLib__(N_Acceleration_dim,N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t1 = allocZ__MatrixLib__(N_Acceleration_dim,N_Nodes_Mesh);
  Matrix Nodal_Velocity = allocZ__MatrixLib__(N_Velocity_dim,N_Nodes_Mesh);
 
  /* GPs values of the fields */
  double Mass_GP; /* Mass of the GP */
  Matrix Vel_GP; /* Velocity of the GP */
  Matrix Acc_GP; /* Stress of the GP */
  
  /* Shape function auxiliar variables */
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_GP; /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double Mass_GP_I; /* Nodal contribution of each GP */  
  int GP_I;
  
  /* 1º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<N_GPs ; i++){

    /* 2º Define element of the GP */
    Nodes_p = get_particle_Set(i, MPM_Mesh.ListNodes[i], MPM_Mesh.NumberNodes[i]);
    
    /* 3º Evaluate the shape function in the coordinates of the GP */
    N_GP = compute_N__ShapeFun__(Nodes_p, MPM_Mesh,FEM_Mesh);
   
    /* 4º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Acc_GP.nV = MPM_Mesh.Phi.acc.nM[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];

    /* 5º Get the nodal kinetics (I) */
    for(int k = 0 ; k<Nodes_p.NumberNodes ; k++){
      /* Get the node for the GP */
      GP_I = Nodes_p.Connectivity[k];
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
    free__MatrixLib__(N_GP);
    free(Nodes_p.Connectivity);
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

  Matrix Nodal_Kinetics;
  /*   MatAssign(N_Kinetics_dim,N_Nodes_Mesh,NAN,NULL, */
  /* 	      (double**)malloc(SizeKinetics)); */
  /* strcpy(Nodal_Kinetics.Info, */
  /* 	 "MASS;ACCELERATION_t0;ACCELERATION_t1;VELOCITY"); */
  
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

static Matrix GetNodalVelocityDisplacement(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
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
  Matrix Nodal_Mass = allocZ__MatrixLib__(1,N_Nodes_Mesh);
  Matrix Nodal_Displacement = allocZ__MatrixLib__(Ndim,N_Nodes_Mesh);
  Matrix Nodal_Velocity = allocZ__MatrixLib__(Ndim,N_Nodes_Mesh);
  Matrix Nodal_Acceleration = allocZ__MatrixLib__(Ndim,N_Nodes_Mesh);
 
  /* GPs values of the fields */
  double Mass_GP; /* Mass of the GP */
  Matrix Disp_GP; /* Stress of the GP */
  Matrix Vel_GP; /* Velocity of the GP */
  Matrix Acel_GP; /* Velocity of the GP */
  
  /* Shape function auxiliar variables */
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_GP; /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double Mass_GP_I; /* Nodal contribution of each GP */  
  int GP_I;
  
  /* 1º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<N_GPs ; i++){

    /* 2º Define element of the GP */
    Nodes_p = get_particle_Set(i, MPM_Mesh.ListNodes[i], MPM_Mesh.NumberNodes[i]);
    
    /* 3º Evaluate the shape function in the coordinates of the GP */
    N_GP = compute_N__ShapeFun__(Nodes_p, MPM_Mesh,FEM_Mesh);
   
    /* 4º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];
    Disp_GP.nV = MPM_Mesh.Phi.dis.nM[i];
    Acel_GP.nV = MPM_Mesh.Phi.acc.nM[i];

    /* 5º Get the nodal kinetics (I) */
    for(int k = 0 ; k<Nodes_p.NumberNodes ; k++){

      /* Get the node for the GP */
      GP_I = Nodes_p.Connectivity[k];

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
    free__MatrixLib__(N_GP);
    free(Nodes_p.Connectivity);
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

  Matrix Nodal_Kinetics;
  /*   MatAssign(N_Kinetics_dim,N_Nodes_Mesh,NAN,NULL, */
  /* 	      (double**)malloc(SizeKinetics)); */
  /* strcpy(Nodal_Kinetics.Info, */
  /* 	 "MASS;DISPLACEMENT;VELOCITY;ACCELERATION"); */
  
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

static void update_Particles(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
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
    N_GP = compute_N__ShapeFun__(GP_Element, MPM_Mesh, FEM_Mesh);


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
    free__MatrixLib__(N_GP);
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
