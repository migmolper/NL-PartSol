#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Matrix compute_Nodal_Mass(GaussPoint, Mesh);
static Matrix compute_Velocity_Predictor(GaussPoint, Mesh, Matrix,Matrix,double,double);
static void   imposse_Velocity(Mesh, Matrix, int);
static Matrix compute_Velocity_Corrector(Mesh, Matrix, Matrix,Matrix,double,double);
static void   update_Particles(GaussPoint,Mesh,Matrix,Matrix,Matrix,double);
static void   update_LocalState(Matrix, GaussPoint,Mesh, double);
static Matrix compute_InternalForces(Matrix, GaussPoint,Mesh);
static Matrix compute_BodyForces(Matrix, GaussPoint, Mesh, int);
static Matrix compute_ContacForces(Matrix, GaussPoint, Mesh, int);
static Matrix compute_Reactions(Mesh, Matrix);

/**************************************************************/

void U_Newmark_Predictor_Corrector(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /*!
    Control parameters of the generalized-alpha algorithm 
    all the parameters are controled by a simple parameter :
    SpectralRadius 
  */
  double gamma = 0.5;

  /*!
    Auxiliar variable for the mass and momentum 
  */
  Matrix Lumped_Mass;
  Matrix Velocity;
  Matrix Forces;
  Matrix Reactions;

  
  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {
      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("First step : Predictor stage ... WORKING",TimeStep);
      Lumped_Mass = compute_Nodal_Mass(MPM_Mesh, FEM_Mesh);      
      Velocity = compute_Velocity_Predictor(MPM_Mesh,FEM_Mesh,Velocity, Lumped_Mass,
					    gamma, DeltaTimeStep);
      imposse_Velocity(FEM_Mesh,Velocity,TimeStep);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute equilibrium ... WORKING",TimeStep);
      update_LocalState(Velocity, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
      Forces = allocZ__MatrixLib__(Nnodes,Ndim);    
      Forces = compute_InternalForces(Forces, MPM_Mesh, FEM_Mesh);    
      Forces = compute_BodyForces(Forces, MPM_Mesh, FEM_Mesh, TimeStep);
      Forces = compute_ContacForces(Forces, MPM_Mesh, FEM_Mesh, TimeStep);
      Reactions = compute_Reactions(FEM_Mesh, Forces);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Third step : Corrector stage ... WORKING",TimeStep);
      Velocity = compute_Velocity_Corrector(FEM_Mesh,Velocity,Forces,
					    Lumped_Mass, gamma,
					    DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update particles lagrangian ... WORKING",TimeStep);
      update_Particles(MPM_Mesh, FEM_Mesh, Lumped_Mass, Velocity, Forces,DeltaTimeStep);
      local_search__Particles__(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
      	{
	       /*!
	           Print Nodal values after appling the BCCs
	       */
	 //      nodal_results_vtk__InOutFun__("Mesh",FEM_Mesh,Reactions,TimeStep,(int)TimeStep/ResultsTimeStep);
	       /*!
	         Print particle results 
	       */
	       particle_results_vtk__InOutFun__("MPM_VALUES",MPM_Mesh,"ALL",(int)TimeStep/ResultsTimeStep,ResultsTimeStep);
	     }

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Reset nodal values ... WORKING",TimeStep);
      free__MatrixLib__(Lumped_Mass);
      free__MatrixLib__(Velocity);
      free__MatrixLib__(Forces);
      free__MatrixLib__(Reactions);
      print_Status("DONE !!!",TimeStep);

    } 
}

/*******************************************************/

static Matrix compute_Nodal_Mass(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  int Np = MPM_Mesh.NumGP;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pI; /* Evaluation of the GP in the node */
  double m_p; /* Mass of the GP */ 
  Element Nodes_p; /* Element for each Gauss-Point */
  int Ip;

  /* Output */
  Matrix Lumped_Mass;

  /* 1º Allocate the output list of fields */
  Lumped_Mass = allocZ__MatrixLib__(Nnodes,1);
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int p = 0 ; p<Np ; p++){

    /* 3º Define element of the GP */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 4º Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
    /* 5º Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* 6º Get the nodal mass and mommentum */
    for(int k = 0 ; k<Nodes_p.NumberNodes ; k++){
      /* Get the node for the GP */
      Ip = Nodes_p.Connectivity[k];
      /* Evaluate the GP function in the node */
      ShapeFunction_pI = ShapeFunction_p.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(ShapeFunction_pI == 0){
	continue;
      }
      /* Nodal mass */
      Lumped_Mass.nV[Ip] += m_p*ShapeFunction_pI;
    }

    /* 7º Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }
 
  return Lumped_Mass;
  
}

/*******************************************************/

static Matrix compute_Velocity_Predictor(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
					 Matrix Velocity,Matrix Lumped_Mass,
					 double gamma,double DeltaTimeStep)
/*
  Get the nodal velocity using : 
  v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
  Initialize nodal velocities 
*/
{
  
  int Np = MPM_Mesh.NumGP;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;

  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pI; /* Evaluation of the GP in the node */
  double m_p; /* Mass of the GP */ 
  Element Nodes_p; /* Element for each Gauss-Point */
  int Ip;
  
  /* Matrix Vel_Mesh */
  Velocity = allocZ__MatrixLib__(Nnodes,Ndim);
  strcpy(Velocity.Info,"VELOCITY");
 
  /* 2º Iterate over the GP to get the nodal values */
  for(int p = 0 ; p<Np ; p++){
    
    /* 3º Define element of the GP */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 4º Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
    /* 5º Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* 6º Get the nodal mass and mommentum */
    for(int I = 0 ; I<Nodes_p.NumberNodes ; I++){
      /* Get the node for the GP */
      Ip = Nodes_p.Connectivity[I];
      /* Evaluate the GP function in the node */
      ShapeFunction_pI = ShapeFunction_p.nV[I];
      /* If this node has a null Value of the SHF continue */
      if(ShapeFunction_pI == 0){
	continue;
      }
      /* Nodal momentum */
      if(Lumped_Mass.nV[Ip] > 0){
	for(int i = 0 ; i<Ndim ; i++){
	  Velocity.nM[Ip][i] +=
	    m_p*ShapeFunction_pI*
	    (MPM_Mesh.Phi.vel.nM[p][i] +
	     MPM_Mesh.Phi.acc.nM[p][i]*(1-gamma)*DeltaTimeStep);
	}
      }
    }

    /* 7º Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /* 1º Get nodal values of the velocity */
  for(int I = 0 ; I<Nnodes ; I++){
    if(Lumped_Mass.nV[I] > 0){
      for(int i = 0 ; i<Ndim ; i++){
	Velocity.nM[I][i] = (double)Velocity.nM[I][i]/Lumped_Mass.nV[I];
      }
    }    
  }
  
  return Velocity;
}

/**********************************************************************/

static void imposse_Velocity(Mesh FEM_Mesh, Matrix V_I, int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for(int j = 0 ; j<NumNodesBound ; j++){
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for(int k = 0 ; k<NumDimBound ; k++){
	/* 7º Apply only if the direction is active (1) */
	if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1){
	  /* 8º Check if the curve it is on time */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num)){
	    printf("%s : %s \n",
		   "Error in imposse_NodalMomentum()",
		   "The time step is out of the curve !!");
	    exit(EXIT_FAILURE);
	  }
	  /* 9º Assign the boundary condition */
	  V_I.nM[Id_BCC][k] =
	    FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
	}
      }
    }    
  }

}

/*******************************************************/

static Matrix compute_Velocity_Corrector(Mesh FEM_Mesh,Matrix Velocity,
					 Matrix Forces,Matrix Lumped_Mass,
					 double gamma,double DeltaTimeStep)
/*
  Get the nodal velocity using : 
  v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
  Initialize nodal velocities 
*/  
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
 
  /* 1º Get nodal values of the velocity */
  for(int I = 0 ; I<Nnodes ; I++)
    {
      if((FEM_Mesh.NumParticles[I] > 0) && (Lumped_Mass.nV[I] > 0))
	{
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      Velocity.nM[I][i] += gamma*DeltaTimeStep*Forces.nM[I][i]/Lumped_Mass.nV[I];
	    }
	}
    }
  
  return Velocity;
}

/*******************************************************/

static void update_Particles(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
			     Matrix Lumped_Mass,Matrix Velocity,
			     Matrix Forces,double DeltaTimeStep)
{
  
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
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nnodes);

    /* 3º Evaluate shape function in the GP i */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

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
      mass_I = Lumped_Mass.nV[Ip];
      /* Update GP cuantities with nodal values */
      if(mass_I>0){
	for(int i = 0 ; i<Ndim ; i++){
	  /* Update the GP velocities */
	  MPM_Mesh.Phi.acc.nM[p][i] +=
	    ShapeFunction_pI*Forces.nM[Ip][i]/mass_I;
	  /* Update the GP velocities */
	  MPM_Mesh.Phi.vel.nM[p][i] +=
	    ShapeFunction_pI*DeltaTimeStep*Forces.nM[Ip][i]/mass_I;
	  /* Update the GP displacement */
	  MPM_Mesh.Phi.x_GC.nM[p][i] +=
	    ShapeFunction_pI*DeltaTimeStep*Velocity.nM[Ip][i] +
	    ShapeFunction_pI*0.5*pow(DeltaTimeStep,2)*Forces.nM[Ip][i]/mass_I;
	}
      } 
    }
    
    /* 5º Free memory */
    free(Nodes_p.Connectivity);
    free__MatrixLib__(ShapeFunction_p);
  }  
}

/*************************************************************/

static void update_LocalState(Matrix V_I, GaussPoint MPM_Mesh,Mesh FEM_Mesh, double TimeStep)
{
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix Gradient_p; /* Shape functions gradients */
  Matrix Nodal_Velocity_p; /* Velocity of the element nodes */
  Material Material_p; /* Properties of the Gauss-Point material */
  Tensor Rate_Strain_p; /* Increment of strain tensor */
  Tensor Strain_p; /*  Strain tensor */
  Tensor Stress_p; /* Stress tensor */
  double rho_p; /* Density of the Gauss-Point */
  int Np = MPM_Mesh.NumGP;
  int Nn;
  int Idx_Mat_p;

  /* Loop in the GPs */
  for(int p = 0 ; p<Np ; p++){

    /* Get the value of the density */
    rho_p = MPM_Mesh.Phi.rho.nV[p];

    /* Asign memory to tensors */
    Strain_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Strain.nM[p], 2);
    Stress_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nn = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

    /* Get the velocity of the nodes of the element */
    Nodal_Velocity_p = get_set_field_old__MeshTools__(V_I, Nodes_p);

    /* Compute gradient of the shape function in each node */
    Gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    
    /* Get the material properties */
    Idx_Mat_p = MPM_Mesh.MatIdx[p];
    Material_p = MPM_Mesh.Mat[Idx_Mat_p];

    /* Update Strain tensor */
    Rate_Strain_p = rate_inifinitesimal_Strain__Particles__(Nodal_Velocity_p,Gradient_p);
    Strain_p = infinitesimal_Strain__Particles__(Strain_p, Rate_Strain_p, TimeStep);

    /* Update density field */
    MPM_Mesh.Phi.rho.nV[p] = update_density__Particles__(rho_p, TimeStep, Rate_Strain_p);
    free__TensorLib__(Rate_Strain_p);

    /* Compute stress tensor */
    Stress_p = explicit_integration_stress__Particles__(Strain_p,Stress_p,Material_p);

    /* Compute deformation energy */
    MPM_Mesh.Phi.W.nV[p] = 0.5*inner_product__TensorLib__(Strain_p, Stress_p);
        
    /* Free the matrix with the nodal velocity of the element */
    free__MatrixLib__(Nodal_Velocity_p);
    
    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(Gradient_p);
    free(Nodes_p.Connectivity);
    
  }
  
  /* Loop in the particles to compute the damage */
  for(int p = 0 ; p<Np ; p++){    
    /* Compute damage of the particles */
    compute_particle_Damage(p, MPM_Mesh, FEM_Mesh);
  }
  
}

/*************************************************************/

static Matrix compute_InternalForces(Matrix F_I, GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{
  int Ndim = NumberDimensions;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix Gradient_p; /* Shape functions gradients */
  Tensor Stress_p; /* Stress tensor */
  Tensor Gradient_pI;
  Tensor InternalForcesDensity_Ip;
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V_p; /* Volume of the Gauss-Point */
  int Np = MPM_Mesh.NumGP;
  int Ip;
  int Nn;

  /* Loop in the GPs */
  for(int p = 0 ; p<Np ; p++){

    /* Get the value of the density */
    rho_p = MPM_Mesh.Phi.rho.nV[p];

    /* Get the value of the mass */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Asign memory to tensors */
    Stress_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nn = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

    /* Compute gradient of the shape function in each node */
    Gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Compute the volume of the Gauss-Point */
    V_p = m_p/rho_p;

    /* Compute nodal forces */
    for(int I = 0 ; I<Nn ; I++){
      /* Pass by reference the nodal gradient to the tensor */
      Gradient_pI = memory_to_tensor__TensorLib__(Gradient_p.nM[I], 1);
      
      /* Compute the nodal forces of the Gauss-Point */
      InternalForcesDensity_Ip =
	vector_linear_mapping__TensorLib__(Stress_p, Gradient_pI);
      
      /* Get the node of the mesh for the contribution */
      Ip = Nodes_p.Connectivity[I];
      
      /* Asign the nodal forces contribution to the node */
      for(int i = 0 ; i<Ndim ; i++){
	F_I.nM[Ip][i] -= InternalForcesDensity_Ip.n[i]*V_p;
      }

      /* Free the internal forces density */
      free__TensorLib__(InternalForcesDensity_Ip);
    }
        
    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(Gradient_p);
    free(Nodes_p.Connectivity);
  }

  return F_I;
  
}

/*********************************************************************/

static Matrix compute_BodyForces(Matrix F_I, GaussPoint MPM_Mesh,
				 Mesh FEM_Mesh, int TimeStep)
{
  int Ndim = NumberDimensions;
  Load * B = MPM_Mesh.B;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pI;
  Tensor b = alloc__TensorLib__(1); /* Body forces vector */
  double m_p; /* Mass of the Gauss-Point */
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  int NumNodesLoad;
  int p;
  int Ip;
  int Nn; /* Number of nodes of each Gauss-Point */

  for(int i = 0 ; i<NumBodyForces; i++){

    NumNodesLoad = B[i].NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++){

      /* Get the index of the Gauss-Point */
      p = B[i].Nodes[j];

      /* Get the value of the mass */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Define element for each GP */
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

      /* Compute shape functions */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Fill vector of body forces */
      for(int k = 0 ; k<Ndim ; k++){
	if(B[i].Dir[k]){
	  if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num)){
	    printf("%s : %s\n",
		   "Error in compute_BodyForces()",
		   "The time step is out of the curve !!");
	    exit(EXIT_FAILURE);
	  }
	  b.n[k] = B[i].Value[k].Fx[TimeStep];
	}
      }

      /* Get the node of the mesh for the contribution */
      for(int I = 0 ; I<Nn ; I++){

	/* Node for the contribution */
	Ip = Nodes_p.Connectivity[I];
	
	/* Pass the value of the nodal shape function to a scalar */
	ShapeFunction_pI = ShapeFunction_p.nV[I];
	
	/* Compute External forces */
	for(int k = 0 ; k<Ndim ; k++){
	  F_I.nM[Ip][k] += ShapeFunction_pI*b.n[k]*m_p;
	}
	
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
	
    }

      
  }

  free__TensorLib__(b);
  
  return F_I;

}

/*********************************************************************/

static Matrix compute_ContacForces(Matrix F_I, GaussPoint MPM_Mesh,
				   Mesh FEM_Mesh, int TimeStep)
{
  int Ndim = NumberDimensions;
  Load * F = MPM_Mesh.F;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pI;
  Tensor t = alloc__TensorLib__(1); /* Body forces vector */
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V_p; /* Volumen of the Gauss-Point */
  double thickness_p; /* Thickness of the Gauss-Point */
  int Mat_p; /* Index of tha material for each Gauss-Point */
  int NumContactForces = MPM_Mesh.NumNeumannBC;
  int NumNodesLoad;
  int p;
  int Ip;
  int Nn; /* Number of nodes of each Gauss-Point */

  for(int i = 0 ; i<NumContactForces; i++){

    NumNodesLoad = F[i].NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++){

      /* Get the index of the Gauss-Point */
      p = F[i].Nodes[j];

      /* Get the value of the mass */
      m_p = MPM_Mesh.Phi.mass.nV[p];
      
      /* Get the value of the density */
      rho_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the value of the volum */
      V_p = m_p/rho_p;

      /* Get the thickness of the material point */
      Mat_p = MPM_Mesh.MatIdx[p];
      thickness_p = MPM_Mesh.Mat[Mat_p].thickness;

      /* Define element for each GP */
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

      /* Compute shape functions */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Fill vector of body forces */
      for(int k = 0 ; k<Ndim ; k++){
	if(F[i].Dir[k]){
	  if( (TimeStep < 0) || (TimeStep > F[i].Value[k].Num)){
	    printf("%s : %s\n",
		   "Error in compute_ContacForces()",
		   "The time step is out of the curve !!");
	  }
	  t.n[k] = F[i].Value[k].Fx[TimeStep];
	}
      }

      /* Get the node of the mesh for the contribution */
      for(int I = 0 ; I<Nn ; I++){

	/* Node for the contribution */
	Ip = Nodes_p.Connectivity[I];
	
	/* Pass the value of the nodal shape function to a scalar */
	ShapeFunction_pI = ShapeFunction_p.nV[I];
	
	/* Compute Contact forces */
	for(int k = 0 ; k<Ndim ; k++){
	  F_I.nM[Ip][k] += ShapeFunction_pI*(t.n[k]/thickness_p)*V_p;
	}
	
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
	
    }

      
  }

  free__TensorLib__(t);
  
  return F_I;

}


/**********************************************************************/

static Matrix compute_Reactions(Mesh FEM_Mesh, Matrix F_I)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Ndim = NumberDimensions;

  Matrix R_I = allocZ__MatrixLib__(FEM_Mesh.NumNodesMesh,Ndim);
  strcpy(R_I.Info,"REACTIONS");

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for(int j = 0 ; j<NumNodesBound ; j++){
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for(int k = 0 ; k<NumDimBound ; k++){
	/* 7º Apply only if the direction is active (1) */
	if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1){
	  /* 8º Set to zero the forces in the nodes where velocity is fixed */
	  R_I.nM[Id_BCC][k] = F_I.nM[Id_BCC][k];
	  F_I.nM[Id_BCC][k] = 0;
	}
      }
    }    
  }

  return R_I;
}

/**********************************************************************/




