#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Matrix compute_Nodal_Mass(GaussPoint, Mesh);
static Matrix compute_Velocity_Predictor(GaussPoint, Mesh, Matrix,Matrix,double,double);
static void imposse_Velocity(Mesh, Matrix, int);
static Matrix compute_Velocity_Corrector(Mesh, Matrix, Matrix,Matrix,double,double);
static void update_Particles(GaussPoint,Mesh,Matrix,Matrix,Matrix,double);

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
      Forces = compute_equilibrium_U(Velocity,MPM_Mesh,FEM_Mesh,TimeStep); 
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
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
	{
	  /*!
	    Print Nodal values after appling the BCCs
	  */
	  WriteVtk_FEM("Mesh",FEM_Mesh,Reactions,TimeStep);
	  /*!
	    Print particle results 
	  */
	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,"ALL",TimeStep,ResultsTimeStep);
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
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 4º Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_N__ShapeFun__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
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
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 4º Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_N__ShapeFun__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
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
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nnodes);

    /* 3º Evaluate shape function in the GP i */
    ShapeFunction_p = compute_N__ShapeFun__(Nodes_p, MPM_Mesh, FEM_Mesh);

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

/*******************************************************/
