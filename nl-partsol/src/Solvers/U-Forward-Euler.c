#include "nl-partsol.h"

/*
  Call global variables
*/
double DeltaTimeStep;

/*
  Auxiliar functions 
*/
static void   update_Particles(GaussPoint, Mesh, Matrix, Matrix, double);
static Matrix compute_NodalMomentumMass(GaussPoint, Mesh);
static void   imposed_Momentum(Mesh, Matrix, int);
static Matrix compute_Nodal_Velocity(Mesh, Matrix);
static void   update_Nodal_Momentum(Mesh, Matrix, Matrix);

/**************************************************************/

void U_Forward_Euler(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /*!
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /*!
    Auxiliar variable for the mass and momentum 
  */
  Matrix Phi_I;
  strcpy(Phi_I.Info,"MOMENTUM;MASS");
  Matrix V_I;
  Matrix F_I;
  Matrix R_I;
  

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status(" First step : Get the nodal fields ... WORKING",TimeStep);
      Phi_I = compute_NodalMomentumMass(MPM_Mesh,FEM_Mesh);
      imposed_Momentum(FEM_Mesh,Phi_I,TimeStep);    
      V_I = compute_Nodal_Velocity(FEM_Mesh, Phi_I);
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
	{
	  /*!
	    Print Nodal values after appling the BCCs
	  */
	  WriteVtk_FEM("Mesh",FEM_Mesh,Phi_I,
		       (int)TimeStep/ResultsTimeStep);
	  /*!
	    Print particle results 
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
      print_Status(" Third step : Update nodal momentum ... WORKING",TimeStep);
      update_Nodal_Momentum(FEM_Mesh,Phi_I,F_I);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Four step : Update lagrangian ... WORKING",TimeStep);
      update_Particles(MPM_Mesh, FEM_Mesh, Phi_I, F_I, DeltaTimeStep);
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Five step : Reset nodal values ... WORKING",TimeStep);
      free__MatrixLib__(Phi_I);
      free__MatrixLib__(V_I);
      free__MatrixLib__(F_I);
      free__MatrixLib__(R_I);
      print_Status("DONE !!!",TimeStep);

    } 

}

/*******************************************************/

static Matrix compute_NodalMomentumMass(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*
  This function performs a information trasference from the lagrangian particles to 
  the nodes of the eulerian mesh.
  Output : Phi_I = {P_I | M_I}
*/
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
  Phi_I = allocZ__MatrixLib__(Nnodes, Ndim + 1);
  
  /* Iterate over the GP to get the nodal values */
  for(int p = 0 ; p<Np ; p++){

    /* Define element of the GP */
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_N__ShapeFun__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
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
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }
 
  return Phi_I;
  
}

/**********************************************************************/

static void imposed_Momentum(Mesh FEM_Mesh, Matrix Phi_I, int TimeStep)
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
		   "Error in imposed_Momentum()",
		   "The time step is out of the curve !!");
	    exit(EXIT_FAILURE);
	  }
	  /* 9º Assign the boundary condition */
	  Phi_I.nM[Id_BCC][k] =
	    FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
	}
      }
    }    
  }

}

/**************************************************************/

static void update_Particles(GaussPoint MPM_Mesh, Mesh FEM_Mesh,
				Matrix Phi_I, Matrix F_I, double Dt)
{

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
    N_p = compute_N__ShapeFun__(Nodes_p, MPM_Mesh, FEM_Mesh);
    
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
    free__MatrixLib__(N_p);
  }  
}

/*******************************************************/

static Matrix compute_Nodal_Velocity(Mesh FEM_Mesh, Matrix Phi_I)
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
  Matrix V_I = allocZ__MatrixLib__(Nnodes,Ndim);
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

static void update_Nodal_Momentum(Mesh FEM_Mesh, Matrix Phi_I, Matrix F_I)
/*!
 * \brief Brief description of UpdateGridNodalMomentum.
 *        Compute the nodal contribution of each GP to the total forces.
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param Phi_I {P_I | M_I}
 *  @param F_I : Nodal value of the total forces.
 *
 */
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  
  /* Update the grid nodal momentum */
  for(int i = 0 ; i<Nnodes ; i++){
    for(int j = 0 ; j<Ndim ; j++){
      if(FEM_Mesh.NumParticles[i] > 0){
	Phi_I.nM[i][j] += DeltaTimeStep*F_I.nM[i][j];
      }
    }
  }  
}

/**************************************************************/


