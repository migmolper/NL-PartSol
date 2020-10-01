#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

/*
  Auxiliar functions 
*/
static Matrix compute_Nodal_Lumped_Mass(GaussPoint, Mesh, Mask, double);
static Matrix compute_Nodal_Momentum(GaussPoint, Mesh, Mask);
static Matrix compute_Nodal_Velocity_Predicted(GaussPoint, Mesh, Mask, Matrix,double ,double);
static void   update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);
static Matrix compute_Nodal_Forces(Mask, GaussPoint, Mesh, int);
static void   compute_Nodal_Internal_Forces(Matrix,Mask,GaussPoint, Mesh);
static void   compute_Nodal_Body_Forces(Matrix,Mask, GaussPoint, Mesh, int);
static void   update_D_Displacement(Matrix, Matrix, Matrix, Matrix, double);
static void   imposed_displacements(Matrix, Mask, Mesh, int);
static Matrix compute_D_Velocity(Matrix, Matrix,double);
static void   update_Particles(Matrix, Matrix, GaussPoint, Mesh, Mask, double);

/**************************************************************/

void U_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /*
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Nactivenodes;
  /*
    Auxiliar variables for the solver
  */
  Matrix Lumped_Mass;
  Matrix Forces;
  Matrix Velocity;
  Mask ActiveNodes;
  double TOL = 0.000000000001;
  double epsilon = 1.0;
  double DeltaTimeStep;
  bool Convergence;
  int Iter = 0;
  int MaxIter = 500;



  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);


      print_Status("*************************************************",TimeStep);
      print_Status("First step : Generate Mask ... WORKING",TimeStep);
      /*
      	With the active set of nodes generate a mask to help the algorithm to compute
	      the equilibrium only in the active nodes
      */
      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute effective mass ... WORKING",TimeStep);
      /*
	     Compute the effective mass matrix as a convex combination of the consistent mass
	     matrix and the lumped mass matrix.
      */
      Lumped_Mass = compute_Nodal_Lumped_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);
      print_Status("DONE !!!",TimeStep);
      
      /*
        Compute the predicted nodal valocity.
      */      
      print_Status("*************************************************",TimeStep);
      print_Status("Third step : Compute predicted nodal velocity ... WORKING",TimeStep);
      Velocity = compute_Nodal_Velocity_Predicted(MPM_Mesh,FEM_Mesh,ActiveNodes
        Lumped_Mass,gamma,DeltaTimeStep);
      /* imposed_velocty(); */
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Compute equilibrium ... WORKING",TimeStep);
      
      /*
        Compute the stress-strain state for each particle
      */
      update_Local_State(Velocity, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
      /*
       Compute the nodal forces
      */
      Forces = compute_Nodal_Forces(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Compute velocity corrector ... WORKING",TimeStep);
      /*
	     Once the equilibrium is reached, obtain the increment of nodal velocity
      */
      D_Velocity = compute_Nodal_Velocity_Corrected(Velocity,D_Displacement,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);

      
      print_Status("*************************************************",TimeStep);
      print_Status("Six step : Update particles lagrangian ... WORKING",TimeStep);
      /*
	      Update Lagrangians with D_Displacement
      */
      update_Particles(D_Displacement,D_Velocity,MPM_Mesh,FEM_Mesh,ActiveNodes,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);
      /*
	     Reload the connectivity information for each particle
      */
      local_search__Particles__(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
      
      /*
	     Outputs
      */
      if(TimeStep % ResultsTimeStep == 0)
      {
	     /* WriteVtk_FEM("Mesh",FEM_Mesh,R_I,TimeStep); */
	     WriteVtk_MPM("MPM_VALUES",MPM_Mesh,"ALL",TimeStep,ResultsTimeStep);
      }

      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Reset nodal values ... WORKING",TimeStep);
      /*
	Free memory.
      */
      free__MatrixLib__(Lumped_Mass); 
      free__MatrixLib__(Velocity);
      free__MatrixLib__(Forces);
      free(ActiveNodes.Mask2Nodes);
      free(ActiveNodes.Nodes2Mask);
      
      print_Status("DONE !!!",TimeStep);

    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Lumped_Mass(GaussPoint MPM_Mesh,
					   Mesh FEM_Mesh,
					   Mask ActiveNodes,
					   double epsilon)
/*
  This function computes the lumped mass matrix.
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndof*Nnodes_mask;
  int Ap;
  int A_mask;
  int idx_AB_mask_i;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;  

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA, ShapeFunction_pB;
  /* Mass of the particle */
  double m_p;
  /* Nodal contribution A of the particle p */
  double m_A_p;
  /* Element for each particle */
  Element Nodes_p;


  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Nnodes_mask*Ndof, 1);

  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {

      /* Define tributary nodes of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
      /* Get the mass of the particle */
      m_p = MPM_Mesh.Phi.mass.nV[p];
      

      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
      {

	     /* Get the node in the mass matrix with the mask */
	     Ap = Nodes_p.Connectivity[A];
	     A_mask = ActiveNodes.Nodes2Mask[Ap];

	     /* Get the value of the shape function */
	     ShapeFunction_pA = ShapeFunction_p.nV[A];

	     /* Compute the nodal A contribution of the particle p */
    	 m_A_p = m_p*ShapeFunction_pA;


       /* Fill the Lumped mass matrix considering the number of dofs */
  	   for(int i = 0 ; i<Ndof ; i++)
	     {
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      Lumped_MassMatrix.nV[idx_A_mask_i] += m_A_p;	      
	     }
	  	  
	    }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);      
      
    }

  /* Add some usefulll info */
  strcpy(Lumped_MassMatrix.Info,"Lumped-Mass-Matrix");

  return Lumped_MassMatrix; 
}

/**************************************************************/

static Matrix compute_Nodal_Velocity_Predicted(GaussPoint MPM_Mesh, 
  Mesh FEM_Mesh, Mask ActiveNodes,Matrix Lumped_Mass,double gamma,double DeltaTimeStep)
/*
  Get the nodal velocity using : 
  v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
  Initialize nodal velocities 
  */
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

 /* Define and allocate the momentum vector */
  Matrix Velocity = allocZ__MatrixLib__(Ndim,Nnodes_mask);
    
  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {

      /* Define element of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Get the mass of the GP */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the nodal Velocity */
      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
      {
      /*
	     Get the node in the nodal Velocity with the mask
       */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];

        /* Evaluate the GP function in the node */
        ShapeFunction_pA = ShapeFunction_p.nV[A];

        /* Nodal momentum */
        for(int i = 0 ; i<Ndim ; i++)
        {
          idx_A_mask_i = A_mask + i*Nnodes_mask;

          Velocity.nV[idx_A_mask_i] += m_p*ShapeFunction_pA*(MPM_Mesh.Phi.vel.nM[p][i]+
            MPM_Mesh.Phi.acc.nM[p][i]*(1-gamma)*DeltaTimeStep);
        }
      }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*
    At this point the effective mass matrix coincides with the consistent mass
    matrix. We can tune it by a convecx combination with the lumped mass matrix
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      for(int i = 0 ; i<Ndof ; i++)
        {

        idx_A_mask_i = A + i*Nnodes_mask;

        Velocity.nV[idx_A_mask_i] = Velocity.nV[idx_A_mask_i]/Lumped_Mass.nV[idx_A_mask_i];
      }
    }

  /*
     Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");
 
  return Velocity;
}



/**************************************************************/

static void update_Local_State(Matrix D_Displacement,
                  			       Mask ActiveNodes,
			                         GaussPoint MPM_Mesh,
	                   		       Mesh FEM_Mesh,
                  			       double TimeStep)
{

  /*
    Auxiliar variables
   */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int MatIndx_p;
  int Nnodes_p;
  double J_p;  
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  Matrix Velocity_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor f_n1_p;
  Tensor S_p;
  
  /*
    Loop in the material point set 
  */
  for(int p = 0 ; p<Np ; p++)
    {
      /*
	       Define tributary nodes of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
      /*
      	Get the nodal increment of displacement using the mask
       */
      D_Displacement_Ap = get_set_field__MeshTools__(D_Displacement, Nodes_p, ActiveNodes);

      /*
      	 Evaluate the shape function gradient in the coordinates of the particle
       */
      gradient_p = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);
	  
      /*
      	Take the values of the deformation gradient from the previous step
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      
      /*
      	Compute the increment of the deformation gradient
      */
      f_n1_p = increment_Deformation_Gradient__Particles__(D_Displacement_Ap,gradient_p);

      /*
      	Update the deformation gradient in t = n + 1 with the information
	      from t = n and the increment of deformation gradient.
      */  
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, f_n1_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	     integration rule.
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      
      S_p = average_strain_integration_Stress__Particles__(S_p,F_n1_p,F_n_p,MatProp_p);
      
      /*
	     Free memory 
      */
      free__TensorLib__(f_n1_p);
      free__MatrixLib__(D_Displacement_Ap);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Forces(Mask ActiveNodes,
                        				   GaussPoint MPM_Mesh,
				                           Mesh FEM_Mesh,
                                   int TimeStep)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Ndim,Nnodes_mask);

  /*
    Add internal forces contribution
  */
  compute_Nodal_Internal_Forces(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh);

  /*
    Add body forces contribution
  */
  compute_Nodal_Body_Forces(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);
  
  return Forces;
}


/**************************************************************/

static void compute_Nodal_Internal_Forces(Matrix Forces,
					  Mask ActiveNodes,
					  GaussPoint MPM_Mesh,
					  Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;
  int idx_A_mask_i;

  Tensor P_p; /* First Piola-Kirchhoff Stress tensor */
  Tensor S_p; /* Second Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ap;

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n1_p;
  Tensor transpose_F_n1_p;
  
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V0_p; /* Volume of the Gauss-Point */
  double J_p; /* Jacobian of the deformation gradient */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /*
      	Get the value of the density 
      */
      rho_p = MPM_Mesh.Phi.rho.nV[p];

      /*
      	Get the value of the mass 
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /*
      	Define nodes for each particle
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

      /*
      	Compute gradient of the shape function in each node 
      */
      gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    	  
      /*
	      Take the values of the deformation gradient ant t = n and t = n + 1. 
	      Later compute the midpoint deformation gradient and 
	      the transpose of the deformation gradient.
      */
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);

      /*
	      Compute the first Piola-Kirchhoff stress tensor
      */
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);
      P_p = matrix_product__TensorLib__(F_n_p, S_p);

      /*
      	Compute the volume of the particle in the reference configuration 
      */
      J_p = I3__TensorLib__(F_n1_p);
      V0_p = (1/J_p)*(m_p/rho_p);
    
      for(int A = 0 ; A<NumNodes_p ; A++)
      {
      
  	  /*
	     Compute the gradient in the reference configuration 
	    */
	     gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
	     GRADIENT_pA = vector_linear_mapping__TensorLib__(transpose_F_n1_p,gradient_pA);
      
  	  /*
	     Compute the nodal forces of the particle 
	    */
	     InternalForcesDensity_Ap = vector_linear_mapping__TensorLib__(P_p, GRADIENT_pA);
      
  	  /*
	     Get the node of the mesh for the contribution 
	    */
	     Ap = Nodes_p.Connectivity[A];
	     A_mask = ActiveNodes.Nodes2Mask[Ap];
      
  	  /*
	     Asign the nodal forces contribution to the node 
	    */
	     for(int i = 0 ; i<Ndim ; i++)
       {
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      Forces.nV[idx_A_mask_i] += InternalForcesDensity_Ap.n[i]*V0_p;
	     }

  	  /*
	      Free memory 
	   */
	   free__TensorLib__(InternalForcesDensity_Ap);
	   free__TensorLib__(GRADIENT_pA);
	}
        
  /* 
	 Free memory 
  */
  free__TensorLib__(transpose_F_n1_p);
  free__TensorLib__(P_p);
  free__MatrixLib__(gradient_p);
  free(Nodes_p.Connectivity);
  }
   
}

/**************************************************************/

static void compute_Nodal_Body_Forces(Matrix Forces,
				      Mask ActiveNodes,
				      GaussPoint MPM_Mesh,
				      Mesh FEM_Mesh,
				      int TimeStep)
{
  /* Define auxilar variables */
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  int NumParticles_i; /* Number of particles with the */
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int idx_A_mask_k; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */
  int p; /* Particle index */
  
  double m_p; /* Mass of the particle */
  Load * B = MPM_Mesh.B; /* List with the load cases */
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA; /* Evaluation in the node I for the particle p */
  
  Tensor b = alloc__TensorLib__(1); /* Body forces vector */
  
  for(int i = 0 ; i<NumBodyForces ; i++)
    {
      /* Get the number of particles with the body load i */
      NumParticles_i = B[i].NumNodes;

      for(int j = 0 ; j<NumParticles_i ; j++){

	/* Get the index of the Gauss-Point */
	p = B[i].Nodes[j];
	
	/* Get the value of the mass */
	m_p = MPM_Mesh.Phi.mass.nV[p];

	/* Define tributary nodes of the particle */
	NumNodes_p = MPM_Mesh.NumberNodes[p];
	Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

	/* Compute shape functions */
	ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);


	/* Fill vector of body forces */
	for(int k = 0 ; k<Ndim ; k++){
	  if(B[i].Dir[k]){
	    if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num)){
	      printf("%s : %s\n",
		     "Error in compute_Nodal_Body_Forces()",
		     "The time step is out of the curve !!");
	      exit(EXIT_FAILURE);
	    }
	    b.n[k] = B[i].Value[k].Fx[TimeStep];
	  }
	}


	/* Get the node of the mesh for the contribution */
	for(int A = 0 ; A<NumNodes_p ; A++)
	  {
	    
	    /* Pass the value of the nodal shape function to a scalar */
	    ShapeFunction_pA = ShapeFunction_p.nV[A];
	    
	    /* Get the node of the mesh for the contribution */
	    Ap = Nodes_p.Connectivity[A];
	    A_mask = ActiveNodes.Nodes2Mask[Ap];
	    
	    /* Compute body forces */
	    for(int k = 0 ; k<Ndim ; k++)
	      {
		idx_A_mask_k = A_mask + k*Nnodes_mask;
		Forces.nV[idx_A_mask_k] += ShapeFunction_pA*b.n[k]*m_p;
	      } 
	    
	  }
	
	/* Free the matrix with the nodal gradient of the element */
	free__MatrixLib__(ShapeFunction_p);
	free(Nodes_p.Connectivity);
	
      }
      
    }

  /*
    Free auxiliar tensor
   */
  free__TensorLib__(b);

  
}

/**************************************************************/

static void imposed_displacements(Matrix D_Displacement,
				  Mask ActiveNodes,
				  Mesh FEM_Mesh,
				  int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* Define auxilar variables */
  int NumBoundaryConditions = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int A_BCC; /* Index of the node where we apply the BCC */
  int A_mask_BCC; /* Index of the node where we apply the BCC */
  
  for(int i_boundary = 0 ; i_boundary<NumBoundaryConditions ; i_boundary++)
    {
      
      /*
	Get the number of nodes of this boundary and 
	the number of dimensions where it is applied
      */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i_boundary].NumNodes;
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i_boundary].Dim;
    
    for(int A = 0 ; A<NumNodesBound ; A++)
      {
	
      /* 
	 Get the index of the node and get the mask node value 
      */
      A_BCC = FEM_Mesh.Bounds.BCC_i[i_boundary].Nodes[A];
      A_mask_BCC = ActiveNodes.Nodes2Mask[A_BCC];

      for(int i_dim = 0 ; i_dim<NumDimBound ; i_dim++)
	{
	  
	if(FEM_Mesh.Bounds.BCC_i[i_boundary].Dir[i_dim] == 1)
	  {
	  /*
	    Check if the curve it is on time 
	  */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i_boundary].Value[i_dim].Num))
	    {
	      printf("%s : %s \n",
		     "Error in imposed_displacements()",
		     "The time step is out of the curve !!");
	      exit(EXIT_FAILURE);
	    }
	  /* 
	     Assign the boundary condition 
	  */
	  D_Displacement.nM[i_dim][A_mask_BCC] =
	    FEM_Mesh.Bounds.BCC_i[i_boundary].Value[i_dim].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i_boundary].Dir[i_dim];
	  
	  }
	}
      }    
    }
  
}

/**************************************************************/


static Matrix compute_Nodal_Velocity_Predicted(Matrix Velocity,
				 Matrix D_Displacement,
				 double Dt)
{
  int Nnodes_mask = Velocity.N_cols;
  int Ndim = NumberDimensions;
  
  Matrix D_Velocity = allocZ__MatrixLib__(Ndim,Nnodes_mask);

  /*
    Compute the velocity in the midd-point 
   */
  for(int idx_A_i = 0 ; idx_A_i < Nnodes_mask*Ndim ; idx_A_i++)
    {	
      D_Velocity.nV[idx_A_i] = 2*(D_Displacement.nV[idx_A_i]/Dt - Velocity.nV[idx_A_i]);
    }


  return D_Velocity;
}


/**************************************************************/

static void update_Particles(Matrix D_Displacement,
			     Matrix D_Velocity,
			     GaussPoint MPM_Mesh,
			     Mesh FEM_Mesh,
			     Mask ActiveNodes,
			     double DeltaTimeStep)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  int idx_ij;
  Matrix D_Displacement_Ap;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  Matrix gradient_p;
  double ShapeFunction_pI; /* Nodal value for the particle */
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor f_n1_p;
  double J_p;
  Element Nodes_p; /* Element for each particle */

  /* iterate over the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /* Define element of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
      /*
	Get the nodal increment of displacement using the mask
      */
      D_Displacement_Ap = get_set_field__MeshTools__(D_Displacement, Nodes_p, ActiveNodes);
      
      /*
	Evaluate the shape function and gradient in the coordinates of the particle 
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      gradient_p      = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);

      /*
      	Take the values of the deformation gradient from the previous step
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

      /*
      	Compute the increment of the deformation gradient
      */
      f_n1_p = increment_Deformation_Gradient__Particles__(D_Displacement_Ap,gradient_p);

      /*
	Update the deformation gradient in t = n + 1
      */  
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, f_n1_p);

      /*
	Update density with the jacobian of the increment deformation gradient
       */
      J_p = I3__TensorLib__(f_n1_p);
      MPM_Mesh.Phi.rho.nV[p] = MPM_Mesh.Phi.rho.nV[p]/J_p;

      /*
	Replace the deformation gradient at t = n with the converged deformation gradient
      */
      for(int i = 0 ; i<Ndim  ; i++)
	{
	  for(int j = 0 ; j<Ndim  ; j++)
	    {
	      F_n_p.N[i][j] = F_n1_p.N[i][j];
	    }
	}      
      
      /* Iterate over the nodes of the particle */
      for(int A = 0; A<Nodes_p.NumberNodes; A++)
	{

	  /*
	    Get the node in the nodal momentum with the mask
	  */
	  Ap = Nodes_p.Connectivity[A];
	  A_mask = ActiveNodes.Nodes2Mask[Ap];
	  
	  /*
	    Evaluate the GP function in the node 
	  */
	  ShapeFunction_pI = ShapeFunction_p.nV[A];
	  
	  /*
	    Update velocity position and deformation gradient of the particles
	   */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      MPM_Mesh.Phi.vel.nM[p][i]  += ShapeFunction_pI*D_Velocity.nV[idx_A_mask_i];
	      MPM_Mesh.Phi.x_GC.nM[p][i] += ShapeFunction_pI*D_Displacement.nV[idx_A_mask_i];
	    } 
	}

      /*
	Free memory
       */
      free(Nodes_p.Connectivity);
      free__MatrixLib__(ShapeFunction_p);
      free__MatrixLib__(gradient_p);
      free__TensorLib__(f_n1_p);
    }  
}

/**************************************************************/
