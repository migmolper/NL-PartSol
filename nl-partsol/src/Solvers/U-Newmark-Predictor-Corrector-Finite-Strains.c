#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

/*
  Call global variables
*/
Event * Out_nodal_path_csv;
Event * Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

/*
  Auxiliar functions 
*/
static Matrix compute_Nodal_Lumped_Mass(Particle,Mesh,Mask);
static Matrix compute_Nodal_Momentum(Particle,Mesh,Mask);
static Matrix compute_Nodal_Velocity_Predicted(Particle,Mesh,Mask,Matrix,double,double);
static void   imposse_Velocity(Mesh,Matrix,Mask,int);
static Matrix compute_Nodal_D_Displacement(Matrix,Mask,double);
static void   update_Local_State(Matrix,Mask,Particle,Mesh,double);
static Matrix compute_Nodal_Forces(Mask,Particle,Mesh,int);
static void   compute_Nodal_Internal_Forces(Matrix,Mask,Particle,Mesh);
static void   compute_Nodal_Body_Forces(Matrix,Mask,Particle,Mesh,int);
static Matrix compute_Reactions(Mesh,Matrix,Mask);
static void   compute_Nodal_Velocity_Corrected(Matrix,Matrix,Matrix,double,double);
static void   update_Particles(Particle,Mesh,Matrix,Matrix,Matrix,Mask,double);
static void   output_selector(Particle, Mesh, Mask, Matrix, Matrix, Matrix, Matrix, int, int);

/**************************************************************/

void U_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, int InitialStep)
{

  /*
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Nactivenodes;

  /*!
    Control parameters of the generalized-alpha algorithm 
    all the parameters are controled by a simple parameter :
    SpectralRadius 
  */
  double gamma = 0.5;
  double DeltaTimeStep;

  /*
    Auxiliar variables for the solver
  */
  Matrix Lumped_Mass;
  Matrix Velocity;
  Matrix D_Displacement;
  Matrix Forces;
  Matrix Reactions;
  Mask ActiveNodes;


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
      Lumped_Mass = compute_Nodal_Lumped_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes);
      print_Status("DONE !!!",TimeStep);
      
      /*
        Compute the predicted nodal valocity.
      */      
      print_Status("*************************************************",TimeStep);
      print_Status("Third step : Compute predicted nodal velocity ... WORKING",TimeStep);
      Velocity = compute_Nodal_Velocity_Predicted(MPM_Mesh,FEM_Mesh,ActiveNodes,
                                                  Lumped_Mass,gamma,DeltaTimeStep);
      /*
        Imposse velocity values in the boundary conditions nodes.
      */
      imposse_Velocity(FEM_Mesh,Velocity,ActiveNodes,TimeStep);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Compute equilibrium ... WORKING",TimeStep);
      /*
        Compute the nodal displacement
      */
      D_Displacement = compute_Nodal_D_Displacement(Velocity, ActiveNodes, DeltaTimeStep);
      /*
        Compute the stress-strain state for each particle
      */
      update_Local_State(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
      /*
	Compute the nodal forces
      */
      Forces = compute_Nodal_Forces(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);
      /*
        Compute reactions
      */
      Reactions = compute_Reactions(FEM_Mesh,Forces,ActiveNodes);

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Compute velocity corrector ... WORKING",TimeStep);
      /*
	Correct the predicted velocity
      */
      compute_Nodal_Velocity_Corrected(Velocity,Forces,Lumped_Mass,gamma,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Six step : Update particles lagrangian ... WORKING",TimeStep);
      /*
	Update Lagrangians with D_Displacement
      */
      update_Particles(MPM_Mesh,FEM_Mesh,Lumped_Mass,Forces,Velocity,ActiveNodes,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);
      /*
	Reload the connectivity information for each particle
      */
      local_search__Particles__(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
      
      /*
	       Outputs
      */
      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, Velocity, D_Displacement,
                      Forces, Reactions, TimeStep, ResultsTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Reset nodal values ... WORKING",TimeStep);
      /*
      	Free memory.
      */
      free__MatrixLib__(Lumped_Mass); 
      free__MatrixLib__(Velocity);
      free__MatrixLib__(D_Displacement);
      free__MatrixLib__(Forces);
      free__MatrixLib__(Reactions);
      free(ActiveNodes.Nodes2Mask);
      
      print_Status("DONE !!!",TimeStep);

    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Lumped_Mass(Particle MPM_Mesh,
					Mesh FEM_Mesh,
					Mask ActiveNodes)
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
	      Lumped_MassMatrix.nV[A_mask*Ndof + i] += m_A_p;	      
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

static Matrix compute_Nodal_Velocity_Predicted(Particle MPM_Mesh, 
                                               Mesh FEM_Mesh,
                                               Mask ActiveNodes,
                                               Matrix Lumped_Mass,
                                               double gamma,
                                               double DeltaTimeStep)
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
  Matrix Velocity = allocZ__MatrixLib__(Nnodes_mask,Ndim);
    
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
	      idx_A_mask_i = A_mask*Ndim + i;
	      Velocity.nV[idx_A_mask_i] += m_p*ShapeFunction_pA*(MPM_Mesh.Phi.vel.nM[p][i] + MPM_Mesh.Phi.acc.nM[p][i]*(1-gamma)*DeltaTimeStep);
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
      for(int i = 0 ; i<Ndim ; i++)
        {
          idx_A_mask_i = A*Ndim + i;
          Velocity.nV[idx_A_mask_i] = Velocity.nV[idx_A_mask_i]/Lumped_Mass.nV[idx_A_mask_i];
	}
    }

  /*
    Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");
 
  return Velocity;
}

/**********************************************************************/

static void imposse_Velocity(Mesh FEM_Mesh,
                             Matrix Velocity,
                             Mask ActiveNodes,
                             int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++)
    {

      /* 
	 Get the number of nodes of this boundarie 
      */
      NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

      /* 
	 Get the number of dimensions where the BCC it is applied 
      */
      NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

      for(int j = 0 ; j<NumNodesBound ; j++)
        {
          /* 
	     Get the index of the node 
          */
          Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
          Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

          /*
            The boundary condition is not affecting any active node,
            continue interating
          */
          if(Id_BCC_mask == -1)
	    {
	      continue;
	    }

          /* 
	     Loop over the dimensions of the boundary condition 
          */
          for(int k = 0 ; k<NumDimBound ; k++)
            {

              /* 
		 Apply only if the direction is active (1) 
              */
              if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
                {
    
                  /* 
		     Check if the curve it is on time 
                  */
                  if( (TimeStep < 0) ||
		      (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num))
                    {
                      printf("%s : %s \n",
                             "Error in imposse_NodalMomentum()",
                             "The time step is out of the curve !!");
                      exit(EXIT_FAILURE);
                    }

                  /* 
		     Assign the boundary condition 
                  */
                  Id_BCC_mask_k = Id_BCC_mask*NumDimBound + k; 
                  Velocity.nV[Id_BCC_mask_k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
		    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
                }
            }
        }    
    }

}

/**************************************************************/

static Matrix compute_Nodal_D_Displacement(Matrix Velocity,
                                           Mask ActiveNodes,
                                           double Dt)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = Velocity.N_rows;
  int Order = Ndim*Nnodes_mask;
  Matrix D_Displacement = allocZ__MatrixLib__(Nnodes_mask,Ndim);

  /*
    Compute nodal displacement
  */
  for(int idx_A = 0 ; idx_A<Order ; idx_A++)
    {
      D_Displacement.nV[idx_A] = Dt*Velocity.nV[idx_A];
    }

  return D_Displacement;
}

/**************************************************************/

static void update_Local_State(Matrix D_Displacement,
			       Mask ActiveNodes,
			       Particle MPM_Mesh,
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
  double Vol_0_p;
  double rho_n_p;
  double J_n1_p;  
  double Delta_J_p;
  Plastic_status Input_Plastic_Parameters;
  Plastic_status Output_Plastic_Parameters;
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;
  Tensor F_plastic_p;
  Tensor P_p;
  
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
      DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);

      
      /*
      	Compute the increment of the deformation gradient
      */
      update_increment_Deformation_Gradient__Particles__(DF_p,D_Displacement_Ap,gradient_p);

      /*
      	Update the deformation gradient in t = n + 1 with the information
	from t = n and the increment of deformation gradient.
      */  
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	integration rule.
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      

      if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
      {
          P_p = compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(P_p, F_n1_p, MatProp_p);
      }
      else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
      {
          J_n1_p = I3__TensorLib__(F_n1_p);
          P_p = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(P_p, F_n1_p, J_n1_p, MatProp_p);
      }
      else if(strcmp(MatProp_p.Type,"Von-Mises") == 0)
      {
          J_n1_p = I3__TensorLib__(F_n1_p);
          F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
          Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
          Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

          /* Run the plastic solver */
          Output_Plastic_Parameters = finite_strains_plasticity_Von_Mises(P_p,F_plastic_p,F_n1_p,Input_Plastic_Parameters,MatProp_p,J_n1_p);

          /* Update variables (cohesion and EPS) */
          MPM_Mesh.Phi.cohesion.nV[p] = Output_Plastic_Parameters.Yield_stress;
          MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;
      }
      else if((strcmp(MatProp_p.Type,"Drucker-Prager-Plane-Strain") == 0) || 
              (strcmp(MatProp_p.Type,"Drucker-Prager-Outer-Cone") == 0))
      {
          J_n1_p = I3__TensorLib__(F_n1_p);
          F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
          Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
          Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

          /* Run the plastic solver */
          Output_Plastic_Parameters = finite_strains_plasticity_Drucker_Prager_Sanavia(P_p,F_plastic_p,F_n1_p,Input_Plastic_Parameters,MatProp_p,J_n1_p);

          /* Update variables (cohesion and EPS) */
          MPM_Mesh.Phi.cohesion.nV[p] = Output_Plastic_Parameters.Cohesion;
          MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;

        }
      else
        {
          fprintf(stderr,"%s : %s %s %s \n",
		  "Error in update_Local_State()",
		  "The material",MatProp_p.Type,"has not been yet implemnented");
          exit(EXIT_FAILURE);
        }
      

      /*
        Update density with the jacobian of the increment deformation gradient
      */
      Delta_J_p = I3__TensorLib__(DF_p);
      rho_n_p = MPM_Mesh.Phi.rho.nV[p];
      MPM_Mesh.Phi.rho.nV[p] = rho_n_p/Delta_J_p;

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

      /* Compute the deformation energy */
  //    Vol_0_p = MPM_Mesh.Phi.Vol_0.nV[p];
  //    MPM_Mesh.Phi.W.nV[p]= finite_strains_internal_energy__Particles__(F_n_p, MatProp_p, Vol_0_p);

      /*
	Free memory 
      */
      free__MatrixLib__(D_Displacement_Ap);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Forces(Mask ActiveNodes,
				   Particle MPM_Mesh,
				   Mesh FEM_Mesh,
                                   int TimeStep)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask,Ndim);

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
					  Particle MPM_Mesh,
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
      transpose_F_n1_p = transpose__TensorLib__(F_n1_p);

      /*
	Compute the first Piola-Kirchhoff stress tensor
      */
      P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

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
	      idx_A_mask_i = A_mask*Ndim + i;
	      Forces.nV[idx_A_mask_i] -= InternalForcesDensity_Ap.n[i]*V0_p;
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
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
    }
   
}

/**************************************************************/

static void compute_Nodal_Body_Forces(Matrix Forces,
				      Mask ActiveNodes,
				      Particle MPM_Mesh,
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
		idx_A_mask_k = A_mask*Ndim + k;
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

/**********************************************************************/

static Matrix compute_Reactions(Mesh FEM_Mesh, Matrix Forces, Mask ActiveNodes)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int Ndim = NumberDimensions;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  Matrix Reactions = allocZ__MatrixLib__(Nnodes_mask ,Ndim);
  strcpy(Reactions.Info,"REACTIONS");

  /*
    Loop over the the boundaries 
  */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++)
    {
      /* 
	 Get the number of nodes of this boundarie 
      */
      NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    
      /* 
	 Get the number of dimensions where the BCC it is applied 
      */
      NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    
      for(int j = 0 ; j<NumNodesBound ; j++)
	{
	  /* 
	     Get the index of the node 
	  */
	  Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
	  Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

	  /*
	    The boundary condition is not affecting any active node,
	    continue interating
	  */
	  if(Id_BCC_mask == -1)
	    {
	      continue;
	    }
      
	  /* 
	     Loop over the dimensions of the boundary condition 
	  */
	  for(int k = 0 ; k<NumDimBound ; k++)
	    {

	      /* 
		 Apply only if the direction is active (1) 
	      */
	      if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
		{
		  /* 
		     Set to zero the forces in the nodes where velocity is fixed 
		  */
		  Id_BCC_mask_k = Id_BCC_mask*NumDimBound + k; 
		  Reactions.nV[Id_BCC_mask_k] = Forces.nV[Id_BCC_mask_k];
		  Forces.nV[Id_BCC_mask_k] = 0;
		}
	    }
	}    
    }

  return Reactions;
}

/**************************************************************/


static void compute_Nodal_Velocity_Corrected(Matrix Velocity,
                                             Matrix Forces,
					     Matrix Lumped_Mass,
                                             double gamma,
                                             double Dt)
{
  int Nnodes_mask = Velocity.N_rows;
  int Ndim = NumberDimensions;
  int Order = Nnodes_mask*Ndim;
  
  /*
    Correct the nodal velocity field
  */
  for(int idx_A_i = 0 ; idx_A_i <Order ; idx_A_i++)
    {	
      Velocity.nV[idx_A_i] += gamma*Dt*Forces.nV[idx_A_i]/Lumped_Mass.nV[idx_A_i];
    }

}


/**************************************************************/

static void update_Particles(Particle MPM_Mesh,
                             Mesh FEM_Mesh,
                             Matrix Lumped_Mass,
                             Matrix Forces,
                             Matrix Velocity,
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
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  Matrix gradient_p;
  double ShapeFunction_pI; /* Nodal value for the particle */
  double mass_I;
  double D_U_pI; /* Increment of displacement */
  Element Nodes_p; /* Element for each particle */

  /* iterate over the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /*
        Define element of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /*
	Evaluate the shape function and gradient in the coordinates of the particle 
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      gradient_p      = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);
      
      /* 
	 Iterate over the nodes of the particle 
      */
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
	      idx_A_mask_i = A_mask*Ndim + i;

	      /* Get the nodal mass */
	      mass_I = Lumped_Mass.nV[idx_A_mask_i];
	      /* Update the particles accelerations */
	      MPM_Mesh.Phi.acc.nM[p][i] += ShapeFunction_pI*Forces.nV[idx_A_mask_i]/mass_I;
	      /* Update the particles velocities */
	      MPM_Mesh.Phi.vel.nM[p][i] += ShapeFunction_pI*DeltaTimeStep*Forces.nV[idx_A_mask_i]/mass_I;

        /* Compute the nodal contribution of the increment of displacement */
        D_U_pI = ShapeFunction_pI*DeltaTimeStep*Velocity.nV[idx_A_mask_i] + 
                 ShapeFunction_pI*0.5*pow(DeltaTimeStep,2)*Forces.nV[idx_A_mask_i]/mass_I;

        /* Update the particle displacement */
        MPM_Mesh.Phi.dis.nM[p][i] += D_U_pI;
	      /* Update the particles position */
	      MPM_Mesh.Phi.x_GC.nM[p][i] += D_U_pI;
	    } 
    	}

      /*
	Free memory
      */
      free(Nodes_p.Connectivity);
      free__MatrixLib__(ShapeFunction_p);
      free__MatrixLib__(gradient_p);
    }  
}

/**************************************************************/

static void output_selector(Particle MPM_Mesh, Mesh FEM_Mesh, Mask ActiveNodes,
                            Matrix Velocity, Matrix D_Displacement, Matrix Forces,
                            Matrix Reactions, int TimeStep, int ResultsTimeStep)
{

  /*
    vtk results
  */
  if(TimeStep % ResultsTimeStep == 0)
  {
    particle_results_vtk__InOutFun__(MPM_Mesh,TimeStep,ResultsTimeStep);

    nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes, Reactions, TimeStep, ResultsTimeStep);
  }

  /* 
    csv results 
  */
  for(int i = 0 ; i<Number_Out_nodal_path_csv ; i++)
  {

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity)
    {
      path_nodes_analysis_csv__InOutFun__(Velocity, FEM_Mesh.Coordinates,"Nodal_path_velocity_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_D_Displacement)
    {
      path_nodes_analysis_csv__InOutFun__(D_Displacement, FEM_Mesh.Coordinates,"Nodal_path_displacement_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Forces)
    {
      path_nodes_analysis_csv__InOutFun__(Forces, FEM_Mesh.Coordinates,"Nodal_path_forces_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Reactions)
    {
      path_nodes_analysis_csv__InOutFun__(Reactions, FEM_Mesh.Coordinates,"Nodal_path_reactions_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

  }


  for(int i = 0 ; i<Number_Out_particles_path_csv ; i++)
  {
    if(Out_particles_path_csv[i].Out_csv_particles_path_Damage)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.chi, MPM_Mesh.Phi.x_GC, "Particles_path_damage_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Velocity)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.vel, MPM_Mesh.Phi.x_GC, "Particles_path_velocity_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Acceleration)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.acc, MPM_Mesh.Phi.x_GC, "Particles_path_acceleration_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Displacement)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.dis, MPM_Mesh.Phi.x_GC, "Particles_path_displacement_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Stress)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Stress, MPM_Mesh.Phi.x_GC, "Particles_path_stress_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Strain)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Strain, MPM_Mesh.Phi.x_GC, "Particles_path_strain_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Deformation_gradient)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.F_n, MPM_Mesh.Phi.x_GC, "Particles_path_deformation_gradient_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

  }

}

/**************************************************************/
