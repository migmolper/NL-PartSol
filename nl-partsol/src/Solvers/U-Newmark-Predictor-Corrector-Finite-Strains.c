#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif


/*
  Call global variables
*/
double Thickness_Plain_Stress;
Event * Out_nodal_path_csv;
Event * Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

/*
  Auxiliar functions 
*/

/* Step 1 */
static Matrix compute_Mass_Matrix(Particle MPM_Mesh, Mesh FEM_Mesh, Mask ActiveNodes);
/* Step 2 */
static  void  compute_Explicit_Newmark_Predictor(Particle,double,double);
/* Step 3 */
static Matrix compute_Nodal_Gravity_field(Mask, Particle, int);
static Matrix compute_Nodal_D_Displacement(Particle,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Velocity(Particle,Mesh,Mask,Matrix);
static void   impose_Dirichlet_Boundary_Conditions(Mesh,Matrix,Matrix,Mask,int);
/* Step 4 */
static void   update_Local_State(Matrix,Mask,Particle,Mesh,double);
/* Step 5 */
static Matrix compute_Nodal_Forces(Mask,Particle,Mesh,int);
static void   compute_Nodal_Internal_Forces(Matrix,Mask,Particle,Mesh);
static void   compute_Nodal_Nominal_traction_Forces(Matrix,Mask,Particle,Mesh,int);
static Matrix solve_Nodal_Equilibrium(Matrix,Matrix,Matrix,Matrix,Particle,Mesh,Mask,Mask);
/* Step 6 */
static  void  compute_Explicit_Newmark_Corrector(Particle,double,double);
/* Step 7 */
static void   output_selector(Particle, Mesh, Mask, Matrix, Matrix, Matrix, Matrix, double, int, int);

/**************************************************************/

void U_Newmark_Predictor_Corrector_Finite_Strains(
  Mesh FEM_Mesh,
  Particle MPM_Mesh,
  Time_Int_Params Parameters_Solver)
{

  /*
    Auxiliar variables for the solver
  */
  int Ndim = NumberDimensions;
  int Nactivenodes;
  int InitialStep = Parameters_Solver.InitialTimeStep;
  int NumTimeStep = Parameters_Solver.NumTimeStep;  

  double gamma = 0.5;
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  Matrix Lumped_Mass;
  Matrix Gravity_field;
  Matrix Velocity;
  Matrix D_Displacement;
  Matrix Forces;
  Matrix Reactions;

  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;


  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = U_DeltaT__SolversLib__(MPM_Mesh, DeltaX, CFL);
      print_step(TimeStep,DeltaTimeStep);
      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      Free_and_Restricted_Dofs = generate_Mask_for_static_condensation__MeshTools__(ActiveNodes,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("First step : Compute effective mass ... WORKING",TimeStep);

      Lumped_Mass = compute_Mass_Matrix(MPM_Mesh,FEM_Mesh,ActiveNodes);

      print_Status("DONE !!!",TimeStep);
         
      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Explicit Newmark predictor",TimeStep);
      print_Status("WORKING ...",TimeStep);
      
      compute_Explicit_Newmark_Predictor(MPM_Mesh, gamma, DeltaTimeStep);

      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Third step : Compute nodal magnitudes",TimeStep);
      print_Status("WORKING ...",TimeStep);

      Gravity_field = compute_Nodal_Gravity_field(ActiveNodes, MPM_Mesh, TimeStep);

      D_Displacement = compute_Nodal_D_Displacement(MPM_Mesh,FEM_Mesh,ActiveNodes,Lumped_Mass);

      Velocity = compute_Nodal_Velocity(MPM_Mesh,FEM_Mesh,ActiveNodes,Lumped_Mass);

      impose_Dirichlet_Boundary_Conditions(FEM_Mesh,D_Displacement,Velocity,ActiveNodes,TimeStep);

      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update local state",TimeStep);
      print_Status("WORKING ...",TimeStep);

      update_Local_State(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Compute equilibrium ... WORKING",TimeStep);

      Forces = compute_Nodal_Forces(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);

      Reactions = solve_Nodal_Equilibrium(Lumped_Mass,Gravity_field,Forces,D_Displacement,MPM_Mesh,FEM_Mesh,ActiveNodes,Free_and_Restricted_Dofs);

      print_Status("*************************************************",TimeStep);
      print_Status("Six step : Compute corrector",TimeStep);
      print_Status("WORKING ...",TimeStep);

      compute_Explicit_Newmark_Corrector(MPM_Mesh,gamma,DeltaTimeStep);
      local_search__Particles__(MPM_Mesh,FEM_Mesh);
      
      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Output variables and reset nodal values",TimeStep);
      print_Status("WORKING ...",TimeStep);

      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, Velocity, D_Displacement,Forces, Reactions, DeltaTimeStep, TimeStep, ResultsTimeStep);

      /*
      	Free memory
      */
      free__MatrixLib__(Lumped_Mass); 
      free__MatrixLib__(Gravity_field);
      free__MatrixLib__(Velocity);
      free__MatrixLib__(D_Displacement);
      free__MatrixLib__(Forces);
      free__MatrixLib__(Reactions);
      free(ActiveNodes.Nodes2Mask); 
      free(Free_and_Restricted_Dofs.Nodes2Mask);
      print_Status("DONE !!!",TimeStep);

    }
  
}

/**************************************************************/

static Matrix compute_Mass_Matrix(
  Particle MPM_Mesh, // Variable with information of the particles 
  Mesh FEM_Mesh, // Variable with information of the nodes
  Mask ActiveNodes) // Variable with information of the active nodes
/*
  This function computes the lumped mass matrix with the displacements degree of freedom
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndim*Nnodes_mask;
  int Ap;
  int A_mask;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;  

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Mass of the particle */
  double m_p;
  /* Nodal contribution A of the particle p */
  double m_A_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Order, Order);

  /*
    Iterate over the particles to get the nodal values 
  */
  for(int p = 0 ; p<Np ; p++)
    {

      /*
        Define tributary nodes of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* 
        Evaluate the shape function in the coordinates of the particle
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
       Get the mass of the particle 
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];


      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
        {

          /* 
             Get the node in the mass matrix with the mask
          */
          Ap = Nodes_p.Connectivity[A];
          A_mask = ActiveNodes.Nodes2Mask[Ap];

          /* 
            Get the value of the shape function 
          */
          ShapeFunction_pA = ShapeFunction_p.nV[A];

          /*
            Compute the nodal A contribution of the particle p
          */
          m_A_p = m_p*ShapeFunction_pA;

          /* 
             Fill the Lumped mass matrix considering the number of dofs
          */
          for(int i = 0 ; i<Ndim ; i++)
            {
              Lumped_MassMatrix.nM[A_mask*Ndim + i][A_mask*Ndim + i] += m_A_p;       
            }

          }

        /* 
          Free the value of the shape functions 
        */
        free__MatrixLib__(ShapeFunction_p);
        free(Nodes_p.Connectivity);      

    }

  /* 
     Add some usefulll info 
  */
  strcpy(Lumped_MassMatrix.Info,"Lumped-Mass-Matrix");

  return Lumped_MassMatrix; 
}


/**************************************************************/

static void compute_Explicit_Newmark_Predictor(
  Particle MPM_Mesh, // Information related with particles
  double gamma, // Newmark integration parameter
  double Dt) // Time step
/*
  The predictor stage is computed in the particles
*/
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;

  for(int p = 0 ; p<Np ; p++)
    {

    /* 
      Compute velocity predictor and increment of displacements 
    */
    for(int i = 0 ; i<Ndim ; i++)
    {

      MPM_Mesh.Phi.D_dis.nM[p][i] = Dt*MPM_Mesh.Phi.vel.nM[p][i] + 0.5*DSQR(Dt)*MPM_Mesh.Phi.acc.nM[p][i];

      MPM_Mesh.Phi.vel.nM[p][i] += (1-gamma)*Dt*MPM_Mesh.Phi.acc.nM[p][i];

    }


  }

}

/**************************************************************/

static Matrix compute_Nodal_Gravity_field(
  Mask ActiveNodes,
  Particle MPM_Mesh,
  int TimeStep)
/*

*/
{
  /* Define auxilar variables */
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  
  double m_p; /* Mass of the particle */
  Load * B = MPM_Mesh.B; /* List with the load cases */
  
  Matrix Gravity_field = allocZ__MatrixLib__(Nnodes_mask, Ndim);

  /*
    Initialise distance accelerations
  */
  Tensor b = alloc__TensorLib__(1);
  

  for(int i = 0 ; i<NumBodyForces ; i++)
    {

     /* Fill vector b of body acclerations */
     for(int k = 0 ; k<Ndim ; k++)
     {
       if(B[i].Dir[k])
       {
         if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num))
         {
            printf("%s : %s\n", "Error in compute_Nodal_Gravity_field()","The time step is out of the curve !!");
            exit(EXIT_FAILURE);
         }
         
         b.n[k] += B[i].Value[k].Fx[TimeStep];

       }
     }
      
    }

  /* Get the node of the mesh for the contribution */
     for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      /* Compute body forces */
      for(int k = 0 ; k<Ndim ; k++)
        {
          Gravity_field.nM[A][k] = b.n[k];
        }  
    }

    free__TensorLib__(b);


  return Gravity_field;
}

/**************************************************************/

static Matrix compute_Nodal_D_Displacement(
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Lumped_Mass)
/*
  Compute the nodal increment of displacement. The operation is linearized and
  all the dof split the increment of displacement array in n components like :
  | M 0 |   |D_u.x|   | M*D_u.x |
  | 0 M | * |D_u.y| = | M*D_u.y |
*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = Nnodes_mask*Ndim;
  int Ap;
  int A_mask;
  int AB;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p; /* Mass of the particle */

  /* 
    Define and allocate the nodal increment of displacement vector 
  */
  Matrix D_Displacement = allocZ__MatrixLib__(Nnodes_mask,Ndim);


  for(int p = 0 ; p<Np ; p++)
    {

      /* 
        Define element of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* 
        Evaluate the shape function in the coordinates of the particle
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* 
        Get the mass of the GP 
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
        {

          /*
            Get the node in the nodal momentum with the mask
          */
          Ap = Nodes_p.Connectivity[A];
          A_mask = ActiveNodes.Nodes2Mask[Ap];

          /*
             Evaluate the GP function in the node 
          */
          ShapeFunction_pA = ShapeFunction_p.nV[A];

          /*
            Nodal velocity and D_Displacement 
          */
          for(int i = 0 ; i<Ndim ; i++)
            {
              D_Displacement.nM[A_mask][i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.D_dis.nM[p][i];
            }
        }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*
    Compute the D_Displacements
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      AB = A*Ndim + i;
      D_Displacement.nM[A][i] = D_Displacement.nM[A][i]/Lumped_Mass.nM[AB][AB];
    }
  }


  /*
    Add some usefulll info
  */
  strcpy(D_Displacement.Info,"Nodal-D-Displacement");

  return D_Displacement;
}


/**************************************************************/

static Matrix compute_Nodal_Velocity(
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Lumped_Mass)
/*
  Compute the nodal velocity. The operation is linearized and
  all the dof split the velocity array in n components like :
  | M 0 |   |V.x|   | p.x |
  | 0 M | * |V.y| = | p.y |

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndim*Nnodes_mask;
  int Ap;
  int A_mask;
  int AB;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p; /* Mass of the particle */

  /* Define and allocate the velocity vector */
  Matrix Velocity = allocZ__MatrixLib__(Nnodes_mask,Ndim);


  for(int p = 0 ; p<Np ; p++)
    {

      /*
         Define element of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* 
        Evaluate the shape function in the coordinates of the particle 
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Get the mass of the GP
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* 
        Get the nodal mommentum 
      */
      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
        {
          /*
            Get the node with the mask
          */
          Ap = Nodes_p.Connectivity[A];
          A_mask = ActiveNodes.Nodes2Mask[Ap];

          /* 
            Evaluate the GP function in the node 
          */
          ShapeFunction_pA = ShapeFunction_p.nV[A];


          for(int i = 0 ; i<Ndim ; i++)
            {
              Velocity.nM[A_mask][i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.vel.nM[p][i];
            }
        }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*
    Compute the nodal velocities
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      AB = A*Ndim + i;
      Velocity.nM[A][i] = Velocity.nM[A][i]/Lumped_Mass.nM[AB][AB];
    }
  }

 
  /*
    Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");

  return Velocity;
}


/**************************************************************/

static void impose_Dirichlet_Boundary_Conditions(
  Mesh FEM_Mesh,
  Matrix D_Displacement,
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
  int Ndim = NumberDimensions; /* Number of dimensions */
  int Ndof = NumberDOF; /* Number of degree of freedom */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;

  /* 
    Loop over the the boundaries to set boundary conditions
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
    Ndof = FEM_Mesh.Bounds.BCC_i[i].Dim;

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
      for(int k = 0 ; k<Ndof ; k++)
      {

        /* 
          Apply only if the direction is active (1) 
        */
        if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
        {
    
          /* 
            Check if the curve it is on time 
          */
          if( (TimeStep < 0) || (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num))
          {
            printf("%s : %s \n","Error in initialise_Nodal_Increments()","The time step is out of the curve !!");
            exit(EXIT_FAILURE);
          }

          /*
            Initialise increments using newmark and the value of the boundary condition
          */
          D_Displacement.nM[Id_BCC_mask][k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];                    

        }
      }
    }    
  }

}

/**************************************************************/

static void update_Local_State(
  Matrix D_Displacement,
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
  int Nnodes_p;
  int MatIndx_p;
  double rho_n_p;
  double Delta_J_p;
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;
  Tensor P_p;

  /*
    Loop in the material point set to update strains
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
      Compute Jacobian of the deformation gradient
    */
    MPM_Mesh.Phi.J.nV[p] = I3__TensorLib__(F_n1_p);
            
    /*
      Update density with the jacobian of the increment deformation gradient
    */
    Delta_J_p = I3__TensorLib__(DF_p);
    rho_n_p = MPM_Mesh.Phi.rho.nV[p];
    MPM_Mesh.Phi.rho.nV[p] = rho_n_p/Delta_J_p;

    /*
      Free memory 
    */
    free__MatrixLib__(D_Displacement_Ap);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
	  
  }

  /*
    Loop in the material point set to update stress
  */
  for(int p = 0 ; p<Np ; p++)
  {
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];

    /*
      Activate locking control technique (F-bar)
    */
    if(MPM_Mesh.Mat[MPM_Mesh.MatIdx[p]].Locking_Control_Fbar)
    {
      get_locking_free_Deformation_Gradient_n1__Particles__(p,MPM_Mesh,FEM_Mesh);
    }

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    P_p = forward_integration_Stress__Particles__(p,MPM_Mesh,MatProp_p); 

  }

  
}

/**************************************************************/

static Matrix compute_Nodal_Forces(
  Mask ActiveNodes,
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
    Add contact forces contribution
  */
  compute_Nodal_Nominal_traction_Forces(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

  
  return Forces;
}


/**************************************************************/

static void compute_Nodal_Internal_Forces(
  Matrix Forces,
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
  Tensor F_n_p;
  Tensor transpose_F_n_p;
  
  double V0_p; /* Volume of the Gauss-Point */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
  {
    /*
      Get the volume of the particle in the reference configuration 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

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
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);

    /*
      Get the first Piola-Kirchhoff stress tensor.
    */
    Tensor P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 

    for(int A = 0 ; A<NumNodes_p ; A++)
    {

      /*
        Compute the gradient in the reference configuration 
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pA);
      
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
	      Forces.nM[A_mask][i] -= InternalForcesDensity_Ap.n[i]*V0_p;
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
    free__TensorLib__(transpose_F_n_p);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }
   
}

/**************************************************************/


static void compute_Nodal_Nominal_traction_Forces(
  Matrix Forces,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{

  int Ndim = NumberDimensions;
  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA;
  Tensor T = alloc__TensorLib__(1); // Nominal traction
  double A0_p; /* Area of the particle in the reference configuration */ 

  int NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  int NumNodesLoad;
  int p;
  int Ap;
  int A_mask;
  int NumNodes_p; /* Number of nodes of each particle */

  for(int i = 0 ; i<NumContactForces; i++)
  {

    /*
      Read load i
    */
    Load_i = MPM_Mesh.Neumann_Contours.BCC_i[i];

    NumNodesLoad = Load_i.NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++)
    {

      /*
        Get the index of the particle
      */
      p = Load_i.Nodes[j];

      /*
        Get the area of each particle
      */      
      if(Ndim == 2)
      {
        A0_p = MPM_Mesh.Phi.Vol_0.nV[p]/Thickness_Plain_Stress;
      }
      else if(Ndim == 3)
      {
        A0_p = MPM_Mesh.Phi.Area_0.nV[p];
      }

      /*
        Define tributary nodes of the particle 
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],NumNodes_p);

      /* 
        Evaluate the shape function in the coordinates of the particle
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Fill vector of contact forces
      */
      for(int k = 0 ; k<Ndim ; k++)
      {
        if(Load_i.Dir[k] == 1)
        {
          if((TimeStep < 0) || (TimeStep > Load_i.Value[k].Num))
          {
            printf("%s : %s",
              "Error in compute_Nodal_Nominal_traction_Forces()",
              "The time step is out of the curve !!");
            exit(EXIT_FAILURE);
          }

          T.n[k] = Load_i.Value[k].Fx[TimeStep];
        }
      }

      /*
        Get the node of the mesh for the contribution
      */
      for(int A = 0 ; A<NumNodes_p ; A++)
      {

        /*
          Pass the value of the nodal shape function to a scalar
        */
        ShapeFunction_pA = ShapeFunction_p.nV[A];

        /*
          Node for the contribution
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];
  
        /*
          Compute Contact forces
        */
        for(int k = 0 ; k<Ndim ; k++)
        {
          Forces.nM[A_mask][k] += ShapeFunction_pA*T.n[k]*A0_p;
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
  
    }

      
  }

  free__TensorLib__(T);

}

/**************************************************************/

static Matrix solve_Nodal_Equilibrium(
  Matrix Lumped_Mass,
  Matrix Gravity_field,
  Matrix Total_Forces,
  Matrix D_Displacement,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Mask Free_and_Restricted_Dofs)
/*
  Call the LAPACK solver to compute the accelerations and velocities
  Solve equilibrium equation to get the nodal values of the aceleration 
  at t = n + 1
*/
{
  
  /*
    General varibles
  */
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
  int Nnodes = ActiveNodes.Nactivenodes;
  int NumNodes_p;
  int Order = Nnodes*Ndim;
  int AB;
  int Ap;
  int A_mask;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p;
  double ShapeFunction_pA;

  /* 
    Solution nodal variable
  */
  Matrix Acceleration = allocZ__MatrixLib__(Nnodes,Ndim);

  /*
    Output
  */
  Matrix Reactions = allocZ__MatrixLib__(Nnodes, Ndim);
  strcpy(Reactions.Info,"REACTIONS");

  /*
    The solution is now stored in the internal forces vector
  */
  for(int A = 0; A<Nnodes; A++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      if(Free_and_Restricted_Dofs.Nodes2Mask[A*Ndof + i] != -1)
      {
        AB = A*Ndim + i;
        Acceleration.nM[A][i] = Gravity_field.nM[A][i] + Total_Forces.nM[A][i]/Lumped_Mass.nM[AB][AB];
      }
      else
      {
        Acceleration.nM[A][i] = 0.0;
        Reactions.nM[A][i] = Total_Forces.nM[A][i];
      }
      
    }
  }

  /*
    Update particle acceleration
  */
  for(int p = 0 ; p<Np ; p++)
  {

    /*
      Set to zero particle acceleration for interpolation
    */
    for(int i = 0 ; i<Ndim ; i++)
    {
      MPM_Mesh.Phi.acc.nM[p][i] = 0.0;
      MPM_Mesh.Phi.D_dis.nM[p][i] = 0.0;
    }

    /*
      Define element of the particle 
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Evaluate the shape function in the coordinates of the particle 
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Interpolate the new value of the acceleration
    */
    for(int A = 0; A<NumNodes_p; A++)
    {
      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];
    
      /*
        Evaluate the GP function in the node 
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];
      
      for(int i = 0 ; i<Ndim ; i++)
       {
          MPM_Mesh.Phi.acc.nM[p][i] += ShapeFunction_pA*Acceleration.nM[A_mask][i];
          MPM_Mesh.Phi.D_dis.nM[p][i] += ShapeFunction_pA*D_Displacement.nM[A_mask][i];
       }
    }

    /*
      Free auxiliar data
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);

  }

  /*
    Free acceleration vector
  */
  free__MatrixLib__(Acceleration);


  return Reactions;
}

/**************************************************************/

static void compute_Explicit_Newmark_Corrector(
  Particle MPM_Mesh,
  double gamma,
  double Dt)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  Tensor F_n_p;
  Tensor F_n1_p;


  for(int p = 0 ; p<Np ; p++)
    {
      
      /*
        Replace the deformation gradient at t = n with the new one
      */
      F_n_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

      /* 
        Update/correct tensor and vector variables
      */
      for(int i = 0 ; i<Ndim  ; i++)
       {

          /* 
            Correct particle velocity
          */
          MPM_Mesh.Phi.vel.nM[p][i] += gamma*Dt*MPM_Mesh.Phi.acc.nM[p][i];

          /*
            Update the particles position and displacement
          */
          MPM_Mesh.Phi.x_GC.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];
          MPM_Mesh.Phi.dis.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];

          /* Update deformation gradient tensor */
         for(int j = 0 ; j<Ndim  ; j++)
           {
             F_n_p.N[i][j] = F_n1_p.N[i][j];
           }
       }

    }  
}

/**************************************************************/

static void output_selector(
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Velocity,
  Matrix D_Displacement,
  Matrix Forces,
  Matrix Reactions,
  double DeltaTimeStep,
  int TimeStep,
  int ResultsTimeStep)
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
