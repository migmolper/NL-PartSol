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

/* Step 1 */
static Matrix compute_Mass_Matrix_Mixture(GaussPoint,Mesh,Mask,double);
static Matrix compute_Compressibility_Matrix_Fluid(GaussPoint, Mesh, Mask);
/* Step 2 */
static  void  compute_Explicit_Newmark_Predictor(GaussPoint,double,double);
/* Step 3 */
static Matrix compute_Nodal_Gravity_field(Mask, GaussPoint, int);
static Matrix compute_Nodal_D_Displacement(GaussPoint,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Acceleration(GaussPoint,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Velocity(GaussPoint,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Pore_water_pressure(GaussPoint,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Rate_Pore_water_pressure(GaussPoint,Mesh,Mask,Matrix);
static  void  Imposse_Velocity(Mesh,Matrix,Mask,int);
/* Step 4 */
static  void  update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);
/* Step 5 */
static Matrix compute_Internal_Forces_Mixture(Mask,GaussPoint,Mesh);
static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor,double,Tensor);
static Matrix compute_Reactions(Mesh,Matrix,Mask);
static  void  solve_Nodal_Equilibrium_Mixture(Matrix,Matrix,Matrix,Matrix);
/* Step 6 */
static Matrix compute_Viscous_Forces_Fluid(Mask,GaussPoint,Mesh,Matrix);
static Matrix compute_C_Viscosity(Mask,GaussPoint,Mesh);
static Matrix compute_Permeability_Forces_Fluid(Mask,GaussPoint,Mesh,Matrix);
static Matrix compute_K_Permeability(Mask,GaussPoint,Mesh);
static Matrix compute_Permeability_Inertial_Forces_Fluid(Mask, GaussPoint, Mesh, Matrix, Matrix);
static Matrix compute_C_Permeability(Mask,GaussPoint,Mesh);
static  void  solve_Nodal_Generalized_Darcy_Law(Matrix,Matrix,Matrix,Matrix,Matrix);
/* Step 7 */
static  void  compute_Explicit_Newmark_Corrector(GaussPoint,Mesh,Matrix,Matrix,Matrix,Mask,double, double);
/* Step 8 */
static  void  output_selector(GaussPoint, Mesh, Mask, Matrix, Matrix, Matrix, Matrix, int, int);

/**************************************************************/

void upw_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
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
  Matrix Mass_Matrix_Mixture;
  Matrix Compressibility_Matrix_Fluid;
  Matrix Velocity;
  Matrix Pore_water_pressure;
  Matrix Acceleration;
  Matrix Rate_Pore_water_pressure;
  Matrix Gravity_field;
  Matrix D_Displacement;
  Matrix Internal_Forces_Mixture;
  Matrix Reactions;
  Mask ActiveNodes;

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("First step : Compute lumped mass of the mixture and compressibility matrix ... WORKING",TimeStep);

      Mass_Matrix_Mixture = compute_Mass_Matrix_Mixture(MPM_Mesh,FEM_Mesh,ActiveNodes,1);

      Compressibility_Matrix_Fluid = compute_Compressibility_Matrix_Fluid(MPM_Mesh,FEM_Mesh,ActiveNodes,1);

      print_Status("DONE !!!",TimeStep);
   
      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Explicit Newmark predictor ... WORKING",TimeStep);
      
      compute_Explicit_Newmark_Predictor(MPM_Mesh, gamma, Dt);

      print_Status("*************************************************",TimeStep);
      print_Status("Third step : Compute nodal magnitudes ... WORKING",TimeStep);

      Gravity_field = compute_Nodal_Gravity_field(ActiveNodes, MPM_Mesh, TimeStep);

      D_Displacement = compute_Nodal_D_Displacement(MPM_Mesh, FEM_Mesh, ActiveNodes, Mass_Matrix_Mixture);

      Velocity = compute_Nodal_Velocity(MPM_Mesh, FEM_Mesh, ActiveNodes, Mass_Matrix_Mixture);

      Acceleration = compute_Nodal_Acceleration(MPM_Mesh, FEM_Mesh, ActiveNodes, Mass_Matrix_Mixture);
      
      Pore_water_pressure = compute_Pore_water_pressure(MPM_Mesh, FEM_Mesh, ActiveNodes, Compressibility_Matrix_Fluid);

      Rate_Pore_water_pressure = compute_Rate_Pore_water_pressure(MPM_Mesh, FEM_Mesh, ActiveNodes, Compressibility_Matrix_Fluid);
 
      Imposse_Velocity(FEM_Mesh,Velocity,ActiveNodes,TimeStep);

      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update local state ... WORKING",TimeStep);

      update_Local_State(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);

      print_Status("DONE !!!",TimeStep);


      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Compute equilibrium mixture ... WORKING",TimeStep);

      Internal_Forces_Mixture = compute_Internal_Forces_Mixture(ActiveNodes, MPM_Mesh, FEM_Mesh);

      Reactions = compute_Reactions(FEM_Mesh,Internal_Forces_Mixture,ActiveNodes);

      solve_Nodal_Equilibrium_Mixture(Acceleration,Mass_Matrix_Mixture,Gravity_field,Internal_Forces_Mixture);

      print_Status("*************************************************",TimeStep);
      print_Status("Six step : Compute equilibrium fluid ... WORKING",TimeStep);

      Viscous_Forces_Fluid = compute_Viscous_Forces_Fluid(ActiveNodes, MPM_Mesh, FEM_Mesh, Velocity);

      Permeability_Forces_Fluid = compute_Permeability_Forces_Fluid(ActiveNodes, MPM_Mesh, FEM_Mesh, Pore_water_pressure);

      Permeability_Inertial_Forces_Fluid = compute_Permeability_Inertial_Forces_Fluid(ActiveNodes, MPM_Mesh, FEM_Mesh, Acceleration, Gravity_field);

      solve_Nodal_Generalized_Darcy_Law(Rate_Pore_water_pressure,Compressibility_Matrix_Fluid,Viscous_Forces_Fluid,Permeability_Forces_Fluid,Permeability_Inertial_Forces_Fluid);

      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Compute corrector ... WORKING",TimeStep);

      compute_Explicit_Newmark_Corrector(MPM_Mesh,FEM_Mesh,D_Displacement,Acceleration,Rate_Pore_water_pressure,ActiveNodes,gamma,Dt);

      local_search__Particles__(MPM_Mesh,FEM_Mesh);

      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Eight step : Output variables and reset nodal values ... WORKING",TimeStep);

      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, Velocity, D_Displacement,
                      Internal_Forces_Mixture, Reactions, TimeStep, ResultsTimeStep);
      /*
      	Free memory.
      */
      free__MatrixLib__(Mass_Matrix_Mixture); 
      free__MatrixLib__(Compressibility_Matrix_Fluid);
      free__MatrixLib__(Velocity);
      free__MatrixLib__(D_Displacement);
      free__MatrixLib__(Internal_Forces_Mixture);
      free__MatrixLib__(Gravity_field);
      free__MatrixLib__(Reactions);
      free(ActiveNodes.Nodes2Mask);
      
      print_Status("DONE !!!",TimeStep);

    }
  
}

/**************************************************************/

static Matrix compute_Mass_Matrix_Mixture(
  GaussPoint MPM_Mesh, // Variable with information of the particles 
  Mesh FEM_Mesh, // Variable with information of the nodes
  Mask ActiveNodes, // Variable with information of the active nodes
  double epsilon) // Combianation varible
/*
  This function computes the effective mass matrix (mixture) as a convex combination
  of the lumped mass matrix and the consistent mass matrix. Later assemble
  a total mass matrix with the contribution of each degree of freedom.

  | M_eff |   0   |              | M_cons |   0    |          | M_lump |   0    |
  -----------------  = (1-eps) * -------------------  + eps * -------------------
  |    0  | M_eff |              |   0    | M_cons |          |   0    | M_lump |
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndim*Nnodes_mask;
  int Ap;
  int Bp;
  int A_mask;
  int B_mask;
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
  /* Nodal contribution A-B of the particle p */
  double m_AB_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Define and allocate the effective mass matrix */
  Matrix Effective_MassMatrix = allocZ__MatrixLib__(Order, Order);

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

          for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
            {       
              /* 
               Get the node in the mass matrix with the mask
              */
              Bp = Nodes_p.Connectivity[B];
              B_mask = ActiveNodes.Nodes2Mask[Bp];

              /* 
                Get the value of the shape function 
              */
              ShapeFunction_pB = ShapeFunction_p.nV[B];

              /*
               Compute the nodal AB contribution of the particle p
              */
              m_AB_p = m_p*ShapeFunction_pA*ShapeFunction_pB;

              /* 
               Fill the effective mass matrix considering the number of dofs
              */
              for(int i = 0 ; i<Ndim ; i++)
              {
                /*
                  Compute the vectorized index
                */
                Effective_MassMatrix.nM[A_mask*Ndim+i][A_mask*Ndim+i] += m_AB_p;
              }

            }
          }

        /* 
          Free the value of the shape functions 
        */
        free__MatrixLib__(ShapeFunction_p);
        free(Nodes_p.Connectivity);      

    }

  /*
    At this point the effective mass matrix coincides with the consistent mass
    matrix. We can tune it by a convex combination with the lumped mass matrix
  */
  for(int A = 0 ; A<Order ; A++)
    {
      for(int B = 0 ; B<Order ; B++)
      {    
        Effective_MassMatrix.nM[A][B] = (1-epsilon)*Effective_MassMatrix.nM[A][B] + (A == B)*epsilon*Lumped_MassMatrix.nM[A][B];
      }
    }

  /*
    Free lumped mass matrix.
  */
  free__MatrixLib__(Lumped_MassMatrix);

  /* 
     Add some usefulll info 
  */
  strcpy(Effective_MassMatrix.Info,"Effective-Mass-Matrix");

  return Effective_MassMatrix; 
}

/**************************************************************/


static Matrix compute_Compressibility_Matrix_Fluid(
  GaussPoint MPM_Mesh, // Variable with information of the particles
  Mesh FEM_Mesh, // Variable with information of the nodes
  Mask ActiveNodes, // Variable with information of the active nodes
  double epsilon) // Combianation varible
/*
  Compute the compresibility matrix
*/
{
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = 1*Nnodes_mask;
  int Mat_idx;
  int Ap;
  int A_mask;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double V0_p; /* Mass of the particle (mixture) */
  double rho_f_p; /* Material density of the fluid */
  double phi_f_p; /* Volume fractions of the fluid */
  double K_f; /* Compressibility (fluid) */
  double compressibility_density_f_p; /* Compressibilidy density fo the fluid */
  double compressibility_density_A_p; /* Nodal contribution A of the particle p */
  double compressibility_density_AB_p; /* Nodal contribution AB of the particle p */  

  /* 
    Define and allocate the effective Compressibility matrix 
  */
  Matrix Compressibility_Matrix = allocZ__MatrixLib__(Order, Order);

  /* 
    Define and allocate the lumped Compressibility matrix 
  */
  Matrix Lumped_Compressibility_Matrix = allocZ__MatrixLib__(Order, Order);


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
        Get the initial volume of the particle 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Get the current material density, volume fraction 
        and compressibility for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
      phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
      Mat_idx = MPM_Mesh.MatIdx[p];
      K_f = MPM_Mesh.Mat[Mat_idx].K_f;

      /*
        Compute relative density 
      */
      relative_rho_f_p = phi_f_p*rho_f_p;

      /* 
        Compute the compressibility density 
      */
      compressibility_density_f_p = (relative_rho_f_p/K_f)*V0_p;
      
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
        compressibility_density_A_p = compressibility_density_f_p*ShapeFunction_pA;

        /* 
          Fill the Lumped Compressibility matrix
        */
        Lumped_Compressibility_Matrix.nM[A_mask][A_mask] += compressibility_density_A_p;      


        for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
        {       
              /* 
               Get the node in the mass matrix with the mask
              */
              Bp = Nodes_p.Connectivity[B];
              B_mask = ActiveNodes.Nodes2Mask[Bp];

              /* 
                Get the value of the shape function 
              */
              ShapeFunction_pB = ShapeFunction_p.nV[B];

              /*
               Compute the nodal AB contribution of the particle p
              */
              compressibility_density_AB_p = compressibility_density_f_p*ShapeFunction_pA*ShapeFunction_pB;

              /* 
                Fill the effective Compressibility matrix
              */
              Compressibility_Matrix.nM[A_mask][B_mask] += compressibility_density_AB_p;
        }       
      }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);      
      
    }

  /*
    At this point the effective Compressibility matrix coincides with the consistent Compressibility
    matrix. We can tune it by a convex combination with the lumped Compressibility matrix
  */
  for(int A = 0 ; A<Order ; A++)
    {
      for(int B = 0 ; B<Order ; B++)
      {    
        Compressibility_Matrix.nM[A][B] = (1-epsilon)*Compressibility_Matrix.nM[A][B] + (A == B)*epsilon*Lumped_Compressibility_Matrix.nM[A][B];
      }
    }

  /*
    Free lumped Compressibility matrix.
  */
  free__MatrixLib__(Lumped_Compressibility_Matrix);

  /* Add some usefulll info */
  strcpy(Compressibility_Matrix.Info,"Compressibility-Matrix");

  return Compressibility_Matrix_Fluid;
}

/**************************************************************/

static  void  compute_Explicit_Newmark_Predictor(
  GaussPoint MPM_Mesh, // Information related with particles
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
      Compute pore water pressure predictor
    */
    MPM_Mesh.Phi.Pw.nV[p] += (1-gamma)*Dt*MPM_Mesh.Phi.d_Pw.nV[p];

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
  GaussPoint MPM_Mesh,
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
  Tensor b = alloc__TensorLib__(1); /* Body forces vector */
  
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
         
         b.n[k] = B[i].Value[k].Fx[TimeStep];

       }
     }

     /* Get the node of the mesh for the contribution */
     for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      /* Compute body forces */
      for(int k = 0 ; k<Ndim ; k++)
        {
          Gravity_field.nV[A*Ndim + k] += b.n[k];
        }  
     }
      
    }

  /*
    Free auxiliar tensor
  */
  free__TensorLib__(b);

  return Gravity_field;
}

/**************************************************************/

static Matrix compute_Nodal_D_Displacement(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Mass_Matrix_Mixture)
/*
  Call the LAPACK solver to compute the nodal increment of displacement. The operation is linearized and
  all the dof split the increment of displacement array in n components like :
  | M 0 |   |D_u.x|   | M*D_u.x |
  | 0 M | * |D_u.y| = | M*D_u.y |
*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
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
    Call the LAPACK solver to compute the D_Displacements
  */
  int Order = Nnodes_mask*Ndim;
  int LDA   = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Mass_Matrix_Mixture.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the D_Displacement
  */
  dgetrs_(&TRANS,&Order,&NRHS,Mass_Matrix_Mixture.nV,&LDA,IPIV,D_Displacement.nV,&LDB,&INFO);
  free(IPIV);

  /*
    Add some usefulll info
  */
  strcpy(D_Displacement.Info,"Nodal-D-Displacement");

  return D_Displacement;
}


/**************************************************************/

static Matrix compute_Nodal_Acceleration(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Mass_Matrix_Mixture)
/*
  Call the LAPACK solver to compute the nodal Acceleration. The operation is linearized and
  all the dof split the acceleration array in n components like :
  | M 0 |   |A.x|   | dt_p.x |
  | 0 M | * |A.y| = | dt_p.y |
*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p; /* Mass of the particle */

  /* Define and allocate the nodal acceleration vector */
  Matrix Acceleration = allocZ__MatrixLib__(Nnodes_mask,Ndim);


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
            Get the node in the nodal momentum with the mask
          */
          Ap = Nodes_p.Connectivity[A];
          A_mask = ActiveNodes.Nodes2Mask[Ap];

          /* 
            Evaluate the GP function in the node 
          */
          ShapeFunction_pA = ShapeFunction_p.nV[A];

          /* 
            Nodal velocity and acceleration 
          */
          for(int i = 0 ; i<Ndim ; i++)
            {
              Acceleration.nM[A_mask][i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.acc.nM[p][i];
            }
        }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }


  /*
    Call the LAPACK solver to compute the nodal accelerations
  */
  int Order = Nnodes_mask*Ndim;
  int LDA   = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Mass_Matrix_Mixture.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the nocal acceleration
  */
  dgetrs_(&TRANS,&Order,&NRHS,Mass_Matrix_Mixture.nV,&LDA,IPIV,Acceleration.nV,&LDB,&INFO);
  free(IPIV);
 
  /*
    Add some usefulll info
  */
  strcpy(Acceleration.Info,"Nodal-Acceleration");

  return Acceleration;
}

/**************************************************************/

static Matrix compute_Nodal_Velocity(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Mass_Matrix_Mixture)
/*
  Call the LAPACK solver to compute the nodal velocity. The operation is linearized and
  all the dof split the velocity array in n components like :
  | M 0 |   |V.x|   | p.x |
  | 0 M | * |V.y| = | p.y |

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
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
    Call the LAPACK solver to compute the velocities
  */
  int Order = Nnodes_mask*Ndim;
  int LDA   = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Mass_Matrix_Mixture.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the velocity
  */
  dgetrs_(&TRANS,&Order,&NRHS,Mass_Matrix_Mixture.nV,&LDA,IPIV,Velocity.nV,&LDB,&INFO);
  free(IPIV);
 
  /*
    Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");

  return Velocity;
}

/**************************************************************/

static Matrix compute_Nodal_Pore_water_pressure(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Compressibility_Matrix_Fluid)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Mat_idx;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double Pw_p; /* Particle pore water pressure */
  double V0_p; /* Volume of the particle */
  double rho_f_p; /* Material density of the particle (fluid phase) */
  double phi_f_p; /* Volume fraction of the particle (fluid phase) */
  double relative_rho_f_p; /* Relative density of the particle (fluid phase) */
  double K_f; /* Compressibility (fluid) */
  double compressibility_density_f_p; /* Compressibily density for the particle (fluid phase) */
  double compressibility_density_A_p; /* Contribution of the particle p to the compressibility in the node A */

  /* 
    Define and allocate the Pore water pressure vector 
  */
  Matrix Pore_water_pressure = allocZ__MatrixLib__(Nnodes_mask,1);

  /* Iterate over the particles to get the nodal values */
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
         Get the pore water pressure 
      */
      Pw_p = MPM_Mesh.Phi.Pw.nV[p];

      /* 
        Get the initial volume of the particle 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Get the current material density, volume fraction 
        and compressibility for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
      phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
      Mat_idx = MPM_Mesh.MatIdx[p];
      K_f = MPM_Mesh.Mat[Mat_idx].K_f;

      /*  
        Compute relative density 
      */
      relative_rho_f_p = phi_f_p*rho_f_p;

      /* 
        Compute the compressibility density 
      */
      compressibility_density_f_p = (relative_rho_f_p/K_f)*V0_p;

      /* Get the nodal mommentum */
      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
        {
          /*
            Get the node in the nodal momentum with the mask
          */
          Ap = Nodes_p.Connectivity[A];
          A_mask = ActiveNodes.Nodes2Mask[Ap];

          /* Evaluate the GP function in the node */
          ShapeFunction_pA = ShapeFunction_p.nV[A];

          /* 
            Compute the nodal A contribution of the particle p 
          */
          compressibility_density_A_p = compressibility_density_f_p*ShapeFunction_pA;

          /* Nodal Pore water pressure */
          Pore_water_pressure.nV[A_mask] += compressibility_density_A_p*pw_p;
        }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }


  /*
    Call the LAPACK solver to compute the nodal pore water pressure
  */
  int Order = Nnodes_mask*Ndim;
  int LDA   = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Compressibility_Matrix_Fluid.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the sistem
  */
  dgetrs_(&TRANS,&Order,&NRHS,Compressibility_Matrix_Fluid.nV,&LDA,IPIV,Pore_water_pressure.nV,&LDB,&INFO);
  free(IPIV);
 
  /*
    Add some usefulll info
  */
  strcpy(Pore_water_pressure.Info,"Nodal-Pore-water-pressure");

  return Pore_water_pressure;
}

/**************************************************************/

static Matrix compute_Nodal_Rate_Pore_water_pressure(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Compressibility_Matrix_Fluid)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Mat_idx;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double Pw_p; /* Particle pore water pressure */
  double V0_p; /* Volume of the particle */
  double rho_f_p; /* Material density of the particle (fluid phase) */
  double phi_f_p; /* Volume fraction of the particle (fluid phase) */
  double K_f; /* Compressibility (fluid) */
  double relative_rho_f_p; /* Relative density of the particle (fluid phase) */
  double compressibility_density_f_p; /* Compressibily density for the particle (fluid phase) */
  double compressibility_density_A_p; /* Contribution of the particle p to the compressibility in the node A */

  /* Define and allocate the Pore water pressure vector */
  Matrix Rate_Pore_water_pressure = allocZ__MatrixLib__(Nnodes_mask,1);

  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {

      /* Define element of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Get the rate of pore water pressure */
      d_Pw_p = MPM_Mesh.Phi.d_Pw.nV[p];

      /* Get the initial volume of the particle */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Get the current material density, volume fraction 
        and compressibility for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
      phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
      Mat_idx = MPM_Mesh.MatIdx[p];
      K_f = MPM_Mesh.Mat[Mat_idx].K_f;

      /* Compute relative density */
      relative_rho_f_p = phi_f_p*rho_f_p;

      /* Get the nodal mommentum */
      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
        {
          /*
            Get the node in the nodal momentum with the mask
          */
          Ap = Nodes_p.Connectivity[A];
          A_mask = ActiveNodes.Nodes2Mask[Ap];

          /* Evaluate the GP function in the node */
          ShapeFunction_pA = ShapeFunction_p.nV[A];

          /* 
            Compute the nodal A contribution of the particle p 
          */
          compressibility_density_A_p = compressibility_density_f_p*ShapeFunction_pA;

          /* Nodal Pore water pressure */
          Rate_Pore_water_pressure.nV[A_mask] += compressibility_density_A_p*d_Pw_p;
        }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }


  /*
    Call the LAPACK solver to compute the rate of pore water pressure
  */
  int Order = Nnodes_mask*Ndim;
  int LDA   = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Compressibility_Matrix_Fluid.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the sistem
  */
  dgetrs_(&TRANS,&Order,&NRHS,Compressibility_Matrix_Fluid.nV,&LDA,IPIV,Rate_Pore_water_pressure.nV,&LDB,&INFO);
  free(IPIV);
 
  /*
    Add some usefulll info
  */
  strcpy(Rate_Pore_water_pressure.Info,"Nodal-Pore-water-pressure");

  return Rate_Pore_water_pressure;
}

/**********************************************************************/

static void Imposse_Velocity(
  Mesh FEM_Mesh,
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
                  if( (TimeStep < 0) || (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num))
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
                  Velocity.nV[Id_BCC_mask_k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
                }
            }
        }    
    }

}

/**************************************************************/

static void update_Local_State(
  Matrix D_Displacement, // Nodal values of the increment of displacement
  Mask ActiveNodes, // Information with the active nodes
  GaussPoint MPM_Mesh, // Information related with the particles
  Mesh FEM_Mesh, // Information related with the nodes
  double TimeStep) // Time step
/*
  
*/
{

  /*
    Auxiliar variables
  */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int MatIndx_p;
  int Nnodes_p;
  double J_n1_p; /* Jacobian of the deformation gradient at t = n + 1 */
  double Delta_J_p; /* Jacobian of the increment of deformation gradient */
  double K_f; /* Compressibility (fluid) */
  double rho_f_0; /* Initial density of the fluid */
  Plastic_status Input_Plastic_Parameters; /* Input parameters for plasticity */
  Plastic_status Output_Plastic_Parameters; /* Output parameters for plasticity */
  Element Nodes_p; /* Element for each particle */
  Material MatProp_p; /* Variable with the material properties of each particle */
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Tensor F_n_p; /* Deformation gradient of the soil skeleton (t = n) */
  Tensor F_n1_p; /* Deformation gradient of the soil skeleton (t = n + 1) */
  Tensor DF_p; /* Increment of the deformation gradient of the soil skeleton */
  Tensor F_plastic_p; /* Plastic deformation gradient of the soil skeleton */
  Tensor S_p; /* Second Piola-Kirchhoff stress tensor */
  Tensor C_n1_p; /* Right Cauchy-Green tensor */
  double Pw_0; /* Cauchy pore water pressure at t = 0 */
  double Pw_n1; /* Cauchy pore water pressure at t = n + 1 */
  
  /*
    Loop in the material point set 
  */
  for(int p = 0 ; p<Np ; p++)
    {
      /*
        Read material properties
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      K_f = MatProp_p.K_f;
      rho_f_0 = MatProp_p.rho_f_0;

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
        Compute the Jacobian of the deformation gradient
      */
      Delta_J_p = I3__TensorLib__(DF_p);
      J_n1_p = I3__TensorLib__(F_n1_p);
      
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
        Compute the right Cauchy Green tensor
      */
      C_n1_p = right_Cauchy_Green__Particles__(F_n1_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	      integration rule.
      */
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      

      if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
        {
          S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, C_n1_p, MatProp_p);
        }
      else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
        {
          S_p = grad_energy_Neo_Hookean_Wriggers(S_p, C_n1_p, J_n1_p, MatProp_p);
        }
      else if(strcmp(MatProp_p.Type,"Von-Mises") == 0)
        {
          F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
          Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
          Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

          Output_Plastic_Parameters = finite_strains_plasticity_Von_Mises(S_p, C_n1_p, F_plastic_p, F_n1_p, 
                                                                          Input_Plastic_Parameters, MatProp_p, J_n1_p);

          MPM_Mesh.Phi.cohesion.nV[p] = Output_Plastic_Parameters.Yield_stress;
          MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;
        }
      else if((strcmp(MatProp_p.Type,"Drucker-Prager-Plane-Strain") == 0) || (strcmp(MatProp_p.Type,"Drucker-Prager-Outer-Cone") == 0))
        {
          F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
          Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
          Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

          Output_Plastic_Parameters = finite_strains_plasticity_Drucker_Prager_Sanavia(S_p, C_n1_p, F_plastic_p, F_n1_p, 
                                                                          Input_Plastic_Parameters, MatProp_p, J_n1_p);

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
        Update state parameters
      */
      Pw_0 = MPM_Mesh.Phi.Pw_0[p]; /* Get the initial pressure */
      Pw_n1 = MPM_Mesh.Phi.Pw.nV[p]/J_n1_p; /* From the Kirchhoff pressure compure the cauchy pore water pressure */

      MPM_Mesh.Phi.rho_f.nV[p] = rho_f_0*exp((Pw_n1-Pw_0)/K_f); /* Update the fluid pressure */

      MPM_Mesh.Phi.phi_s.nV[p] *= (1/Delta_J_p); /* Update the volume fraction of the solid phase */
      MPM_Mesh.Phi.phi_f.nV[p] = 1 - MPM_Mesh.Phi.phi_s.nV[p]; /* Update the volume fraction of the fluid phase */

      MPM_Mesh.Phi.rho.nV[p] = MPM_Mesh.Phi.rho_s.nV[p]*MPM_Mesh.Phi.phi_s.nV[p] +
                               MPM_Mesh.Phi.rho_f.nV[p]*MPM_Mesh.Phi.phi_f.nV[p];  /* Update density of the mixture */

      /*
	       Free memory 
      */
      free__TensorLib__(C_n1_p);
      free__MatrixLib__(D_Displacement_Ap);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Internal_Forces_Mixture(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh, Mesh
  FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;
  int idx_A_mask_i;

  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask,Ndim);

  Tensor P_p; /* Total First Piola-Kirchhoff Stress tensor */
  Tensor P_effective_p; /* Effective First Piola-Kirchhoff Stress tensor */
  Tensor S_effective_p; /* Effective Second Piola-Kirchhoff Stress tensor */
  double theta_p; /* Kirchhoff pore fluid pressure */

  Tensor InternalForcesDensity_Ap;

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor transpose_F_n_p;
  double V0_p; /* Volume of the Gauss-Point */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
    {

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
	       Take the values of the deformation gradient at t = n. And transpose it
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

      /*
        Get the volume of the particle in the reference configuration 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
	       Get the first Piola-Kirchhoff stress tensor
      */
      S_effective_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);
      P_effective_p = matrix_product__TensorLib__(F_n1_p, S_effective_p);

      /*
        Get the Kirchhoff pressure
      */
      theta_p = MPM_Mesh.Phi.Pw.nV[p];

      /*
        Following Terzaghi's idea, the effective stress tensor in the reference configuration
        is computed:
      */
      P_p = compute_total_first_Piola_Kirchhoff_stress(P_effective_p,theta_p,F_n1_p);

    
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
      free__TensorLib__(transpose_F_n_p);
      free__TensorLib__(P_effective_p);
      free__TensorLib__(P_p);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
    }
   

   return Forces;
}

/**************************************************************/

static Tensor compute_total_first_Piola_Kirchhoff_stress(
  Tensor P_effective_p,
  double theta_p,
  Tensor F_n1_p)
/*
  This function returns : P = P' - theta*F^{-T}
*/
{
  int Ndim = NumberDimensions;

  Tensor P_p = alloc__TensorLib__(2);
  Tensor inverse_Fm1_n1_p = Tensor Inverse__TensorLib__(F_n1_p);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      P_p.N[i][j] = P_effective_p.N[i][j] - theta_p*Fm1_n1_p.N[j][i];
    }
  }

  free__TensorLib__(Fm1_n1_p);

  return P_p;
}

/**********************************************************************/

static Matrix compute_Reactions(
  Mesh FEM_Mesh,
  Matrix Forces,
  Mask ActiveNodes)
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

  Matrix Reactions = allocZ__MatrixLib__(Nnodes_mask, Ndim);
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

static void solve_Nodal_Equilibrium_Mixture(
  Matrix Acceleration,
  Matrix Mass_Matrix_Mixture,
  Matrix Gravity_field,
  Matrix Internal_Forces_Mixture)
/*
  Solve equilibrium equation to get the nodal values of the aceleration 
  at t = n + 1
*/
{
  
  /*
    Call the LAPACK solver to compute the accelerations and velocities
  */
  int Ndim = NumberDimensions;
  int Nnodes = Acceleration.N_rows;
  int Order = Nnodes*Ndim;
  int LDA = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO = 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Mass_Matrix_Mixture.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the sistem
  */
  dgetrs_(&TRANS,&Order,&NRHS,Mass_Matrix_Mixture.nV,&LDA,IPIV,Internal_Forces_Mixture.nV,&LDB,&INFO);
  free(IPIV);

  /*
    The solution is now stored in the internal forces vector
  */
  for(int A_indx = 0; A_indx<Order; A_indx++)
  {
    Acceleration.nV[A_indx] = Gravity_field.nV[A_indx] + Internal_Forces_Mixture.nV[A_indx];
  }

}

/**************************************************************/

static Matrix compute_Viscous_Forces_Fluid(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Matrix Velocity)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask*Ndim;
  int idx_B, idx_AB;

  Matrix Viscous_Forces = allocZ__MatrixLib__(Nnodes_mask,1);

  Matrix C_Viscosity = compute_C_Viscosity(ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Compute the permabily forces (Vectorized)
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      for(int B = 0 ; B<Nnodes_mask ; B++)
      {
        for(int i = 0 ; i<Ndim ; i++)
        {
          idx_B = B*Ndim + i; 
          idx_AB = A*Order + B*Ndim + i;
          Viscous_Forces.nV[A] += K_Permeability.nV[idx_AB]*Velocity.nV[idx_B];
        }
      }
    }

  free__MatrixLib__(C_Viscosity);  
   

  return Viscous_Forces;
}

/**************************************************************/

static Matrix compute_C_Viscosity(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */
  int B_mask; /* Index of the node where we apply the body force */
  int Bp; /* Tributary node B of particle p */

  Element Nodes_p; /* Element for each particle p */
  Matrix ShapeFunction_p; /* Matrix with the value of the shape function in the particle p */ 
  Matrix gradient_p; /* Matrix with the value of the shape function gradient in the particle p */ 
  Tensor F_n1_p; /* Deformation gradient of the solid skeleton in t = n + 1 */
  Tensor gradient_pB; /* Vector with the gradient of the shape function in node B for the particle p */
  double ShapeFunction_pA; /* Value of the shape funtion in node A for the particle p */
  double rho_f_p; /* Material density of the fluid phase for particle p */
  double V0_p; /* Volume of the particle */
  double J_n1_p; /* Determinant of the solid skeleton deformation gradient in t = n + 1 */

  Matrix C_Viscosity = allocZ__MatrixLib__(Nnodes_mask,Nnodes_mask*Ndim); /* Nodal viscosity */

  for(int p = 0 ; p<Np ; p++)
    {
      
      /* 
        Get the current material density for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];

      /*
        Get the reference volume for each material point (mixture) 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Define nodes for each particle
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

      /*
        Evaluate the shape function and its gradient in the coordinates of the particle 
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
        
      /*
         Get the value of the deformation gradient at t = n + 1. 
      */
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

      /*
        Compute the jacobian for each material point
      */
      J_n1_p = I3__TensorLib__(F_n1_p);

    
      for(int A = 0 ; A<NumNodes_p ; A++)
       {
      
         /*
           Compute the gradient in the reference configuration 
         */
        ShapeFunction_pA = ShapeFunction_p.nV[A];   

        /*
          Get the node of the mesh for the contribution 
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];

        for(int B = 0 ; B<NumNodes_p ; B++)
        {
          /*
            Compute the gradient in the deformed configuration
          */
          gradient_pB = memory_to_tensor__TensorLib__(gradient_p.nM[B],1);

          /*
            Get the node of the mesh for the contribution 
          */
          Bp = Nodes_p.Connectivity[B];
          B_mask = ActiveNodes.Nodes2Mask[Bp];

          /*
            Fill the viscosity matrix
          */
          for(int i = 0 ; i<Ndim ; i++)
          {
            C_Viscosity.nM[A_mask][B_mask*Ndim + i] += ShapeFunction_pA*J_n1_p*rho_f_p*gradient_pB.n[i]*V0_p;
          }


        }
    

      }
        
      /* 
   Free memory 
      */
      free__MatrixLib__(ShapeFunction_p);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
    }

    return C_Viscosity;
}

/**************************************************************/

static Matrix compute_Permeability_Forces_Fluid(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Matrix Pore_water_pressure)
/*

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;

  Matrix Permeability_Forces = allocZ__MatrixLib__(Nnodes_mask, 1);
  Matrix K_Permeability = compute_K_Permeability(ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Compute the permabily forces (Vectorized)
  */
  for(int idx_A = 0 ; idx_A<Nnodes_mask ; idx_A++)
    {
      for(int idx_B = 0 ; idx_B<Nnodes_mask ; idx_B++)
      {
        idx_AB = idx_A*Nnodes_mask + idx_B;
        Permeability_Forces.nV[idx_A] += K_Permeability.nV[idx_AB]*Pore_water_pressure.nV[idx_B];
      }
    }

  free__MatrixLib__(K_Permeability);

  return Permeability_Forces;
}

/**************************************************************/

static Matrix compute_K_Permeability(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh)
/*

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;

  Matrix K_Permeability = allocZ__MatrixLib__(Nnodes_mask, Nnodes_mask);

  Element Nodes_p; /* Element for each particle */
  Matrix gradient_p;  /* Shape functions gradients */
  Tensor gradient_pA; /* Shape functions gradients (Node A), def config */
  Tensor gradient_pB; /* Shape functions gradients (Node B), def config */
  Tensor GRADIENT_pA; /* Shape functions gradients (Node A), ref config */
  Tensor GRADIENT_pB; /* Shape functions gradients (Node B), ref config */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor transpose_F_n_p; /*  /* Transpose of the deformation gradient t = n */
  Tensor k_p; /* Spatial permebility tensor */
  Tensor K_p; /* Material permeability tensor */
  double g = 9.81;
  double rho_f_p;
  double V0_p; /* Volume of the particle at the reference configuration */
  double J_n1_p; /* Determinant of the deformation gradient at t = n + 1 */
  double aux1_K_Permeability; /* Intermediate result */
  double aux2_K_Permeability; /* Intermediate result */

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
        Compute gradient of the shape function in each node 
      */
      gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Take the values of the deformation gradient at t = n and t = n + 1. 
        Later compute the transpose of the deformation gradient and the determinant.
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);
      J_n1_p = I3__TensorLib__(F_n1_p);
        
      /* 
        Get the current material density for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];

      /*
        Get the reference volume for each material point (mixture) 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Describe the permeability tensor using the reference configuration 
        of the soil
      */
      k_p = MPM_Mesh.Mat.Permeability;
      K_p = alloc__TensorLib__(2);
      contravariant_pull_back_tensor__TensorLib__(K_p, k_p, F_n1_p);


      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
      {

        /* 
          Get the node in the mass matrix with the mask
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];
    
        /*
          Compute the gradient in the reference configuration for the node A
        */
        gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
        GRADIENT_pA = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pA);
  

          for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
          {       
            /* 
              Get the node in the mass matrix with the mask
            */
            Bp = Nodes_p.Connectivity[B];
            B_mask = ActiveNodes.Nodes2Mask[Bp];

            /*
              Compute the gradient in the reference configuration 
            */
            gradient_pB = memory_to_tensor__TensorLib__(gradient_p.nM[B], 1);
            GRADIENT_pB = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pA);

            /*
              Compute (1/gJ)*GRAD(N_A)*K*GRAD(N_B)
            */
            aux1_K_Permeability = 0;

            for(int i = 0 ; i<Ndim ; i++)
            {

              aux2_K_Permeability = 0;

              for(int j = 0 ; j<Ndim ; j++)
              {
                aux2_K_Permeability += K_p.N[i][j]*GRADIENT_pB.n[j];
              }

              aux1_K_Permeability += GRADIENT_pA.n[i]*aux2_K_Permeability;

            }

            /*
              Assign the calculated value to the permeability matrix
            */
            K_Permeability.nM[A_mask][B_mask] += (V0_p/(g*J_n1_p))*aux1_K_Permeability;


           free__TensorLib__(GRADIENT_pB);

          }

          free__TensorLib__(GRADIENT_pA);
      }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(gradient_p);
      free__TensorLib__(transpose_F_n_p);
      free__TensorLib__(K_p);
      free(Nodes_p.Connectivity);      

    }

  return K_Permeability;
}

/**************************************************************/

static Matrix compute_Permeability_Inertial_Forces_Fluid(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Matrix Acceleration,
  Matrix Gravity)
/*

*/
{
  
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask*Ndim;
  int idx_B, idx_AB;

  Matrix Permeability_Inertial_Forces = allocZ__MatrixLib__(Nnodes_mask,1);

  Matrix C_Permeability = compute_C_Permeability(ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Compute the permabily inertial forces (Vectorized)
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      for(int B = 0 ; B<Nnodes_mask ; B++)
      {
        for(int i = 0 ; i<Ndim ; i++)
        {
          idx_B = B*Ndim + i; 
          idx_AB = A*Order + B*Ndim + i;
          Permeability_Inertial_Forces.nV[A] += C_Permeability.nV[idx_AB]*(Acceleration.nV[idx_B] - Gravity.nV[idx_B]);
        }
      }
    }

  free__MatrixLib__(C_Permeability);  
   

  return Permeability_Inertial_Forces;
}

/**************************************************************/

static Matrix compute_C_Permeability(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh)
/*

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */
  int B_mask; /* Index of the node where we apply the body force */
  int Bp; /* Tributary node B of particle p */

  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p;
  Matrix gradient_p;  /* Shape functions gradients */
  Tensor gradient_pA; /* Shape functions gradients (Node A), def config */
  Tensor GRADIENT_pA; /* Shape functions gradients (Node A), ref config */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor transpose_F_n_p; /* Transpose of the deformation gradient t = n */
  Tensor inverse_F_n1_p; /* Inverse of the deformation gradient t = n + 1 */
  Tensor k_p; /* Spatial permebility tensor */
  Tensor k_prim_p; /* Two-point permeability tensor */
  Tensor transpose_k_prim_p; /* Transpose of the two-point permeability tensor */
  Tensor GRADIENT_pA_k_prim; /* Intermediate result */
  double ShapeFunction_pB; /* Nodal value of the shape function in node B */
  double rho_f_p; /* Intrinsic or material density (fluid phase) */
  double g = 9.81; /* Gravity constant */
  double V0_p; /* Inital volume of the particle p */
  double J_n1_p; /* Determianant of the soil skeleton deformation gradient at t = n + 1 */

  Matrix C_Permeability = allocZ__MatrixLib__(Nnodes_mask,Nnodes_mask*Ndim);

  for(int p = 0 ; p<Np ; p++)
    {
      
      /* 
        Get the current material density for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];

      /*
        Get the reference volume for each material point (mixture) 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Define nodes for each particle
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

      /*
        Evaluate the shape function in the coordinates of the particle 
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Compute gradient of the shape function in each node 
      */
      gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
        
      /*
        Take the value of the deformation gradient at t = n + 1, and t = n. 
        Get intermediate results
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);
      inverse_F_n1_p = Inverse__TensorLib__(F_n1_p);
      J_n1_p = I3__TensorLib__(F_n1_p);

      /*
        Compute the tensor k' and its transpose
      */
      k_prim_p = matrix_product__TensorLib__(inverse_F_n1_p,k_prim_p);
      for(int i = 0 ; i<Ndim ; i++)
      {
        for(int j = 0 ; j<Ndim ; j++)
        {
          k_prim_p.N[i][j] *= J_n1_p;
        }
      }
      transpose_k_prim_p = transpose__TensorLib__(k_prim_p);

    
      for(int A = 0 ; A<NumNodes_p ; A++)
       {

        /*
          Get the gradient of the shape function and multiply it by the permeability tensor
        */
        gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A],1);
        GRADIENT_pA = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pA);
        GRADIENT_pA_k_prim = vector_linear_mapping__TensorLib__(transpose_k_prim_p,GRADIENT_pA);
     
        /*
          Get the node of the mesh for the contribution 
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];

        for(int B = 0 ; B<NumNodes_p ; B++)
        {
          /*
            Get the shape function in the node B
          */
          ShapeFunction_pB = ShapeFunction_p.nV[B];  

          /*
            Get the node of the mesh for the contribution 
          */
          Bp = Nodes_p.Connectivity[B];
          B_mask = ActiveNodes.Nodes2Mask[Bp];

          /*
            Fill the viscosity matrix
          */
          for(int i = 0 ; i<Ndim ; i++)
          {
            C_Permeability.nM[A_mask][B_mask*Ndim + i] += GRADIENT_pA_k_prim.n[i]*(rho_f_p/g)*ShapeFunction_pB*V0_p;
          }


        }
  
        free__TensorLib__(GRADIENT_pA);
        free__TensorLib__(GRADIENT_pA_k_prim);

      }
        
      /* 
       Free memory 
      */
      free__TensorLib__(transpose_F_n_p);
      free__TensorLib__(inverse_F_n1_p);
      free__TensorLib__(k_prim_p);
      free__TensorLib__(transpose_k_prim_p);
      free__MatrixLib__(ShapeFunction_p);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
    }

    return C_Permeability;
}

/**************************************************************/

static void solve_Nodal_Generalized_Darcy_Law(
  Matrix Rate_Pore_water_pressure,
  Matrix Compressibility_Matrix_Fluid,
  Matrix Viscous_Forces_Fluid,
  Matrix Permeability_Forces_Fluid,
  Matrix Permeability_Inertial_Forces_Fluid)
/*

*/
{
  int Nnodes_mask = Viscous_Forces_Fluid.N_rows;
  int Order = Nnodes_mask;
  int LDA   = Nnodes_mask;
  int LDB   = Nnodes_mask;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO = 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  Matrix Fluid_Forces = allocZ__MatrixLib__(Nnodes_mask,1);

  for(int A = 0 ; A<Order ; A++)
  {
    Fluid_Forces.nV[A] = 
    - Viscous_Forces_Fluid.nV[A]
    + Permeability_Forces_Fluid.nV[A] 
    + Permeability_Inertial_Forces_Fluid.nV[A];
  }

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,Compressibility_Matrix_Fluid.nV,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descompistion  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
        "Error in solve_non_reducted_system",
        "The function",
        "dgetrf_",
        "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Compressibility_Matrix_Fluid.nV,&LDA,IPIV,Fluid_Forces.nV,&LDB,&INFO);
  free(IPIV);
  
  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
        "Error in solve_non_reducted_system",
        "The function",
        "dgetrs_",
        "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    The solution is now stored in the fluid forces vector
  */
  for(int A_indx = 0; A_indx<Order; A_indx++)
  {
    Rate_Pore_water_pressure.nV[A_indx] = Fluid_Forces.nV[A_indx];
  }

  /*
    Free auxiliar vector
  */
  free__MatrixLib__(Fluid_Forces);

  return Fluid_Forces;
}

/**************************************************************/

static void compute_Explicit_Newmark_Corrector(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Matrix D_Displacement,
  Matrix Acceleration,
  Matrix Rate_Pore_water_pressure,
  Mask ActiveNodes,
  double gamma,
  double Dt)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  int idx_ij;
  Tensor F_n_p;
  Tensor F_n1_p;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  double ShapeFunction_pI; /* Nodal value for the particle */
  double mass_I;
  double D_U_pI; /* Increment of displacement */
  Element Nodes_p; /* Element for each particle */

  /* iterate over the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /*
        Replace the deformation gradient at t = n with the new one
      */
      F_n_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

      for(int i = 0 ; i<Ndim  ; i++)
       {
         for(int j = 0 ; j<Ndim  ; j++)
           {
             F_n_p.N[i][j] = F_n1_p.N[i][j];
           }
       }

       /*
        Set to zero the particle acceleration to get later it vale from  
          direct interpolation using nodal values
        */
       for(int i = 0 ; i<Ndim ; i++)
       {
          MPM_Mesh.Phi.acc.nM[p][i] = 0.0;
       }

      /*
        Define element of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /*
	     Evaluate the shape function in the coordinates of the particle 
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      
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

    	      /* Update the particles accelerations */
    	      MPM_Mesh.Phi.acc.nM[p][i] += ShapeFunction_pI*Acceleration.nV[idx_A_mask_i];
    	      /* Update the particles velocities */
    	      MPM_Mesh.Phi.vel.nM[p][i] += gamma*ShapeFunction_pI*Dt*Acceleration.nV[idx_A_mask_i];
            /* Compute the nodal contribution of the increment of displacement */
            D_U_pI = ShapeFunction_pI*D_Displacement.nV[idx_A_mask_i];
            /* Update the particle displacement */
            MPM_Mesh.Phi.dis.nM[p][i] += D_U_pI;
    	      /* Update the particles position */
    	      MPM_Mesh.Phi.x_GC.nM[p][i] += D_U_pI;
	       } 

         /*
          Update pressure field
         */
          MPM_Mesh.Phi.Pw.nM[p] += gamma*ShapeFunction_pI*Dt*Rate_Pore_water_pressure.nV[A_mask];

      	}

      /*
	Free memory
      */
      free(Nodes_p.Connectivity);
      free__MatrixLib__(gradient_p);
    }  
}

/**************************************************************/

static void output_selector(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Velocity,
  Matrix D_Displacement,
  Matrix Forces,
  Matrix Reactions,
  int TimeStep,
  int ResultsTimeStep)
/*

*/
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
