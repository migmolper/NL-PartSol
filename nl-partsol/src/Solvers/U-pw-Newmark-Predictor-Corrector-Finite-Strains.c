#include "nl-partsol.h"

/*
  Call global variables
*/
Mixture * Soil_Water_Mixtures; // Structure with the properties of the sample
Event * Out_nodal_path_csv;
Event * Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

/*
  Auxiliar functions and variables
*/

static char Error_message[MAXW];
static void standard_error();

/* Step 1 */
static Matrix compute_Mass_Matrix_Mixture(GaussPoint,Mesh,Mask);
static Matrix compute_Compressibility_Matrix_Fluid(GaussPoint, Mesh, Mask);
/* Step 2 */
static  void  compute_Explicit_Newmark_Predictor(GaussPoint,double,double);
/* Step 3 */
static Matrix compute_Nodal_Gravity_field(Mask, GaussPoint, int);
static Matrix compute_Nodal_D_Displacement(GaussPoint,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Velocity(GaussPoint,Mesh,Mask,Matrix);
static Matrix compute_Nodal_Pore_water_pressure(GaussPoint,Mesh,Mask,Matrix);
static  void  impose_Dirichlet_Boundary_Conditions(Mesh,Matrix,Matrix,Mask,int);
/* Step 4 */
static  void  update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);
/* Step 5 */
static Matrix compute_Total_Forces_Mixture(Mask,GaussPoint,Mesh,int);
static  void  compute_Internal_Forces_Mixture(Matrix,Mask,GaussPoint,Mesh);
static  void  compute_Contact_Forces_Mixture(Matrix,Mask,GaussPoint,Mesh,int);
static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor,double,Tensor);
static Matrix compute_Reactions_Mixture(Mesh,Matrix,Mask);
static Matrix solve_Nodal_Equilibrium_Mixture(Matrix,Matrix,Matrix,Mask);
/* Step 6 */
static Matrix compute_Total_Forces_Fluid(Matrix,Matrix,Matrix,Matrix,Mask,GaussPoint,Mesh,int);
static  void  compute_Viscous_Forces_Fluid(Matrix,Mask,GaussPoint,Mesh,Matrix);
static Matrix compute_C_Viscosity(Mask,GaussPoint,Mesh);
static  void  compute_Permeability_Forces_Fluid(Matrix,Mask,GaussPoint,Mesh,Matrix);
static Matrix compute_K_Permeability(Mask,GaussPoint,Mesh);
static  void  compute_Permeability_Inertial_Forces_Fluid(Matrix,Mask,GaussPoint,Mesh,Matrix,Matrix);
static Matrix compute_C_Permeability(Mask,GaussPoint,Mesh);
static Matrix compute_Reactions_Fluid(Mesh,Matrix,Mask);
static Matrix solve_Nodal_Equilibrium_Fluid(Matrix,Matrix);
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
  Matrix Total_Forces_Mixture;
  Matrix Reactions_Mixture;
  Matrix Total_Forces_Fluid;
  Matrix Reactions_Fluid;
  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_Coussy__SolversLib__(MPM_Mesh, FEM_Mesh.DeltaX, 1.0);

      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      Free_and_Restricted_Dofs = generate_Mask_for_static_condensation__MeshTools__(ActiveNodes,FEM_Mesh);
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("First step : Compute mass and compressibility matrix",TimeStep);
      print_Status("WORKING ...",TimeStep);

      Mass_Matrix_Mixture = compute_Mass_Matrix_Mixture(MPM_Mesh,FEM_Mesh,ActiveNodes);

      Compressibility_Matrix_Fluid = compute_Compressibility_Matrix_Fluid(MPM_Mesh,FEM_Mesh,ActiveNodes);

      print_Status("DONE !!!",TimeStep);
   
      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Explicit Newmark predictor",TimeStep);
      print_Status("WORKING ...",TimeStep);
      
      compute_Explicit_Newmark_Predictor(MPM_Mesh, gamma, TimeStep);

      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Third step : Compute nodal magnitudes",TimeStep);
      print_Status("WORKING ...",TimeStep);

      Gravity_field = compute_Nodal_Gravity_field(ActiveNodes, MPM_Mesh, TimeStep);

      D_Displacement = compute_Nodal_D_Displacement(MPM_Mesh, FEM_Mesh, ActiveNodes, Mass_Matrix_Mixture);

      Velocity = compute_Nodal_Velocity(MPM_Mesh, FEM_Mesh, ActiveNodes, Mass_Matrix_Mixture);
      
      Pore_water_pressure = compute_Nodal_Pore_water_pressure(MPM_Mesh, FEM_Mesh, ActiveNodes, Compressibility_Matrix_Fluid);

      puts("Pore_water_pressure");
      print__MatrixLib__(Pore_water_pressure,6,1);

      impose_Dirichlet_Boundary_Conditions(FEM_Mesh,Velocity,Pore_water_pressure,ActiveNodes,TimeStep);

      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update local state",TimeStep);
      print_Status("WORKING ...",TimeStep);

      update_Local_State(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);

      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Compute equilibrium mixture",TimeStep);
      print_Status("WORKING ...",TimeStep);

      Total_Forces_Mixture = compute_Total_Forces_Mixture(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);

      Reactions_Mixture = compute_Reactions_Mixture(FEM_Mesh,Total_Forces_Mixture,ActiveNodes);

      Acceleration = solve_Nodal_Equilibrium_Mixture(Mass_Matrix_Mixture,Gravity_field,Total_Forces_Mixture,Free_and_Restricted_Dofs);

      puts("Forces");
      print__MatrixLib__(Total_Forces_Mixture,6,2);
      puts("Mass");
      print__MatrixLib__(Mass_Matrix_Mixture,12,12);
      puts("Acceleration");
      print__MatrixLib__(Acceleration,6,2);


      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Six step : Compute equilibrium fluid",TimeStep);
      print_Status("WORKING ...",TimeStep);

      Total_Forces_Fluid = compute_Total_Forces_Fluid(Velocity,Acceleration,Pore_water_pressure,Gravity_field,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

      Reactions_Fluid = compute_Reactions_Fluid(FEM_Mesh,Total_Forces_Mixture,ActiveNodes);

      Rate_Pore_water_pressure = solve_Nodal_Equilibrium_Fluid(Compressibility_Matrix_Fluid,Total_Forces_Fluid);

      puts("Rate_Pore_water_pressure");
      print__MatrixLib__(Rate_Pore_water_pressure,6,1);

      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Compute corrector",TimeStep);
      print_Status("WORKING ...",TimeStep);

      compute_Explicit_Newmark_Corrector(MPM_Mesh,FEM_Mesh,D_Displacement,Acceleration,Rate_Pore_water_pressure,ActiveNodes,gamma,TimeStep);

      local_search__Particles__(MPM_Mesh,FEM_Mesh);

      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Eight step : Output variables and reset nodal values",TimeStep);
      print_Status("WORKING ...",TimeStep);

      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, Velocity, D_Displacement,
                    Total_Forces_Mixture, Reactions_Mixture, TimeStep, ResultsTimeStep);

      /*
      	Free memory.
      */
      free__MatrixLib__(Mass_Matrix_Mixture); 
      free__MatrixLib__(Compressibility_Matrix_Fluid);
      free__MatrixLib__(Gravity_field);
      free__MatrixLib__(Velocity);
      free__MatrixLib__(D_Displacement);
      free__MatrixLib__(Total_Forces_Mixture);
      free__MatrixLib__(Total_Forces_Fluid);
      free__MatrixLib__(Reactions_Mixture);
      free__MatrixLib__(Reactions_Fluid);
      free__MatrixLib__(Acceleration);
      free__MatrixLib__(Rate_Pore_water_pressure);
      free(ActiveNodes.Nodes2Mask);
      
      print_Status("DONE !!!",TimeStep);

      if(TimeStep == 4)
      {
        exit(0);
      }

    }
  
}

/**************************************************************/

static Matrix compute_Mass_Matrix_Mixture(
  GaussPoint MPM_Mesh, // Variable with information of the particles 
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


static Matrix compute_Compressibility_Matrix_Fluid(
  GaussPoint MPM_Mesh, // Variable with information of the particles
  Mesh FEM_Mesh, // Variable with information of the nodes
  Mask ActiveNodes) // Variable with information of the active nodes
/*
  Compute the lumped compresibility matrix
*/
{
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = 1*Nnodes_mask;
  int Mixture_idx;
  int Material_Water_idx;
  int Ap;
  int A_mask;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
  double ShapeFunction_pA; /* Evaluation of the particle in the node A */
  double V0_p; /* Mass of the particle (mixture) */
  double rho_f_p; /* Material density of the fluid */
  double relative_rho_f_p; /* Relative density of the fluid */
  double phi_f_p; /* Volume fractions of the fluid */
  double K_f; /* Compressibility (fluid) */
  double compressibility_density_f_p; /* Compressibilidy density fo the fluid */
  double compressibility_density_A_p; /* Nodal contribution A of the particle p */
  /* 
    Define and allocate the lumped Compressibility matrix 
  */
  Matrix Lumped_Compressibility_Matrix_Fluid = allocZ__MatrixLib__(Order, Order);


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

      /*
        Load intrinsic properties for the fluid phase to get the compressibility
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
      MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];
      K_f = MatProp_Water_p.Compressibility;

      /* 
        Load intrinsic properties for the fluid phase to get the compressibility
      */
      K_f = MatProp_Water_p.Compressibility;

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
        Lumped_Compressibility_Matrix_Fluid.nM[A_mask][A_mask] += compressibility_density_A_p;      
     
      }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);      
      
    }


  /* Add some usefulll info */
  strcpy(Lumped_Compressibility_Matrix_Fluid.Info,"Lumped-Matrix");

  return Lumped_Compressibility_Matrix_Fluid;
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
  for(int A = 0 ; A<Order ; A++)
  {
    AB = A*Order + A;
    D_Displacement.nV[A] = D_Displacement.nV[A]/Mass_Matrix_Mixture.nV[AB];
  }


  /*
    Add some usefulll info
  */
  strcpy(D_Displacement.Info,"Nodal-D-Displacement");

  return D_Displacement;
}


/**************************************************************/

static Matrix compute_Nodal_Velocity(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Matrix Mass_Matrix_Mixture)
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
  for(int A = 0 ; A<Order ; A++)
  {
    AB = A*Order + A;
    Velocity.nV[A] = Velocity.nV[A]/Mass_Matrix_Mixture.nV[AB];
  }
 
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
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Mixture_idx;
  int Material_Water_idx;
  int Ap;
  int A_mask;
  int AB;
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
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

      /*
        Get intrinsic material properties for the fluid phase (Compressibility)
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
      MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];
      K_f = MatProp_Water_p.Compressibility;

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
          Pore_water_pressure.nV[A_mask] += compressibility_density_A_p*Pw_p;
        }

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*
    Compute the nodal velocities
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
  {
    AB = A*Nnodes_mask + A;
    Pore_water_pressure.nV[A] = Pore_water_pressure.nV[A]/Compressibility_Matrix_Fluid.nV[AB];
  }

  /*
    Add some usefulll info
  */
  strcpy(Pore_water_pressure.Info,"Nodal-Pore-water-pressure");

  return Pore_water_pressure;
}

/**************************************************************/

static void impose_Dirichlet_Boundary_Conditions(
  Mesh FEM_Mesh,
  Matrix Velocity,
  Matrix Pore_water_pressure,
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
            printf("%s : %s \n","Error in imposse_NodalMomentum()","The time step is out of the curve !!");
            exit(EXIT_FAILURE);
          }

          /* 
            Assign the boundary condition 
          */
          if(k<Ndim)
          {
            Id_BCC_mask_k = Id_BCC_mask*Ndim + k;
            Velocity.nV[Id_BCC_mask_k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];                    
          }
          else
          {
            Pore_water_pressure.nV[Id_BCC_mask] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
          }                
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
  int Mixture_idx;
  int Material_Soil_idx;
  int Material_Water_idx;
  int Nnodes_p;
  double J_n1_p; /* Jacobian of the deformation gradient at t = n + 1 */
  double K_f; /* Compressibility (fluid) */
  double rho_f_0; /* Initial density of the fluid */
  double phi_s_0; /* Initial volume fraction (solid) */
  double phi_f_0; /* Initial volume fraction (fluid) */
  Plastic_status Input_Plastic_Parameters; /* Input parameters for plasticity */
  Plastic_status Output_Plastic_Parameters; /* Output parameters for plasticity */
  Element Nodes_p; /* Element for each particle */
  Material MatProp_Soil_p; /* Variable with the material properties of the solid phase */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
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
  

  for(int p = 0 ; p<Np ; p++)
    {
      /* 
        Load material properties for each phase
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
      Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
      MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
      MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];

      /*
        Read intrinsic material properties (fluid)
      */
      K_f = MatProp_Water_p.Compressibility;
      rho_f_0 = MatProp_Water_p.rho;

      /*
        Read reference volume fraction for each phase
      */
      phi_f_0 = Soil_Water_Mixtures[Mixture_idx].phi_f_0;
      phi_s_0 = Soil_Water_Mixtures[Mixture_idx].phi_s_0;

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
        Compute the Jacobian of the deformation gradient
      */
      J_n1_p = I3__TensorLib__(F_n1_p);

      /*
        Compute the right Cauchy Green tensor
      */
      C_n1_p = right_Cauchy_Green__Particles__(F_n1_p);

      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	      integration rule.
      */
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      

      if(strcmp(MatProp_Soil_p.Type,"Saint-Venant-Kirchhoff") == 0)
        {
          S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, C_n1_p, MatProp_Soil_p);
        }
      else if(strcmp(MatProp_Soil_p.Type,"Neo-Hookean-Wriggers") == 0)
        {
          S_p = grad_energy_Neo_Hookean_Wriggers(S_p, C_n1_p, J_n1_p, MatProp_Soil_p);
        }
      else if(strcmp(MatProp_Soil_p.Type,"Von-Mises") == 0)
        {
          F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
          Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
          Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

          Output_Plastic_Parameters = finite_strains_plasticity_Von_Mises(S_p, C_n1_p, F_plastic_p, F_n1_p, 
                                                                          Input_Plastic_Parameters, MatProp_Soil_p, J_n1_p);

          MPM_Mesh.Phi.cohesion.nV[p] = Output_Plastic_Parameters.Yield_stress;
          MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;
        }
      else if((strcmp(MatProp_Soil_p.Type,"Drucker-Prager-Plane-Strain") == 0) || (strcmp(MatProp_Soil_p.Type,"Drucker-Prager-Outer-Cone") == 0))
        {
          F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
          Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
          Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

          Output_Plastic_Parameters = finite_strains_plasticity_Drucker_Prager_Sanavia(S_p, C_n1_p, F_plastic_p, F_n1_p, 
                                                                          Input_Plastic_Parameters, MatProp_Soil_p, J_n1_p);

          MPM_Mesh.Phi.cohesion.nV[p] = Output_Plastic_Parameters.Cohesion;
          MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;

        }
      else
        {
          fprintf(stderr,"%s : %s %s %s \n",
		  "Error in update_Local_State()",
		  "The material",MatProp_Soil_p.Type,"has not been yet implemnented");
          exit(EXIT_FAILURE);
        }

      /*
        Update state parameters
      */
      Pw_0 = MPM_Mesh.Phi.Pw_0.nV[p]; /* Get the initial pressure */
      Pw_n1 = MPM_Mesh.Phi.Pw.nV[p]/J_n1_p; /* From the Kirchhoff pressure compure the cauchy pore water pressure */

      MPM_Mesh.Phi.rho_f.nV[p] = rho_f_0*exp((Pw_n1-Pw_0)/K_f); /* Update the fluid density */

      MPM_Mesh.Phi.phi_s.nV[p] = phi_s_0/J_n1_p; /* Update the volume fraction of the solid phase */
      MPM_Mesh.Phi.phi_f.nV[p] = 1 - (1 - phi_f_0)/J_n1_p; /* Update the volume fraction of the fluid phase */

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

static Matrix compute_Total_Forces_Mixture(
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{
  
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask,Ndim);

  /*
    Add internal forces contribution
  */
  compute_Internal_Forces_Mixture(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh);

  /*
    Add contact forces contribution
  */
  compute_Contact_Forces_Mixture(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);
  
  return Forces;
}

/**************************************************************/


static void compute_Internal_Forces_Mixture(
  Matrix Forces,
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh)
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
  Tensor inverse_F_n1_p = Inverse__TensorLib__(F_n1_p);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      P_p.N[i][j] = P_effective_p.N[i][j] - theta_p*inverse_F_n1_p.N[j][i];
    }
  }

  free__TensorLib__(inverse_F_n1_p);

  return P_p;
}

/*********************************************************************/

static void compute_Contact_Forces_Mixture(
  Matrix Forces,
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{
  int Ndim = NumberDimensions;
  Load T_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Material MatProp_Soil_p; /* Information with properties of the soil phase */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA;
  Tensor t = alloc__TensorLib__(1); /* Body forces vector */
  double V0_p; /* Volumen of the particle in the reference configuration */
  double thickness_p; /* Thickness of the particle */
  double A0_p; /* Area of the particle in the reference configuration */ 

  int Mixture_idx; /* Index of the mixture for each particle */
  int Material_Soil_idx; /* Index of the soil properties for each particle */
  int NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  int NumNodesLoad;
  int p;
  int Ap;
  int A_mask;
  int idx_A_mask_k;
  int NumNodes_p; /* Number of nodes of each particle */

  for(int i = 0 ; i<NumContactForces; i++)
  {

    /*
      Read load i
    */
    T_i = MPM_Mesh.Neumann_Contours.BCC_i[i];

    NumNodesLoad = T_i.NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++)
    {

      /*
        Get the index of the particle
      */
      p = T_i.Nodes[j];

      /*
        Get the volume of the particle in the reference configuration 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Get the material properties of the soil phase for each particle
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
      MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];

      /*
        Get the thickness of each particle
      */
      thickness_p = MatProp_Soil_p.thickness;

      /*
        Get the area of each particle
      */
      A0_p = V0_p/thickness_p;

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
        if(T_i.Dir[k])
        {
          if( (TimeStep < 0) || (TimeStep > T_i.Value[k].Num))
          {
            sprintf(Error_message,"%s : %s",
              "Error in compute_Contact_Forces_Mixture()",
              "The time step is out of the curve !!");
            standard_error();
          }
          t.n[k] = T_i.Value[k].Fx[TimeStep];
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
          idx_A_mask_k = A_mask*Ndim + k;
          Forces.nV[idx_A_mask_k] += ShapeFunction_pA*t.n[k]*A0_p;
        }
  
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
  
    }

      
  }

  free__TensorLib__(t);

}

/**********************************************************************/

static Matrix compute_Reactions_Mixture(
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
      for(int k = 0 ; k<Ndim ; k++)
	    {

        /* 
		      Apply only if the direction is active (1) 
        */
	      if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
        {
          /* 
            Set to zero the forces in the nodes where velocity is fixed 
          */
          Id_BCC_mask_k = Id_BCC_mask*Ndim + k; 
          Reactions.nV[Id_BCC_mask_k] = Forces.nV[Id_BCC_mask_k];
          Forces.nV[Id_BCC_mask_k] = 0;
        }
	    }
    }    
  }

  return Reactions;
}

/**************************************************************/

static Matrix solve_Nodal_Equilibrium_Mixture(
  Matrix Mass_Matrix_Mixture,
  Matrix Gravity_field,
  Matrix Total_Forces_Mixture,
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
  int Nnodes = Total_Forces_Mixture.N_rows;
  int Order = Nnodes*Ndim;
  int AB_indx;

  /* 
    Output variable
  */
  Matrix Acceleration = allocZ__MatrixLib__(Nnodes,Ndim);


  /*
    The solution is now stored in the internal forces vector
  */
  for(int A_indx = 0; A_indx<Order; A_indx++)
  {
    AB_indx = A_indx*Order + A_indx;
    Acceleration.nV[A_indx] = Gravity_field.nV[A_indx] + Total_Forces_Mixture.nV[A_indx]/Mass_Matrix_Mixture.nV[AB_indx];
  }

  return Acceleration;

}


/**************************************************************/

static Matrix compute_Total_Forces_Fluid(
  Matrix Velocity,
  Matrix Acceleration,
  Matrix Pore_water_pressure,
  Matrix Gravity_field,
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{
  
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask,Ndim);

  /*
    Compute the diferrent contributions to the total forces
  */

  compute_Viscous_Forces_Fluid(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh,Velocity);

  compute_Permeability_Forces_Fluid(Forces,ActiveNodes, MPM_Mesh, FEM_Mesh, Pore_water_pressure);
  
  compute_Permeability_Inertial_Forces_Fluid(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh, Acceleration, Gravity_field);

  
  return Forces;
}

/**************************************************************/

static void compute_Viscous_Forces_Fluid(
  Matrix Forces,
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

  Matrix C_Viscosity = compute_C_Viscosity(ActiveNodes, MPM_Mesh, FEM_Mesh);

  print__MatrixLib__(C_Viscosity,C_Viscosity.N_rows,C_Viscosity.N_cols);
  
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
          Forces.nV[A] -= C_Viscosity.nV[idx_AB]*Velocity.nV[idx_B];
        }
      }
    }

  free__MatrixLib__(C_Viscosity);  
   
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

  /*
    Add some usefulll info
  */
  strcpy(C_Viscosity.Info,"C-Viscosity");

  return C_Viscosity;
}

/**************************************************************/

static void compute_Permeability_Forces_Fluid(
  Matrix Forces,
  Mask ActiveNodes,
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh,
  Matrix Pore_water_pressure)
/*

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int idx_AB;

  Matrix K_Permeability = compute_K_Permeability(ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Compute the permabily forces (Vectorized)
  */
  for(int idx_A = 0 ; idx_A<Nnodes_mask ; idx_A++)
    {
      for(int idx_B = 0 ; idx_B<Nnodes_mask ; idx_B++)
      {
        idx_AB = idx_A*Nnodes_mask + idx_B;
        Forces.nV[idx_A] += K_Permeability.nV[idx_AB]*Pore_water_pressure.nV[idx_B];
      }
    }

  free__MatrixLib__(K_Permeability);

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
  int Mixture_idx;
  int Material_Soil_idx;
  int Ap;
  int Bp;
  int A_mask;
  int B_mask;

  Matrix K_Permeability = allocZ__MatrixLib__(Nnodes_mask, Nnodes_mask);

  Material MatProp_Soil_p; /* Variable with the material propertie of solid phase for each particle */
  Element Nodes_p; /* Element for each particle */
  Matrix gradient_p;  /* Shape functions gradients */
  Tensor gradient_pA; /* Shape functions gradients (Node A), def config */
  Tensor gradient_pB; /* Shape functions gradients (Node B), def config */
  Tensor GRADIENT_pA; /* Shape functions gradients (Node A), ref config */
  Tensor GRADIENT_pB; /* Shape functions gradients (Node B), ref config */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor transpose_F_n_p; /* Transpose of the deformation gradient t = n */
  Tensor inverse_F_n1_p; /* */
  Tensor k_p; /* Spatial permebility tensor */
  Tensor Fk_p; /* Product of the inverse of the defomration gradient and the permeability tensor */
  Tensor GRADIENT_pA__x__Fk; 
  double GRADIENT_pA__x__Fk__x__gradient_pB;
  double g = 9.81;
  double rho_f_p;
  double V0_p; /* Volume of the particle at the reference configuration */

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
      inverse_F_n1_p = Inverse__TensorLib__(F_n1_p);
        
      /* 
        Get the current material density for each material point (fluid) 
      */
      rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];

      /*
        Get the reference volume for each material point (mixture) 
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /* 
        Load intrinsic properties for the solid phase to compute the permeability
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
      MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
      k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;

      /*
        Intermediate result 1
      */
      Fk_p = matrix_product__TensorLib__(inverse_F_n1_p,k_p);


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
  
      /*
        Intermediate result 2
      */
        GRADIENT_pA__x__Fk = vector_linear_mapping__TensorLib__(GRADIENT_pA,Fk_p);

          for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
          {       
            /* 
              Get the node in the mass matrix with the mask
            */
            Bp = Nodes_p.Connectivity[B];
            B_mask = ActiveNodes.Nodes2Mask[Bp];

            /*
              Compute the gradient in the current configuration 
            */
            gradient_pB = memory_to_tensor__TensorLib__(gradient_p.nM[B], 1);

        /*
            Intermediate result 3
        */
            GRADIENT_pA__x__Fk__x__gradient_pB = inner_product__TensorLib__(GRADIENT_pA__x__Fk, gradient_pB);

            /*
              Assign the calculated value to the permeability matrix
            */
            K_Permeability.nM[A_mask][B_mask] += (V0_p/g)*GRADIENT_pA__x__Fk__x__gradient_pB;

          }

          free__TensorLib__(GRADIENT_pA);
          free__TensorLib__(GRADIENT_pA__x__Fk);
      }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(gradient_p);
      free__TensorLib__(transpose_F_n_p);
      free__TensorLib__(inverse_F_n1_p);
      free__TensorLib__(Fk_p);
      free(Nodes_p.Connectivity);      
    }

  /*
    Add some usefulll info
  */
  strcpy(K_Permeability.Info,"K-Permeability");

  return K_Permeability;
}

/**************************************************************/

static void compute_Permeability_Inertial_Forces_Fluid(
  Matrix Forces,
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
          Forces.nV[A] += C_Permeability.nV[idx_AB]*(Acceleration.nV[idx_B] - Gravity.nV[idx_B]);
        }
      }
    }

  free__MatrixLib__(C_Permeability);  
   
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
  int Mixture_idx; /* Index for the material point mixture parameters */

  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p;
  Matrix gradient_p;  /* Shape functions gradients */
  Tensor gradient_pA; /* Shape functions gradients (Node A), def config */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor transpose_F_n_p; /* Transpose of the deformation gradient t = n */
  Tensor inverse_F_n1_p; /* Inverse of the deformation gradient t = n + 1 */
  Tensor k_p; /* Spatial permebility tensor */
  Tensor k_prim_p; /* Two-point permeability tensor */
  Tensor transpose_k_prim_p; /* Transpose of the two-point permeability tensor */
  Tensor gradient_pA_k_prim; /* Intermediate result */
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
        Get the permeability tensor
      */
      Mixture_idx = MPM_Mesh.MixtIdx[p];
      k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;

      /*
        Compute the tensor k' and its transpose
      */
      k_prim_p = matrix_product__TensorLib__(inverse_F_n1_p,k_p);
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
        gradient_pA_k_prim = vector_linear_mapping__TensorLib__(transpose_k_prim_p,gradient_pA);
     
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
            C_Permeability.nM[A_mask][B_mask*Ndim + i] += gradient_pA_k_prim.n[i]*(rho_f_p/g)*ShapeFunction_pB*V0_p;
          }

        }
  
        free__TensorLib__(gradient_pA_k_prim);

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


  /*
    Add some usefulll info
  */
  strcpy(C_Permeability.Info,"C-Permeability");

  return C_Permeability;
}

/**************************************************************/

static Matrix compute_Reactions_Fluid(
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
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  Matrix Reactions = allocZ__MatrixLib__(Nnodes_mask, 1);
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
        Apply only if the direction is active (1) 
      */
      if(FEM_Mesh.Bounds.BCC_i[i].Dir[Ndim] == 1)
      {
        /* 
          Set to zero the forces in the nodes where pressure is fixed 
        */
        Reactions.nV[Id_BCC_mask] = Forces.nV[Id_BCC_mask_k];
        Forces.nV[Id_BCC_mask] = 0;
      }
      
    }    
  }

  return Reactions;
}


/**************************************************************/

static Matrix solve_Nodal_Equilibrium_Fluid(
  Matrix Compressibility_Matrix_Fluid,
  Matrix Total_Forces_Fluid)
/*

*/
{
  int Nnodes_mask = Total_Forces_Fluid.N_rows;
  int Order = Nnodes_mask;
  int AB_indx;

  Matrix Rate_Pore_water_pressure = allocZ__MatrixLib__(Nnodes_mask,1);

  /*
    The solution is now stored in the fluid forces vector
  */
  for(int A_indx = 0; A_indx<Order; A_indx++)
  {
    AB_indx = A_indx*Order + A_indx;
    Rate_Pore_water_pressure.nV[A_indx] = Total_Forces_Fluid.nV[A_indx]/Compressibility_Matrix_Fluid.nV[AB_indx];
  }


  return Rate_Pore_water_pressure;
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
          MPM_Mesh.Phi.Pw.nV[p] += gamma*ShapeFunction_pI*Dt*Rate_Pore_water_pressure.nV[A_mask];

      	}

      /*
	Free memory
      */
      free(Nodes_p.Connectivity);
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

static void standard_error()
{
  fprintf(stderr,"%s !!! \n",Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/
