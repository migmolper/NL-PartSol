#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*************************************************************/

void update_LocalState(Matrix V_I, GaussPoint MPM_Mesh,
		       Mesh FEM_Mesh, double TimeStep)
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
    Strain_p = memory_to_Tensor(MPM_Mesh.Phi.Strain.nM[p], 2);
    Stress_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nn = MPM_Mesh.NumberNodes[p];
    Nodes_p = get_Element(p, MPM_Mesh.ListNodes[p], Nn);

    /* Get the velocity of the nodes of the element */
    Nodal_Velocity_p = get_Element_Field(V_I, Nodes_p);

    /* Compute gradient of the shape function in each node */
    Gradient_p = compute_ShapeFunction_Gradient(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the material properties */
    Idx_Mat_p = MPM_Mesh.MatIdx[p];
    Material_p = MPM_Mesh.Mat[Idx_Mat_p];

    /* Update Strain tensor */
    Rate_Strain_p = compute_RateOfStrain(Nodal_Velocity_p,Gradient_p);
    Strain_p = update_Strain(Strain_p, Rate_Strain_p, TimeStep);

    /* Update density field */
    MPM_Mesh.Phi.rho.nV[p] = update_Density(rho_p, TimeStep, Rate_Strain_p);

    /* Compute stress tensor */
    Stress_p = compute_Stress(Strain_p,Stress_p,Material_p);

    /* Compute deformation energy */
    MPM_Mesh.Phi.W.nV[p] = 0.5*get_innerProduct_Of(Strain_p, Stress_p);
        
    /* Free the matrix with the nodal velocity of the element */
    FreeMat(Nodal_Velocity_p);
    
    /* Free the matrix with the nodal gradient of the element */
    FreeMat(Gradient_p);
    
  }

  /* Compute damage of the particles */
  ComputeDamage(MPM_Mesh, FEM_Mesh);

}

/*******************************************************/

void ComputeDamage(GaussPoint MPM_Mesh, Mesh FEM_Mesh){

  double DeltaX = FEM_Mesh.DeltaX;
  Matrix ji = MPM_Mesh.Phi.ji;
  Matrix W = MPM_Mesh.Phi.W;
  Matrix Mass = MPM_Mesh.Phi.mass;
  Matrix Stress = MPM_Mesh.Phi.Stress;
  Matrix Strain = MPM_Mesh.Phi.Strain;
  Matrix StrainF = MPM_Mesh.Phi.StrainF;
  int * MatIdx = MPM_Mesh.MatIdx;
  Material * MatProp = MPM_Mesh.Mat;
  ChainPtr * Beps = MPM_Mesh.Beps;

  if(MPM_Mesh.Mat[0].Eigenerosion){
    UpdateBeps(MPM_Mesh,FEM_Mesh);
    EigenerosionAlgorithm(ji, W, Mass, Stress,
			  MatIdx, MatProp,
			  Beps, DeltaX);
  }

  if(MPM_Mesh.Mat[0].Eigensoftening){
    UpdateBeps(MPM_Mesh,FEM_Mesh);
    EigensofteningAlgorithm(ji, Strain, StrainF,
			    Mass, Stress, MatIdx,
			    MatProp, Beps);
  }
}

/*******************************************************/
