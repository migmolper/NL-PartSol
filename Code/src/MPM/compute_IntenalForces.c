#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*************************************************************/

Matrix compute_InternalForces(Matrix F_I, Matrix V_I,
			      GaussPoint MPM_Mesh,
			      Mesh FEM_Mesh, double TimeStep)
{
  int Ndim = 3;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix Gradient_p; /* Shape functions gradients */
  Matrix Nodal_Velocity_p; /* Velocity of the element nodes */
  Material Material_p; /* Properties of the Gauss-Point material */
  Tensor Rate_Strain_p; /* Increment of strain tensor */
  Tensor Strain_p; /*  Strain tensor */
  Tensor Stress_p; /* Stress tensor */
  Tensor Gradient_pI;
  Tensor InternalForcesDensity_Ip;
  double W_p; /* Internal energy of the Gauss-Point */
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V_p; /* Volumen of the Gauss-Point */
  int Np = MPM_Mesh.NumGP;
  int Ip;
  int Nn;
  int Idx_Mat_p;

  /* Loop in the GPs */
  for(int p = 0 ; p<Np ; p++){

    /* Get the value of the density */
    rho_p = MPM_Mesh.Phi.rho.nV[p];

    /* Get the value of the mass */
    m_p = MPM_Mesh.Phi.mass.nV[p];

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
    rho_p = update_Density(rho_p, TimeStep, Rate_Strain_p);

    /* Compute stress tensor */
    Stress_p = compute_Stress(Strain_p,Stress_p,Material_p);

    /* Compute deformation energy */
    W_p = 0.5*get_innerProduct_Of(Strain_p, Stress_p);

    /* Compute the volume of the Gauss-Point */
    V_p = m_p/rho_p;

    /* Compute nodal forces */
    for(int I = 0 ; I<Nn ; I++){
      /* Pass by reference the nodal gradient to the tensor */
      Gradient_pI = memory_to_Tensor(Gradient_p.nM[I], 1);
      
      /* Compute the nodal forces of the Gauss-Point */
      InternalForcesDensity_Ip =
	get_firstOrderContraction_Of(Stress_p, Gradient_pI);
      
      /* Get the node of the mesh for the contribution */
      Ip = Nodes_p.Connectivity[I];
      
      /* Asign the nodal forces contribution to the node */
      for(int i = 0 ; i<Ndim ; i++){
	F_I.nM[Ip][i] -= InternalForcesDensity_Ip.n[i]*V_p;
      }

      /* Free the internal forces density */
      free_Tensor(InternalForcesDensity_Ip);
    }

    /* Update memory */
    MPM_Mesh.Phi.rho.nV[p] = rho_p;
    MPM_Mesh.Phi.W.nV[p] = W_p;
    
    
    /* Free the matrix with the nodal velocity of the element */
    FreeMat(Nodal_Velocity_p);
    
    /* Free the matrix with the nodal gradient of the element */
    FreeMat(Gradient_p);
    
  }

  return F_I;
  
}

/*******************************************************/


