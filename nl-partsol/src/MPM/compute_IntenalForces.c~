#include "nl-partsol.h"

/*************************************************************/

Matrix compute_InternalForces(Matrix F_I, GaussPoint MPM_Mesh, Mesh FEM_Mesh)
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
    Stress_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nn = MPM_Mesh.NumberNodes[p];
    Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nn);

    /* Compute gradient of the shape function in each node */
    Gradient_p = compute_ShapeFunction_gradient(Nodes_p, MPM_Mesh, FEM_Mesh);

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
        
    /* Free the matrix with the nodal gradient of the element */
    FreeMat(Gradient_p);
    free(Nodes_p.Connectivity);
  }

  return F_I;
  
}

/*******************************************************/


