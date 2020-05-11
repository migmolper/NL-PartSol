#include "nl-partsol.h"

/*********************************************************************/

Matrix compute_BodyForces(Matrix F_I, GaussPoint MPM_Mesh,
			  Mesh FEM_Mesh, int TimeStep)
{
  int Ndim = NumberDimensions;
  Load * B = MPM_Mesh.B;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pI;
  Tensor b = alloc_Tensor(1); /* Body forces vector */
  double m_p; /* Mass of the Gauss-Point */
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  int NumNodesLoad;
  int p;
  int Ip;
  int Nn; /* Number of nodes of each Gauss-Point */

  for(int i = 0 ; i<NumBodyForces; i++){

    NumNodesLoad = B[i].NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++){

      /* Get the index of the Gauss-Point */
      p = B[i].Nodes[j];

      /* Get the value of the mass */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Define element for each GP */
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_Element(p, MPM_Mesh.ListNodes[p], Nn);

      /* Compute shape functions */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Fill vector of body forces */
      for(int k = 0 ; k<Ndim ; k++){
	if(B[i].Dir[k]){
	  if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num)){
	    printf("%s : %s\n",
		   "Error in compute_BodyForces()",
		   "The time step is out of the curve !!");
	    exit(0);
	  }
	  b.n[k] = B[i].Value[k].Fx[TimeStep];
	}
      }

      /* Get the node of the mesh for the contribution */
      for(int I = 0 ; I<Nn ; I++){

	/* Node for the contribution */
	Ip = Nodes_p.Connectivity[I];
	
	/* Pass the value of the nodal shape function to a scalar */
	ShapeFunction_pI = ShapeFunction_p.nV[I];
	
	/* Compute External forces */
	for(int k = 0 ; k<Ndim ; k++){
	  F_I.nM[Ip][k] += ShapeFunction_pI*b.n[k]*m_p;
	}
	
      }

      /* Free the matrix with the nodal gradient of the element */
      FreeMat(ShapeFunction_p);
      free(Nodes_p.Connectivity);
	
    }

      
  }

  free_Tensor(b);
  
  return F_I;

}

/*********************************************************************/

Matrix compute_ContacForces(Matrix F_I, GaussPoint MPM_Mesh,
			    Mesh FEM_Mesh, int TimeStep)
{
  int Ndim = NumberDimensions;
  Load * F = MPM_Mesh.F;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pI;
  Tensor t = alloc_Tensor(1); /* Body forces vector */
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V_p; /* Volumen of the Gauss-Point */
  double thickness_p; /* Thickness of the Gauss-Point */
  int Mat_p; /* Index of tha material for each Gauss-Point */
  int NumContactForces = MPM_Mesh.NumNeumannBC;
  int NumNodesLoad;
  int p;
  int Ip;
  int Nn; /* Number of nodes of each Gauss-Point */

  for(int i = 0 ; i<NumContactForces; i++){

    NumNodesLoad = F[i].NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++){

      /* Get the index of the Gauss-Point */
      p = F[i].Nodes[j];

      /* Get the value of the mass */
      m_p = MPM_Mesh.Phi.mass.nV[p];
      
      /* Get the value of the density */
      rho_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the value of the volum */
      V_p = m_p/rho_p;

      /* Get the thickness of the material point */
      Mat_p = MPM_Mesh.MatIdx[p];
      thickness_p = MPM_Mesh.Mat[Mat_p].thickness;

      /* Define element for each GP */
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_Element(p, MPM_Mesh.ListNodes[p], Nn);

      /* Compute shape functions */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Fill vector of body forces */
      for(int k = 0 ; k<Ndim ; k++){
	if(F[i].Dir[k]){
	  if( (TimeStep < 0) || (TimeStep > F[i].Value[k].Num)){
	    printf("%s : %s\n",
		   "Error in compute_ContacForces()",
		   "The time step is out of the curve !!");
	  }
	  t.n[k] = F[i].Value[k].Fx[TimeStep];
	}
      }

      /* Get the node of the mesh for the contribution */
      for(int I = 0 ; I<Nn ; I++){

	/* Node for the contribution */
	Ip = Nodes_p.Connectivity[I];
	
	/* Pass the value of the nodal shape function to a scalar */
	ShapeFunction_pI = ShapeFunction_p.nV[I];
	
	/* Compute Contact forces */
	for(int k = 0 ; k<Ndim ; k++){
	  F_I.nM[Ip][k] += ShapeFunction_pI*(t.n[k]/thickness_p)*V_p;
	}
	
      }

      /* Free the matrix with the nodal gradient of the element */
      FreeMat(ShapeFunction_p);
      free(Nodes_p.Connectivity);
	
    }

      
  }

  free_Tensor(t);
  
  return F_I;

}

