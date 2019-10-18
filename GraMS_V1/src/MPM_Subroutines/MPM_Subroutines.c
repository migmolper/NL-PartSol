#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MeshTools/MeshTools.h"
#include "../InOutFun/InOutFun.h"
#include "MPM_Subroutines.h"


/*********************************************************************/

void GlobalSearchGaussPoints(GaussPoint MPM_Mesh, Mesh FEM_Mesh){

  /* 0º Variable declaration */
  Matrix X_GC_GP;
  X_GC_GP.N_rows = NumberDimensions;
  X_GC_GP.N_cols = 1;  
  X_GC_GP.n = NAN;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;  
  X_EC_GP.n = NAN;
  Matrix Poligon = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  int * Poligon_Connectivity;

  /* 1º Set to zero the active/non-active elements */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }
  
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      Poligon_Connectivity = FEM_Mesh.Connectivity[j];
      /* 4º Fill the poligon Matrix */
      for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	for(int l = 0 ; l<NumberDimensions ; l++){
	  Poligon.nM[k][l] =
	    FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][l];
	}
      }
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon) == 1){

	/* 6º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.Element_id[i] = j;

	/* 7º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon);

	/* 8º Create a connectivity for the GP */
	for(int k = 0 ; k<MPM_Mesh.Nodes.nV[i] ; k++){
	  MPM_Mesh.Nodes.nM[i][k] = Poligon_Connectivity[k];
	}
	
	/* 9º Active those nodes that interact with the GP */
	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	  FEM_Mesh.ActiveNode[Poligon_Connectivity[k]] += 1;
	}
	
      }      
    } /* Loop over the elements */

  } /* Loop over the GP */

  /* 10º Free memory */
  FreeMat(Poligon);
  
}

/*********************************************************************/

void LocalSearchGaussPoints(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*
  Local search algorithm based on the velocity of the particle
*/
{

  /* 0º Variable declaration */
  Matrix X_GC_GP;
  X_GC_GP.N_rows = NumberDimensions;
  X_GC_GP.N_cols = 1;
  X_GC_GP.n = NAN;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  X_EC_GP.n = NAN;
  Matrix Element_GP_Coordinates =
    MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  int * Element_GP_Connectivity;
  int Element_GP_i; /* Index of the initial element */
  Matrix V_GP; /* Velocity array */
  V_GP.N_rows = NumberDimensions;
  V_GP.N_cols = 1;
  V_GP.nM = NULL;
  V_GP.n = NAN;
  strcpy(V_GP.Info,"V_GP");
  Matrix Search_Direction ; /* Auxiliar array for the local search */
  Search_Direction = MatAllocZ(1,NumberDimensions);
  Matrix V_GP_n;
  V_GP_n = MatAllocZ(1,NumberDimensions);
  int SearchVertex; /* Index to start the search */
  int * SearchList; /* Pointer to store the search list */
  
  /* 1º Set to zero the active/non-active elements */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }

  /* 2º Loop over the GP */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Get the global coordinate of the GP */
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    /* 4º Get the index of the initial element */
    Element_GP_i = MPM_Mesh.Element_id[i];

    /* 5º Get the connectivity of the initial element  */
    Element_GP_Connectivity = MPM_Mesh.Nodes.nM[i];

    /* 6º Fill the matrix with the nodal coordinates of 
       the initial element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Element_GP_Coordinates.nM[j][k] =
	  FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[j]][k];
      }
    }

    /* 6º Check if the GP is in the same element */
    if(InOut_Poligon(X_GC_GP,Element_GP_Coordinates) == 1){
      
      /* 6bº If the GP is in the element, get its natural coordinates */
      X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];     
      Get_X_EC_Q4(X_EC_GP,X_GC_GP,Element_GP_Coordinates);
      
      /* 6aº Active those nodes that interact with the GP */
      for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
	FEM_Mesh.ActiveNode[Element_GP_Connectivity[j]] += 1;
      }
      
    }
    /* 7º If the GP is not in the same element, search in the neighbour */
    else{
      
      /* 7aº As we are in a new element we set to zero the element coordinates */
      for(int j = 0 ; j<NumberDimensions ; j++){
	MPM_Mesh.Phi.x_EC.nM[i][j] = 0;
      }
      
      /* 7bº Get the velocity vector of the GP */
      V_GP.nV = MPM_Mesh.Phi.vel.nM[i];

      /* 7bº Set to a NAN the SearchVertex in order to avoid bugs */
      SearchVertex = -999;

      /* 7dº Get the search direction for the vertex of the element
       and check the search direction */
      for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){

	if(j == 0){ /* First vertex */
	  for(int k = 0 ; k<NumberDimensions ; k++){
	    Search_Direction.nV[k] =
	      2*FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[0]][k] -
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[1]][k] -
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[FEM_Mesh.NumNodesElem-1]][k];
	  }
	}
	else if(j == FEM_Mesh.NumNodesElem-1){ /* Last vertex */
	  for(int k = 0 ; k<NumberDimensions ; k++){
	    Search_Direction.nV[k] =
	      2*FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[FEM_Mesh.NumNodesElem-1]][k] -
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[0]][k] -
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[FEM_Mesh.NumNodesElem-2]][k];
	  }
	}
	else{ /* The rest of the elements */
	  for(int k = 0 ; k<NumberDimensions ; k++){
	    Search_Direction.nV[k] =
	      2*FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[j]][k] -
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[j+1]][k] -
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[j-1]][k];
	  }
	}

	for(int k = 0 ; k<NumberDimensions ; k++){
	  V_GP_n.nV[k] = V_GP.nV[k]*Search_Direction.nV[k];
	}

	if( (V_GP_n.nV[0] >= 0) && (V_GP_n.nV[1] >= 0)){
	  SearchVertex = Element_GP_Connectivity[j];
	  break;
	}

      }
     
      /* 7eº Check for errors */
      if (SearchVertex<0 || SearchVertex>=FEM_Mesh.NumNodesMesh){
	puts("Error in LocalSearchGaussPoints() : Search algorithm fails !!! ");
	exit(0);
      }

      /* 7fº Create the search list of this vertex */
      SearchList = FEM_Mesh.NodeNeighbour[SearchVertex];
     
      /* 7gº Search in the search list */
      for(int j  = 1 ; j<(FEM_Mesh.NodeNeighbour[SearchVertex][0]+1) ; j++){

	/* Discard the initial element for the search */
	if(SearchList[j] == Element_GP_i) continue; 

	/* Connectivity of the Element in the list */
	Element_GP_Connectivity = FEM_Mesh.Connectivity[SearchList[j]];

	/* Fill the matrix with the nodal coordinates */
	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Element_GP_Coordinates.nM[k][l] =
	      FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[k]][l];
	  }
	}

	if(InOut_Poligon(X_GC_GP,Element_GP_Coordinates) == 1){

	  /* Asign to the GP a element in the background mesh, just for 
	     searching porpuses */
	  MPM_Mesh.Element_id[i] = SearchList[j];

	  /* Get its natural coordinates */
	  X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	  Get_X_EC_Q4(X_EC_GP,X_GC_GP,Element_GP_Coordinates);

	  /* Generate a GP connectivity */
	  for(int k = 0 ; k<MPM_Mesh.Nodes.nV[i] ; k++){
	    MPM_Mesh.Nodes.nM[i][k] = Element_GP_Connectivity[k];
	  }

	  /* Active those nodes that interact with the GP */
	  for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	    FEM_Mesh.ActiveNode[Element_GP_Connectivity[k]] += 1;
	  }

	  /* If this is true, stop the search */
	  break;
	}
	
      }

      if(Element_GP_i == MPM_Mesh.Element_id[i]){
	printf(" %s %i %s %i !!! \n",
	       "Error in LocalSearchGaussPoints() : GP",i,
	       "is not in the neighbours of",SearchVertex);
	exit(0);
      }
      
    }

  } /* Loop over the GP */

  /* 8º Free memory */
  FreeMat(V_GP_n);
  FreeMat(Search_Direction);
  FreeMat(Element_GP_Coordinates);

 
}

/*********************************************************************/

/* int SearchGaussPoint(int idx0, Matrix X_GP, Matrix V_GP, Mesh FEM_Mesh){ */

/*   int idx; */
  
/*   switch( idx0 >= 0 ? 0 : 1){ */
/*   case 0:  */
/*     idx = LocalSearchGaussPoints(idx0, X_GP, V_GP, FEM_Mesh); */
/*     break; */
/*   case 1: */
/*     idx = GlobalSearchGaussPoints(X_GP, FEM_Mesh); */
/*     break; */
/*   default : */
/*     puts("Error in SearchGaussPoint() : Logical error"); */
/*     exit(0); */
/*   } */
  
/*   return idx; */
/* } */

/*********************************************************************/


Matrix GetNodalMassMomentum(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  /* 0º Variable declaration */

  /* Output */
  Matrix Nodal_FIELDS;
  
  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Particle voxel */
  Matrix N_GP; /* Value of the shape-function */

  /* Element for each Gauss-Point */
  Matrix Elem_Coor;
  int Num_Nodes; /* Number of nodes */
  int * Elem_Nods; /* Connectivity of the element */

  /* 1º Allocate the output list of fields */
  Nodal_FIELDS = MatAllocZ(NumberDimensions + 1,FEM_Mesh.NumNodesMesh);
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 1º Define element of the GP */
    Num_Nodes = MPM_Mesh.Nodes.nV[i];
    Elem_Nods = MPM_Mesh.Nodes.nM[i];
    Elem_Coor = MatAllocZ(Num_Nodes,NumberDimensions); 
    for(int j = 0 ; j<Num_Nodes ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Elem_Coor.nM[j][k] =
	  FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 3º Evaluate the shape function in the coordinates of the GP */
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
    }
    else if(strcmp(FEM_Mesh.TypeElem,"GIMP2D") == 0){      
      X_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      N_GP = GIMP_2D(X_GP,lp,Elem_Coor,FEM_Mesh.DeltaX);      
    }
    /* Free memory */
    FreeMat(Elem_Coor);

    /* 8º Get the nodal mass and mommentum */
    for(int k = 0 ; k<Num_Nodes ; k++){
      /* 8aº Nodal mass */
      Nodal_FIELDS.nM[0][Elem_Nods[k]] +=
	MPM_Mesh.Phi.mass.nV[i]*
	N_GP.nV[k];
      /* 8bº Nodal momentum */
      for(int l = 0 ; l<NumberDimensions ; l++){
	Nodal_FIELDS.nM[l+1][Elem_Nods[k]] +=
	  MPM_Mesh.Phi.mass.nV[i]*
	  MPM_Mesh.Phi.vel.nM[i][l]*
	  N_GP.nV[k];
      }
    }

    /* 9º Free the value of the shape functions */
    FreeMat(N_GP);
    
  }
    
  return Nodal_FIELDS;
  
}

/*******************************************************/

Matrix GetNodalVelocity(Mesh FEM_Mesh,
			Matrix Nodal_MOMENTUM,
			Matrix Nodal_MASS){
/*
  Get the nodal velocity using : 
  v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
  Initialize nodal velocities 
*/
  Matrix Vel_Mesh;
  Vel_Mesh = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Vel_Mesh.Info,"VELOCITY");
  
/* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      if(FEM_Mesh.ActiveNode[i] > 0){
	Vel_Mesh.nM[j][i] = Nodal_MOMENTUM.nM[j][i] / Nodal_MASS.nV[i];
      }
    }    
  }
  
return Vel_Mesh;
}

/*******************************************************/

void UpdateGaussPointStrain(GaussPoint MPM_Mesh,
			    Mesh FEM_Mesh,
			    Matrix Mesh_Vel)
/*
  Calcule the particle stress increment :

  \Delta\Epsilon_{ij,p}^{k-1/2} = 
  \frac{\Delta t}{2} \cdot
  \sum_{I=0}^{Nn}(N^{k}_{Ip,j} \cdot
  v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot 
  v_{jI}^{k-1/2})
*/
{
  /* Gauss-Point properties */
  Matrix X_GP; /* Element coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Particle voxel */
  
  /* Element for each Gauss-Point */
  int Num_Nodes; /* Number of nodes */
  int * Elem_Nods; /* Connectivity of the element */
  Matrix Elem_Coor; /* Coordinates of the nodes */

  /* Mesh variables */
  Matrix Elem_Vel; /* Array with the nodal velocities */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix B; /* B marix to get the deformation */
  Matrix Increment_Strain_GP; /* Vectoriced Strain tensor */
  double Incr_TraceStrain; /* Increment of the trace of the Stress tensor */
 
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
   
    /* 1º Define element of the GP */
    Num_Nodes = MPM_Mesh.Nodes.nV[i];    
    Elem_Nods = MPM_Mesh.Nodes.nM[i];
    Elem_Coor = MatAllocZ(Num_Nodes,NumberDimensions);
    
    for(int j = 0 ; j<Num_Nodes ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Elem_Coor.nM[j][k] =
	  FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }
    
    /* 2º Get the B matrix */
    /* Get the element gradient */
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      dNdx_GP = Get_dNdX_Q4(X_GP,Elem_Coor);
    }
    else if(strcmp(FEM_Mesh.TypeElem,"GIMP2D") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      dNdx_GP = dGIMP_2D(X_GP,lp,Elem_Coor,FEM_Mesh.DeltaX);
    }
    /* Free memory */
    FreeMat(Elem_Coor);
    /* Calcule the B matrix */
    B = Get_B_GP(dNdx_GP);
    /* Free shape-function derivatives */
    FreeMat(dNdx_GP);

    /* 3º Allocate and array with the velocities of the element */
    Elem_Vel = MatAllocZ(Num_Nodes*NumberDimensions,1);

    /* 4º Get the nodal velocities in the element */
    for(int j = 0 ; j<Num_Nodes ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Elem_Vel.nV[j*NumberDimensions + k] =
	  Mesh_Vel.nM[k][Elem_Nods[j]];
      }
    }

    /* 5º Multiply B by the velocity array and by the time step to get
       the increment stress tensor */
    Increment_Strain_GP = Scalar_prod(B,Elem_Vel);
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      Increment_Strain_GP.nV[j] *= DeltaTimeStep;
    }
    /* Free memory */
    FreeMat(Elem_Vel);
    FreeMat(B);

    /* 6º Set to zero the trace of the stress tensor */
    Incr_TraceStrain = 0;

    /* 7º Udate the Gauss-Point strain tensor */
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      MPM_Mesh.Phi.Strain.nM[i][j] += Increment_Strain_GP.nV[j];
      /* 7º Get the trace of the stress tensor */
      if(j<NumberDimensions)
	Incr_TraceStrain += Increment_Strain_GP.nV[j];
    }
    FreeMat(Increment_Strain_GP);

    /* 8º Update the density of the GP */
    MPM_Mesh.Phi.rho.nV[i] =
      UpdateGaussPointDensity(MPM_Mesh.Phi.rho.nV[i],
			      Incr_TraceStrain);

    
  }

}


/*******************************************************/


double UpdateGaussPointDensity(double rho_n,
			       double Incr_TraceStrain){

  /* 1º Density for the next step */
  double rho_n1;

  /* 2º Update the density */
  rho_n1 = (double)rho_n/(1 + Incr_TraceStrain);

  /* 3º Return density updated */
  return rho_n1;  
}


/*******************************************************/

void UpdateGaussPointStress(GaussPoint MPM_Mesh){

  /* 0º Variable declaration  */
  Matrix StrainTensor_GP;
  Matrix StressTensor_GP;
  Matrix D = MPM_Mesh.D;

  /* 1º Switch the dimensions of the aulixiar strain tensor */
  switch(NumberDimensions){
  case 1:
    StrainTensor_GP.N_rows = 1;
    StrainTensor_GP.N_cols = 1;
    StrainTensor_GP.nM = NULL;
    break;
  case 2:
    StrainTensor_GP.N_rows = 3;
    StrainTensor_GP.N_cols = 1;
    StrainTensor_GP.nM = NULL;
    break;
  case 3:
    StrainTensor_GP.N_rows = 6;
    StrainTensor_GP.N_cols = 1;
    StrainTensor_GP.nM = NULL;
    break;
  default :
    puts("Error in UpdateGaussPointStress() : Wrong number of dimensions !!! ");
    exit(0);
  }

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    /* 3º Store in an auxiliar variable the strain tensor in the GP */
    StrainTensor_GP.nV = MPM_Mesh.Phi.Strain.nM[i];
    /* 4º Get the new stress tensor */
    StressTensor_GP = Scalar_prod(D,StrainTensor_GP);
    /* 5º Update the stress tensor with the new-one */
    for(int j = 0 ; j<StrainTensor_GP.N_rows ; j++){
      MPM_Mesh.Phi.Stress.nM[i][j] = StressTensor_GP.nV[j];
    }
    /* 6º Free memory */
    FreeMat(StressTensor_GP);
  }
  
}

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep)
{
  
  /* 0º Auxiliar variable declaration */

  /* Properties of each Gauss-Point */
  Matrix X_GP; /* Coordinate for each Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Particle voxel */
  double Vol_GP; /* Gauss-Point volumen */

  /* Mesh properties evaluated in Gauss-Point coords */
  Matrix N_GP; /* Matrix with the nodal shape functions */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix B, B_T;

  /* Element for each Gauss-Point */
  int Num_Nodes; /* Number of nodes */
  int * Elem_Nods; /* Connectivity of the element */
  Matrix Elem_Coor; /* Coordinates of the nodes */

  /* Stress tensor of a Gauss-Point */
  Matrix StressTensor_GP; 
  StressTensor_GP.N_rows = MPM_Mesh.Phi.Stress.N_cols;
  StressTensor_GP.N_cols = 1;
  StressTensor_GP.nM = NULL;

  /* Internal forces for each node in a element (by a GP)*/  
  Matrix Div_Stress_Tensor;

  /* Total forces */
  Matrix Nodal_TOT_FORCES; 
  Nodal_TOT_FORCES = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES");

  /* Contact forces */
  Matrix Contact_Forces_t;
  Contact_Forces_t = MatAllocZ(NumberDimensions,MPM_Mesh.NumGP);

  /* Body Forces */
  Matrix Body_Forces_t;
  Body_Forces_t = MatAllocZ(NumberDimensions,MPM_Mesh.NumGP);
  int GP_Force;
  
  /* 1º Fill matrix with the body forces for TimeStep */ 
  if(MPM_Mesh.B.NumLoads>0){
    for(int i = 0 ; i<MPM_Mesh.B.NumLoads; i++){
      for(int j = 0 ; j<MPM_Mesh.B.Load_i[i].NumNodes ; j++){
	GP_Force = MPM_Mesh.B.Load_i[i].Nodes[j];
	for(int k = 0 ; k<MPM_Mesh.B.Load_i[i].Dim ; k++){
	  if( (MPM_Mesh.B.Load_i[i].Dir[k] == 1) ||
	      (MPM_Mesh.B.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > MPM_Mesh.B.Load_i[i].Value[k].Num)){
	      puts("Error in GetNodalForces() : The time step is out of the curve !!");
	      exit(0);
	    }
	    Body_Forces_t.nM[k][GP_Force] +=
	      MPM_Mesh.B.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)MPM_Mesh.B.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  
  /* 2º Fill matrix with the contact forces for TimeStep */
  if(MPM_Mesh.F.NumLoads>0){
    for(int i = 0 ; i<MPM_Mesh.F.NumLoads; i++){
      for(int j = 0 ; j<MPM_Mesh.F.Load_i[i].NumNodes ; j++){
	GP_Force = MPM_Mesh.F.Load_i[i].Nodes[j];
	for(int k = 0 ; k<MPM_Mesh.F.Load_i[i].Dim ; k++){
	  if( (MPM_Mesh.F.Load_i[i].Dir[k] == 1) ||
	      (MPM_Mesh.F.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > MPM_Mesh.F.Load_i[i].Value[k].Num)){
	      puts("Error in GetNodalForces() : The time step is out of the curve !!");
	      exit(0);
	    }
	    Contact_Forces_t.nM[k][GP_Force] +=
	      MPM_Mesh.F.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)MPM_Mesh.F.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  /* 3º Iterate over all the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 4º Define element of the GP */
    Num_Nodes = MPM_Mesh.Nodes.nV[i];
    Elem_Nods = MPM_Mesh.Nodes.nM[i];
    Elem_Coor = MatAllocZ(Num_Nodes,NumberDimensions); 
    for(int j = 0 ; j<Num_Nodes ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Elem_Coor.nM[j][k] =
	  FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }
    
    /* 5º Evaluate the shape function in the GP */	
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
      dNdx_GP = Get_dNdX_Q4(X_GP,Elem_Coor);
    }
    else if(strcmp(FEM_Mesh.TypeElem,"GIMP2D") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      N_GP = GIMP_2D(X_GP,lp,Elem_Coor,FEM_Mesh.DeltaX);
      dNdx_GP = dGIMP_2D(X_GP,lp,Elem_Coor,FEM_Mesh.DeltaX);
    }
    /* Free memory */
    FreeMat(Elem_Coor);
       
    /* 6º Get the B_T matrix for the derivates */
    B = Get_B_GP(dNdx_GP);
    FreeMat(dNdx_GP);
    B_T = Transpose_Mat(B);
    FreeMat(B);

    /* 7º Asign to an auxiliar variable the value of the stress tensor */
    StressTensor_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

    /* 8º Get the divergence stress tensor evaluates in the Gauss-Point 
     and free the B_T matrix */
    Div_Stress_Tensor = Scalar_prod(B_T,StressTensor_GP);
    FreeMat(B_T);

    /* 9º Calcule the volumen of the Gauss-Point */
    Vol_GP = MPM_Mesh.Phi.mass.nV[i]/MPM_Mesh.Phi.rho.nV[i];

    /* 10º Acumulate this forces to the total array with the internal forces */
    for(int k = 0 ; k<Num_Nodes ; k++){  
      for(int l = 0 ; l<NumberDimensions ; l++){

	/* 10aº Add the internal forces */
	Nodal_TOT_FORCES.nM[l][Elem_Nods[k]] -=
	  Div_Stress_Tensor.nV[k*NumberDimensions + l]*
	  Vol_GP;

	/* 10bº Add the body forces */
	Nodal_TOT_FORCES.nM[l][Elem_Nods[k]] +=
	  MPM_Mesh.Phi.mass.nV[i]*
	  N_GP.nV[k]*
	  Body_Forces_t.nM[l][i];
	
	/* 10cº Add the contact forces */
	Nodal_TOT_FORCES.nM[l][Elem_Nods[k]] +=
	  N_GP.nV[k]*
	  Contact_Forces_t.nM[l][i]*
	  Vol_GP;
	
      }
    }
	
    /* 10dº Free memory */
    FreeMat(Div_Stress_Tensor);
    FreeMat(N_GP);

  }

  /* 11º Free memory */
  FreeMat(Contact_Forces_t);
  FreeMat(Body_Forces_t);
  
  return Nodal_TOT_FORCES;
  
}

/*******************************************************/

void UpdateGridNodalMomentum(Mesh FEM_Mesh,
			     Matrix Nodal_MOMENTUM,
			     Matrix Nodal_TOT_FORCES)
{
  /* Update the grid nodal momentum */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      if(FEM_Mesh.ActiveNode[i] > 0){
	Nodal_MOMENTUM.nM[j][i] +=
	  DeltaTimeStep*Nodal_TOT_FORCES.nM[j][i];
      }
    }
  } 
}

/*******************************************************/

void UpdateVelocityAndPositionGP(GaussPoint MPM_Mesh,
				 Mesh FEM_Mesh,
				 Matrix Nodal_MASS,
				 Matrix Nodal_MOMENTUM,
				 Matrix Nodal_TOT_FORCES){

  /* 0º Variable declaration */

  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Particle voxel */
  Matrix N_GP; /* Value of the shape-function in the GP */

  /* Element for each GP */
  Matrix Elem_Coor;
  int Num_Nodes; /* Number of nodes for the GP */
  int * Elem_Nods; /* Connectivity of the element */

  /* Mesh properties */

  /* 1º iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Define element of the GP */
    Num_Nodes = MPM_Mesh.Nodes.nV[i];
    Elem_Nods = MPM_Mesh.Nodes.nM[i];
    Elem_Coor = MatAllocZ(Num_Nodes,NumberDimensions); 
    for(int j = 0 ; j<Num_Nodes ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Elem_Coor.nM[j][k] =
	  FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 3º Evaluate the shape function in the coordinates of the GP */
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
    }
    else if(strcmp(FEM_Mesh.TypeElem,"GIMP2D") == 0){      
      X_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      N_GP = GIMP_2D(X_GP,lp,Elem_Coor,FEM_Mesh.DeltaX);      
    }

    /* 4º Iterate over the nodes of the element */
    for(int j = 0 ; j<Num_Nodes ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){

	/* 4aº Update the GP velocities */
	MPM_Mesh.Phi.vel.nM[i][k] +=
	  DeltaTimeStep*N_GP.nV[j]*
	  Nodal_TOT_FORCES.nM[k][Elem_Nods[j]]/
	  Nodal_MASS.nV[Elem_Nods[j]];
	
	/* 4bº Update the GP position */
	MPM_Mesh.Phi.x_GC.nM[i][k] +=
	  DeltaTimeStep*N_GP.nV[j]*
	  Nodal_MOMENTUM.nM[k][Elem_Nods[j]]/
	  Nodal_MASS.nV[Elem_Nods[j]];
      }      
    }

    /* 5º Free memory */
    FreeMat(N_GP);
    FreeMat(Elem_Coor);
    
  }  
}


/*******************************************************/


  
