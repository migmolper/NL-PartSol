#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../InOutFun/InOutFun.h"


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
  free(FEM_Mesh.ActiveNode);
  FEM_Mesh.ActiveNode = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
   
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      Poligon_Connectivity = FEM_Mesh.Connectivity[j];
      /* 4º Fill the poligon Matrix */
      for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	for(int l = 0 ; l<NumberDimensions ; l++){
	  Poligon.nM[k][l] = FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][l];
	}
      }
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon) == 1){
	/* 6º If the GP is in the element, set the index of the position and 
	 update the array to set if an element is active or not */

	/* 6º If the GP is in the element, set the index of the position and 
	   update the array to set if an element is active or not */
	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	  FEM_Mesh.ActiveNode[Poligon_Connectivity[k]] += 1;
	}	
	MPM_Mesh.Element_id[i] = j;

	/* 7º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	X_EC_GP = GetNaturalCoordinates(X_EC_GP,X_GC_GP,Poligon);
      }      
    } /* Loop over the elements */

  } /* Loop over the GP */

  /* 8º Free memory */
  free(Poligon.nM);
  
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
  Matrix Element_GP_Coordinates = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
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
  free(FEM_Mesh.ActiveNode);
  FEM_Mesh.ActiveNode = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));

  /* 2º Loop over the GP */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Get the global coordinate of the GP */
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    /* 4º Get the index of the initial element */
    Element_GP_i = MPM_Mesh.Element_id[i];

    /* 5º Get the connectivity of the initial element  */
    Element_GP_Connectivity = FEM_Mesh.Connectivity[Element_GP_i];

    /* 6º Fill the matrix with the nodal coordinates of the initial element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Element_GP_Coordinates.nM[j][k] =
	  FEM_Mesh.Coordinates.nM[Element_GP_Connectivity[j]][k];
      }
    }

    /* 6º Check if the GP is in the same element */
    if(InOut_Poligon(X_GC_GP,Element_GP_Coordinates) == 1){
      /* 6aº If the GP is in the element, set the index of the position and
	 update the array to set if an element is active or not,
	 it is not necessary to update the index of the element */
      for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
	FEM_Mesh.ActiveNode[Element_GP_Connectivity[j]] += 1;
      }
      /* 6bº If the GP is in the element, get its natural coordinates */
      X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      X_EC_GP = GetNaturalCoordinates(X_EC_GP,X_GC_GP,Element_GP_Coordinates);
    }
    /* 7º If the GP is not in the same element, search in the neighbour */
    else{

      /* 7aº Get the velocity vector of the GP */
      V_GP.nV = MPM_Mesh.Phi.vel.nM[i];

      /* 7bº Set to a NAN the SearchVertex in order to avoid bugs */
      SearchVertex = -999;

      /* 7cº Get the search direction for the vertex of the element
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
      
      /* 7dº Check for errors */
      if (SearchVertex<0 || SearchVertex>=FEM_Mesh.NumNodesMesh){
	puts("Error in LocalSearchGaussPoints() : Search algorithm fails !!! ");
	exit(0);
      }

      /* 7eº Create the search list of this vertex */
      SearchList = FEM_Mesh.NodeNeighbour[SearchVertex];
      
      /* 7fº Search in the search list */
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
	  /* If the GP is in the element, set the index of the position and
	     update the array to set if an element is active or not */
	  for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	    FEM_Mesh.ActiveNode[Element_GP_Connectivity[k]] += 1;
	  }
	  MPM_Mesh.Element_id[i] = SearchList[j];
	  /* If the GP is in the element, get its natural coordinates */
	  X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	  X_EC_GP = GetNaturalCoordinates(X_EC_GP,X_GC_GP,Element_GP_Coordinates);
	  /* If this is true, stop the search */
	  break;
	}
	
      }

      if(Element_GP_i == MPM_Mesh.Element_id[i]){
	printf("Error in LocalSearchGaussPoints() : GP %i is not in the neighbours of %i !!! \n",
	       i,SearchVertex);
	exit(0);
      }
      
    }

  } /* Loop over the GP */

  /* 8º Free memory */
  free(V_GP_n.nV);
  free(Search_Direction.nV);
  free(Element_GP_Coordinates.nM);

 
}
/*********************************************************************/


Matrix GetNodalMassMomentum(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  /* 0º Variable declaration */
  Matrix Nodal_FIELDS; /* Output */
  int Elem_GP; /* Index of the element */
  int * Elem_Nods;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  Matrix N_Ref_GP;
  Matrix Elem_Coords = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  strcpy(Elem_Coords.Info,FEM_Mesh.TypeElem);

  /* 1º Allocate the output list of fields */
  Nodal_FIELDS = MatAllocZ(NumberDimensions + 1,FEM_Mesh.NumNodesMesh);
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Get the coordinates of the GP in the Element */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 4º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* 5º Nodes of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 6º Coordinates of the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<FEM_Mesh.Dimension ; k++){
	Elem_Coords.nM[j][k] = FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 7º Evaluate the shape function in the GP */
    N_Ref_GP = FEM_Mesh.N_ref(X_EC_GP);

    /* 8º Get the nodal mass and mommentum */
    for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
      /* 8aº Nodal mass */
      Nodal_FIELDS.nM[0][Elem_Nods[k]] +=
	MPM_Mesh.Phi.mass.nV[i]*
	N_Ref_GP.nV[k];
      /* 8bº Nodal momentum */
      for(int l = 0 ; l<NumberDimensions ; l++){
	Nodal_FIELDS.nM[l+1][Elem_Nods[k]] +=
	  MPM_Mesh.Phi.mass.nV[i]*
	  MPM_Mesh.Phi.vel.nM[i][l]*
	  N_Ref_GP.nV[k];
      }
    }

    /* 9º Free the value of the shape functions */
    free(N_Ref_GP.nV);
    
  }
    
  return Nodal_FIELDS;
  
}

/*******************************************************/

Matrix GetNodalVelocity(Mesh FEM_Mesh,
			Matrix Nodal_MOMENTUM,
			Matrix Nodal_MASS)
/*
  Get the nodal velocity using : 
  v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
*/
{

  /* 0º Initialize nodal velocities */
  Matrix Nodal_VELOCITY;
  Nodal_VELOCITY = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_VELOCITY.Info,"VELOCITY");

  /* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      if(FEM_Mesh.ActiveNode[i] > 0){
	Nodal_VELOCITY.nM[j][i] = Nodal_MOMENTUM.nM[j][i] / Nodal_MASS.nV[i];
      }
    }
  }

  return Nodal_VELOCITY;
}

/*******************************************************/

void UpdateGaussPointStrain(GaussPoint MPM_Mesh,
			    Mesh FEM_Mesh,
			    Matrix Nodal_VELOCITY)
/*
  Calcule the particle stress increment :

  \Delta\Epsilon_{ij,p}^{k-1/2} = 
  \frac{\Delta t}{2} \cdot
  \sum_{I=0}^{Nn}(N^{k}_{Ip,j} \cdot
  v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot 
  v_{jI}^{k-1/2})
*/
{
  /* 0º Define variables */
  Matrix X_EC_GP; /* Element coordinates of the Gauss-Point */
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  int Elem_GP; /* Index of the element of the Gauss-Point */
  int * Elem_Nods; /* Connectivity of the element */
  Matrix Elem_Coords; /* Coordinates of the nodes of the element */
  Matrix Nodal_VELOCITY_Elem; /* Array with the nodal velocities */
  Matrix B; /* B marix to get the deformation */
  Matrix Increment_Strain_GP; /* Vectoriced strain tensor */

  /* 1º Allocate auxiliar variables */
  Elem_Coords = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  strcpy(Elem_Coords.Info,FEM_Mesh.TypeElem);
  Nodal_VELOCITY_Elem = MatAllocZ(FEM_Mesh.NumNodesElem*NumberDimensions,1);
  
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Get the coordinates of the GP in the Element */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 3º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* 4º Nodes of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 5º Coordinates of each node for the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<FEM_Mesh.Dimension ; k++){
	Elem_Coords.nM[j][k] = FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 6º Get the B matrix */
    B = Get_B_GP(X_EC_GP,Elem_Coords);

    /* 7º Get the nodal velocities in the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Nodal_VELOCITY_Elem.nV[j*NumberDimensions + k] = Nodal_VELOCITY.nM[k][Elem_Nods[j]];
      }
    }

    /* 8º Multiply B by the velocity array */
    Increment_Strain_GP = Scalar_prod(B,Nodal_VELOCITY_Elem);

    /* 9º Udate the Gauss-Point strain tensor */
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      MPM_Mesh.Phi.Strain.nM[i][j] += Increment_Strain_GP.nV[j]*DeltaTimeStep;
    }

    /* 10º Free memory for the next step */
    free(B.nM);
    free(Increment_Strain_GP.nV);
    
  }

  /* 11º Free memory for the global variables */
  free(Elem_Coords.nM);
  free(Nodal_VELOCITY_Elem.nV);

}


/*******************************************************/


void UpdateGaussPointDensity(GaussPoint MPM_Mesh){

  /* 0º Define auxiliar variables */
  double TraceStrainTensor;

  /* 1º Loop over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Set to zero the trace of the strain tensor */
    TraceStrainTensor = 0;

    /* 3º Calcule the trace of the strain tensor */
    for(int j = 0 ; j<NumberDimensions ; j++){
      TraceStrainTensor += MPM_Mesh.Phi.Strain.nM[i][j];      
    }

    /* 4º Update the density */
    MPM_Mesh.Phi.rho.nV[i] *= (double)1/(1 + TraceStrainTensor);
  }
  
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
    free(StressTensor_GP.nV);
  }
  
}

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{
  
  /* 0bº Auxiliar variable declaration */
  int Elem_GP; /* Index of the element */
  int * Elem_Nods;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  Matrix N_Ref_GP;
  Matrix B, B_T;
  Matrix Elem_Coords = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  strcpy(Elem_Coords.Info,FEM_Mesh.TypeElem);
  Matrix StressTensor_GP; /* Stress tensor of a Gauss Point */
  StressTensor_GP.N_rows = MPM_Mesh.Phi.Stress.N_cols;
  StressTensor_GP.N_cols = 1;
  StressTensor_GP.nM = NULL;
  Matrix Element_INT_FORCES; /* Internal forces for each node in a element (by a GP)*/  
  Matrix Nodal_TOT_FORCES; /* Total forces */
  Nodal_TOT_FORCES = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES");

  /* 1º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Get the coordinates of the GP in the Element */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 3º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* 4º Nodes of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 5º Coordinates of the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<FEM_Mesh.Dimension ; k++){
	Elem_Coords.nM[j][k] = FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 6º Evaluate the shape function in the GP */
    N_Ref_GP = FEM_Mesh.N_ref(X_EC_GP);

    /* 7º Asign to an auxiliar variable the value of the stress tensor */
    StressTensor_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

    /* 8º Multiply the stress tensor by the mass and the inverse of the density */
    for(int k = 0 ; k<StressTensor_GP.N_rows ; k++){
      StressTensor_GP.nV[k] *= (MPM_Mesh.Phi.mass.nV[i]/MPM_Mesh.Phi.rho.nV[i]);	  
    }
	
    /* 9º Get the B_T matrix for the derivates */
    B = Get_B_GP(X_EC_GP,Elem_Coords);	
    B_T = Transpose_Mat(B), free(B.nM);

    /* 10º Get forces in the nodes of the element created by the Gauss-Point 
     and free the B_T matrix */
    Element_INT_FORCES = Scalar_prod(B_T,StressTensor_GP);
    free(B_T.nV);

    /* 11º Acumulate this forces to the total array with the internal forces */
    for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){  
      for(int l = 0 ; l<NumberDimensions ; l++){
	/* 11aº Add the internal forces */
	Nodal_TOT_FORCES.nM[l][Elem_Nods[k]] -=
	  Element_INT_FORCES.nV[k*NumberDimensions + l];
	/* 11bº Add the body forces */
	Nodal_TOT_FORCES.nM[l][Elem_Nods[k]] +=
	  MPM_Mesh.Phi.mass.nV[i]*N_Ref_GP.nV[k]*g.nV[l];
	/* 11cº Add the contact forces */
	Nodal_TOT_FORCES.nM[l][Elem_Nods[k]] +=
	  N_Ref_GP.nV[k]*MPM_Mesh.Phi.F.nM[i][l]*
	  (MPM_Mesh.Phi.mass.nV[i]/MPM_Mesh.Phi.rho.nV[i]);
      }
    }
	
    /* 12º Free memory */
    free(Element_INT_FORCES.nV);
    free(N_Ref_GP.nV);

  }  

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
      Nodal_MOMENTUM.nM[j][i] +=
	DeltaTimeStep*Nodal_TOT_FORCES.nM[j][i];
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
  Matrix X_EC_GP; /* Local coordinates of the Gauss-Point */
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  int Elem_GP; /* Index of the element */
  int * Elem_Nods; /* Connectivity of the element */
  Matrix N_Ref_GP; /* Value of the shape-function in the GP */

  /* 1º iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Get the coordinates of the GP (in the element) */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 3º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* Get the connectivity of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 4º Evaluate the shape function in the coordinates of the GP */
    N_Ref_GP = FEM_Mesh.N_ref(X_EC_GP);

    /* 5º Iterate over the nodes of the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){

	/* 6aº Update the GP velocities */
	MPM_Mesh.Phi.vel.nM[i][k] +=
	  DeltaTimeStep*N_Ref_GP.nV[j]*
	  Nodal_TOT_FORCES.nM[k][Elem_Nods[j]]/
	  Nodal_MASS.nV[Elem_Nods[j]];
	
	/* 6bº Update the GP position */
	MPM_Mesh.Phi.x_GC.nM[i][k] +=
	  DeltaTimeStep*N_Ref_GP.nV[j]*
	  Nodal_MOMENTUM.nM[k][Elem_Nods[j]]/
	  Nodal_MASS.nV[Elem_Nods[j]];
      }      
    }

    /* 7º Free The value of the shape functions */
    free(N_Ref_GP.nV);
    
  }  
}


/*******************************************************/


  