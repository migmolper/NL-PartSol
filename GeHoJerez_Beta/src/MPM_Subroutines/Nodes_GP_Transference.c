#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../InOutFun/InOutFun.h"

Matrix GetNodalValuesFromGP(GaussPoint MPM_Mesh, Mesh FEM_Mesh, char LisOfFields[MAXW])
{

  /* 0º Variable declaration */
  Matrix Nod_Fields; /* Output list of fields */
  Matrix GP_Fields; /* Auxiliar structure to do the selection */
  int Elem_GP; /* Index of the element */
  int * Elem_Nods;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  Matrix N_Ref_GP;
  Matrix B, B_T;
  char * FieldsList[MAXW] = {NULL};
  int NumberFields;
  int Size_Nod_Fields = 0; /* Number of fields ouput */
  int i_Field;
  Matrix Elem_Coords = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  strcpy(Elem_Coords.Info,FEM_Mesh.TypeElem);
  Matrix StressTensor_GP; /* Stress tensor of a Gauss Point */
  Matrix F_INT_ELEM; /* Internal forces for each node in a element (by a GP)*/

  /* 1º Get those fields to transfeer */
  strcpy(Nod_Fields.Info,LisOfFields); /* Transfeer this information first */
  NumberFields = parse (FieldsList,LisOfFields,";\n");

  for(int i = 0 ; i<NumberFields ; i++){
    if(strcmp(FieldsList[i],"MASS") == 0){
      Size_Nod_Fields += 1;
    }
    if(strcmp(FieldsList[i],"MOMENTUM") == 0){
      Size_Nod_Fields += NumberDimensions;
    }
    if(strcmp(FieldsList[i],"F_INT") == 0){
      Size_Nod_Fields += NumberDimensions;
      StressTensor_GP.N_rows = MPM_Mesh.Phi.Stress.N_cols;
      StressTensor_GP.N_cols = 1;
      StressTensor_GP.nM = NULL;
    }
  }
  
  if(NumberFields<=0){
    puts("Error in GaussPointToMesh() : Wrong number of fields !!!");
    exit(0);
  }
  
  /* 2º Allocate the output list of fields */
  Nod_Fields.N_rows = Size_Nod_Fields;
  Nod_Fields.N_cols = FEM_Mesh.NumNodesMesh;
  Nod_Fields.n = NAN;
  Nod_Fields.nV = NULL;
  Nod_Fields.nM = (double **)Allocate_MatrixZ(Size_Nod_Fields,
					      FEM_Mesh.NumNodesMesh,
					      sizeof(double));
  
  /* 3º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 4º Get the coordinates of the GP in the Element */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 5º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* 6º Nodes of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 7º Coordinates of the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<FEM_Mesh.Dimension ; k++){
	Elem_Coords.nM[j][k] = FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 8º Evaluate the shape function and it derivarive in the GP */
    N_Ref_GP = FEM_Mesh.N_ref(X_EC_GP);
   
    /* 10º Set to zero the auxiliar index for the input fields */
    i_Field = 0;
    
    /* 11º Loop for each nodal field that we want to calc : */
    for(int j = 0 ; j<NumberFields ; j++){

      /* 12aº Asign GP mass to the nodes */
      if(strcmp(FieldsList[j],"MASS") == 0){
      	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
      	  Nod_Fields.nM[i_Field][Elem_Nods[k]] += MPM_Mesh.Phi.mass.nV[i]*N_Ref_GP.nV[k];
      	}
	/* Update the index of the field */
	i_Field += 1;
      }

      /* 12bº Asign GP momentum to the nodes */
      if( strcmp(FieldsList[j],"MOMENTUM") == 0 ){
      	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
      	  for(int l = 0 ; l<NumberDimensions ; l++){
      	    Nod_Fields.nM[i_Field+l][Elem_Nods[k]] +=
      	      MPM_Mesh.Phi.mass.nV[i]*
      	      MPM_Mesh.Phi.vel.nM[i][l]*
      	      N_Ref_GP.nV[j];
      	  }
      	}
	/* Update the index of the field */
	i_Field += NumberDimensions; /* Number of components of the Velocity field */
      }

      /* 12cº Asign GP internal forces to the nodes */
      if(strcmp(FieldsList[j],"F_INT") == 0){

	/* Asign to an auxiliar variable the value of the stress tensor */
	StressTensor_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

	/* Multiply the stress tensor by the mass and the inverse of the density */
	for(int k = 0 ; k<StressTensor_GP.N_rows ; k++){
	  StressTensor_GP.nV[k] *= (MPM_Mesh.Phi.mass.nV[i]/MPM_Mesh.Phi.rho.nV[i]);
	}
	
	/* Get the B_T matrix to for the derivates */
	B = Get_B_GP(X_EC_GP,Elem_Coords);
	
	B_T = Transpose_Mat(B), free(B.nM);

	/* Get forces in the nodes of the element created by the Gauss-Point */
	F_INT_ELEM = Scalar_prod(B_T,StressTensor_GP);

	/* Acumulate this forces to the total array with the internal forces */
	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){  
	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Nod_Fields.nM[i_Field+l][Elem_Nods[k]] -=  F_INT_ELEM.nV[k*NumberDimensions + l];
	  }
	}

	/* Free memory */
	free(F_INT_ELEM.nV);

	/* Update the index of the field */
	if(NumberDimensions == 1)
	  i_Field += 1; /* Number of components of the 2D stress tensor */
	if(NumberDimensions == 2)
	  i_Field += 3; /* Number of components of the 2D stress tensor */
	if(NumberDimensions == 3)
	  i_Field += 6; /* Number of components of the 3D stress tensor */


      }

      /* 12dº Asign GP body forces to the nodes */
      if(strcmp(FieldsList[j],"F_BODY") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }

      /* 12eº Asign GP local forces to the nodes */
      if(strcmp(FieldsList[j],"F_LOC") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }

      /* 12fº Asign GP total forces to the nodes */
      if(strcmp(FieldsList[j],"F_TOT") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }
    
    }
    
  }

  /* 13º Free The value of the shape functions */
  free(N_Ref_GP.nV);

  
  return Nod_Fields;
  
}

/*******************************************************/

Matrix GetNodalVelocity(Matrix Nodal_MOMENTUM,
			Matrix Nodal_MASS)
/*
  Get the nodal velocity using : 
  v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
*/
{

  Matrix Nodal_VELOCITY;

  Nodal_VELOCITY = MatAllocZ(Nodal_MOMENTUM.N_rows,Nodal_MOMENTUM.N_cols);

  for(int i = 0 ; i<Nodal_MOMENTUM.N_cols ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      Nodal_VELOCITY.nM[j][i] = Nodal_MOMENTUM.nM[j][i] / Nodal_MASS.nV[i];
    }
  }

  return Nodal_VELOCITY;
}

/*******************************************************/

void GetGaussPointStrainIncrement(GaussPoint MPM_Mesh,
				  Mesh FEM_Mesh,
				  Matrix Nodal_VELOCITY)
/*
  Calcule the particle stress increment :

  \Delta\Epsilon_{ij,p}^{k-1/2} = \frac{\Delta t}{2} \cdot \sum_{I=0}^{Nn}(N^{k}_{Ip,j} \cdot v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot v_{jI}^{k-1/2})

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
  Matrix Strain_GP; /* Vectoriced strain tensor */

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
    Strain_GP = Scalar_prod(B,Nodal_VELOCITY_Elem);

    /* 9º Set the Gauss-Point increment strain tensor */
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      MPM_Mesh.Phi.Strain.nM[i][j] = Strain_GP.nV[j];
    }

    /* 10º Free memory for the next step */
    free(B.nM);
    free(Strain_GP.nV);
    
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
    MPM_Mesh.Phi.rho.nV[i] /= 1 + TraceStrainTensor;
  }
  
}


/*******************************************************/

void UpdateGaussPointStressTensor(GaussPoint MPM_Mesh, Matrix D){

  /* 0º Variable declaration  */
  Matrix DeltaStrainTensor_GP;
  Matrix DeltaStressTensor_GP;

  /* 1º Switch the dimensions of the aulixiar strain tensor */
  switch(NumberDimensions){
  case 1:
    DeltaStrainTensor_GP.N_rows = 1;
    DeltaStrainTensor_GP.N_cols = 1;
    DeltaStrainTensor_GP.nM = NULL;
    break;
  case 2:
    DeltaStrainTensor_GP.N_rows = 3;
    DeltaStrainTensor_GP.N_cols = 1;
    DeltaStrainTensor_GP.nM = NULL;
    break;
  case 3:
    DeltaStrainTensor_GP.N_rows = 6;
    DeltaStrainTensor_GP.N_cols = 1;
    DeltaStrainTensor_GP.nM = NULL;
    break;
  default :
    puts("Error in UpdateGaussPointStressTensor() : Wrong number of dimensions !!! ");
    exit(0);
  }

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    /* 3º Store in an auxiliar variable the increment of the strain tensor in the GP */
    DeltaStrainTensor_GP.nV = MPM_Mesh.Phi.Strain.nM[i];
    /* 4º Get the increment of the stress tensor */
    DeltaStressTensor_GP = Scalar_prod(D,DeltaStrainTensor_GP);
    /* 5º Update the stress tensor with the increment */
    for(int j = 0 ; j<DeltaStrainTensor_GP.N_rows ; j++){
      MPM_Mesh.Phi.Stress.nM[i][j] += DeltaStressTensor_GP.nV[j];
    }
    /* 6º Free memory */
    free(DeltaStressTensor_GP.nV);
  }
  
}
