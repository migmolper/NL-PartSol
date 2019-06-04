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
  
  /* 0aº Nodal variable declaration */
  Matrix Nodal_FIELDS; /* Output list of fields */
  Matrix Nodal_MASS; /* Nodal values of the mass */
  Matrix Nodal_MOMENTUM; /* Nodal values of the momentum */
  Matrix Nodal_INT_FORCES; /* Nodal values of the internal forces */
  Matrix Nodal_GRAVITY_FORCES; /* Nodal values of the gravity forces */

  /* 0bº Auxiliar variable declaration */
  int Elem_GP; /* Index of the element */
  int * Elem_Nods;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  Matrix N_Ref_GP;
  Matrix B, B_T;
  char * FieldsList[MAXW] = {NULL};
  int NumberFields;
  int Size_Nodal_FIELDS = 0; /* Number of fields ouput */
  int i_Field;
  Matrix Elem_Coords = MatAllocZ(FEM_Mesh.NumNodesElem,NumberDimensions);
  strcpy(Elem_Coords.Info,FEM_Mesh.TypeElem);
  Matrix StressTensor_GP; /* Stress tensor of a Gauss Point */
  Matrix F_INT_ELEM; /* Internal forces for each node in a element (by a GP)*/

  /* 1º Get those fields to transfeer */
  strcpy(Nodal_FIELDS.Info,LisOfFields); /* Transfeer this information first */
  NumberFields = parse (FieldsList,LisOfFields,";\n");
  if(NumberFields<=0){ /* Check */
    puts("Error in GaussPointToMesh() : Wrong number of fields !!!");
    exit(0);
  }

  for(int i = 0 ; i<NumberFields ; i++){
    /* Nodal mass */
    if(strcmp(FieldsList[i],"MASS") == 0){
      Size_Nodal_FIELDS += 1;
      Nodal_MASS = MatAllocZ(1,FEM_Mesh.NumNodesMesh);
    }
    /* Nodal momentum */
    if(strcmp(FieldsList[i],"MOMENTUM") == 0){
      Size_Nodal_FIELDS += NumberDimensions;
      Nodal_MOMENTUM = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
    }
    /* Nodal internal forces */
    if(strcmp(FieldsList[i],"F_INT") == 0){
      Size_Nodal_FIELDS += NumberDimensions;
      Nodal_INT_FORCES = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
      StressTensor_GP.N_rows = MPM_Mesh.Phi.Stress.N_cols;
      StressTensor_GP.N_cols = 1;
      StressTensor_GP.nM = NULL;
    }
    /* Nodal value of the gravity forces */
    if(strcmp(FieldsList[i],"F_GRAV") == 0){
      Size_Nodal_FIELDS += NumberDimensions;
      Nodal_GRAVITY_FORCES = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
    }
  }
  
  /* 2º Allocate the output list of fields */
  Nodal_FIELDS.N_rows = Size_Nodal_FIELDS;
  Nodal_FIELDS.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_FIELDS.n = NAN;
  Nodal_FIELDS.nV = NULL;
  Nodal_FIELDS.nM = (double **)malloc((unsigned)Size_Nodal_FIELDS*sizeof(double*));

  /* 3º Set to zero the auxiliar index for the input fields */
  i_Field = 0;

  /* 4º Asign memory */
  for(int i = 0 ; i<NumberFields ; i++){
    /* Nodal mass */
    if(strcmp(FieldsList[i],"MASS") == 0){
      Nodal_FIELDS.nM[i_Field] = Nodal_MASS.nV;
      i_Field++; /* Update the field index */
    }
    /* Nodal momentum */
    if(strcmp(FieldsList[i],"MOMENTUM") == 0){
      for(int j = 0 ; j<NumberDimensions ; j++){
	Nodal_FIELDS.nM[i_Field] = Nodal_MOMENTUM.nM[j];
	i_Field++; /* Update the field index */
      }
    }
    /* Nodal internal forces */
    if(strcmp(FieldsList[i],"F_INT") == 0){
      for(int j = 0 ; j<NumberDimensions ; j++){
	Nodal_FIELDS.nM[i_Field] = Nodal_INT_FORCES.nM[j];
	i_Field++; /* Update the field index */
      }
    }
    /* Nodal value of the gravity forces */
    if(strcmp(FieldsList[i],"F_GRAV") == 0){
      for(int j = 0 ; j<NumberDimensions ; j++){
	Nodal_FIELDS.nM[i_Field] = Nodal_GRAVITY_FORCES.nM[j];
	i_Field++; /* Update the field index */
      }
    }    
  }
  
  /* 5º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 6º Get the coordinates of the GP in the Element */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 7º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* 8º Nodes of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 9º Coordinates of the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<FEM_Mesh.Dimension ; k++){
	Elem_Coords.nM[j][k] = FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 10º Evaluate the shape function and it derivarive in the GP */
    N_Ref_GP = FEM_Mesh.N_ref(X_EC_GP);
    
    /* 11º Loop for each nodal field that we want to calc : */
    for(int j = 0 ; j<NumberFields ; j++){

      /******* 12aº Asign GP mass to the nodes ****/
      if(strcmp(FieldsList[j],"MASS") == 0){
      	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	  Nodal_MASS.nV[Elem_Nods[k]] +=
	    MPM_Mesh.Phi.mass.nV[i]*
	    N_Ref_GP.nV[k];
      	}
      }
      /********************************************/

      /**** 12bº Asign GP momentum to the nodes ***/
      if( strcmp(FieldsList[j],"MOMENTUM") == 0 ){
      	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
      	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Nodal_MOMENTUM.nM[l][Elem_Nods[k]] +=
      	      MPM_Mesh.Phi.mass.nV[i]*
      	      MPM_Mesh.Phi.vel.nM[i][l]*
      	      N_Ref_GP.nV[j];
      	  }
	}
      }
      /********************************************/

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
	    Nodal_INT_FORCES.nM[l][Elem_Nods[k]] -=
	      F_INT_ELEM.nV[k*NumberDimensions + l];
	  }
	}

	
	/* Free memory */
	free(F_INT_ELEM.nV);
      }
      /********************************************/

      /** 12dº Asign GP body forces to the nodes **/
      if(strcmp(FieldsList[j],"F_GRAV") == 0){
	/* Set the gravity force */
	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){  
	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Nodal_GRAVITY_FORCES.nM[l][Elem_Nods[k]] +=
	      MPM_Mesh.Phi.mass.nV[i]*
	      N_Ref_GP.nV[j]*g.nV[l];
	  }
	}	
      }
      /********************************************/

      /** 12eº Asign GP local forces to the nodes */
      if(strcmp(FieldsList[j],"F_LOC") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }
      /********************************************/
   
    }
    
  }

  /* 13º Free The value of the shape functions */
  free(N_Ref_GP.nV);

  
  return Nodal_FIELDS;
  
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
      if(Nodal_MASS.nV[i] != 0){
	Nodal_VELOCITY.nM[j][i] = Nodal_MOMENTUM.nM[j][i] / Nodal_MASS.nV[i];
      }
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

/*******************************************************/

void UpdateGridNodalMomentum(Mesh FEM_Mesh,
			     Matrix Nodal_MOMENTUM,
			     Matrix Nodal_TOT_FORCES)
{
  /* Auxiliar variables */
  Matrix Nodal_INT_FORCES;
  Matrix Nodal_EXT_FORCES;
  
  /* Select first wich forces are over the Gauss-Point */
  Nodal_INT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  Nodal_EXT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  Nodal_INT_FORCES.nM[0] = Nodal_TOT_FORCES.nM[0];
  Nodal_INT_FORCES.nM[1] = Nodal_TOT_FORCES.nM[1];
  Nodal_EXT_FORCES.nM[0] = Nodal_TOT_FORCES.nM[2];
  Nodal_EXT_FORCES.nM[1] = Nodal_TOT_FORCES.nM[3];
  
  /* Update the grid nodal momentum */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      Nodal_MOMENTUM.nM[j][i] += DeltaTimeStep*(Nodal_INT_FORCES.nM[j][i] +
						Nodal_EXT_FORCES.nM[j][i]);
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
  Matrix Nodal_INT_FORCES;
  Nodal_INT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  Nodal_INT_FORCES.nM[0] = Nodal_TOT_FORCES.nM[0];
  Nodal_INT_FORCES.nM[1] = Nodal_TOT_FORCES.nM[1];
  Matrix Nodal_EXT_FORCES;
  Nodal_EXT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  Nodal_EXT_FORCES.nM[0] = Nodal_TOT_FORCES.nM[2];
  Nodal_EXT_FORCES.nM[1] = Nodal_TOT_FORCES.nM[3];
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
	MPM_Mesh.Phi.vel.nM[i][k] += DeltaTimeStep*N_Ref_GP.nV[j]*
	  (Nodal_INT_FORCES.nM[k][Elem_Nods[j]] + Nodal_EXT_FORCES.nM[k][Elem_Nods[j]])/
	  Nodal_MASS.nV[Elem_Nods[j]];
	
	/* 6bº Update the GP position */
	MPM_Mesh.Phi.x_GC.nM[i][k] += DeltaTimeStep*N_Ref_GP.nV[j]*
	  Nodal_MOMENTUM.nM[k][Elem_Nods[j]]/
	  Nodal_MASS.nV[Elem_Nods[j]];
      }      
    }

    /* 7º Free The value of the shape functions */
    free(N_Ref_GP.nV);
    
  }  
}
