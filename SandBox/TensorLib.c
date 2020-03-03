#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*************************************************************/

/* Matrix definition */
typedef struct{
  int N_rows; /* Number of rows */
  int N_cols; /* Number of columns */
  double n; /* Value if is an scalar */
  double * nV; /* Pointer for a vector */
  double ** nM; /* Table of pointers for a matrix */
  char Info [100]; /* Aditional information */
} Matrix;

/*************************************************************/

/* Tensor definition */
typedef struct{
  int Order; /* Order of the tensor */
  double *n; /* First order tensor */
  double *N[3]; /* Second order tensor */
  char Info [100]; /* Aditional information */
} Tensor;

/*************************************************************/

Tensor alloc_Tensor(int Order){
  /* Define output */
  Tensor A;
  /* Swith cases */
    switch(Order){
    case 0:
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in alloc_Tensor()",
	      "Not posible to allocate a scalar");
      exit(EXIT_FAILURE);       
    case 1:
      A.Order = 1;
      A.n = malloc(sizeof(double[3]));  
      break;
    case 2:
      A.Order = 2;
      A.N[0] = malloc(sizeof(double[3]));
      A.N[1] = malloc(sizeof(double[3]));
      A.N[2] = malloc(sizeof(double[3]));  
      break;
    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in alloc_Tensor()",
	      "Invalird order of the tensor");
      exit(EXIT_FAILURE); 
    }
  return A;
}

/*************************************************************/

Tensor memory_to_Tensor(double * A_mem, int Order){
  /* Define output */
  Tensor A_tens;
  /* Swith cases */
  switch(Order){
  case 0:
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in memory_to_Tensor()",
	    "Tensor type not include scalar properties");
    exit(EXIT_FAILURE);       
  case 1:
    A_tens.Order = 1;
    A_tens.n = A_mem;  
    break;
  case 2:
    A_tens.Order = 2;
    A_tens.N[0] = A_mem;
    A_tens.N[1] = A_mem+3;
    A_tens.N[2] = A_mem+6;
    break;
  default :
    fprintf(stderr,"%s : %s \n",
	    "Error in memory_to_Tensor()",
	    "Invalird order of the tensor");
    exit(EXIT_FAILURE); 
  }

  return A_tens;
}

/*************************************************************/

void free_Tensor(Tensor A){
  switch(A.Order){
    case 0:
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in free_Tensor()",
	      "Not posible to free a scalar");
      exit(EXIT_FAILURE);       
    case 1:
      free(A.n);
      break;
    case 2:
      free(A.N[0]);
      free(A.N[1]);
      free(A.N[2]);
      break;
    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in free_Tensor()",
	      "Invalird order of the tensor");
      exit(EXIT_FAILURE); 
    }
}

/*************************************************************/

double get_I1_Of(Tensor A)
{
  /* Define output */
  double I1;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    I1 =
      A.N[0][0] + A.N[1][1] + A.N[2][2];
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_I1_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I1;  
}

/*************************************************************/

double get_I2_Of(Tensor A)
{
  /* Define output */
  double I2;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    I2 =
      A.N[0][0]*A.N[1][1] +
      A.N[1][1]*A.N[2][2] +
      A.N[0][0]*A.N[2][2] -
      A.N[0][1]*A.N[1][0] -
      A.N[1][2]*A.N[2][3] -
      A.N[0][2]*A.N[2][0];
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_I2_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I2;  
}

/*************************************************************/

double get_I3_Of(Tensor A)
{
  /* Define output */
  double I3;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    I3 =
      A.N[0][0]*A.N[1][1]*A.N[2][2] - 
      A.N[0][0]*A.N[1][2]*A.N[2][1] + 
      A.N[0][1]*A.N[1][2]*A.N[2][0] - 
      A.N[0][1]*A.N[1][0]*A.N[2][2] + 
      A.N[0][2]*A.N[1][0]*A.N[2][1] - 
      A.N[0][2]*A.N[1][1]*A.N[2][0] ; 
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_I3_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I3;
}

/*************************************************************/

double get_J1_Of(Tensor A)
{
  /* Define output */
  double J1;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    J1 = get_I1_Of(A);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_J1_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return J1;
}

/*************************************************************/

double get_J2_Of(Tensor A)
{
  /* Define output */
  double J2;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    double I1 = get_I1_Of(A);
    double I2 = get_I2_Of(A);
    J2 = pow(I1,2) - 2*I2;
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_J2_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return J2;
}

/*************************************************************/

double get_J3_Of(Tensor A)
{
  /* Define output */
  double J3;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    double I1 = get_I1_Of(A);
    double I2 = get_I2_Of(A);
    double I3 = get_I3_Of(A);
    J3 = pow(I1,3) - 3*I1*I2 + 3*I3;
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_J3_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return J3;
}

/*************************************************************/

double get_EuclideanNorm_Of(Tensor A)
{
  /* Define output */
  double Out;
  /* Check if the input is a first order tensor */
  if(A.Order == 1){
    /* Compute norm */    
    double Aux = 0;
    for(int i = 0 ; i<3 ; i++){
      Aux += A.n[i]*A.n[i] ;
    }
    Out = pow(Aux,0.5);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_EuclideanNorm_Of()",
	    "The input should be of order 1");
    exit(EXIT_FAILURE);    
  }
  return Out;
}

/*************************************************************/

Tensor get_I(){
  Tensor I = alloc_Tensor(2);
  for(int i = 0 ; i<3 ; i++){
    I.N[i][i] = 1;
  }
  return I;
}

/*************************************************************/


Tensor get_Inverse_Of(Tensor A)
{
  /* Allocate the output */
  Tensor Am1 = alloc_Tensor(2);
  /* Check if the input is a second order tensor */
  if (A.Order == 2){  
    /* Get the determinant of the matrix */
    double detA = get_I3_Of(A);
    if(detA == 0){
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in get_Inverse_Of()",
	      "Determinant null");
      exit(EXIT_FAILURE);    
    }
    /* Compute each component  */
    Am1.N[0][0] =
      (double)1/detA*(A.N[1][1]*A.N[2][2] - A.N[1][2]*A.N[2][1]);
    Am1.N[0][1] =
      -(double)1/detA*(A.N[0][1]*A.N[2][2] - A.N[0][2]*A.N[2][1]);
    Am1.N[0][2] =
      (double)1/detA*(A.N[0][1]*A.N[1][2] - A.N[0][2]*A.N[1][1]);
    Am1.N[1][0] =
      -(double)1/detA*(A.N[1][0]*A.N[2][2] - A.N[1][2]*A.N[2][0]);
    Am1.N[1][1] =
      (double)1/detA*(A.N[0][0]*A.N[2][2] - A.N[0][2]*A.N[2][0]);
    Am1.N[1][2] =
      -(double)1/detA*(A.N[0][0]*A.N[1][2] - A.N[0][2]*A.N[1][0]);
    Am1.N[2][0] =
      (double)1/detA*(A.N[1][0]*A.N[2][1] - A.N[1][1]*A.N[2][0]);
    Am1.N[2][1] =
      -(double)1/detA*(A.N[0][0]*A.N[2][1] - A.N[0][1]*A.N[2][0]);
    Am1.N[2][2] =
      (double)1/detA*(A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0]);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_Inverse_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the inverse matrix */
  return Am1;  
}

/*************************************************************/

Tensor get_Transpose_Of(Tensor A)
{
  /* Allocate the output */
  Tensor AT = alloc_Tensor(2);  
  /* Check if the input is a second order tensor */
  if (A.Order == 2){
    /* Get the transpose */
    for(int i = 0 ; i < 3 ; i++){
      for(int j = 0 ; j < 3 ; j++){
	AT.N[i][j] = A.N[j][i];
      }
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_Transpose_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the transpose */
  return AT;
}

/*************************************************************/

double get_innerProduct_Of(Tensor A, Tensor B)
{

  /* Variable declaration output matrix */
  double AdotB = 0;

  /* Inner product of second order tensors */
  if ( (A.Order == 2) && (B.Order == 2) ) { 
    for(int i = 0 ; i < 3 ; i++){
      for(int j = 0 ; j < 3 ; j++){
	  AdotB += A.N[i][j]*B.N[i][j];
      }
    }
  }
  /* Inner product of first order tensors */
  else if( (A.Order == 1) && (B.Order == 1) ){
    for(int i = 0 ; i < 3 ; i++){
      AdotB += A.n[i]*B.n[i];
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_innerProduct_Of()",
	    "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);        
  }

  return AdotB;
    
}

/*************************************************************/

Tensor get_vectorProduct_Of(Tensor a, Tensor b){
  /* Allocate output */
  Tensor axb = alloc_Tensor(1);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)){
    /* Operate vetor product */
    axb.n[0] = a.n[1]*b.n[2] - a.n[2]*b.n[1];
    axb.n[1] = a.n[2]*b.n[0] - a.n[0]*b.n[2];
    axb.n[2] = a.n[0]*b.n[1] - a.n[1]*b.n[0];
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_vectorProduct_Of()",
	    "The input should be two tensors of first order");
    exit(EXIT_FAILURE);    
  }
  /* Return vector product */
  return axb; 
}

/*************************************************************/

Tensor get_dyadicProduct_Of(Tensor a, Tensor b)
{
  /* Tensor declaration */
  Tensor aob = alloc_Tensor(2);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)){
    /* Operate tensor product */
    for(int i = 0 ; i<3 ; i++){
      for(int j = 0 ; j<3 ; j++){
	aob.N[i][j] = a.n[i]*b.n[j];
      } 
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_dyadicProduct_Of()",
	    "The input should be two tensors of first order");
    exit(EXIT_FAILURE);    
  }
  /* Return tensorial product */
  return aob;
}

/*************************************************************/

Tensor get_firstOrderContraction_Of(Tensor A, Tensor b)
{
  /* Tensor declaration */
  Tensor Adotb = alloc_Tensor(1);
  /* Check in the input its is ok */
  if ((A.Order == 2) && (b.Order == 1)){
    /* Auxiliar variable */
    double Aux;
    /* Operate tensors */
    for(int i = 0 ; i<3 ; i++){
      Aux = 0;
      for(int j = 0 ; j<3 ; j++){
	Aux += A.N[i][j]*b.n[j];
      }
      Adotb.n[i] = Aux;
    }    
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_firstOrderContraction_Of()",
	    "The input should be 2ord tensor and a 1rd tensor");
    exit(EXIT_FAILURE);
  }
  /* Return tensor */
  return Adotb;
}

/*************************************************************/

Tensor compute_IncrementStrain(Matrix Velocity,
			       Matrix Gradient,
			       double DeltaTimeStep)
{
  Tensor Increment_Strain = alloc_Tensor(2);
  Tensor Velocity_I;
  Tensor Gradient_I;
  Tensor vog_I;

  int NodesElem = Gradient.N_rows;

  /* Compute rate of strain */
  for(int I = 0 ; I<NodesElem ; I++){
    /* Assign from matrix to tensor */
    Velocity_I = memory_to_Tensor(Velocity.nM[I], 1);
    Gradient_I = memory_to_Tensor(Gradient.nM[I], 1);
    /* Compute the dyadic product of the nodal velocity and the
       gradient of the shape functions */
    vog_I = get_dyadicProduct_Of(Velocity_I, Gradient_I);
    /* Ad the nodal contribution to the train tensor */
    for(int i = 0 ; i<3 ; i++){
      for(int j = 0 ; j<3 ; j++){
	Increment_Strain.N[i][j] +=
	  0.5*(vog_I.N[i][j] + vog_I.N[j][i]);
      }
    }
    /* Free memory */
    free_Tensor(vog_I);
  }

  /* Compute increment of strain */
  for(int i = 0 ; i<3 ; i++){
    for(int j = 0 ; j<3 ; j++){
      Increment_Strain.N[i][j] =
	Increment_Strain.N[i][j]*DeltaTimeStep;
    }
  }
  
  return Increment_Strain;
}

/* /\*************************************************************\/ */

/* Tensor compute_Stress(Tensor Strain, Tensor Stress, Material Mat) */
/* { */
/*   /\* Variable definition  *\/ */
/*   Tensor Strain_n1;  */
    
/*   /\* Select the constitutive model *\/ */
/*   if(strcmp(Mat.Type,"LE") == 0){ */
/*     Stress = LinearElastic(Strain,Stress,Mat); */
/*   } */
/*   else{ */
/*     exit(0); */
/*   } */
  
/*   /\* Return the stress tensor *\/ */
/*   return Stress; */
/* } */


/* /\*************************************************************\/ */

/* void compute_InternalForces(Matrix F_I, Matrix V_I, */
/* 			    GaussPoint MPM_Mesh, */
/* 			    Mesh FEM_Mesh){ */

/*   Element Nodes_p; /\* Element for each Gauss-Point *\/ */
/*   Matrix Gradient_p; /\* Shape functions gradients *\/ */
/*   Matrix Nodal_Velocity_p; /\* Velocity of the element nodes *\/ */
/*   Material Material_p; /\* Properties of the Gauss-Point material *\/ */
/*   Tensor Increment_Strain_p; /\* Increment of strain tensor *\/ */
/*   Tensor Strain_p; /\*  Strain tensor *\/ */
/*   Tensor Stress_p; /\* Stress tensor *\/ */
/*   Tensor Gradient_pI; */
/*   Tensor InternalForcesDensity_Ip; */
/*   double W_p; /\* Internal energy of the Gauss-Point *\/ */
/*   double m_p; /\* Mass of the Gauss-Point *\/ */
/*   double rho_p; /\* Density of the Gauss-Point *\/ */
/*   double V_p; /\* Volumen of the Gauss-Point *\/ */
/*   int Ip; */
/*   int Nn, Np; */

/*   /\* Loop in the GPs *\/ */
/*   for(int p = 0 ; p<Np ; p++){ */

/*     /\* Get the value of the density *\/ */
/*     rho_p = MPM_Mesh.Phi.rho.nV[i]; */

/*     /\* Get the value of the mass *\/ */
/*     m_p = MPM_Mesh.Phi.mass.nV[i]; */

/*     /\* Asign memory to tensors *\/ */
/*     Strain_p = memory_to_Tensor(MPM_Mesh.Phi.Strain.nM[I], 2); */
/*     Stress_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[I], 2); */

/*     /\* Define element for each GP *\/ */
/*     Nodes_p = */
/*       get_Element(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]); */

/*     /\* Get the velocity of the nodes of the element *\/ */
/*     Nodal_Velocity_p = get_Element_velocity(Nodes_p, V_I); */

/*     /\* Compute gradient of the shape function in each node *\/ */
/*     Gradient_p = */
/*       compute_ShapeFunction_Gradient(Nodes_p, MPM_Mesh, FEM_Mesh); */

/*     /\* Get the material properties *\/ */
/*     Idx_Mat_p = MPM_Mesh.MatIdx[p]; */
/*     Material_p = MPM_Mesh.Mat[Idx_Mat_p];  */

/*     /\* Compute Strain tensor *\/ */
/*     Increment_Strain_p = */
/*       compute_IncrementStrain(Nodal_Velocity_p,dNdx_p,DeltaTimeStep); */
/*     for(int i = 0 ; i<3 ; i++){ */
/*       for(int j = 0 ; j<3 ; j++){ */
/* 	Strain_p.N[i][j] += Increment_Strain_p.N[i][j]; */
/*       } */
/*     } */

/*     /\* Update density field *\/ */
/*     rho_p = rho_p/(1 + get_I1_Of(Increment_Strain_p)); */

/*     /\* Compute stress tensor *\/ */
/*     Stress_p = compute_Stress(Strain_p,Stress_p,Material_p); */

/*     /\* Compute deformation energy *\/ */
/*     W_p = 0.5*get_innerProduct_Of(Strain_p, Stress_p); */

/*     /\* Compute the volume of the Gauss-Point *\/ */
/*     V_p = m_p/rho_p; */

/*     /\* Compute nodal forces *\/ */
/*     for(int I = 0 ; I<Nn ; I++){ */
/*       /\* Pass by reference the nodal gradient to the tensor *\/ */
/*       Gradient_pI = memory_to_Tensor(Gradient_p.nM[I], 1); */
/*       /\* Compute the nodal forces of the Gauss-Point *\/ */
/*       InternalForcesDensity_Ip = */
/* 	get_firstOrderContraction_Of(Stress_p, Gradient_pI); */
/*       /\* Get the node of the mesh for the contribution *\/       */
/*       Ip = Nodes_p.Connectivity[I]; */
/*       /\* Asign the nodal forces contribution to the node *\/ */
/*       for(int i = 0 ; i<3 ; i++){ */
/* 	F_I.nM[Ip][i] += InternalForcesDensity_Ip.n[i]*V_p; */
/*       } */
/*       /\* Free the internal forces density *\/ */
/*       free_Tensor(InternalForcesDensity_Ip); */
/*     } */

/*     /\* Update memory *\/ */
/*     MPM_Mesh.Phi.rho.nV[i] = rho_p; */
    
    
/*     /\* Free the matrix with the nodal velocity of the element *\/ */
/*     FreeMat(Nodal_Velocity_p); */
    
/*     /\* Free the matrix with the nodal gradient of the element *\/ */
/*     FreeMat(Gradient_p); */
    
/*   }   */
  
/* } */

/*************************************************************/

void main(){
  int X = 10;
  int Y = 10;
  Tensor A = alloc_Tensor(2);
  Tensor v = alloc_Tensor(1);
  Tensor I = get_I();
  Tensor Im1 = get_Inverse_Of(I);
  double * B_mem = malloc(9*sizeof(double));
  B_mem[0] = 1;
  B_mem[4] = 1;
  B_mem[8] = 1;
  Tensor B_tens = memory_to_Tensor(B_mem, 2);
  printf("%f %f \n",get_I3_Of(A),get_I3_Of(I));
  printf("Im1\n");
  for(int i = 0 ; i<3 ; i++){
    for(int j = 0 ; j<3 ; j++){
      printf("%f ",Im1.N[i][j]);
    }
    printf("\n");
  }
  printf("B\n");
  for(int i = 0 ; i<3 ; i++){
    for(int j = 0 ; j<3 ; j++){
      printf("%f ",B_tens.N[i][j]);
    }
    printf("\n");
  }
  free(B_mem);
  free_Tensor(A);
  free_Tensor(v);
  free_Tensor(I);
  free_Tensor(Im1);  
}
