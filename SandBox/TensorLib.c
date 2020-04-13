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

int InOut_Poligon(Matrix X_Point, Matrix Poligon)
/*! 
 * Check if a point is or not (1/0) inside of a Poligon.
 * Inputs :
 * - \a X_Point : Coordinates of the point 
 * - \a Poligon : Coordinates of the vertex 0,1,....,n,0
 */
{
  /* By default, we suppose that the point is in the poligon */
  int InOut = 1;

  Matrix a = MatAllocZ(3,1);
  Matrix b = MatAllocZ(3,1);
  Matrix c;
  Matrix n;
  Matrix nxc;

  /* Get the normal vector */
  a.nV[0] = Poligon.nM[1][0] - Poligon.nM[0][0];
  a.nV[1] = Poligon.nM[1][1] - Poligon.nM[0][1];
  a.nV[2] = Poligon.nM[1][2] - Poligon.nM[0][2];  
  b.nV[0] = Poligon.nM[Poligon.N_rows-1][0] - Poligon.nM[0][0];
  b.nV[1] = Poligon.nM[Poligon.N_rows-1][1] - Poligon.nM[0][1];
  b.nV[2] = Poligon.nM[Poligon.N_rows-1][2] - Poligon.nM[0][2];
  n = Vectorial_prod(a,b);
  n.N_rows = 1;
  n.N_cols = 3;

  /* Fill a and b for the First search */
  a.nV[0] = Poligon.nM[0][0] - Poligon.nM[Poligon.N_rows-1][0];
  a.nV[1] = Poligon.nM[0][1] - Poligon.nM[Poligon.N_rows-1][1];
  a.nV[2] = Poligon.nM[0][2] - Poligon.nM[Poligon.N_rows-1][2];

  b.nV[0] = X_Point.nV[0] - Poligon.nM[Poligon.N_rows-1][0];
  b.nV[1] = X_Point.nV[1] - Poligon.nM[Poligon.N_rows-1][1];
  b.nV[2] = X_Point.nV[2] - Poligon.nM[Poligon.N_rows-1][2];
  
  for(int i = 0 ; i<Poligon.N_rows-1 ; i++){

    c = Vectorial_prod(a,b);
    nxc = Scalar_prod(n,c);
    FreeMat(c);

    if(nxc.n < 0){
      InOut = 0;
      break;
    }
    
    a.nV[0] = Poligon.nM[i+1][0] - Poligon.nM[i][0];
    a.nV[1] = Poligon.nM[i+1][1] - Poligon.nM[i][1];
    a.nV[2] = Poligon.nM[i+1][2] - Poligon.nM[i][2];

    b.nV[0] = X_Point.nV[0] - Poligon.nM[i][0];
    b.nV[1] = X_Point.nV[1] - Poligon.nM[i][1];
    b.nV[2] = X_Point.nV[2] - Poligon.nM[i][2];
    
  }

  /* Last check */
  c = Vectorial_prod(a,b);
  nxc = Scalar_prod(n,c);
  if(nxc.n < 0){
    InOut = 0;
  }

  FreeMat(a);
  FreeMat(b);
  FreeMat(c);
  FreeMat(n);

  return InOut;
}

/*************************************************************/

void main(){
  int X = 10;
  int Y = 10;
  Tensor A = alloc_Tensor(2);
  Tensor v = alloc_Tensor(1);
  v.n[0] = 1;
  v.n[1] = 1;
  v.n[2] = 1;
  Tensor w = alloc_Tensor(1);
  w.n[0] = 1;
  w.n[1] = 1;
  w.n[2] = 1;
  Tensor I = get_I();
  Tensor Im1 = get_Inverse_Of(I);
  double * B_mem = malloc(9*sizeof(double));
  B_mem[0] = 1;
  B_mem[4] = 1;
  B_mem[8] = 1;
  Tensor B_tens = memory_to_Tensor(B_mem, 2);
  Tensor vow = get_dyadicProduct_Of(v,w);
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
  printf("vow \n");
  for(int i = 0 ; i<3 ; i++){
    for(int j = 0 ; j<3 ; j++){
      printf("%f ",vow.N[i][j]);
    }
    printf("\n");
  }
  free(B_mem);
  free_Tensor(A);
  free_Tensor(v);
  free_Tensor(w);
  free_Tensor(vow);
  free_Tensor(I);
  free_Tensor(Im1);  
}
