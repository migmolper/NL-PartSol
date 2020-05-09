#include "grams.h"

/*************************************************************/

Tensor alloc_Tensor(int Order)
{
  int Ndim = NumberDimensions;
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
    A.n = (double *)malloc(sizeof(double[Ndim]));
    for(int i = 0 ; i<Ndim ; i++){
      A.n[i] = 0.0;
    }
    break;
  case 2:
    A.Order = 2;
    for(int i = 0 ; i<Ndim ; i++){
      A.N[i] = (double *)malloc(sizeof(double[Ndim]));
      for(int j = 0 ; j<Ndim ; j++){
	A.N[i][j] = 0.0;
      }
    }
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

Tensor memory_to_Tensor(double * A_mem, int Order)
{
  long Ndim = NumberDimensions;
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
    for(long i = 0 ; i<Ndim ; i++){
      A_tens.N[i] = A_mem+i*Ndim;
    }
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

void free_Tensor(Tensor A)
{
  int Ndim = NumberDimensions;
  
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
      for(int i = 0 ; i<Ndim ; i++){
	free(A.N[i]);
      }
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
  int Ndim = NumberDimensions;
  /* Define output */
  double I1 = 0;
  /* Check if is the order is order 2 */
  if(A.Order == 2){
    for(int i = 0 ; i<Ndim ; i++){
      I1 += A.N[i][i];
    }
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
  int Ndim = NumberDimensions;  
  /* Define output */
  double I2 = 0;
  /* Check if is the order is order 2 */
  if(A.Order == 2){
    if(Ndim == 3){
      I2 =
	A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0] +
	A.N[1][1]*A.N[2][2] - A.N[1][2]*A.N[2][1] +
	A.N[0][0]*A.N[2][2] - A.N[2][0]*A.N[0][2];
    }
    if(Ndim == 2){
      I2 =
	A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0];
    }
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
  int Ndim = NumberDimensions;  
  /* Define output */
  double I3 = 0;
  /* Check if is the order is order 2 */
  if(A.Order == 2){
    if(Ndim == 3){
      I3 =
	A.N[0][0]*A.N[1][1]*A.N[2][2] - 
	A.N[0][0]*A.N[1][2]*A.N[2][1] + 
	A.N[0][1]*A.N[1][2]*A.N[2][0] - 
	A.N[0][1]*A.N[1][0]*A.N[2][2] + 
	A.N[0][2]*A.N[1][0]*A.N[2][1] - 
	A.N[0][2]*A.N[1][1]*A.N[2][0] ; 
    }
    if(Ndim == 2){
      I3 =
	A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0];
    }
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

Tensor get_Eigenvalues_Of(Tensor A)
{
  /* Auxiliar variables */
  int Ndim = NumberDimensions;
  double I1, I2, I3;
  /* Define output */
  Tensor w = alloc_Tensor(1);
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    if(Ndim == 2){      
      I1 = get_I1_Of(A);
      I2 = get_I2_Of(A);
      w.n[0] = 0.5*(I1+sqrt(DMAX(0,I1*I1-4*I2)));
      w.n[1] = 0.5*(I1-sqrt(DMAX(0,I1*I1-4*I2)));
    }
    if(Ndim == 3){      
      I1 = get_I1_Of(A);
      I2 = get_I2_Of(A);
      I3 = get_I3_Of(A);
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in get_Eigenvalues_Of()",
	      "3D cases are not implemented");
      exit(EXIT_FAILURE);          
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in get_Eigenvalues_Of()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }

  return w;
}

/*************************************************************/

double get_EuclideanNorm_Of(Tensor A)
{
  int Ndim = NumberDimensions;
  /* Define output */
  double Out;
  /* Check if the input is a first order tensor */
  if(A.Order == 1){
    /* Compute norm */    
    double Aux = 0;
    for(int i = 0 ; i<Ndim ; i++){
      Aux += DSQR(A.n[i]);
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

Tensor get_I()
{
  int Ndim = NumberDimensions;
  Tensor I = alloc_Tensor(2);
  for(int i = 0 ; i<Ndim ; i++){
    I.N[i][i] = 1;
  }
  return I;
}

/*************************************************************/

Tensor get_Inverse_Of(Tensor A)
{
  int Ndim = NumberDimensions;
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
    if(Ndim == 3){
      /* Compute each component  */
      Am1.N[0][0] = + (double)1/detA*(A.N[1][1]*A.N[2][2] - A.N[1][2]*A.N[2][1]);
      Am1.N[0][1] = - (double)1/detA*(A.N[0][1]*A.N[2][2] - A.N[0][2]*A.N[2][1]);
      Am1.N[0][2] = + (double)1/detA*(A.N[0][1]*A.N[1][2] - A.N[0][2]*A.N[1][1]);
      Am1.N[1][0] = - (double)1/detA*(A.N[1][0]*A.N[2][2] - A.N[1][2]*A.N[2][0]);
      Am1.N[1][1] = + (double)1/detA*(A.N[0][0]*A.N[2][2] - A.N[0][2]*A.N[2][0]);
      Am1.N[1][2] = - (double)1/detA*(A.N[0][0]*A.N[1][2] - A.N[0][2]*A.N[1][0]);
      Am1.N[2][0] = + (double)1/detA*(A.N[1][0]*A.N[2][1] - A.N[1][1]*A.N[2][0]);
      Am1.N[2][1] = - (double)1/detA*(A.N[0][0]*A.N[2][1] - A.N[0][1]*A.N[2][0]);
      Am1.N[2][2] = + (double)1/detA*(A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0]);
    }
    if(Ndim == 2){
      Am1.N[0][0] = + (double)1/detA*A.N[1][1];
      Am1.N[0][1] = - (double)1/detA*A.N[1][0];
      Am1.N[1][0] = - (double)1/detA*A.N[0][1];
      Am1.N[1][1] = + (double)1/detA*A.N[0][0]; 
    }    
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
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor AT = alloc_Tensor(2);  
  /* Check if the input is a second order tensor */
  if (A.Order == 2){
    /* Get the transpose */
    for(int i = 0 ; i < Ndim ; i++){
      for(int j = 0 ; j < Ndim ; j++){
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

  int Ndim = NumberDimensions;
  /* Variable declaration output matrix */
  double AdotB = 0;
  /* Inner product of second order tensors */
  if ( (A.Order == 2) && (B.Order == 2) ) { 
    for(int i = 0 ; i < Ndim ; i++){
      for(int j = 0 ; j < Ndim ; j++){
	  AdotB += A.N[i][j]*B.N[i][j];
      }
    }
  }
  /* Inner product of first order tensors */
  else if( (A.Order == 1) && (B.Order == 1) ){
    for(int i = 0 ; i < Ndim ; i++){
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

Tensor get_vectorProduct_Of(Tensor a, Tensor b)
{
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
  int Ndim = NumberDimensions;
  /* Tensor declaration */
  Tensor aob = alloc_Tensor(2);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)){
    /* Operate tensor product */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
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
  int Ndim = NumberDimensions;
  /* Tensor declaration */
  Tensor Adotb = alloc_Tensor(1);
  /* Check in the input its is ok */
  if ((A.Order == 2) && (b.Order == 1)){
    /* Auxiliar variable */
    double Aux;
    /* Operate tensors */
    for(int i = 0 ; i<Ndim ; i++){
      Aux = 0;
      for(int j = 0 ; j<Ndim ; j++){
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



