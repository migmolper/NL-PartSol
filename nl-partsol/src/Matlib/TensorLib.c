#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

/*************************************************************/

Tensor alloc__TensorLib__(int Order)
{
  int Ndim = NumberDimensions;
  /* Define output */
  Tensor A;
  /* Swith cases */
  switch(Order)
  {
  case 0:
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in alloc__TensorLib__()",
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
	      "Error in alloc__TensorLib__()",
	      "Invalird order of the tensor");
      exit(EXIT_FAILURE); 
    }
  return A;
}

/*************************************************************/

Tensor memory_to_tensor__TensorLib__(double * A_mem, int Order)
{
  long Ndim = NumberDimensions;
  /* Define output */
  Tensor A_tens;
  /* Swith cases */
  switch(Order)
  {
  case 0:
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in memory_to_tensor__TensorLib__()",
	    "Tensor type not include scalar properties");
    exit(EXIT_FAILURE);       
  case 1:
    A_tens.Order = 1;
    A_tens.n = A_mem;  
    break;
  case 2:
    A_tens.Order = 2;
    for(long i = 0 ; i<Ndim ; i++)
    {
      A_tens.N[i] = A_mem+i*Ndim;
    }
    break;
  default :
    fprintf(stderr,"%s : %s \n",
	    "Error in memory_to_tensor__TensorLib__()",
	    "Invalird order of the tensor");
    exit(EXIT_FAILURE); 
  }

  return A_tens;
}

/*************************************************************/

void free__TensorLib__(Tensor A)
{
  int Ndim = NumberDimensions;
  
  switch(A.Order){
    case 0:
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in free__TensorLib__()",
	      "Not posible to free a scalar");
      exit(EXIT_FAILURE);       
    case 1:
      free(A.n);
      break;
    case 2:
      for(int i = 0 ; i<Ndim ; i++)
      {
	     free(A.N[i]);
      }
      break;
    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in free__TensorLib__()",
	      "Invalird order of the tensor");
      exit(EXIT_FAILURE); 
    }
}

/*************************************************************/

double I1__TensorLib__(const Tensor A)
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
	    "Error in I1__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I1;  
}

/*************************************************************/

double I2__TensorLib__(const Tensor A)
{
  int Ndim = NumberDimensions;  
  /* Define output */
  double I2 = 0;
  /* Check if is the order is order 2 */
  if(A.Order == 2)
  {

#if NumberDimensions == 3 
    I2 =A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0] +
    A.N[1][1]*A.N[2][2] - A.N[1][2]*A.N[2][1] +
    A.N[0][0]*A.N[2][2] - A.N[2][0]*A.N[0][2];
#endif

#if NumberDimensions == 2
    I2 = A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0];
#endif

  }
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in I2__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I2;  
}

/*************************************************************/

double I3__TensorLib__(const Tensor A)
{
  int Ndim = NumberDimensions;  
  /* Define output */
  double I3 = 0;
  /* Check if is the order is order 2 */
  if(A.Order == 2)
  {

#if NumberDimensions == 3 
      I3 =
      A.N[0][0]*A.N[1][1]*A.N[2][2] - 
      A.N[0][0]*A.N[1][2]*A.N[2][1] + 
      A.N[0][1]*A.N[1][2]*A.N[2][0] - 
      A.N[0][1]*A.N[1][0]*A.N[2][2] + 
      A.N[0][2]*A.N[1][0]*A.N[2][1] - 
      A.N[0][2]*A.N[1][1]*A.N[2][0] ; 
#endif

#if NumberDimensions == 2
      I3 = A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0];
#endif

  }
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in I3__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I3;
}

/*************************************************************/

double J1__TensorLib__(const Tensor A)
{
  /* Define output */
  double J1;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    J1 = I1__TensorLib__(A);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in J1__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return J1;
}

/*************************************************************/

double J2__TensorLib__(const Tensor A)
{
  /* Define output */
  double J2;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    double I1 = I1__TensorLib__(A);
    double I2 = I2__TensorLib__(A);
    J2 = pow(I1,2) - 2*I2;
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in J2__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return J2;
}

/*************************************************************/

double J3__TensorLib__(const Tensor A)
{
  /* Define output */
  double J3;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    double I1 = I1__TensorLib__(A);
    double I2 = I2__TensorLib__(A);
    double I3 = I3__TensorLib__(A);
    J3 = pow(I1,3) - 3*I1*I2 + 3*I3;
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in J3__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return J3;
}

/*************************************************************/

EigenTensor Eigen_analysis__TensorLib__(const Tensor A)
{

  /* Output variable */
  EigenTensor EigenA;

  /* Locals */
  int Ndim = NumberDimensions;
  int n = Ndim;
  int lda = Ndim;
  int ldvl = Ndim;
  int ldvr = Ndim;
  int info; 
  int lwork;
  double wkopt;
  double * work;

  /* Local arrays */
  double wr[NumberDimensions];
  double wi[NumberDimensions];
  double vl[NumberDimensions*NumberDimensions];
  double vr[NumberDimensions*NumberDimensions];

#if NumberDimensions == 3 
  double a[NumberDimensions*NumberDimensions] = 
  {
    A.N[0][0],A.N[0][1],A.N[0][2],
    A.N[1][0],A.N[1][1],A.N[1][2],
    A.N[2][0],A.N[2][1],A.N[2][2]
  }; 
#endif

#if NumberDimensions == 2
  double a[NumberDimensions*NumberDimensions] = 
  {
    A.N[0][0],A.N[0][1],
    A.N[1][0],A.N[1][1]
  }; 
#endif

  
  /* 
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dgeev_("N", "V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc(lwork*sizeof(double));
        
  /* Solve eigenproblem */
  dgeev_("N","V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,work, &lwork, &info );
  
  /* Check for convergence */
  if(info > 0)
  {
    printf("Error in Eigen_analysis__TensorLib__() : The algorithm failed to compute eigenvalues.\n" );
    exit(EXIT_FAILURE); 
  }

  free(work);

  /*
    Fill output
  */
  EigenA.Value = alloc__TensorLib__(1);
  EigenA.Vector = alloc__TensorLib__(2);
  
  for(int i = 0 ; i<Ndim ; i++)
  {
    EigenA.Value.n[i] = wr[i];

    for(int j = 0 ; j<Ndim ; j++)
    {
      EigenA.Vector.N[i][j] = vr[i+j*Ndim];
    }
  }

  return EigenA;
}


/*************************************************************/

double EuclideanNorm__TensorLib__(const Tensor A)
{
  int Ndim = NumberDimensions;
  double Aux = 0; 
  /* Define output */
  double Out;
  /* Check if the input is a first order tensor */
  if(A.Order == 1)
  {
    /* Compute norm */    
    for(int i = 0 ; i<Ndim ; i++){
      Aux += DSQR(A.n[i]);
    }
    Out = pow(Aux,0.5);
  }
  else if (A.Order == 2)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        Aux += DSQR(A.N[i][j]);
      }  
    }
    Out = pow(Aux,0.5);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in EuclideanNorm__TensorLib__()",
	    "The input should be of order 1");
    exit(EXIT_FAILURE);    
  }
  return Out;
}

/*********************************************************************/

double Generalised_norm__TensorLib__(const Tensor a, const Tensor G)
{
  double norm;
  Tensor G_dot_a;

  G_dot_a = vector_linear_mapping__TensorLib__(G, a);
  norm = inner_product__TensorLib__(a,G_dot_a);

  free__TensorLib__(G_dot_a);

  return norm;
}

/*************************************************************/

Tensor Identity__TensorLib__()
{
  int Ndim = NumberDimensions;
  Tensor I = alloc__TensorLib__(2);
  for(int i = 0 ; i<Ndim ; i++){
    I.N[i][i] = 1;
  }
  return I;
}

/*************************************************************/

Tensor Inverse__TensorLib__(const Tensor A)
{
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor Am1 = alloc__TensorLib__(2);
  /* Check if the input is a second order tensor */
  if (A.Order == 2){  
    /* Get the determinant of the matrix */
    double detA = I3__TensorLib__(A);
    if(fabs(detA) < TOL_zero)
    {
      fprintf(stderr,"%s\n","Error in Inverse__TensorLib__(A)");
      printf("%s\n","Input tensor A should be invertible");
      printf("det(A) : %e\n",detA);
      print__TensorLib__(A);
      exit(EXIT_FAILURE);   
    }

#if NumberDimensions == 3 
      Am1.N[0][0] = + (double)1/detA*(A.N[1][1]*A.N[2][2] - A.N[1][2]*A.N[2][1]);
      Am1.N[0][1] = - (double)1/detA*(A.N[0][1]*A.N[2][2] - A.N[0][2]*A.N[2][1]);
      Am1.N[0][2] = + (double)1/detA*(A.N[0][1]*A.N[1][2] - A.N[0][2]*A.N[1][1]);
      Am1.N[1][0] = - (double)1/detA*(A.N[1][0]*A.N[2][2] - A.N[1][2]*A.N[2][0]);
      Am1.N[1][1] = + (double)1/detA*(A.N[0][0]*A.N[2][2] - A.N[0][2]*A.N[2][0]);
      Am1.N[1][2] = - (double)1/detA*(A.N[0][0]*A.N[1][2] - A.N[0][2]*A.N[1][0]);
      Am1.N[2][0] = + (double)1/detA*(A.N[1][0]*A.N[2][1] - A.N[1][1]*A.N[2][0]);
      Am1.N[2][1] = - (double)1/detA*(A.N[0][0]*A.N[2][1] - A.N[0][1]*A.N[2][0]);
      Am1.N[2][2] = + (double)1/detA*(A.N[0][0]*A.N[1][1] - A.N[0][1]*A.N[1][0]);
#endif

#if NumberDimensions == 2
      Am1.N[0][0] = + (double)1/detA*A.N[1][1];
      Am1.N[0][1] = - (double)1/detA*A.N[0][1];
      Am1.N[1][0] = - (double)1/detA*A.N[1][0];
      Am1.N[1][1] = + (double)1/detA*A.N[0][0]; 
#endif

 
  }
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in Inverse__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the inverse matrix */
  return Am1;  
}

/*************************************************************/

Tensor Solve_system__TensorLib__(Tensor A, Tensor b)
{
  Tensor Am1;
  Tensor x;

  if ((A.Order == 2) && (b.Order == 1))
  {

    Am1 = Inverse__TensorLib__(A);

    x = vector_linear_mapping__TensorLib__(Am1, b);
 
  }
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
      "Error in Solve_system__TensorLib__()",
      "The input should be 2ord tensor and a 1rd tensor");
    exit(EXIT_FAILURE);
  }

  return x;
}

/*************************************************************/

Tensor transpose__TensorLib__(const Tensor A)
{
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor AT = alloc__TensorLib__(2);  
  /* Check if the input is a second order tensor */
  if (A.Order == 2)
  {
    /* Get the transpose */
    for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
      {
        AT.N[i][j] = A.N[j][i];
      }
    }
  }
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in transpose__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the transpose */
  return AT;
}

/*************************************************************/

Tensor subtraction__TensorLib__(Tensor A, Tensor B)
{
  int Ndim = NumberDimensions;

  /* Variable declaration output matrix */
  Tensor A_minus_B;

  /* Subtraction of second order tensors */
  if ((A.Order == 2) && (B.Order == 2))
  { 

    A_minus_B = alloc__TensorLib__(2); 

    for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
      {
        A_minus_B.N[i][j] = A.N[i][j] - B.N[i][j];
      }
    }
  }
  /* Subtraction of first order tensors */
  else if((A.Order == 1) && (B.Order == 1))
  {

    A_minus_B = alloc__TensorLib__(1);

    for(int i = 0 ; i < Ndim ; i++)
    {
      A_minus_B.n[i] = A.n[i] - B.n[i];
    }

  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
      "Error in subtraction__TensorLib__()",
      "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);        
  }

  return A_minus_B;
}

/*************************************************************/

Tensor addition__TensorLib__(Tensor A, Tensor B)
{
  int Ndim = NumberDimensions;

  /* Variable declaration output matrix */
  Tensor A_plus_B;

  /* Addition of second order tensors */
  if ((A.Order == 2) && (B.Order == 2))
  { 

    A_plus_B = alloc__TensorLib__(2); 

    for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
      {
        A_plus_B.N[i][j] = A.N[i][j] + B.N[i][j];
      }
    }
  }
  /* Addition of first order tensors */
  else if((A.Order == 1) && (B.Order == 1))
  {

    A_plus_B = alloc__TensorLib__(1);

    for(int i = 0 ; i < Ndim ; i++)
    {
      A_plus_B.n[i] = A.n[i] + B.n[i];
    }

  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
      "Error in addition__TensorLib__()",
      "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);        
  }

  return A_plus_B;
}

/*************************************************************/

double inner_product__TensorLib__(Tensor A, Tensor B)
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
	    "Error in inner_product__TensorLib__()",
	    "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);        
  }

  return AdotB;
    
}

/*************************************************************/

Tensor vector_product__TensorLib__(Tensor a, Tensor b)
{
  /* Allocate output */
  Tensor axb = alloc__TensorLib__(1);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)){
    /* Operate vetor product */
    axb.n[0] = a.n[1]*b.n[2] - a.n[2]*b.n[1];
    axb.n[1] = a.n[2]*b.n[0] - a.n[0]*b.n[2];
    axb.n[2] = a.n[0]*b.n[1] - a.n[1]*b.n[0];
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in vector_product__TensorLib__()",
	    "The input should be two tensors of first order");
    exit(EXIT_FAILURE);    
  }
  /* Return vector product */
  return axb; 
}

/*************************************************************/

Tensor dyadic_Product__TensorLib__(Tensor a, Tensor b)
{
  int Ndim = NumberDimensions;
  /* Tensor declaration */
  Tensor aob = alloc__TensorLib__(2);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1))
  {
    /* Operate tensor product */
    for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        aob.N[i][j] = a.n[i]*b.n[j];
      } 
    }
  }
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in dyadic_Product__TensorLib__()",
	    "The input should be two tensors of first order");
    exit(EXIT_FAILURE);    
  }
  /* Return tensorial product */
  return aob;
}

/*************************************************************/

Tensor vector_linear_mapping__TensorLib__(Tensor A, Tensor b)
{
  int Ndim = NumberDimensions;
  /* Tensor declaration */
  Tensor Adotb = alloc__TensorLib__(1);
  /* Check in the input its is ok */
  if ((A.Order == 2) && (b.Order == 1))
  {
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
  else
  {
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in vector_linear_mapping__TensorLib__()",
	    "The input should be 2ord tensor and a 1rd tensor");
    exit(EXIT_FAILURE);
  }
  /* Return tensor */
  return Adotb;
}

/*************************************************************/

Tensor matrix_product__TensorLib__(Tensor A, Tensor B)
{
  int Ndim = NumberDimensions;  
  Tensor A_x_B = alloc__TensorLib__(2);

  if ( (A.Order == 2) && (B.Order == 2) )
    { 
    for(int i = 0 ; i < Ndim  ; i++)
      {
	for(int j = 0 ; j < Ndim  ; j++)
	  {
	    for(int k = 0 ; k < Ndim  ; k++)
	      {
		A_x_B.N[i][j] += A.N[i][k]*B.N[k][j];
	      }
	  }
      }
  }
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in matrix_product__TensorLib__()",
	      "The input should be two second order tensors");
      exit(EXIT_FAILURE);        
    }
  
  return A_x_B;
}


/*************************************************************/


Tensor Convex_combination__TensorLib__(Tensor F_n1,Tensor F_n,double alpha)
{
  /* Define output */
  Tensor F_alpha = alloc__TensorLib__(2);
  /* Define the number of dimensions */
  int Ndim = NumberDimensions;

  for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
	{
	  F_alpha.N[i][j] = alpha*F_n1.N[i][j] + (1-alpha)*F_n.N[i][j];
	}
    }

  return F_alpha;
}

/*************************************************************/

Tensor volumetric_component__TensorLib__(Tensor A)
{
  int Ndim = NumberDimensions;
  Tensor A_vol = alloc__TensorLib__(2); 
  double trA = I1__TensorLib__(A);

  for(int i = 0 ; i<Ndim ; i++)
  {
    A_vol.N[i][i] = trA/(double)Ndim;
  }

  return A_vol;
}

/*************************************************************/

Tensor deviatoric_component__TensorLib__(Tensor A, Tensor A_vol)
{
  int Ndim = NumberDimensions;
  Tensor A_dev = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      A_dev.N[i][j] =  A.N[i][j] - A_vol.N[i][j];
    }
  }

  return A_dev;
}

/*************************************************************/

Tensor rotate__TensorLib__(Tensor In, Tensor R)
/* 
  Return the rotated tensor of the imput using the
  rotation matrix R.
*/
{
  int Ndim = NumberDimensions;

  Tensor Rm1 = Inverse__TensorLib__(R);

  Tensor In__x__Rm1 = matrix_product__TensorLib__(In,Rm1);

  Tensor Out = matrix_product__TensorLib__(R,In__x__Rm1);

  free__TensorLib__(Rm1);
  free__TensorLib__(In__x__Rm1);

  return Out;

}

/*************************************************************/

Tensor symmetrise__TensorLib__(Tensor A)
{
  Tensor symA = alloc__TensorLib__(2);

  int Ndim = NumberDimensions;

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      if(i != j)
      {
        symA.N[i][j] = 0.5*(A.N[i][j] + A.N[j][i]);
      }
      else
      {
        symA.N[i][j] = A.N[i][j];
      }
    }
  }

  return symA;
}


/*************************************************************/

Tensor covariant_push_forward_tensor__TensorLib__(Tensor A, Tensor F)
/* 
  Covariant push forward operation for any tensor. 
  a = F^-T A F^-1
  From literature, this operation moves a tensor from the material description (A) to the
  spatial description (a). 
  a (out)
  A (in)
  F (in)
*/
{
  int Ndim = NumberDimensions;
  Tensor F_m1 = Inverse__TensorLib__(F);
  Tensor F_mT = transpose__TensorLib__(F_m1);

  Tensor A__x__F_m1 = matrix_product__TensorLib__(A,F_m1);
  Tensor F_mT__x__A__x__F_m1 = matrix_product__TensorLib__(F_mT,A__x__F_m1);

  free__TensorLib__(F_m1);
  free__TensorLib__(F_mT);
  free__TensorLib__(A__x__F_m1);

  return F_mT__x__A__x__F_m1;

}

/*********************************************************************/

void contravariant_push_forward_tensor__TensorLib__(Tensor a, Tensor A, Tensor F)
/* 
  Contravariant push forward operation for any tensor. 
  a = F A F^T
  From literature, this operation moves a tensor in the dual space from the material description (A) to the
  spatial description (a). 
  a (out)
  A (in)
  F (in)
*/
{
  int Ndim = NumberDimensions;

  if (Ndim == 2)
  {
    double aux_1 = A.N[0][0]*F.N[0][0] + A.N[0][1]*F.N[0][1];
    double aux_2 = A.N[0][0]*F.N[1][0] + A.N[0][1]*F.N[1][1];
    double aux_3 = A.N[1][0]*F.N[0][0] + A.N[1][1]*F.N[0][1];
    double aux_4 = A.N[1][0]*F.N[1][0] + A.N[1][1]*F.N[1][1];

    a.N[0][0] = F.N[0][0]*aux_1 + F.N[0][1]*aux_3;
    a.N[0][1] = F.N[0][0]*aux_2 + F.N[0][1]*aux_4;
    a.N[1][0] = F.N[1][0]*aux_1 + F.N[1][1]*aux_3;
    a.N[1][1] = F.N[1][0]*aux_2 + F.N[1][1]*aux_4;
  }
    else if(Ndim == 3)
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in contravariant_push_forward_tensor__TensorLib__()",
        "This operation it is not implemented for 3D");
      exit(EXIT_FAILURE);  
  }
  else
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in contravariant_push_forward_tensor__TensorLib__()",
        "Wrong number of dimensions");
      exit(EXIT_FAILURE);  
  }

}


/*********************************************************************/

void covariant_pull_back_tensor__TensorLib__(Tensor A, Tensor a, Tensor F)
/* 
  Covariant pull back operation for any tensor. 
  A = F^T a F
  From literature, this operation moves a tensor from the spatial description (a) to the
  material description (A).
  A (out)
  a (in)
  F (int)
*/
{
  int Ndim = NumberDimensions;

  if (Ndim == 2)
  {

    double aux_1 = a.N[0][0]*F.N[0][0] + a.N[0][1]*F.N[1][0];
    double aux_2 = a.N[0][0]*F.N[0][1] + a.N[0][1]*F.N[1][1];
    double aux_3 = a.N[1][0]*F.N[0][0] + a.N[1][1]*F.N[1][0];
    double aux_4 = a.N[1][0]*F.N[0][1] + a.N[1][1]*F.N[1][1];

    A.N[0][0] = F.N[0][0]*aux_1 + F.N[1][0]*aux_3;
    A.N[0][1] = F.N[0][0]*aux_2 + F.N[1][0]*aux_4;
    A.N[1][0] = F.N[0][1]*aux_1 + F.N[1][1]*aux_3;
    A.N[1][1] = F.N[0][1]*aux_2 + F.N[1][1]*aux_4;
  }
    else if(Ndim == 3)
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in covariant_pull_back_tensor__TensorLib__()",
        "This operation it is not implemented for 3D");
      exit(EXIT_FAILURE);  
  }
  else
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in covariant_pull_back_tensor__TensorLib__()",
        "Wrong number of dimensions");
      exit(EXIT_FAILURE);  
  }

}


/*********************************************************************/

Tensor contravariant_pull_back_tensor__TensorLib__(Tensor a, Tensor F)
/* 
  Contravariant pull back operation for any tensor. 
  A = F^-1 a F^-T
  From literature, this operation moves a tensor in the dual space from the spatial description (a) to the
  material description (A).
  A (out)
  a (in)
  F (int)
*/
{
  int Ndim = NumberDimensions;
  Tensor F_m1 = Inverse__TensorLib__(F);
  Tensor F_mT = transpose__TensorLib__(F_m1);

  Tensor a__x__F_mT = matrix_product__TensorLib__(a,F_mT);
  Tensor F_m1__x__a__x__F_mT = matrix_product__TensorLib__(F_m1,a__x__F_mT);

  free__TensorLib__(F_m1);
  free__TensorLib__(F_mT);
  free__TensorLib__(a__x__F_mT);

  return F_m1__x__a__x__F_mT;
}

/*********************************************************************/

void print__TensorLib__(Tensor A)
{
  int Ndim = NumberDimensions;
  
  if(A.Order == 2)
    {  
      for(int i = 0 ; i < Ndim  ; i++)
	{
	  for(int j = 0 ; j < Ndim  ; j++)
	    {
	      printf("%e ",A.N[i][j]);
	    }
	  printf("\n");
	}
    }
  else if(A.Order == 1)
    {
      for(int i = 0 ; i < Ndim  ; i++)
	{
	  printf("%e ",A.n[i]);
	}
      printf("\n");      
    }
}

/*************************************************************/









