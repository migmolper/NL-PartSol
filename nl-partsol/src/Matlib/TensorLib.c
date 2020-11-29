#include "nl-partsol.h"

/*************************************************************/

Tensor alloc__TensorLib__(int Order)
{
  int Ndim = NumberDimensions;
  /* Define output */
  Tensor A;
  /* Swith cases */
  switch(Order){
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
  switch(Order){
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
    for(long i = 0 ; i<Ndim ; i++){
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
      for(int i = 0 ; i<Ndim ; i++){
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

double I1__TensorLib__(Tensor A)
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

double I2__TensorLib__(Tensor A)
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
	    "Error in I2__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I2;  
}

/*************************************************************/

double I3__TensorLib__(Tensor A)
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
	    "Error in I3__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return I3;
}

/*************************************************************/

double J1__TensorLib__(Tensor A)
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

double J2__TensorLib__(Tensor A)
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

double J3__TensorLib__(Tensor A)
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


Tensor Eigenvalues__TensorLib__(Tensor A)
{
  /* Auxiliar variables */
  int Ndim = NumberDimensions;
  double I1, I2, I3;
  double b, c, m, n, t, x;

  /* Define output */
  Tensor lambda = alloc__TensorLib__(1);

  /* Check the order of the input tensor */
  if(A.Order == 2){

    /*
      Solve the second order equation
      x² - I1*x + I2 = 0 
     */
    if(Ndim == 2)
      {
	I1 = I1__TensorLib__(A);
	I2 = I2__TensorLib__(A);

	b = I1*I1 - 4*I2;

	if(b > TOL_zero)
	  {
	    
	    lambda.n[0] = 0.5*(I1 + sqrt(b));
	    lambda.n[1] = 0.5*(I1 - sqrt(b));

	  }
	else if(fabs(b) <= TOL_zero)
	  {
	    lambda.n[0] = 0.5*I1;
	    lambda.n[1] = 0.5*I1;
	  }

	else
	  {

	    printf("%s : %s -> %f \n",
		   "Error in Eigenvalues__TensorLib__()",
		   "Input tensor should be Hermitian",b);
	    exit(EXIT_FAILURE);

	  }
	
           
      }

    /*
      Solve the third order equation via Cardano
    */
    if(Ndim == 3)
      {
	I1 = I1__TensorLib__(A);
	I2 = I2__TensorLib__(A);
	I3 = I3__TensorLib__(A);

	b = I2 - I1*I1/3;
	c = -2/27*I1*I1*I1 + I1*I2/3 - I3;

	if(fabs(b) <= TOL_zero)
	  {
	    
	    lambda.n[0] = - pow(c,1/3.);
	    lambda.n[1] = - pow(c,1/3.);
	    lambda.n[2] = - pow(c,1/3.);
	    
	  }
	else if(b > TOL_zero)
	  {

	    m = 2*sqrt(-b/3);
	    
	  }
	else
	  {

	    printf("%s : %s -> %f \n",
		   "Error in Eigenvalues__TensorLib__()",
		   "Input tensor should be Hermitian",b);
	    exit(EXIT_FAILURE);
	    
	  }

	
      }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in Eigenvalues__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }

  return lambda;
}

/*************************************************************/

Tensor Eigenvectors__TensorLib__(Tensor A,Tensor Lambda)
{

  int Ndim = NumberDimensions;
  Tensor EigenVect = alloc__TensorLib__(2);

  if (Ndim == 2)
  {

    if (A.N[0][1]*A.N[1][0] < 1e-15)
    {
      EigenVect.N[0][0] = EigenVect.N[1][1] = 1.;
      EigenVect.N[0][1] = EigenVect.N[1][0] = 0.;
    }
    else
    {
      for(int i = 0 ; i<Ndim ; i++)
      {
        double aux1 = (Lambda.n[i] - A.N[0][0])/A.N[0][1];
        double aux2 = sqrt(1 + aux1*aux1);
        EigenVect.N[0][i] = 1/aux2;
        EigenVect.N[1][i] = aux1/aux2;
      }
    }


  }
  else if(Ndim == 3)
    {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in spectral_descomposition_symmetric__TensorLib__()",
        "This operation it is not implemented for 3D");
      exit(EXIT_FAILURE);  
    }
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in spectral_descomposition_symmetric__TensorLib__()",
        "Wrong number of dimensions");
      exit(EXIT_FAILURE);  
    }


  return EigenVect;
}

/*************************************************************/

double EuclideanNorm__TensorLib__(Tensor A)
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

Tensor Inverse__TensorLib__(Tensor A)
{
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor Am1 = alloc__TensorLib__(2);
  /* Check if the input is a second order tensor */
  if (A.Order == 2){  
    /* Get the determinant of the matrix */
    double detA = I3__TensorLib__(A);
    if(detA == 0){
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in Inverse__TensorLib__()",
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
	    "Error in Inverse__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the inverse matrix */
  return Am1;  
}

/*************************************************************/

Tensor transpose__TensorLib__(Tensor A)
{
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor AT = alloc__TensorLib__(2);  
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
	    "Error in transpose__TensorLib__()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the transpose */
  return AT;
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

double volumetric_component__TensorLib__(Tensor A)
{
  return I1__TensorLib__(A)/3.0;
}

/*************************************************************/

Tensor deviatoric_component__TensorLib__(Tensor A, double A_vol)
{
  int Ndim = NumberDimensions;
  Tensor A_dev = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      A_dev.N[i][j] =  A.N[i][j] - (i == j)*A_vol;
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
  Tensor Out = alloc__TensorLib__(2);

  if (Ndim == 2)
  {

    double aux_1 = R.N[0][0]*In.N[0][0] + R.N[0][1]*In.N[0][1];
    double aux_2 = R.N[0][0]*In.N[1][0] + R.N[0][1]*In.N[1][1];
    double aux_3 = R.N[1][0]*In.N[0][0] + R.N[1][1]*In.N[0][1];
    double aux_4 = R.N[1][0]*In.N[1][0] + R.N[1][1]*In.N[1][1];

    Out.N[0][0] = R.N[0][0]*aux_1 + R.N[0][1]*aux_2;
    Out.N[0][1] = R.N[0][0]*aux_3 + R.N[0][1]*aux_4;
    Out.N[1][0] = R.N[1][0]*aux_1 + R.N[1][1]*aux_2;
    Out.N[1][1] = R.N[1][0]*aux_3 + R.N[1][1]*aux_4;

  }
  else if(Ndim == 3)
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in rotate__TensorLib__()",
        "This operation it is not implemented for 3D");
      exit(EXIT_FAILURE);  
  }
  else
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in rotate__TensorLib__()",
        "Wrong number of dimensions");
      exit(EXIT_FAILURE);  
  }

  return Out;

}

/*************************************************************/

void covariant_push_forward_tensor__TensorLib__(Tensor a, Tensor A, Tensor F)
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

  if (Ndim == 2)
  {

    double aux_1 = A.N[0][0]*F_m1.N[0][0] + A.N[0][1]*F_m1.N[1][0];
    double aux_2 = A.N[0][0]*F_m1.N[0][1] + A.N[0][1]*F_m1.N[1][1];
    double aux_3 = A.N[1][0]*F_m1.N[0][0] + A.N[1][1]*F_m1.N[1][0];
    double aux_4 = A.N[1][0]*F_m1.N[0][1] + A.N[1][1]*F_m1.N[1][1];

    a.N[0][0] = F_m1.N[0][0]*aux_1 + F_m1.N[1][0]*aux_3;
    a.N[0][1] = F_m1.N[0][0]*aux_2 + F_m1.N[1][0]*aux_4;
    a.N[1][0] = F_m1.N[0][1]*aux_1 + F_m1.N[1][1]*aux_3;
    a.N[1][1] = F_m1.N[0][1]*aux_2 + F_m1.N[1][1]*aux_4;

  }
  else if(Ndim == 3)
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in covariant_push_forward_tensor__TensorLib__()",
        "This operation it is not implemented for 3D");
      exit(EXIT_FAILURE);  
  }
  else
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in covariant_push_forward_tensor__TensorLib__()",
        "Wrong number of dimensions");
      exit(EXIT_FAILURE);  
  }

  free__TensorLib__(F_m1);

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

void contravariant_pull_back_tensor__TensorLib__(Tensor A, Tensor a, Tensor F)
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

  if (Ndim == 2)
  {

    double aux_1 = a.N[0][0]*F_m1.N[0][0] + a.N[0][1]*F_m1.N[0][1];
    double aux_2 = a.N[0][0]*F_m1.N[1][0] + a.N[0][1]*F_m1.N[1][1];
    double aux_3 = a.N[1][0]*F_m1.N[0][0] + a.N[1][1]*F_m1.N[0][1];
    double aux_4 = a.N[1][0]*F_m1.N[1][0] + a.N[1][1]*F_m1.N[1][1];

    A.N[0][0] = F_m1.N[0][0]*aux_1 + F_m1.N[0][1]*aux_3;
    A.N[0][1] = F_m1.N[0][0]*aux_2 + F_m1.N[0][1]*aux_4;
    A.N[1][0] = F_m1.N[1][0]*aux_1 + F_m1.N[1][1]*aux_3;
    A.N[1][1] = F_m1.N[1][0]*aux_2 + F_m1.N[1][1]*aux_4;
  }
    else if(Ndim == 3)
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in contravariant_pull_back_tensor__TensorLib__()",
        "This operation it is not implemented for 3D");
      exit(EXIT_FAILURE);  
  }
  else
  {
      fprintf(stderr,"%s : %s !!! \n",
        "Error in contravariant_pull_back_tensor__TensorLib__()",
        "Wrong number of dimensions");
      exit(EXIT_FAILURE);  
  }

  free__TensorLib__(F_m1);
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
	      printf("%f ",A.N[i][j]);
	    }
	  printf("\n");
	}
    }
  else if(A.Order == 1)
    {
      for(int i = 0 ; i < Ndim  ; i++)
	{
	  printf("%f ",A.n[i]);
	}
      printf("\n");      
    }
}

/*************************************************************/









