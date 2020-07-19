#include "nl-partsol.h"

/*********************************************************************/

double area__MatrixLib__(Matrix Poligon)
/*! 
 * This function returns the the area of a \a Poligon using the 
 * Shoelace formula :
 * https://en.wikipedia.org/wiki/Shoelace_formula
 */
{
  double A_Poligon;
  int N_vertex;

  /* Asign the number of vertex */
  N_vertex = Poligon.N_rows;

  /* Check the number of vertex */
  if(N_vertex<3){
    puts("Error in area__MatrixLib__ : Wrong number of vertex !!!");
    exit(EXIT_FAILURE);
  }
  
  /* Initialize the area value with the initial therm */
  A_Poligon = Poligon.nM[0][1]*(Poligon.nM[N_vertex-1][0] - Poligon.nM[1][0]);

  /* Addd the middle therm */
  for(int i = 1 ; i<N_vertex-1 ; i++){
    A_Poligon += Poligon.nM[i][1]*(Poligon.nM[i-1][0] - Poligon.nM[i+1][0]); 
  }

  /* Addd the final therm */
  A_Poligon += Poligon.nM[N_vertex-1][1]*(Poligon.nM[N_vertex-2][0] -
					  Poligon.nM[0][0]);

  /* Get the area */
  A_Poligon = 0.5*fabs(A_Poligon);
    
  return A_Poligon;
}

/*********************************************************************/

Matrix centroid__MatrixLib__(Matrix Poligon)
/*! 
 * This function returns the centroid of a \a Poligon and the area.
 */
{

  Matrix C_Poligon;
  double A_aux;
  int N_vertex;
  int N_dimensions;

  /* Asign the number of vertex and check it  */
  N_vertex = Poligon.N_rows; 
  if(N_vertex<3){
    puts("Error in centroid__MatrixLib__() : Wrong number of vertex !!!");
    exit(EXIT_FAILURE);
  }

  /* Asign the number of dimensions and check it */
  N_dimensions = Poligon.N_cols;
  if(N_dimensions != 2){
    puts("Error in centroid__MatrixLib__() : Wrong number of dimensions !!!");
    exit(EXIT_FAILURE);
  }
  
  /* Allocate the array for the centroid */
  C_Poligon = allocZ__MatrixLib__(1,N_dimensions);
  
  /* Initialize the area of the poligon */
  C_Poligon.n = 0;
  
  /* Calcule the area and the centroid */
  for(int i = 0 ; i<N_vertex-1 ; i++){

    A_aux = Poligon.nM[i][0]*Poligon.nM[i+1][1] -
      Poligon.nM[i+1][0]*Poligon.nM[i][1];

    C_Poligon.nV[0] += (Poligon.nM[i][0] + Poligon.nM[i+1][0])*A_aux;
    C_Poligon.nV[1] += (Poligon.nM[i][1] + Poligon.nM[i+1][1])*A_aux;
    C_Poligon.n += A_aux;
   
  }

  /* Addd the final therm */
  A_aux = Poligon.nM[N_vertex-1][0]*Poligon.nM[0][1] -
    Poligon.nM[0][0]*Poligon.nM[N_vertex-1][1];

  C_Poligon.nV[0] += (Poligon.nM[N_vertex-1][0] + Poligon.nM[0][0])*A_aux;
  C_Poligon.nV[1] += (Poligon.nM[N_vertex-1][1] + Poligon.nM[0][1])*A_aux;
  C_Poligon.n += A_aux;

  /* Fix the area */
  C_Poligon.n *= 0.5;

  /* Fix the centroid */
  C_Poligon.nV[0] /= 6*C_Poligon.n;
  C_Poligon.nV[1] /= 6*C_Poligon.n;

  /* Fix the sign of the area */
  C_Poligon.n = fabs(C_Poligon.n);

  return C_Poligon;
}


/*********************************************************************/

int inout__MatrixLib__(Matrix X_Point, Matrix Poligon)
/*! 
 * Check if a point is or not (1/0) inside of a Poligon.
 * Inputs :
 * - \a X_Point : Coordinates of the point 
 * - \a Poligon : Coordinates of the vertex 0,1,....,n,0
 */
{
  /* By default, we suppose that the point is in the poligon */
  int InOut = 1;

  Matrix a = allocZ__MatrixLib__(3,1);
  Matrix b = allocZ__MatrixLib__(3,1);
  Matrix c;
  Matrix n;
  Matrix nxc;

  /* Get the normal vector */
  a.nV[0] = Poligon.nM[1][0] - Poligon.nM[0][0];
  a.nV[1] = Poligon.nM[1][1] - Poligon.nM[0][1];
  a.nV[2] = 0;  
  b.nV[0] = Poligon.nM[Poligon.N_rows-1][0] - Poligon.nM[0][0];
  b.nV[1] = Poligon.nM[Poligon.N_rows-1][1] - Poligon.nM[0][1];
  b.nV[2] = 0;
  n = vectorial_product__MatrixLib__(a,b);
  n.N_rows = 1;
  n.N_cols = 3;

  /* Fill a and b for the First search */
  a.nV[0] = Poligon.nM[0][0] - Poligon.nM[Poligon.N_rows-1][0];
  a.nV[1] = Poligon.nM[0][1] - Poligon.nM[Poligon.N_rows-1][1];
  a.nV[2] = 0;

  b.nV[0] = X_Point.nV[0] - Poligon.nM[Poligon.N_rows-1][0];
  b.nV[1] = X_Point.nV[1] - Poligon.nM[Poligon.N_rows-1][1];
  b.nV[2] = 0;
  
  for(int i = 0 ; i<Poligon.N_rows-1 ; i++){

    c = vectorial_product__MatrixLib__(a,b);
    nxc = scalar_product__MatrixLib__(n,c);
    free__MatrixLib__(c);

    if(nxc.n < -TOL_InOut){
      InOut = 0;
      break;
    }
    
    a.nV[0] = Poligon.nM[i+1][0] - Poligon.nM[i][0];
    a.nV[1] = Poligon.nM[i+1][1] - Poligon.nM[i][1];
    a.nV[2] = 0;

    b.nV[0] = X_Point.nV[0] - Poligon.nM[i][0];
    b.nV[1] = X_Point.nV[1] - Poligon.nM[i][1];
    b.nV[2] = 0;
    
  }

  /* Last check */
  c = vectorial_product__MatrixLib__(a,b);
  nxc = scalar_product__MatrixLib__(n,c);
  if(nxc.n < -TOL_InOut){
    InOut = 0;
  }

  free__MatrixLib__(a);
  free__MatrixLib__(b);
  free__MatrixLib__(c);
  free__MatrixLib__(n);

  return InOut;
}

/*********************************************************************/

Matrix solve_polynomial__MatrixLib__(Matrix Coeffs)
/*
  Solve this kind of polinom : 
  a*x^2  + b*x + c = 0

  Inputs :
  Coeffs.nV[0] ->  a
  Coeffs.nV[1] ->  b
  Coeffs.nV[2] ->  c 
*/
{
  if( (Coeffs.N_rows > 1) && (Coeffs.N_cols > 1) ){
    printf("%s : %s \n","Error in solve_polynomial__MatrixLib__",
	   "I am not able to solve a system !!!");
    exit(EXIT_FAILURE);
  }
  
  Matrix Solution;
  int N_Coeffs = Coeffs.N_rows*Coeffs.N_cols;

  if(N_Coeffs == 1){
    printf(" %s : %s \n ",
	   "Error in solve_polynomial__MatrixLib__()",
	   "Dummy polynomial of order 0");
    exit(EXIT_FAILURE);
  }
  if(N_Coeffs == 2){
    printf(" %s : %s \n ",
	   "Error in solve_polynomial__MatrixLib__()",
	   "Dummy polynomial of order 1");
    exit(EXIT_FAILURE);
  }

  double a, b, c, d;
  double aux;

  switch(N_Coeffs){
  case 3 :
    /* Coeficients */
    a = Coeffs.nV[0];
    b = Coeffs.nV[1];
    c = Coeffs.nV[2];
    
    Solution = alloc__MatrixLib__(2,1);
    aux = b*b - 4*a*c; 
    if(fabs(aux)<TOL_zero){
      aux = 0.0;
    }
    else if(aux<TOL_zero){
      printf("%s : %s -> %f \n",
	     "Error in solve_polynomial__MatrixLib__()",
	     "Imaginary solutions not implemented",
	     aux);
      exit(EXIT_FAILURE);
    }
    Solution.nV[0] =
      0.5*(- b + sqrt(aux))/a;
    Solution.nV[1] =
      0.5*(- b - sqrt(aux))/a;
    break;

  case 4 :
    /* Coeficients */
    a = Coeffs.nV[0];
    b = Coeffs.nV[1];
    c = Coeffs.nV[2];
    d = Coeffs.nV[3];
    break;
  default :
    printf(" %s : %s \n ",
	   "Error in solve_polynomial__MatrixLib__() ",
	   "I am only able to solve 2 order polynomials !");
    exit(EXIT_FAILURE);
  }

  return Solution;
  
}

/*********************************************************************/

Matrix nurbs_distance__MatrixLib__(Matrix NURBS){

  int Ndim = NumberDimensions;
  int Nnodes = NURBS.N_rows;  
  Matrix Distance = allocZ__MatrixLib__(Nnodes,1);
  double aux;

  for(int i = 0 ; i<Nnodes ; i++){
    aux = 0;
    for(int j = 0 ; j<Ndim ; j++){
      aux += DSQR(NURBS.nM[i][j]);
    }
    Distance.nV[i] = pow(aux,0.5);
  }
  return Distance;
}

/*********************************************************************/

double point_distance__MatrixLib__(Matrix End, Matrix Init)
/*
  Distance between two points.
*/
{
  int Ndim = NumberDimensions;

  if((End.N_rows != Init.N_rows) ||
     (End.N_cols != Init.N_cols)){
    printf("%s : %s \n",
	   "Error in point_distance__MatrixLib__",
	   "Inputs arrays are not equal");
    exit(EXIT_FAILURE);
  }

  double DIST = 0;

  for(int i = 0 ; i<Ndim ; i++){

    DIST += pow(End.nV[i] - Init.nV[i],2);

  }

  DIST = sqrt(DIST);
  

  return DIST;
}


/*********************************************************************/
