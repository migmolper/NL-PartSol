#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


/*********************************************************************/

double Area_Poligon(Matrix Poligon)
/*
  Get the area of a poligon on n vertex using the 
  Shoelace formula :
  https://en.wikipedia.org/wiki/Shoelace_formula
*/
{
  double A_Poligon;
  int N_vertex;

  /* Asign the number of vertex */
  N_vertex = Poligon.N_rows;

  /* Check the number of vertex */
  if(N_vertex<3){
    puts("Error in Area_Poligon : Wrong number of vertex !!!");
    exit(0);
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

Matrix Centroid_Poligon(Matrix Poligon)
/*
  Get the position of the centroid of a poligon, also return the value of the 
  area
*/
{

  Matrix C_Poligon;
  double A_aux;
  int N_vertex;
  int N_dimensions;

  /* Asign the number of vertex and check it  */
  N_vertex = Poligon.N_rows; 
  if(N_vertex<3){
    puts("Error in Centroid_Poligon() : Wrong number of vertex !!!");
    exit(0);
  }

  /* Asign the number of dimensions and check it */
  N_dimensions = Poligon.N_cols;
  if(N_dimensions != 2){
    puts("Error in Centroid_Poligon() : Wrong number of dimensions !!!");
    exit(0);
  }
  
  /* Allocate the array for the centroid */
  C_Poligon = MatAllocZ(1,N_dimensions);
  
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

int InOut_Poligon(Matrix X_Point, Matrix Poligon)
/*
  Check if a point is or not (1/0) inside of a poligon.
  Inputs :
  - X_Point : Coordinates of the point 
  - Poligon : Coordinates of the vertex 0,1,....,n,0
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

/*********************************************************************/

/* double MinLenght_Poligon(Matrix Poligon){ */

  
/* } */

/*********************************************************************/

double SignumFunct(double x)
/*
  
*/
{
 
  if(x > 0){
    return 1;
  }
  else if (x < 0) {
    return -1;
  }
  else {
    return 0;
  }
  
}

/*********************************************************************/

Matrix SolvePolynomial(Matrix Coeffs)
/*
  
*/
{
  if( (Coeffs.N_rows > 1) && (Coeffs.N_cols > 1) ){
    printf("%s : %s \n","Error in SolvePolynomial",
	   "I am not able to solve a system !!!");
    exit(0);
  }
  
  Matrix Solution;
  int N_Coeffs = Coeffs.N_rows*Coeffs.N_cols;
  double a = Coeffs.nV[0];
  double b = Coeffs.nV[1];
  double c = Coeffs.nV[2];

  switch(N_Coeffs){

  case 1 :
    printf(" %s : %s \n ",
	   "Error in SolvePolynomial()",
	   "Dummy polynomial of order 0");
    exit(0);
    break;
    
  case 2 :
    printf(" %s : %s \n ",
	   "Error in SolvePolynomial()",
	   "Dummy polynomial of order 1");
    exit(0);
    break;
    
  case 3 :
    Solution = MatAlloc(2,1);
    if((b*b - 4*a*c)<0){
      printf("%s : %s \n",
	     "Error in SolvePolynomial()",
	     "Imaginary solutions not implemented");
      exit(0);
    }
    Solution.nV[0] =
      (- b + sqrt(b*b - 4*a*c))/(2*a);
    Solution.nV[1] =
      (- b - sqrt(b*b - 4*a*c))/(2*a);
    break;
       
  default :
    printf(" %s : %s \n ",
	   "Error in SolvePolynomial() ",
	   "I am only able to solve 2 order polynomials !");
    exit(0); 
  }

  return Solution;
  
}

/*********************************************************************/