#include <math.h>
#include "nl-partsol.h"

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
  if ((Coeffs.N_rows > 1) && (Coeffs.N_cols > 1)) {
    printf("%s : %s \n", "Error in solve_polynomial__MatrixLib__",
           "I am not able to solve a system !!!");
    exit(EXIT_FAILURE);
  }

  Matrix Solution;
  int N_Coeffs = Coeffs.N_rows * Coeffs.N_cols;

  if (N_Coeffs == 1) {
    printf(" %s : %s \n ", "Error in solve_polynomial__MatrixLib__()",
           "Dummy polynomial of order 0");
    exit(EXIT_FAILURE);
  }
  if (N_Coeffs == 2) {
    printf(" %s : %s \n ", "Error in solve_polynomial__MatrixLib__()",
           "Dummy polynomial of order 1");
    exit(EXIT_FAILURE);
  }

  double a, b, c, d;
  double aux;

  switch (N_Coeffs) {
  case 3:
    /* Coeficients */
    a = Coeffs.nV[0];
    b = Coeffs.nV[1];
    c = Coeffs.nV[2];

    Solution = alloc__MatrixLib__(2, 1);
    aux = b * b - 4 * a * c;
    if (fabs(aux) < TOL_zero) {
      aux = 0.0;
    } else if (aux < TOL_zero) {
      printf("%s : %s -> %f \n", "Error in solve_polynomial__MatrixLib__()",
             "Imaginary solutions not implemented", aux);
      exit(EXIT_FAILURE);
    }
    Solution.nV[0] = 0.5 * (-b + sqrt(aux)) / a;
    Solution.nV[1] = 0.5 * (-b - sqrt(aux)) / a;
    break;

  case 4:
    /* Coeficients */
    a = Coeffs.nV[0];
    b = Coeffs.nV[1];
    c = Coeffs.nV[2];
    d = Coeffs.nV[3];
    break;
  default:
    printf(" %s : %s \n ", "Error in solve_polynomial__MatrixLib__() ",
           "I am only able to solve 2 order polynomials !");
    exit(EXIT_FAILURE);
  }

  return Solution;
}

/*********************************************************************/
