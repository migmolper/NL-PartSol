#include <stdio.h>

/* Fields declaration */
#define KindProblem 2

#if KindProblem == 1
typedef struct {
  double xx;
} Tensor;

typedef struct {
  double x;  
} Vector;

#elif KindProblem == 2
typedef struct {
  double xx;
  double xy;
  double yx;
  double yy;
} Tensor;

typedef struct {
  double x;
  double y;
} Vector;

#elif KindProblem == 3
typedef struct{
  double xx;
  double xy;
  double yx;
  double yy;
  double xz;
  double yz;
  double zx;
  double zy;
  double zz;
} Tensor;

typedef struct {
  double x;
  double y;
  double z;  
} Vector;

#else
put("Check number of dimensions");
#endif

/* /\* Array definition *\/ */
/* typedef struct { */
/*   int Shape; */
/*   void * n;  */
/* } Array; */


/* /\* Array definition *\/ */
/* typedef struct { */
/*   int Columns; */
/*   int Rows; */
/*   void ** a;  */
/* } Matrix; */

/* Gauss points definitions */
typedef struct {
  int id; /* Identification number of the GP */
  int Element_id; /* Pointer to the element where it is */
  int Material; /* Id of the Material */
  Vector * x_GC; /* Position of the GP with global coordiantes */
  Vector * x_EC; /* Position of the GP with element coordiantes */
  Vector * v; /* Velocity field */
  Vector * a; /* Acceleration field */
  Tensor * Stress; /* Stress field */
  Tensor * Strain; /* Strain field */
} GaussPoint;


/* Element type definition */
typedef struct {
  int id; /* Identification number of the element */
  int NumberNodes; /* Number of nodes of the element */
  int NumberDOF; /* Degrees of freedom for each node*/
  int NumberGP; /* Number of Gauss points inside of the element */
  int * N_id; /* Global index of the nodes (connectivity) */
  double ** X_g; /* Global coordiantes of the element nodes */
  double * (* n)(Vector * ); /* Shape function */
  double ** (* dn)(Vector * ); /* Derivative shape function */
  GaussPoint * GP_e; /* List of GP inside of the element */
} Element;

