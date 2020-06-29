//  Non-Linear Particle Solver
//
//
//  Main author:    Miguel Molinos Perez
//

/*! \file nl-partsol.h
    \brief File that controls the nl-partsol
*/

/***************************************/
/************* Word Parser *************/
/***************************************/
#define MAXW 100
#define MAXC 1000

/***************************************/
/******** Name of inputs files *********/
/***************************************/
extern char * SimulationFile;
extern char * RestartFile;
extern char * FEM_MeshFileName;
extern char * MPM_MeshFileName;
extern char OutputDir[MAXC];

/***************************************/
/*********** Kind of analysis **********/
/***************************************/
extern char * TimeIntegrationScheme;
extern char * ShapeFunctionGP;
extern char * Formulation;

#define NumberDimensions 2

extern int NumberDOF;

/***************************************/
/********* Time integration ************/
/***************************************/
extern double CFL; /* Courant number (0-1) */
extern double DeltaTimeStep;
extern double SpectralRadius;
extern int NumTimeStep;
extern int ResultsTimeStep;


/***************************************/
/********* Numeric Tolerances **********/
/***************************************/
extern double gamma_LME;
extern double TOL_lambda;
#define TOL_InOut 10E-23
#define TOL_NR 10E-6
#define TOL_zero 10E-10

/***************************************/
/******** Boundary conditions **********/
/***************************************/
extern char Field[10];
extern int * Node;
extern double * Value;


/***************************************/
/********** External libraries *********/
/***************************************/
#ifdef linux
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>
#endif

#ifdef _WIN32
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>
#endif

/***************************************/
/********** GRAMS's libraries **********/
/***************************************/
#include "Variables.h"
#include "Matlib.h"
#include "Fields.h"
#include "Constitutive.h"
#include "MPM.h"
#include "ShapeFun.h"
#include "InOutFun.h"
#include "Formulations.h"
