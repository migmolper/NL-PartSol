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
char * SimulationFile;
char * RestartFile;
char * FEM_MeshFileName;
char * MPM_MeshFileName;
char OutputDir[MAXC];

/***************************************/
/*********** Kind of analysis **********/
/***************************************/
char * TimeIntegrationScheme;
char * ShapeFunctionGP;
char * Formulation;

#define NumberDimensions 2

int NumberDOF;

/***************************************/
/********* Time integration ************/
/***************************************/
double CFL; /* Courant number (0-1) */
double DeltaTimeStep;
double SpectralRadius;
int NumTimeStep;
int ResultsTimeStep;


/***************************************/
/********* Numeric Tolerances **********/
/***************************************/
double gamma_LME;
double TOL_lambda;
#define TOL_InOut 10E-10
#define TOL_NR 10E-6
#define TOL_zero 10E-10

/***************************************/
/******** Boundary conditions **********/
/***************************************/
char Field[10];
int * Node;
double * Value;


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
