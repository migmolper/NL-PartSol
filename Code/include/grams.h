/*

  Library with all the structures adopted in the code to extend the
  basic functionalities of C with them : 

  -> Matrix : Usefull structure to deal with algebraic operations.
  -> Chain and ChainPtr : Chain of nodes and pointer to a chain of nodes.
  -> Table : Table of integers.
  -> Curve : Structure created for dealing with complex boundary
  conditions and loads.
  -> BoundaryConditions : Structure that store all the information
  of a boundary condition and allow to me and others programers to
  programe/imposse easily new ones.
  -> Load : Structure that the store all the information about 
  one kind of load over the GPs*.
  -> LoadCase : Structure created to join multiple kind of loads in
  one single structure to deal with complex load cases.
  -> Fields : Structure that store all the physical variables of
  the MPM** simulation related with the GPs*.
  -> GaussPoint : Structure created for store general properties of
  a set of GPs*, this allows to generate several groups of them, with
  diferent properties. Ie : water, soil, air.
  -> Mesh : Structure created for store general properties of a FEM*** 
  mesh, this allows to generate differents mesh for diferents sets of 
  GPs*.

  Note* : GPs -> Gauss Points
  Note** : MPM -> Material Point Method
  Note*** : FEM -> Finite Element Method

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
char * FEM_MeshFileName;
char * MPM_MeshFileName;
char OutputDir[MAXC];

/***************************************/
/*********** Kind of analysis **********/
/***************************************/
char * TimeIntegration;
char * ShapeFunctionGP;
char * Formulation;
int NumberDimensions;
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
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>

/***************************************/
/********** GRAMS's libraries **********/
/***************************************/
#include "Matlib.h"
#include "Domain.h"
#include "Constitutive.h"
#include "MPM.h"
#include "ShapeFun.h"
#include "InOutFun.h"
#include "Formulations.h"
