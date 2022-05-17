/**
 * @file U-pw-Newmark-beta.h
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _UPW_NEWMARK_BETA_H_
#define _UPW_NEWMARK_BETA_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Global libs
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"

// Shape functions and auxilar tools
#include "Nodes/Nodes-Tools.h"
#include "Nodes/Shape-Functions.h"

// Courant
#include "Formulations/Courant.h"

// Material lib
#include "Constitutive/Constitutive.h"

// Functions to evaluate residual/jacobian
#include "Formulations/Displacements-WaterPressure/U-pw-Functions.h"

// 
#include "InOutFun.h"

//  Call global variables
double Thickness_Plain_Stress;
Event *Out_nodal_path_csv;
Event *Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

// Gravity field 
Load gravity_field;

// Global variuables
unsigned InitialStep;
unsigned TimeStep;
unsigned NumTimeStep;


/**
 * @brief u-pw formulation with one single set of material points. 
 * Implicit solver for finite strains based in the  
 * Newmark-beta, see \cite Molinos_et_al_2021_CMAME.
 * 
 * @param FEM_Mesh Variable with the nodal information
 * @param MPM_Mesh  Variable with the particle information
 * @param Parameters_Solver 
 * @return int 
 */
int Newmark_beta__upw__(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver);

#endif