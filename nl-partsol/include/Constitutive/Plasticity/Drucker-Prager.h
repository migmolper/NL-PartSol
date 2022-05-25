/**
 * @file Drucker-Prager.h
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _DRUCKER_PRAGER_CONSTITUTIVE_H_
#define _DRUCKER_PRAGER_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/**
 * @brief Compute the behaviour of a Henky-Hyperelastic material \n
 *  using the Drucker-Prager yield criterium.
 * 
 * @param IO_State State parameters of the particle
 * @param MatProp Material properties of the model
 * @return int STATUS
 */
int compute_1PK_Drucker_Prager(
    State_Parameters IO_State, 
    Material MatProp);

#endif
