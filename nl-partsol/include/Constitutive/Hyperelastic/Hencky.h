/**
 * @file Hencky.h
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _HENCKY_CONSTITUTIVE_H_
#define _HENCKY_CONSTITUTIVE_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
#include "Particles.h"
// clang-format on

/**
 * @brief 
 * 
 * @param IO_State 
 * @param MatProp 
 * @return int 
 */
int compute_Kirchhoff_Stress_Hencky(State_Parameters IO_State,
                                    Material MatProp);

#endif
