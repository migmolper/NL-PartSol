//  Non-Linear Particle Solver
//
//
//  Main author:    Miguel Molinos Perez
//

/*! \file nl-partsol.h
    \brief File that controls the nl-partsol
*/


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
#include "Globals.h"
#include "Types.h"
#include "Matlib.h"
#include "Fields.h"
#include "Constitutive.h"
#include "MPM.h"
#include "ShapeFun.h"
#include "InOutFun.h"
#include "Formulations.h"
