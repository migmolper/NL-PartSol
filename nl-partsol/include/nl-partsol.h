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
#ifdef __linux__
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>


#elif __APPLE__
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stddef.h>

#elif _WIN32
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#endif

/***************************************/
/******** nl-partsol libraries *********/
/***************************************/
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "InOutFun.h"