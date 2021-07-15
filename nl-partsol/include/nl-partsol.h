//  Non-Linear Particle Solver
//
//
//  Main author:    Miguel Molinos Perez
//

/*! \file nl-partsol.h
    \brief File that controls the nl-partsol
*/


#define MAXW 100
#define MAXC 1000
#define NumberDimensions 2
#define TOL_InOut 10E-23
#define TOL_NR 10E-6
#define TOL_zero 10E-23

/***************************************/
/********** External libraries *********/
/***************************************/
#ifdef __linux__
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>

#elif __APPLE__
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>

#elif _WIN32
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
/******** nl-partsol libraries *********/
/***************************************/
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Fields.h"
#include "Constitutive.h"
#include "Particles.h"
#include "T3.h"
#include "Q4.h"
#include "T4.h"
#include "H8.h"
#include "GIMP.h"
#include "LME.h"
#include "Nodes.h"
#include "InOutFun.h"
#include "Solvers.h"
