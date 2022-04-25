#ifndef _READ_GID_MESH_H_
#define _READ_GID_MESH_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"


#include "Nodes/H8.h"
#include "Nodes/Q4.h"
#include "Nodes/T3.h"
#include "Nodes/T4.h"

// Global libs
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"

Mesh     ReadGidMesh__MeshTools__(char *);
/********************************************************************/
#endif