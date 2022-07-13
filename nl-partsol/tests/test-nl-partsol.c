#include <check.h>
#include <petscksp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// clang-format off
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "InOutFun.h"
#include "Nodes/LME.h"
// clang-format on

#include "../tests/Nodes/LME.h"

/*
  Call global variables
*/
char ShapeFunctionGP[MAXC];
char SimulationFile[MAXC];
char Static_conditons[MAXC];
char Formulation[MAXC];
char *TimeIntegrationScheme;
bool Flag_Print_Convergence;
Load gravity_field;
bool Driver_EigenErosion;
bool Driver_EigenSoftening;

int main(int argc, char *argv[]) {

  PetscInitialize(&argc, &argv, 0, 0);

  int number_failed;
  Suite *s;
  SRunner *sr;

  s = LME_suite();
  sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  PetscFinalize();

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

  return 0;
}