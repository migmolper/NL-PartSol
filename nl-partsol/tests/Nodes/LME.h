#ifndef _LME_TEST_H_
#define _LME_TEST_H_

#include <check.h>
#include <petscksp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// clang-format off
#include "Nodes/LME.h"
// clang-format on

START_TEST(test_LME_2D_shape_function) {

  const double Tolerance = 1.E-8;
  PetscErrorCode ierr;
  PetscInt Ndim = NumberDimensions;
  PetscInt N_a = 16;
  PetscScalar X_p[NumberDimensions] = {0.0, 0.0};
  PetscScalar lambda_tr[NumberDimensions] = {0.0, 0.0};
  PetscScalar h_avg = 2.0;
  PetscScalar Gamma = 20.0;
  PetscScalar Beta = Gamma / (h_avg * h_avg);
  PetscScalar X_a[32] = {-1., -1., 1.,  -1., 1.,  1.,  -1., 1., -3., -3., -1.,
                         -3., 1.,  -3., 3.,  -3., -3., -1., 3., -1., -3., 1.,
                         3.,  1.,  -3., 3.,  -1., 3.,  1.,  3., 3.,  3.};
  PetscScalar sf_ans[16] = {0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00,
                            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

  PetscScalar gradsf_ans[32] = {
      -0.25, -0.25, 0.25, -0.25, 0.25, 0.25, -0.25, 0.25, 0,  0, -0,
      0,     0,     0,    0,     0,    0,    -0,    0,    -0, 0, 0,
      0,     0,     0,    0,     -0,   0,    0,     0,    0,  0};

  PetscScalar l_a[32]; 
  for (PetscInt a = 0; a < N_a; a++) {
    for (PetscInt i = 0; i < Ndim; i++) {
      l_a[a * Ndim + i] = X_p[i] - X_a[a * Ndim + i];
    }
  }

  ierr = lambda__LME__(l_a, lambda_tr, Beta, N_a);

  PetscScalar *p_a = p__LME__(l_a, lambda_tr, Beta, N_a);
  PetscScalar *dp_a = dp__LME__(l_a, lambda_tr, Beta, N_a);
  PetscScalar PartitionUnity = 0;

  for (PetscInt a = 0; a < N_a; a++) {
    PartitionUnity += p_a[a];
  }
  
  ck_assert_double_eq_tol(PartitionUnity, 1.0, Tolerance);

  for (PetscInt a = 0; a < N_a; a++) {
    ck_assert_double_eq_tol(p_a[a], sf_ans[a], Tolerance);
    for (PetscInt i = 0; i < Ndim; i++) {
      ck_assert_double_eq_tol(gradsf_ans[a * Ndim + i], dp_a[a * Ndim + i],
                              Tolerance);
    }
  }

  PetscFree(p_a);
  PetscFree(dp_a);
}
END_TEST

Suite *LME_suite(void) {

  Suite *s;
  TCase *tc_core;

  s = suite_create("LME-2D-testing");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_LME_2D_shape_function);
  suite_add_tcase(s, tc_core);

  return s;
}

#endif