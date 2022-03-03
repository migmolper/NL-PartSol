#include <string.h>
#include "nl-partsol.h"

/*
  Call global variables
*/
char Output_Backup_File[MAXC];
bool Backup_damage;
bool Backup_plastic_deformation_gradient;
bool Backup_EPS;

/*
  Auxiliar functions
*/
static void vtk_Out_mass(FILE *, Matrix, int);
static void vtk_Out_density(FILE *, Matrix, int);
static void vtk_Out_damage(FILE *, Matrix, int);
static void vtk_Out_nodal_idx(FILE *, int *, int);
static void vtk_Out_material_idx(FILE *, int *, int);
static void vtk_Out_vel(FILE *, Matrix, int);
static void vtk_Out_acc(FILE *, Matrix, int);
static void vtk_Out_dis(FILE *, Matrix, int);
static void vtk_Out_Pw(FILE *, Matrix, Matrix, int);
static void vtk_Out_dPw_dt(FILE *, Matrix, int);
static void vtk_Out_Deformation_Gradient(FILE *, Matrix, int);
static void vtk_Out_plastic_deformation_gradient(FILE *, Matrix, int);
static void vtk_Out_Equiv_Plastic_Strain(FILE *, Matrix, int);
/*********************************************************************/

void particle_backup_vtk__InOutFun__(Particle MPM_Mesh, int TimeStep_i,
                                     int Backup_TimeStep) {

  /* Number of dimensions */
  int Ndim = NumberDimensions;

  int NumParticles = MPM_Mesh.NumGP;

  FILE *Vtk_file;
  char Name_file_t[80];

  sprintf(Name_file_t, "%s/%s_%i.vtk", OutputDir, Output_Backup_File,
          TimeStep_i);

  Vtk_file = fopen(Name_file_t, "w");

  /* Header */
  fprintf(Vtk_file, "# vtk DataFile Version 3.0 \n");
  fprintf(Vtk_file, "Results time step %i \n", Backup_TimeStep);
  fprintf(Vtk_file, "ASCII \n");
  fprintf(Vtk_file, "DATASET UNSTRUCTURED_GRID \n");

  /* Coordinates */
  fprintf(Vtk_file, "POINTS %i double \n", MPM_Mesh.NumGP);
  for (int i = 0; i < MPM_Mesh.NumGP; i++) {
    for (int j = 0; j < 3; j++) {
      if (j < Ndim) {
        fprintf(Vtk_file, "%.20g ", MPM_Mesh.Phi.x_GC.nM[i][j]);
      } else {
        fprintf(Vtk_file, "%.20g ", 0.0);
      }
    }
    fprintf(Vtk_file, "\n");
  }

  /* Connectivity */
  fprintf(Vtk_file, "CELLS %i %i \n", MPM_Mesh.NumGP, 2 * MPM_Mesh.NumGP);
  for (int i = 0; i < MPM_Mesh.NumGP; i++) {
    fprintf(Vtk_file, "%i %i \n", 1, i);
  }

  /* Type of element */
  fprintf(Vtk_file, "CELL_TYPES %i \n", MPM_Mesh.NumGP);
  for (int i = 0; i < MPM_Mesh.NumGP; i++) {
    fprintf(Vtk_file, "%i \n", 1);
  }

  /* Point data */
  fprintf(Vtk_file, "POINT_DATA %i \n", MPM_Mesh.NumGP);

  /* Cell data */
  fprintf(Vtk_file, "CELL_DATA %i \n", MPM_Mesh.NumGP);

  /* Print particle velocity */
  vtk_Out_vel(Vtk_file, MPM_Mesh.Phi.vel, NumParticles);

  /* Print particle acceleration */
  vtk_Out_acc(Vtk_file, MPM_Mesh.Phi.acc, NumParticles);

  /* Print particle displacement */
  vtk_Out_dis(Vtk_file, MPM_Mesh.Phi.dis, NumParticles);

  /* Print particle deformation gradient */
  vtk_Out_Deformation_Gradient(Vtk_file, MPM_Mesh.Phi.F_n, NumParticles);

  /* Print particle mass */
  vtk_Out_mass(Vtk_file, MPM_Mesh.Phi.mass, NumParticles);

  /* Print particle density */
  vtk_Out_density(Vtk_file, MPM_Mesh.Phi.rho, NumParticles);

  /* Print particle nodal index */
  vtk_Out_nodal_idx(Vtk_file, MPM_Mesh.I0, NumParticles);

  /* Print particle material index */
  vtk_Out_material_idx(Vtk_file, MPM_Mesh.MatIdx, NumParticles);

  /* Print elastic left Cauchy-Green tensor */
  if (Backup_plastic_deformation_gradient) {
    vtk_Out_plastic_deformation_gradient(Vtk_file, MPM_Mesh.Phi.b_e_n,
                                         NumParticles);
  }

  /* Print particle damage */
  if (Backup_damage) {
    vtk_Out_damage(Vtk_file, MPM_Mesh.Phi.chi, NumParticles);
  }

  /* print equivalent plastic strain */
  if (Backup_EPS) {
    vtk_Out_Equiv_Plastic_Strain(Vtk_file, MPM_Mesh.Phi.Equiv_Plast_Str,
                                 NumParticles);
  }

  if (strcmp(Formulation, "-upw") == 0) {
    vtk_Out_Pw(Vtk_file, MPM_Mesh.Phi.Pw, MPM_Mesh.Phi.J_n1, NumParticles);

    vtk_Out_dPw_dt(Vtk_file, MPM_Mesh.Phi.d_Pw_dt_n1, NumParticles);
  }

  /* Close the file */
  fclose(Vtk_file);
}

/*********************************************************************/

static void vtk_Out_mass(FILE *Vtk_file, Matrix mass, int NumParticles) {
  fprintf(Vtk_file, "SCALARS MASS double \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%.20g \n", mass.nV[i]);
  }
}

/*********************************************************************/

static void vtk_Out_density(FILE *Vtk_file, Matrix rho, int NumParticles) {
  fprintf(Vtk_file, "SCALARS DENSITY double \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%.20g \n", rho.nV[i]);
  }
}

/*********************************************************************/

static void vtk_Out_damage(FILE *Vtk_file, Matrix chi, int NumParticles) {
  fprintf(Vtk_file, "SCALARS Damage double \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%.20g \n", chi.nV[i]);
  }
}

/*********************************************************************/

static void vtk_Out_nodal_idx(FILE *Vtk_file, int *I0, int NumParticles) {
  fprintf(Vtk_file, "SCALARS ELEM_i integer \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%i \n", I0[i]);
  }
}

/*********************************************************************/
static void vtk_Out_material_idx(FILE *Vtk_file, int *MatIdx,
                                 int NumParticles) {
  fprintf(Vtk_file, "SCALARS MatIdx integer \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%i \n", MatIdx[i]);
  }
}

/*********************************************************************/

static void vtk_Out_vel(FILE *Vtk_file, Matrix vel, int NumParticles) {
  int Ndim = NumberDimensions;

  fprintf(Vtk_file, "VECTORS VELOCITY double \n");
  for (int i = 0; i < NumParticles; i++) {
    for (int j = 0; j < 3; j++) {
      if (j < Ndim) {
        fprintf(Vtk_file, "%.20g ", vel.nM[i][j]);
      } else {
        fprintf(Vtk_file, "%.20g ", 0.0);
      }
    }
    fprintf(Vtk_file, "\n");
  }
}

/*********************************************************************/

static void vtk_Out_acc(FILE *Vtk_file, Matrix acc, int NumParticles) {
  int Ndim = NumberDimensions;

  fprintf(Vtk_file, "VECTORS ACCELERATION double \n");
  for (int i = 0; i < NumParticles; i++) {
    for (int j = 0; j < 3; j++) {
      if (j < Ndim) {
        fprintf(Vtk_file, "%.20g ", acc.nM[i][j]);
      } else {
        fprintf(Vtk_file, "%.20g ", 0.0);
      }
    }
    fprintf(Vtk_file, "\n");
  }
}

/*********************************************************************/

static void vtk_Out_dis(FILE *Vtk_file, Matrix dis, int NumParticles) {
  int Ndim = NumberDimensions;

  fprintf(Vtk_file, "VECTORS DISPLACEMENT double \n");
  for (int i = 0; i < NumParticles; i++) {
    for (int j = 0; j < 3; j++) {
      if (j < Ndim) {
        fprintf(Vtk_file, "%.20g ", dis.nM[i][j]);
      } else {
        fprintf(Vtk_file, "%.20g ", 0.0);
      }
    }
    fprintf(Vtk_file, "\n");
  }
}

/*********************************************************************/

static void vtk_Out_Pw(FILE *Vtk_file, Matrix Pw, Matrix J, int NumParticles) {
  int Ndim = NumberDimensions;
  Tensor Stress_p;
  double P_p; /* Trace of the stress tensor (volumetric) */
  fprintf(Vtk_file, "SCALARS Pw double \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%.20g \n", Pw.nV[i] / J.nV[i]);
  }
}

/*********************************************************************/

static void vtk_Out_dPw_dt(FILE *Vtk_file, Matrix d_Pw, int NumParticles) {
  int Ndim = NumberDimensions;
  Tensor Stress_p;
  double P_p; /* Trace of the stress tensor (volumetric) */
  fprintf(Vtk_file, "SCALARS dPw_dt double \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%.20g \n", d_Pw.nV[i]);
  }
}

/*********************************************************************/

static void vtk_Out_Deformation_Gradient(FILE *Vtk_file, Matrix F_n,
                                         int NumParticles) {
  int Ndim = NumberDimensions;

  fprintf(Vtk_file, "TENSORS DEFORMATION-GRADIENT double \n");
  for (int i = 0; i < NumParticles; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if ((j < Ndim) && (k < Ndim)) {
          fprintf(Vtk_file, "%.20g ", F_n.nM[i][j * Ndim + k]);
        } else {
          fprintf(Vtk_file, "%.20g ", 0.0);
        }
      }
      fprintf(Vtk_file, "\n");
    }
    fprintf(Vtk_file, "\n");
  }
}

/*********************************************************************/

static void vtk_Out_plastic_deformation_gradient(FILE *Vtk_file,
                                                 Matrix F_m1_plastic,
                                                 int NumParticles) {
  int Ndim = NumberDimensions;

  Matrix F_plastic = inverse__MatrixLib__(F_m1_plastic);

  fprintf(Vtk_file, "TENSORS PLASTIC-DEFORMATION-GRADIENT double \n");
  for (int i = 0; i < NumParticles; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if ((j < Ndim) && (k < Ndim)) {
          fprintf(Vtk_file, "%.20g ", F_plastic.nM[i][j * Ndim + k]);
        } else {
          fprintf(Vtk_file, "%.20g ", 0.0);
        }
      }
      fprintf(Vtk_file, "\n");
    }
    fprintf(Vtk_file, "\n");
  }

  free__MatrixLib__(F_plastic);
}

/*********************************************************************/

static void vtk_Out_Equiv_Plastic_Strain(FILE *Vtk_file, Matrix EPS,
                                         int NumParticles) {

  fprintf(Vtk_file, "SCALARS EPS double \n");
  fprintf(Vtk_file, "LOOKUP_TABLE default \n");
  for (int i = 0; i < NumParticles; i++) {
    fprintf(Vtk_file, "%.20g \n", EPS.nV[i]);
  }
}

/*********************************************************************/