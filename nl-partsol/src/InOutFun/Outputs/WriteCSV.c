#include "nl-partsol.h"

/*
  Call global variables
*/

/*
  Auxiliar functions
*/
static void standard_error(char *);
static void print_nodes_scalar_variable_to_csv(FILE *, Matrix, Matrix, Mask,
                                               Event);
static void print_nodes_vectorial_variable_to_csv(FILE *, Matrix, Matrix, Mask,
                                                  Event);
static void print_nodes_tensorial_variable_to_csv(FILE *, Matrix, Matrix, Mask,
                                                  Event);
static void print_particles_scalar_variable_to_csv(FILE *, Matrix, Matrix,
                                                   Event);
static void print_particles_vectorial_variable_to_csv(FILE *, Matrix, Matrix,
                                                      Event);
static void print_particles_tensorial_variable_to_csv(FILE *, Matrix, Matrix,
                                                      Event);
static void print_Gauss_Point_scalar_variable_to_csv(FILE *, double);
static void print_Gauss_Point_vectorial_variable_to_csv(FILE *, Tensor);
static void print_Gauss_Point_tensorial_variable_to_csv(FILE *, Tensor);
/*****************************************************************/

void path_nodes_analysis_csv__InOutFun__(Matrix Nodal_Field,
                                         Matrix Nodes_Coordinates,
                                         char *Name_File, Mask ActiveNodes,
                                         Event Output_Commands, int idx,
                                         int TimeStep, double DeltaTimeStep) {
  int Ndim = NumberDimensions;
  int NumCols = Nodal_Field.N_cols;
  int k_start = Output_Commands.i_start;
  int k_end = Output_Commands.i_end;
  int k_step = Output_Commands.i_step;
  char Error_message[MAXW];
  char Name_file_Ts[MAXC];
  FILE *csv_file;

  if (k_start <= TimeStep && TimeStep <= k_end && TimeStep % k_step == 0) {

    sprintf(Name_file_Ts, "%s/%s_idx_%i_step_%i.csv", Output_Commands.Directory,
            Name_File, idx, TimeStep);
    csv_file = fopen(Name_file_Ts, "w");

    /* Scalar variable */
    if (NumCols == 1) {
      print_nodes_scalar_variable_to_csv(csv_file, Nodal_Field,
                                         Nodes_Coordinates, ActiveNodes,
                                         Output_Commands);
    }
    /* Vectorial variable */
    else if (NumCols == Ndim) {
      print_nodes_vectorial_variable_to_csv(csv_file, Nodal_Field,
                                            Nodes_Coordinates, ActiveNodes,
                                            Output_Commands);
    }
    /* Tensorial variable */
    else if (NumCols == Ndim * Ndim) {
      print_nodes_tensorial_variable_to_csv(csv_file, Nodal_Field,
                                            Nodes_Coordinates, ActiveNodes,
                                            Output_Commands);
    }

    /* Close the file */
    fclose(csv_file);
  }
}

/*****************************************************************/

void path_particles_analysis_csv__InOutFun__(
    Matrix Particle_Field, Matrix Particles_Coordinates, char *Name_File,
    Event Output_Commands, int idx, int TimeStep, double DeltaTimeStep) {
  int Ndim = NumberDimensions;
  int NumCols = Particle_Field.N_cols;
  int k_start = Output_Commands.i_start;
  int k_end = Output_Commands.i_end;
  int k_step = Output_Commands.i_step;
  char Error_message[MAXW];
  char Name_file_Ts[MAXC];
  FILE *csv_file;

  if (k_start <= TimeStep && TimeStep <= k_end && TimeStep % k_step == 0) {

    sprintf(Name_file_Ts, "%s/%s_idx_%i_step_%i.csv", Output_Commands.Directory,
            Name_File, idx, TimeStep);
    csv_file = fopen(Name_file_Ts, "w");

    /* Scalar variable */
    if (NumCols == 1) {
      print_particles_scalar_variable_to_csv(
          csv_file, Particle_Field, Particles_Coordinates, Output_Commands);
    }
    /* Vectorial variable */
    else if (NumCols == Ndim) {
      print_particles_vectorial_variable_to_csv(
          csv_file, Particle_Field, Particles_Coordinates, Output_Commands);
    }
    /* Tensorial variable */
    else if (NumCols == Ndim * Ndim) {
      print_particles_tensorial_variable_to_csv(
          csv_file, Particle_Field, Particles_Coordinates, Output_Commands);
    }

    /* Close the file */
    fclose(csv_file);
  }
}

/*****************************************************************/

void Gauss_Point_evolution__InOutFun__(Particle PointAnalysis,
                                       Event Output_Commands, char *Name_File,
                                       int k, int idx) {
  int Ndim = NumberDimensions;

  char Name_file_Ts[MAXC];
  FILE *csv_file;

  sprintf(Name_file_Ts, "%s/%s_idx_%i.csv", Output_Commands.Directory,
          Name_File, idx);
  csv_file = fopen(Name_file_Ts, "a+");

  if (Output_Commands.Out_csv_Gauss_Point_evolution_Stress) {
    Tensor Stress_k =
        memory_to_tensor__TensorLib__(PointAnalysis.Phi.Stress.nM[k], 2);
    print_Gauss_Point_tensorial_variable_to_csv(csv_file, Stress_k);
  }
  if (Output_Commands.Out_csv_Gauss_Point_evolution_Strain) {
    Tensor F_k = memory_to_tensor__TensorLib__(PointAnalysis.Phi.F_n.nM[k], 2);
    Tensor C_k = right_Cauchy_Green__Particles__(F_k);
    Tensor E_k = strain_Green_Lagrange__Particles__(C_k);
    print_Gauss_Point_tensorial_variable_to_csv(csv_file, E_k);

    free__TensorLib__(C_k);
    free__TensorLib__(E_k);
  }
  if (Output_Commands.Out_csv_Gauss_Point_evolution_Deformation_gradient) {
    Tensor F_k = memory_to_tensor__TensorLib__(PointAnalysis.Phi.F_n.nM[k], 2);
    print_Gauss_Point_tensorial_variable_to_csv(csv_file, F_k);
  }
  if (Output_Commands.Out_csv_Gauss_Point_evolution_EPS) {
    double EPS_k = PointAnalysis.Phi.Equiv_Plast_Str.nV[k];
    print_Gauss_Point_scalar_variable_to_csv(csv_file, EPS_k);
  }
  if (Output_Commands.Out_csv_Gauss_Point_evolution_Cohesion) {
    double cohesion_k = PointAnalysis.Phi.cohesion.nV[k];
    print_Gauss_Point_scalar_variable_to_csv(csv_file, cohesion_k);
  }

  /* New line */
  fprintf(csv_file, "\n");

  /* Close the file */
  fclose(csv_file);
}

/*****************************************************************/

static void standard_error(char *Error_message) {
  fprintf(stderr, "%s : %s !!! \n", "Error in __InOutFun__()", Error_message);
  exit(EXIT_FAILURE);
}

/*****************************************************************/

static void print_nodes_scalar_variable_to_csv(FILE *csv_file, Matrix Variable,
                                               Matrix Coordinates,
                                               Mask ActiveNodes,
                                               Event Output_Commands) {
  int Idx;
  int Idx_mask;

  fprintf(csv_file, "x,y,z,var \n");

  for (int i = 0; i < Output_Commands.Lenght_Path; i++) {
    Idx = Output_Commands.Idx_Path[i];
    Idx_mask = ActiveNodes.Nodes2Mask[Idx];

    if (Idx_mask != -1) {
      fprintf(csv_file, "%lf,%lf,%lf,%lf \n", Coordinates.nM[Idx_mask][0],
              Coordinates.nM[Idx_mask][1], Coordinates.nM[Idx_mask][2],
              Variable.nV[Idx_mask]);
    }
  }
}

/*****************************************************************/

static void print_nodes_vectorial_variable_to_csv(FILE *csv_file,
                                                  Matrix Variable,
                                                  Matrix Coordinates,
                                                  Mask ActiveNodes,
                                                  Event Output_Commands) {
  int Idx;
  int Idx_mask;
  int Ndim = NumberDimensions;

  if (Ndim == 2) {
    fprintf(csv_file, "x,y,z,var.x,var.y \n");
  } else if (Ndim == 3) {
    fprintf(csv_file, "x,y,z,var.x,var.y,var.z \n");
  }

  for (int i = 0; i < Output_Commands.Lenght_Path; i++) {
    Idx = Output_Commands.Idx_Path[i];
    Idx_mask = ActiveNodes.Nodes2Mask[Idx];

    if (Idx_mask != -1) {
      if (Ndim == 2) {
        fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf \n", Coordinates.nM[Idx_mask][0],
                Coordinates.nM[Idx_mask][1], Coordinates.nM[Idx_mask][2],
                Variable.nM[Idx_mask][0], Variable.nM[Idx_mask][1]);
      } else if (Ndim == 3) {
        fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf \n",
                Coordinates.nM[Idx_mask][0], Coordinates.nM[Idx_mask][1],
                Coordinates.nM[Idx_mask][2], Variable.nM[Idx_mask][0],
                Variable.nM[Idx_mask][1], Variable.nM[Idx_mask][2]);
      }
    }
  }
}

/*****************************************************************/

static void print_nodes_tensorial_variable_to_csv(FILE *csv_file,
                                                  Matrix Variable,
                                                  Matrix Coordinates,
                                                  Mask ActiveNodes,
                                                  Event Output_Commands) {
  int Idx;
  int Idx_mask;
  int Ndim = NumberDimensions;

  if (Ndim == 2) {
    fprintf(csv_file, "x,y,z,var.xx,var.xy,var.yx,var.yy \n");
  } else if (Ndim == 3) {
    fprintf(csv_file, "x,y,z,var.xx,var.xy,var.xz,var.yx,var.yy,var.yz,var.zx,"
                      "var.zy,var.zz \n");
  }

  for (int i = 0; i < Output_Commands.Lenght_Path; i++) {
    Idx = Output_Commands.Idx_Path[i];
    Idx_mask = ActiveNodes.Nodes2Mask[Idx];

    if (Idx_mask != -1) {
      if (Ndim == 2) {
        fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf \n",
                Coordinates.nM[Idx_mask][0], Coordinates.nM[Idx_mask][1],
                Coordinates.nM[Idx_mask][2], Variable.nM[Idx_mask][0],
                Variable.nM[Idx_mask][1], Variable.nM[Idx_mask][2],
                Variable.nM[Idx_mask][3]);
      } else if (Ndim == 3) {
        fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf \n",
                Coordinates.nM[Idx_mask][0], Coordinates.nM[Idx_mask][1],
                Coordinates.nM[Idx_mask][2], Variable.nM[Idx_mask][0],
                Variable.nM[Idx_mask][1], Variable.nM[Idx_mask][2],
                Variable.nM[Idx_mask][3], Variable.nM[Idx_mask][4],
                Variable.nM[Idx_mask][5], Variable.nM[Idx_mask][6],
                Variable.nM[Idx_mask][7], Variable.nM[Idx_mask][8]);
      }
    }
  }
}

/*****************************************************************/

static void print_particles_scalar_variable_to_csv(FILE *csv_file,
                                                   Matrix Variable,
                                                   Matrix Coordinates,
                                                   Event Output_Commands) {
  int Idx;
  int Ndim = NumberDimensions;

  fprintf(csv_file, "x,y,z,var \n");

  for (int i = 0; i < Output_Commands.Lenght_Path; i++) {
    Idx = Output_Commands.Idx_Path[i];

    fprintf(csv_file, "%lf,%lf,%lf,%lf \n", Coordinates.nM[Idx][0],
            Coordinates.nM[Idx][1], Coordinates.nM[Idx][2], Variable.nV[Idx]);
  }
}

/*****************************************************************/

static void print_particles_vectorial_variable_to_csv(FILE *csv_file,
                                                      Matrix Variable,
                                                      Matrix Coordinates,
                                                      Event Output_Commands) {
  int Idx;
  int Ndim = NumberDimensions;

  if (Ndim == 2) {
    fprintf(csv_file, "x,y,z,var.x,var.y \n");
  } else if (Ndim == 3) {
    fprintf(csv_file, "x,y,z,var.x,var.y,var.z \n");
  }

  for (int i = 0; i < Output_Commands.Lenght_Path; i++) {
    Idx = Output_Commands.Idx_Path[i];

    if (Ndim == 2) {
      fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf \n", Coordinates.nM[Idx][0],
              Coordinates.nM[Idx][1], Coordinates.nM[Idx][2],
              Variable.nM[Idx][0], Variable.nM[Idx][1]);
    } else if (Ndim == 3) {
      fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf \n", Coordinates.nM[Idx][0],
              Coordinates.nM[Idx][1], Coordinates.nM[Idx][2],
              Variable.nM[Idx][0], Variable.nM[Idx][1], Variable.nM[Idx][2]);
    }
  }
}

/*****************************************************************/

static void print_particles_tensorial_variable_to_csv(FILE *csv_file,
                                                      Matrix Variable,
                                                      Matrix Coordinates,
                                                      Event Output_Commands) {
  int Idx;
  int Ndim = NumberDimensions;

  if (Ndim == 2) {
    fprintf(csv_file, "x,y,z,var.xx,var.xy,var.yx,var.yy \n");
  } else if (Ndim == 3) {
    fprintf(csv_file, "x,y,z,var.xx,var.xy,var.xz,var.yx,var.yy,var.yz,var.zx,"
                      "var.zy,var.zz \n");
  }

  for (int i = 0; i < Output_Commands.Lenght_Path; i++) {
    Idx = Output_Commands.Idx_Path[i];

    if (Ndim == 2) {
      fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf \n",
              Coordinates.nM[Idx][0], Coordinates.nM[Idx][1],
              Coordinates.nM[Idx][2], Variable.nM[Idx][0], Variable.nM[Idx][1],
              Variable.nM[Idx][2], Variable.nM[Idx][3]);
    } else if (Ndim == 3) {
      fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf \n",
              Coordinates.nM[Idx][0], Coordinates.nM[Idx][1],
              Coordinates.nM[Idx][2], Variable.nM[Idx][0], Variable.nM[Idx][1],
              Variable.nM[Idx][2], Variable.nM[Idx][3], Variable.nM[Idx][4],
              Variable.nM[Idx][5], Variable.nM[Idx][6], Variable.nM[Idx][7],
              Variable.nM[Idx][8]);
    }
  }
}

/*****************************************************************/

static void print_Gauss_Point_scalar_variable_to_csv(FILE *csv_file,
                                                     double Variable) {
  fprintf(csv_file, "%lf,", Variable);
}

/*****************************************************************/

static void print_Gauss_Point_vectorial_variable_to_csv(FILE *csv_file,
                                                        Tensor Variable) {
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  fprintf(csv_file, "%lf,%lf,", Variable.n[0], Variable.n[1]);
#endif
#if NumberDimensions == 3
  fprintf(csv_file, "%lf,%lf,%lf,", Variable.n[0], Variable.n[1],
          Variable.n[2]);
#endif
}

/*****************************************************************/

static void print_Gauss_Point_tensorial_variable_to_csv(FILE *csv_file,
                                                        Tensor Variable) {
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  fprintf(csv_file, "%lf,%lf,%lf,%lf,", Variable.N[0][0], Variable.N[0][1],
          Variable.N[1][0], Variable.N[1][1]);
#endif
#if NumberDimensions == 3
  fprintf(csv_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,", Variable.N[0][0],
          Variable.N[0][1], Variable.N[0][2], Variable.N[1][0],
          Variable.N[1][1], Variable.N[1][2], Variable.N[2][0],
          Variable.N[2][1], Variable.N[2][2]);
#endif
}

/*****************************************************************/
