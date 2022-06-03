#include "nl-partsol.h"

// Global variables
int ResultsTimeStep;

/*********************************************************************/

void print_Status(char *Message, int Time) { puts(Message); }

/*********************************************************************/

void print_step(int Time, int NumTimeStep, double DeltaTimeStep) {
  printf("*************************************************************"
         "***********************\n");
  fprintf(stdout,
          "" GREEN "Step" RESET ": [%i/%i] | " GREEN "DeltaT" RESET
          ": %1.2e \n",
          Time, NumTimeStep, DeltaTimeStep);
}

/*********************************************************************/

void DoProgress(char label[], int step, int total) {
  // progress width
  const int pwidth = 72;

  // minus label len
  int width = pwidth - strlen(label);
  int pos = (step * width) / total;

  int percent = (step * 100) / total;

  // set green text color, only on Windows
#ifdef _WIN32
  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_GREEN);
#endif
  printf("%s[", label);

  // fill progress bar with =
  for (int i = 0; i < pos; i++)
    printf("%c", '=');

  // fill progress bar with spaces
  printf("% *c", width - pos + 1, ']');
  printf(" %3d%%\r", percent);

  // reset text color, only on Windows
#ifdef _WIN32
  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 0x08);
#endif
}

/*********************************************************************/

void print_convergence_stats(int Time, int NumTimeStep, int Iter, int MaxIter,
                             double Error0, double Error_total,
                             double Error_relative) {

  if (NumTimeStep < 10) {
    fprintf(stdout, "" GREEN "Step" RESET ": [%01d/%01d] | ", Time,
            NumTimeStep);
  } else if (NumTimeStep < 100) {
    fprintf(stdout, "" GREEN "Step" RESET ": [%02d/%02d] | ", Time,
            NumTimeStep);
  } else if (NumTimeStep < 1000) {
    fprintf(stdout, "" GREEN "Step" RESET ": [%03d/%03d] | ", Time,
            NumTimeStep);
  } else if (NumTimeStep < 10000) {
    fprintf(stdout, "" GREEN "Step" RESET ": [%04d/%04d] | ", Time,
            NumTimeStep);
  } else if (NumTimeStep < 100000) {
    fprintf(stdout, "" GREEN "Step" RESET ": [%05d/%05d] | ", Time,
            NumTimeStep);
  } else if (NumTimeStep < 1000000) {
    fprintf(stdout, "" GREEN "Step" RESET ": [%i/%i] | ", Time, NumTimeStep);
  }

  if (MaxIter < 10) {
    fprintf(stdout, "" GREEN "Iterations" RESET ": [%01d/%01d] | ", Iter,
            MaxIter);
  } else if (MaxIter < 100) {
    fprintf(stdout, "" GREEN "Iterations" RESET ": [%02d/%02d] | ", Iter,
            MaxIter);
  }

  fprintf(stdout,
          "" GREEN "Initial E" RESET ": %1.2e | " GREEN "Total E" RESET
          ": %1.2e | " GREEN "Relative E" RESET ": %1.2e \n",
          Error0, Error_total, Error_relative);

  FILE *Stats_Solver;
  char Name_file_t[10000];
  sprintf(Name_file_t, "%s/Stats_Solver.csv", OutputDir);
  Stats_Solver = fopen(Name_file_t, "a");

  if (Time == 0) {
    fprintf(Stats_Solver, "%s,%s,%s\n", "Iter", "Total Error",
            "Relative Error");
  }
  fprintf(Stats_Solver, "%i,%1.4e,%1.4e\n", Iter, Error_total, Error_relative);

  fclose(Stats_Solver);

  if (Flag_Print_Convergence) {
    if (Time == 1) {
      FILE *gnuplot;
      gnuplot = popen("gnuplot -persistent", "w");
      fprintf(gnuplot, "set datafile separator ','\n");
      fprintf(gnuplot, "set logscale y 10\n");
      fprintf(gnuplot, "set key autotitle columnhead\n");
      fprintf(gnuplot, "set xlabel 'Steps'\n");
      fprintf(gnuplot, "set y2tics\n");
      fprintf(gnuplot, "set ytics nomirror\n");
      fprintf(gnuplot, "set y2label 'Number of iterations'\n");
      fprintf(gnuplot, "set y2range [0:10]\n");
      fprintf(gnuplot, "while (1) {\n");
      fprintf(gnuplot,
              "\t plot './%s/Stats_Solver.csv' using 2 with lines , "
              "'./%s/Stats_Solver.csv' using 3 with lines , "
              "'./%s/Stats_Solver.csv' using 1 with lines axis x1y2 lt -1 \n",
              OutputDir, OutputDir, OutputDir);
      fprintf(gnuplot, "\t pause 1\n");
      fprintf(gnuplot, "}\n");
      fflush(gnuplot);
    }
  }
}

/*********************************************************************/
