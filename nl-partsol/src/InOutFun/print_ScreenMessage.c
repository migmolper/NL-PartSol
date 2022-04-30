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

void print_convergence_stats(int Time, int NumTimeStep, int Iter, double Error0,
                             double Error_total, double Error_relative) {

  fprintf(stdout,
          "" GREEN "Step" RESET ": [%04d/%04d] | " GREEN "Iterations" RESET ": [%02d/%02d] | " GREEN "Initial E" RESET
          ": %1.2e | " GREEN "Total E" RESET ": %1.2e | " GREEN
          "Relative E" RESET ": %1.2e \n",
          Time, NumTimeStep, Iter, 10, Error0, Error_total, Error_relative);

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
}

/*********************************************************************/
