#include "nl-partsol.h"

// Global variables 
int ResultsTimeStep;

/*********************************************************************/

void print_Status(char * Message,int Time)
{
  if(Time%ResultsTimeStep == 0)
    {
      puts(Message);
    }
}

/*********************************************************************/

void print_step(int Time,double DeltaTimeStep)
{
  if(Time%ResultsTimeStep == 0)
    {
      printf("************* STEP : %i , DeltaT : %e \n",
	     Time,DeltaTimeStep);
    }
}

/*********************************************************************/

void print_convergence_stats(int Time, int Iter, double Error0, double Error_total, double Error_relative)
{
  if(Time%ResultsTimeStep == 0)
    {
      printf("\t Convergence reached after : %i interations \n",Iter);
      printf("\t Initial error : %1.4e \n",Error0);
      printf("\t Total error : %1.4e \n",Error_total);
      printf("\t Relative error : %1.4e \n",Error_relative);
    }

    FILE * Stats_Solver;
    char Name_file_t[80];
    sprintf(Name_file_t,"%s/Stats_Solver.csv",OutputDir);  
    Stats_Solver = fopen(Name_file_t,"a");

    if(Time == 0)
    {
      fprintf(Stats_Solver,"%s,%s,%s\n","Iter","Total Error","Relative Error");
    }
    fprintf(Stats_Solver,"%i,%1.4e,%1.4e\n",Iter,Error_total,Error_relative);

    fclose(Stats_Solver);
}

/*********************************************************************/
