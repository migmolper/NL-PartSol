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
}

/*********************************************************************/
