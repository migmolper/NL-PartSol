#include "nl-partsol.h"

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
      printf("************* STEP : %i , DeltaT : %f \n",
	     Time,DeltaTimeStep);
    }
}

/*********************************************************************/
