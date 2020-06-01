#include "nl-partsol.h"

/*********************************************************************/

void OutputSimulation(GaussPoint Set_Particles,
		      int Numerical_T,
		      double Physical_T,
		      double DeltaT,
		      Event Parameters)
{

  if((Numerical_T >= Parameters.k_start) &&
     (Physical_T  >= Parameters.start)   &&
     (Numerical_T <= Parameters.k_end)   &&
     (Physical_T  <= Parameters.end))
    {
      if(Physical_T%Parameters.step <= DeltaT)
	{
	  WriteVtk_MPM(Parameters.File,Set_Particles,"ALL",
		       (int)Physical_T/Parameters.step);
	}
      else if(Numerical_T%Parameters.k_step == 0)
	{	
	  WriteVtk_MPM(Parameters.File,Set_Particles,"ALL",
		       (int)Numerical_T/Parameters.k_step);  
	}
    }

}

/*********************************************************************/
