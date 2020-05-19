#include "nl-partsol.h"


Curve ReadCurve(char * Name_File)
{
  /*!
    Define load curve 
  */
  Curve DatCurve;
  
  /*!
    Define simulation file 
  */
  FILE * Sim_dat;

  /*!
    Auxiliar variable for reading the lines in the files 
  */
  char line[MAXC] = {0};
  
  /*!
    Number of element in the line , just for check 
  */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /*!
    Axiliar variable for defining default curves 
  */
  int Tc, T0, T1;

  /*!
    Open and check .load file 
  */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL)
    {
      fprintf(stderr,"%s : %s %s \n",
	      "Error in ReadCurve()","during the lecture of",Name_File);
      exit(EXIT_FAILURE);
    }

  while(fgets(line, sizeof line, Sim_dat) != NULL)
    {    
      /*!
	Read the line with the space as separators 
      */
      nkwords = parse (kwords, line," \n\t");

      /*!
	Read general data of the curve 
      */
      if (strcmp(kwords[0],"DAT_CURVE") == 0)
	{
	  /*!
	    Loop over the words 
	  */
	  for(int i  = 1 ; i<nkwords ; i++)
	    { 
	      nparam = parse (param,kwords[i],"#\n");
	      if(nparam == 2)
		{
		  /*!
		    Number of values, it should be the same 
		    as the number of timesteps 
		  */
		  if(strcmp(param[0],"NUM") == 0)
		    {
		      DatCurve.Num = atoi(param[1]);
		    }
		}
	    }
	  /*!
	    ALLOCATE FX MATRIX 
	  */
	  DatCurve.Fx =
	    (double *)Allocate_ArrayZ(DatCurve.Num,sizeof(double));
	}
      
      /*!
	Fill the values of the custom curve 
      */
      if (strcmp(kwords[0],"CUSTOM_CURVE") == 0)
	{
	  for(int i = 0 ; i<DatCurve.Num ; i++)
	    {
	      fgets(line, sizeof line, Sim_dat);
	      nkwords = parse (kwords, line," \n\t");
	      DatCurve.Fx[i] = atof(kwords[0]);
	    }
	}
        
      /*!
	Fill curve constant
      */
      else if (strcmp(kwords[0],"CONSTANT_CURVE") == 0)
	{
	  fill_ConstantCurve(DatCurve,kwords[1]);
	}

      /*!
	Fill Ramp curve 
      */
      else if (strcmp(kwords[0],"RAMP_CURVE") == 0)
	{
	  fill_RampCurve(DatCurve,kwords[1]);
	}

      /*!
	Fill the Heaviside curve 
      */
      else if (strcmp(kwords[0],"HEAVISIDE_CURVE") == 0)
      	{
      	  fill_HeavisideCurve(DatCurve,kwords[1],kwords[2]);
      	}

      /*!
	Fill the delta curve
       */
      else if (strcmp(kwords[0],"DELTA_CURVE") == 0)
	{
	  fill_DeltaCurve(DatCurve,kwords[1],kwords[2]);
	}

      /*!
	Fill the hat curve
      */
      else if (strcmp(kwords[0],"HAT_CURVE") == 0)
	{
	  fill_HatCurve(DatCurve,kwords[1],kwords[2],kwords[3]);
	}

      /* /\*! */
      /* 	Error */
      /* *\/ */
      /* else */
      /* 	{ */
      /* 	  fprintf(stderr,"%s : %s %s \n", */
      /* 		  "Error in ReadCurve()","Invalid keyword",kwords[0]); */
      /* 	  exit(EXIT_FAILURE); */
      /* 	} */
    
    }
  
  fclose(Sim_dat);

  return DatCurve;
}

/*********************************************************************/

void fill_ConstantCurve(Curve DatCurve,
			char * Properties)
{
  /*! 
    Number of key words
  */
  int nparam;
  
  /*! 
    String to store the properties of the curve
  */
  char * param[MAXW] = {NULL};

  /*
    Scale of the constant curve
  */
  double Scale;
  
  /*!
    Scale parameter of the curve 
  */
  nparam = parse (param,Properties,"#\n");

  /*!
    Get the value of the scale 
  */
  if(strcmp(param[0],"SCALE") == 0)
    {
      Scale = atof(param[1]);

      /*!
	Fill the values of the curve 
      */
      for(int i = 0 ; i<DatCurve.Num ; i++)
	{
	  DatCurve.Fx[i] += Scale;
	}
    }

  /*!
    Fail 
  */
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in fill_ConstantCurve()","Wrong parameters");
      exit(EXIT_FAILURE);

    }
  
}

/*********************************************************************/

void fill_RampCurve(Curve DatCurve,
		    char * Properties){

  /*! 
    Number of key words
  */
  int nparam;

  /*! 
    String to store the properties of the curve
  */
  char * param[MAXW] = {NULL};

  /*
    Scale of the ramp curve
  */
  double Scale;
  
  /*!
    Read parameter of the curve 
  */
  nparam = parse (param,Properties,"#\n");
  if(strcmp(param[0],"SCALE") == 0)
    {
      Scale = atof(param[1]);

      /*!
	Fill the values of the curve 
      */
      for(int i = 0 ; i<DatCurve.Num ; i++)
	{
	  DatCurve.Fx[i] = Scale*(double)i/DatCurve.Num;
	}
      
    }

  /*!
    Fail 
  */
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in fill_RampCurve()","Wrong parameters");
      exit(EXIT_FAILURE);

    }  
}


/*********************************************************************/

void fill_HeavisideCurve(Curve DatCurve,
			 char * Property_1,
			 char * Property_2)
{
  /*! 
    Number of key words
  */
  int nparam_1;
  int nparam_2;

  /*! 
    String to store the properties of the curve
  */
  char * param_1[MAXW] = {NULL};
  char * param_2[MAXW] = {NULL};

  /*
    Scale of the heaviside curve
  */
  double Scale;

  /*
    Critical time of the heaviside curve
  */
  double Tc;
  
  /* Read parameter of the curve */
  nparam_1 = parse (param_1,Property_1,"#\n");
  nparam_2 = parse (param_2,Property_2,"#\n");
  
  if((strcmp(param_1[0],"SCALE") == 0) &&
     (strcmp(param_2[0],"Tc") == 0))
    {
      /*! 
	Get the scale parameter
      */
      Scale = atof(param_1[1]);

      /*!
	Get the critical time parameter
      */
      Tc = atoi(param_2[1]);

      /* Check the critical time */
      if((Tc > DatCurve.Num) || (Tc < 0))
	{	  
	  fprintf(stderr,"%s : %s !!! \n",
		  "Error in fill_HeavisideCurve()","Tc is not in [0,Sizecurve]");
	  exit(EXIT_FAILURE); 
	}

      /*!
	Fill the values of the curve 
      */
      for(int i = 0 ; i< DatCurve.Num; i++)
	{
	  if(i<=Tc)
	    {
	      DatCurve.Fx[i] += 0;
	    }
	  else
	    {
	      DatCurve.Fx[i] += Scale;
	    }
	}
 
    }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in fill_HeavisideCurve()","Wrong parameters");
    exit(EXIT_FAILURE);
  }

}

/*********************************************************************/

void fill_DeltaCurve(Curve DatCurve,
		     char * Property_1,
		     char * Property_2)
{

  /*! 
    Number of key words
  */
  int nparam_1;
  int nparam_2;

  /*! 
    String to store the properties of the curve
  */
  char * param_1[MAXW] = {NULL};
  char * param_2[MAXW] = {NULL};

  /*!
    Scale of the heaviside curve
  */
  double Scale;

  /*!
    Critical time of the heaviside curve
  */
  double Tc;

  /*!
    Read parameter of the curve 
  */
  nparam_1 = parse (param_1,Property_1,"#\n");
  nparam_2 = parse (param_2,Property_2,"#\n");

  /*!
    Read parameters and fill the curve
   */
  if((strcmp(param_1[0],"SCALE") == 0) &&
     (strcmp(param_2[0],"Tc") == 0))
    {
      /*! 
	Get the scale parameter
      */
      Scale = atof(param_1[1]);

      /*!
	Get the critical time parameter
      */
      Tc = atoi(param_2[1]);

      /* Check the critical time */
      if((Tc > DatCurve.Num) || (Tc < 0))
	{	  
	  fprintf(stderr,"%s : %s !!! \n",
		  "Error in fill_DeltaCurve()","Tc is not in [0,Sizecurve]");
	  exit(EXIT_FAILURE);	  
	}

      /*!
	Fill the values of the curve 
      */
      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++)
	{
	  if(i == Tc){
	    DatCurve.Fx[i] += Scale;
	  }
	  else
	    {
	      DatCurve.Fx[i] += 0;
	    }
	} 
    }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in fill_DeltaCurve()","Wrong parameters");
    exit(EXIT_FAILURE);
  }
      
  
}

/*********************************************************************/


void fill_HatCurve(Curve DatCurve,
		   char * Property_1,
		   char * Property_2,
		   char * Property_3)
{

  /*! 
    Number of key words
  */
  int nparam_1;
  int nparam_2;
  int nparam_3;

  /*! 
    String to store the properties of the curve
  */
  char * param_1[MAXW] = {NULL};
  char * param_2[MAXW] = {NULL};
  char * param_3[MAXW] = {NULL};

  /*
    Scale of the heaviside curve
  */
  double Scale;

  /*!
    Critical time of the heaviside curve
  */
  double T0,T1;

  /*!
    Read parameter of the curve 
  */
  nparam_1 = parse (param_1,Property_1,"#\n");
  nparam_2 = parse (param_2,Property_2,"#\n");
  nparam_3 = parse (param_3,Property_3,"#\n");
  
  
  /*!
    Read parameters and fill the curve
  */
  if((strcmp(param_1[0],"SCALE") == 0) &&
     ((strcmp(param_2[0],"T0") == 0) && (nparam_2 == 2)) &&
     ((strcmp(param_3[0],"T1") == 0) && (nparam_2 == 2)))
    {

      /*!
	Scale parameter of the curve 
      */
      Scale = atof(param_1[1]);

      /*!
	Read critical times 
      */
      T0 = atoi(param_2[1]);
      T1 = atoi(param_3[1]);

      /*!
	Check the time T0 and T1 of the hat 
      */
      if(T0 > T1)
	{
	  fprintf(stderr,"%s : %s !!! \n",
		  "Error in fill_HatCurve()","T0 > T1");
	  exit(EXIT_FAILURE);
	}
      if(T1 > DatCurve.Num)
	{
	  fprintf(stderr,"%s : %s !!! \n",
		  "Error in fill_HatCurve()","T1 > DatCurve.Num");
	  exit(EXIT_FAILURE);
	}
      if(T0 < 0)
	{
	  fprintf(stderr,"%s : %s !!! \n",
		  "Error in fill_HatCurve()","T0 < 0");
	  exit(EXIT_FAILURE);
	}

      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++)
	{
	  if(i < T0)
	    {
	      DatCurve.Fx[i] += 0;
	    }
	  else if ((i >= T0) || (i <= T1))
	    {
	      DatCurve.Fx[i] += Scale;
	    }
	  else
	    {
	      DatCurve.Fx[i] += 0;
	    }
	}
    
    }
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in fill_HatCurve()","Wrong parameters");
      exit(EXIT_FAILURE);
    }
  
}

/*********************************************************************/

void free_Curve(Curve A){
  free(A.Fx);  
}

/*********************************************************************/
