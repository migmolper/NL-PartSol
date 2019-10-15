#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"


double uGIMP_1D(double L, double lp, double Xp, double Xi){

  double S_Ip;
  
  if (fabs(Xp-Xi) >= L+lp){
    S_Ip = 0;
  }
  else if (-L-lp < Xp-Xi <= -L+lp){
    S_Ip = (1/(4*L*lp))*(L+lp+Xp-Xi)**2;
  }
  else if (-L+lp < Xp-Xi <= -lp){
    S_Ip = 1 + (Xp-Xi)/(L);
  }
  else if (-lp < Xp-Xi <= lp){
    S_Ip = 1 - ((Xp - Xi)**2 + (lp)**2)/(2*L*lp);
  }
  else if (lp < Xp-Xi <= L-lp){
    S_Ip = 1 - (Xp-Xi)/(L);
  }
  else if (L-lp < Xp-Xi <= L+lp){
    S_Ip = (1/(4*L*lp))*(L+lp-Xp+Xi)**2;
  }
  
  return S_Ip;
}

double uGIMP_2D(){
  
}
