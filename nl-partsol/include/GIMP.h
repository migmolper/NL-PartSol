/*! \file GIMP.h
  \brief Generalized interpolated material point shape functions

  Shape functions based in : 
  "" The Generalized Interpolation Material Point Method ""
  by S.G.Bardenhagen and E.M.Kober, 2004

            ^           
          __|__         
        _/  |  \_       
      _/    |    \_     
   __/      |      \__  
  --o-------o-------o-- 
   (-1)    (0)     (1)  
*/

#ifndef _GIMP_H_
#define _GIMP_H_

void     initialize__GIMP__(Particle, Mesh);
Matrix   N__GIMP__(Matrix, Matrix, double);
Matrix   dN__GIMP__(Matrix, Matrix, double);
ChainPtr tributary__GIMP__(Matrix, int,Matrix, Mesh);

#endif
