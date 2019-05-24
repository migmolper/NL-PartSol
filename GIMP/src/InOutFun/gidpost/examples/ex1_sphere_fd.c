#include <stdlib.h>
#include "gidpost.h"

/*
  #define ASCII_FORMAT
*/

int main()
{
  GiD_FILE fdm, fdr;

  GiD_PostInit();
  
#if defined(ASCII_FORMAT)
  fdm = GiD_fOpenPostMeshFile( "ex1_sphere_fd.post.msh", GiD_PostAscii);
  fdr = GiD_fOpenPostResultFile( "ex1_sphere_fd.post.res", GiD_PostAscii);
#else
  fdm = fdr = GiD_fOpenPostResultFile( "ex1_sphere_fd.post.bin", GiD_PostBinary);
#endif
  
  GiD_fBeginMesh(fdm, "esferitas", GiD_3D, GiD_Sphere, 1);
  {
    GiD_fBeginCoordinates(fdm);
    {
      GiD_fWriteCoordinates(fdm, 1, 0, 0, 0); 
      GiD_fWriteCoordinates(fdm, 2, 100, 0, 0); 
    }
    GiD_fEndCoordinates(fdm);
    
    GiD_fBeginElements(fdm);
    {
      GiD_fWriteSphere(fdm, 1, 1, 1.0);
      GiD_fWriteSphere(fdm, 2, 2, 10.0);
    }
    GiD_fEndElements(fdm);
  }
  GiD_fEndMesh(fdm);
  
  GiD_fBeginMesh(fdm, "circulitos", GiD_3D, GiD_Circle, 1);
  {
    GiD_fBeginCoordinates(fdm);
    {
      GiD_fWriteCoordinates(fdm, 3, 30, 0, 0); 
      GiD_fWriteCoordinates(fdm, 4, 70, 0, 0 ); 
    }
    GiD_fEndCoordinates(fdm);
    
    GiD_fBeginElements(fdm);
    {
      GiD_fWriteCircle(fdm, 11, 3, 3.0, 1.0, 0.0, 0.0);
      GiD_fWriteCircle(fdm, 12, 4, 2.0, 0.77, 0.77, 0.0);
    }
    GiD_fEndElements(fdm);
  }
  GiD_fEndMesh(fdm);

#ifdef ASCII_FORMAT
  GiD_fClosePostMeshFile(fdm);
#endif
  GiD_fClosePostResultFile(fdr);

  GiD_PostDone();
  
  return 0;
}
