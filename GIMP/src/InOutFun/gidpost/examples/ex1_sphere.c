#include <stdlib.h>
#include <math.h>
#include "gidpost.h"

/* 
   #define ASCII_FORMAT
*/

int main()
{
#if defined(ASCII_FORMAT)
  GiD_OpenPostMeshFile( "ex1_sphere.post.msh", GiD_PostAscii);
  GiD_OpenPostResultFile( "ex1_sphere.post.res", GiD_PostAscii);
#else
  GiD_OpenPostResultFile( "ex1_sphere.post.bin", GiD_PostBinary);
#endif
  
  GiD_BeginMesh("esferitas", GiD_3D, GiD_Sphere, 1);
  {
    GiD_MeshUnit( "m");
    GiD_BeginCoordinates();
    {
      GiD_WriteCoordinates( 1, 0, 0, 0 ); 
      GiD_WriteCoordinates( 2, 100, 0, 0 ); 
    }
    GiD_EndCoordinates();
    
    GiD_BeginElements();
    {
      GiD_WriteSphere(1, 1, 1.0);
      GiD_WriteSphere(2, 2, 10.0);
    }
    GiD_EndElements();
  }
  GiD_EndMesh();
  
  GiD_BeginMesh("circulitos", GiD_3D, GiD_Circle, 1);
  {
    GiD_MeshUnit( "m");
    GiD_BeginCoordinates();
    {
      GiD_WriteCoordinates( 3, 30, 0, 0 ); 
      GiD_WriteCoordinates( 4, 70, 0, 0 ); 
    }
    GiD_EndCoordinates();
    
    GiD_BeginElements();
    {
      GiD_WriteCircle(11, 3, 3.0, 1.0, 0.0, 0.0);
      GiD_WriteCircle(12, 4, 2.0, 0.77, 0.77, 0.0);
    }
    GiD_EndElements();
  }
  GiD_EndMesh();

#ifdef ASCII_FORMAT
  GiD_ClosePostMeshFile();
#endif

  /* some dummy results */
  {
    double step = 1.0;
    for ( step = 1.0; step < 11.0; step += 1.0) {
      /* a scalar */
      GiD_BeginResultHeader( "a scalar", "test", step, GiD_Scalar, GiD_OnNodes, NULL);
      {
	int in = 0;
	GiD_ResultUnit( "m");
	for ( in = 0; in < 4; in++) {
	  GiD_WriteScalar( in + 1, ( double)in * step);
	}
      }
      GiD_EndResult();
      /* a vector */
      GiD_BeginResultHeader( "a vector", "test", step, GiD_Vector, GiD_OnNodes, NULL);
      {
	/*
	  const char *comp_names[] = { "sin()", "cos()", "sin() * cos()"};
	*/
	int in = 0;
	const char *comp_names[] = { "in*step", "2in*step", "x+y"};
	GiD_ResultUnit( "m");
	GiD_ResultComponents( 3, comp_names);
	for ( in = 0; in < 4; in++) {
	  /*
	    double ang = ( 1.0 + in * step) / 3.14159265355897932384626;
	    double sa = sin( ang);
	    double ca = cos( ang);
	  */
	  double sa = in * step;
	  double ca = in * 2 * step;
	  GiD_WriteVector( in + 1, sa, ca, sa + ca);
	}
      }
      GiD_EndResult();

      /* complex scalar */
      GiD_BeginResultHeader( "complex scalar", "test", step, GiD_ComplexScalar, GiD_OnNodes, NULL);
      {
	int in = 0;
	const char *comp_names[] = { "scalar-real", "scalar-imag"};
	GiD_ResultUnit( "m");
	GiD_ResultComponents( 2, comp_names);
	for ( in = 0; in < 4; in++) {
	  GiD_WriteComplexScalar( in + 1, ( double)in * step, ( double)( in * step + 0.5));
	}
      }
      GiD_EndResult();
      /* comples vector */
      GiD_BeginResultHeader( "complex vector", "test", step, GiD_ComplexVector, GiD_OnNodes, NULL);
      {
	/* 
	   const char *comp_names[] = { "sin().r", "sin().i", "cos().r", "cos().i", "sin*cos.r", "sin*cos.i"};
	*/
	int in = 0;
	const char *comp_names[] = { "in*step", "in*step+0.5", "2in*step", "2in*step+0.5", "x_r+y_r", "x_i+y_i"};
	GiD_ResultUnit( "m");
	GiD_ResultComponents( 6, comp_names);
	for ( in = 0; in < 4; in++) {
	  /*
	    double ang_r = ( 1.0 + in * step) / 3.14159265355897932384626;
	    double sa_r = sin( ang_r);
	    double ca_r = cos( ang_r);
	    double ang_i = ( 1.0 - in * step) / 3.14159265355897932384626;
	    double sa_i = sin( ang_i);
	    double ca_i = cos( ang_i);
	  */
	  double sa_r = in * step;
	  double sa_i = in * step + 0.5;
	  double ca_r = in * 2 * step;
	  double ca_i = in * 2 * step + 0.5;
	  GiD_WriteComplexVector( in + 1, sa_r, sa_i, ca_r, ca_i, sa_r + ca_r, sa_i + ca_i);
	}
      }
      GiD_EndResult();
    }
  }

  GiD_ClosePostResultFile();
  
  return 0;
}
