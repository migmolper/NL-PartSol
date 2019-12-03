#ifndef _GIDPOSTFILES_H_
#define _GIDPOSTFILES_H_

#include "gidpost.h"
#include "gidpostInt.h"

CPostFile *_GiDfiles_GetMeshFile();
CPostFile * _GiDfiles_NewFile(GiD_PostMode Mode);
int _GiDfiles_CheckState( post_state s_req, CPostFile *file );

int _GiDfiles_BeginMesh( CPostFile *File,
                         GP_CONST char * MeshName, GiD_Dimension Dim,
                         GiD_ElementType EType, int NNode );
int _GiDfiles_BeginMeshColor(CPostFile *File,
                             GP_CONST char * MeshName, GiD_Dimension Dim,
                             GiD_ElementType EType, int NNode,
                             double Red, double Green, double Blue);
int _GiDfiles_MeshUnit(CPostFile *File,GP_CONST char * UnitName);
int _GiDfiles_BeginCoordinates( CPostFile *File );
int _GiDfiles_EndCoordinates( CPostFile *File );
int _GiDfiles_BeginMeshGroup( CPostFile *File, GP_CONST char* Name );
int _GiDfiles_EndMeshGroup( CPostFile *fileMesh );
int _GiDfiles_WriteCoordinates(CPostFile *File, int id, 
                               double x, double y, double z);
int _GiDfiles_WriteCoordinates2D(CPostFile *File, int id, double x, double y);
int _GiDfiles_BeginElements( CPostFile *File );
int _GiDfiles_EndElements(CPostFile *File);
int _GiDfiles_WriteElement(CPostFile *File, int id, int nid[]);
int _GiDfiles_WriteElementMat(CPostFile *File, int id, int nid[]);
int _GiDfiles_WriteSphere(CPostFile *File, int id, int nid, double r);
int _GiDfiles_WriteSphereMat(CPostFile * File, int id, int nid, double r, int mat);
int _GiDfiles_WriteCircle(CPostFile *File,
                          int id, int nid, double r,
                          double nx, double ny, double nz);
int _GiDfiles_WriteCircleMat(CPostFile *File, int id, int nid, double r,
                             double nx, double ny, double nz, int mat);
int _GiDfiles_BeginGaussPoint(CPostFile *File,
                              GP_CONST char * name, GiD_ElementType EType,
                              GP_CONST char * MeshName,
                              int GP_number, int NodesIncluded, int InternalCoord);
int _GiDfiles_EndGaussPoint( CPostFile *File );
int _GiDfiles_WriteGaussPoint2D( CPostFile *File, double x, double y );
int _GiDfiles_WriteGaussPoint3D(CPostFile *File, double x, double y, double z);
int _GiDfiles_BeginRangeTable( CPostFile *File, GP_CONST char * name );
int _GiDfiles_EndRangeTable( CPostFile *File );
int _GiDfiles_WriteMinRange( CPostFile* File, double max, GP_CONST char * name );
int _GiDfiles_WriteRange( CPostFile* File,
                          double min, double max, GP_CONST char * name );
int _GiDfiles_WriteMaxRange(CPostFile *File, double min, GP_CONST char * name);
int _GiDfiles_BeginResult(CPostFile *File,
                          GP_CONST char     *Result,
                          GP_CONST char     *Analysis,
                          double             step,
                          GiD_ResultType     Type,
                          GiD_ResultLocation Where,
                          GP_CONST char     *GaussPointsName,
                          GP_CONST char     *RangeTable, 
                          int                compc,
                          GP_CONST char      *compv[]);
int _GiDfiles_BeginResultHeader( CPostFile         *File,
                                 GP_CONST char     *Result,
                                 GP_CONST char     *Analysis,
                                 double             step,
                                 GiD_ResultType     Type,
                                 GiD_ResultLocation Where,
                                 GP_CONST char     *GaussPointsName);
int _GiDfiles_ResultRange( CPostFile *File, GP_CONST char * RangeTable );
int _GiDfiles_ResultComponents( CPostFile *File, int compc, GP_CONST char * compv[] );
int _GiDfiles_ResultUnit( CPostFile *File, GP_CONST char * UnitName );
int _GiDfiles_ResultUserDefined( CPostFile *File, GP_CONST char * Name,GP_CONST char * Value );
int _GiDfiles_BeginResultGroup( CPostFile         *File,
                                GP_CONST char     *Analysis,
                                double             step,
                                GiD_ResultLocation Where,
                                GP_CONST char     *GaussPointsName );
int _GiDfiles_ResultDescription_( CPostFile     *File,
                                  GP_CONST char *Result,
                                  GiD_ResultType Type,
                                  size_t         s );
int _GiDfiles_EndResult(CPostFile *File);
int _GiDfiles_BeginOnMeshGroup( CPostFile *File, char *Name );
int _GiDfiles_EndOnMeshGroup( CPostFile *File );
int _GiDfiles_WriteScalar(CPostFile *File, int id, double v );
int _GiDfiles_Write2DVector(CPostFile *File, int id, double x, double y);
int _GiDfiles_WriteVector( CPostFile *File, int id, double x, double y, double z );
int _GiDfiles_WriteVectorModule( CPostFile *File,
                                 int id, double x, double y, double z, double mod );
int _GiDfiles_Write2DMatrix(CPostFile *File,
                            int id, double Sxx, double Syy, double Sxy);
int _GiDfiles_Write3DMatrix(CPostFile *File,
                            int id, double Sxx, double Syy, double Szz,
                            double Sxy, double Syz, double Sxz);
int _GiDfiles_WritePlainDefMatrix(CPostFile *File, int id,
                                  double Sxx, double Syy, double Sxy, double Szz );
int _GiDfiles_WriteMainMatrix(CPostFile *File, int id,
                              double Si, double Sii, double Siii,
                              double Vix, double Viy, double Viz,
                              double Viix, double Viiy, double Viiz,
                              double Viiix, double Viiiy, double Viiiz);
int _GiDfiles_WriteLocalAxes(CPostFile *File,
                             int id, double euler_1, double euler_2, double euler_3);
int _GiDfiles_WriteComplexScalar( CPostFile *File,
                                  int id, double complex_real, double complex_imag);
int _GiDfiles_Write2DComplexVector( CPostFile *File, int id,
                                    double x_real, double x_imag,
                                    double y_real, double y_imag);
int _GiDfiles_WriteComplexVector( CPostFile *File, int id,
                                  double x_real, double x_imag,
                                  double y_real, double y_imag,
                                  double z_real, double z_imag);
int _GiDfiles_Write2DComplexMatrix(CPostFile *File, int id,
                                   double Sxx_real, double Syy_real, double Sxy_real,
                                   double Sxx_imag, double Syy_imag, double Sxy_imag);
int _GiDfiles_WriteComplexMatrix(CPostFile *File, int id,
                                 double Sxx_real, double Syy_real, double Szz_real,
                                 double Sxy_real, double Syz_real, double Sxz_real,
                                 double Sxx_imag, double Syy_imag, double Szz_imag,
                                 double Sxy_imag, double Syz_imag, double Sxz_imag);
int _GiDfiles_WriteNurbsSurface( CPostFile *File, int id, int n, double* v );
int _GiDfiles_WriteNurbsSurfaceVector( CPostFile *File, int id, int n, int num_comp, double* v );

#endif /* _GIDPOSTFILES_H_ */