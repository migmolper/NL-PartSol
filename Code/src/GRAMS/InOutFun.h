#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef Matlib
#define Matlib
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* GnuPlotOutput.c */
void WriteGnuplot(Matrix, Matrix, double, double, int, int, char[20]);

/* Parse.c */
int parse(char **words, char *str, char *delims);

/* ReadCSV.c */
Matrix Read_CSV(char *, int);

/* ReadCurve.c */
Curve ReadCurve(char *);

/* ReadGidMesh.c */
Mesh ReadGidMesh(char *);

/* ReadDatFile.c */
void Read_GeneralParameters(char *);
Boundaries Set_FEM_BCC(char *, Mesh);
LoadCase Read_MPM_LoadCase_ExtForces(char *,GaussPoint);
LoadCase Read_MPM_LoadCase_BodyForces(char *,GaussPoint);
void Read_MPM_InitVal(char *, GaussPoint);

/* GraMS Interface */
void GramsTime(char * );
Material * GramsMaterials(char *, GaussPoint);
void GramsShapeFun(char * );
void GramsOutputs(char * );
Mesh GramsBox(char *);

/* WriteVtk.c */
void WriteVtk_MPM(char *, GaussPoint, Matrix, int);
void WriteVtk_FEM(char *, Mesh, Matrix, int);
void WriteVtk_Float_Scalar(char *, Matrix);
void WriteVtk_Float_Vector(char *, Matrix);
void WriteVtk_Float_Tensor(char *, Matrix);

