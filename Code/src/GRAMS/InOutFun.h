#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef Matlib
#define Matlib
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* Auxiliar functions */
int parse(char **, char *, char *);
Matrix Read_CSV(char *, int);
Curve ReadCurve(char *);
Mesh ReadGidMesh(char *);
ChainPtr File2Chain(char *);

/* Read .gfd format */
Mesh GramsBox(char *);
Boundaries GramsBoundary(char *,int);
GaussPoint GramsSolid2D(char *,Mesh);
void GramsTime(char * );
Material * GramsMaterials(char *, GaussPoint);
void GramsInitials(char *, GaussPoint);
void GramsShapeFun(char * );
void GramsOutputs(char * );
Load * GramsNeumannBC(char *, int);
Load * GramsBodyForces(char *,GaussPoint);

/* Print .vtk format */
void WriteVtk_MPM(char *, GaussPoint, Matrix, int);
void WriteVtk_FEM(char *, Mesh, Matrix, int);
void WriteVtk_Float_Scalar(char *, Matrix);
void WriteVtk_Float_Vector(char *, Matrix);
void WriteVtk_Float_Tensor(char *, Matrix);

/* Print Gnuplot */
void WriteGnuplot(Matrix, Matrix, double, double, int, int, char[20]);
