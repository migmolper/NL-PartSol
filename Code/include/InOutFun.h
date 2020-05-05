#ifndef _INOUTFUN_H_
#define _INOUTFUN_H_


/* Auxiliar functions */
int parse(char **, char *, char *);
void generate_route(char *, char *);
Matrix Read_CSV(char *, int);
Curve ReadCurve(char *);
Mesh ReadGidMesh(char *);
ChainPtr File2Chain(char *);

/* Read .gfd format */
Mesh GramsBox(char *);
Boundaries GramsBoundary(char *,int);
GaussPoint GramsSolid2D(char *,Mesh);
void GramsTime(char * );
Material * GramsMaterials(char *, GaussPoint, int);
void GramsInitials(char *, GaussPoint, int);
void GramsShapeFun(char * );
void GramsOutputs(char * );
Load * GramsNeumannBC(char *, int, int);
Load * GramsBodyForces(char *, int, int);

/* Print .vtk format */
void WriteVtk_MPM(char *, GaussPoint, Matrix, int);
void WriteVtk_FEM(char *, Mesh, Matrix, int);
void WriteVtk_Float_Scalar(char *, Matrix);
void WriteVtk_Float_Vector(char *, Matrix);
void WriteVtk_Float_Tensor(char *, Matrix);

/* Print Gnuplot */
void WriteGnuplot(Matrix, Matrix, double, double, int, int, char[20]);


#endif
