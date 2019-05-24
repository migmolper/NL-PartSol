GiDPost 2.7
===========

gidpost is a set of functions (library) for writing postprocess results for GiD in ASCII or binary compressed format.

This software is copyrighted by CIMNE www.gidhome.com. The software can be used freely under the terms described in license.terms, all terms described there apply to all files associated with the software unless explicitly disclaimed in individual files. Particular terms apply to the third party code, "cfortran.h", which has its own distribution policy (please read the "cfortran.doc" for this code).
This description asumes that the readear is familiar with the postprocess terminology.
For futher details please check the online help available in GiD ( Postprocess data files chapter).

The library was implemented taking into account two of the must widely used development environments: C/C++ and FORTRAN.

**Look into the *doc* folder for the library manual**

ChangeLog:
==========

*From version 2.6 to 2.7*

Added support for MeshGroups ( BeginMeshGroup(), End..., BeginOnMeshGroup(), End...) in HDF5 as used in GiD, useful for dynamic meshes, refinement, multi-stage simulation, etc.
All non implemented functions for HDF5 now return -1, i.e. all GiD_fXXX() functions.

*From version 2.5 to 2.6*

Added support for ComplexMatrix Result types.
More exhaustive example and some corrections.

*From version 2.4 to 2.5*

Default format to print real numbers changed from "%g" to "%.9g" increasing the amount of significant digits to be printed in ASCII.
New funcion GiD_PostSetFormatReal to allow change the default format of real numbers.

*From version 2.3 to 2.4*

Allow write OnNurbsSurface results for isogeometrical analysis

*From version 2.1 to 2.3*

 Allow use HD5 from FORTRAN interface
 Fixed bug with mesh group

*From version 2.0 to 2.1*

Allow write complex scalar and complex vector results
Write units of mesh or results
Enhanced FORTRAN 77 and 90 interface to be more compatible.

*From version 1.7.2 to 2.0*

New set of functions "GiD_fxxx" that allow specify the file to write, in case of have multiple files of mesh or results.
Library rewritten in C avoiding C++ to be more compatible linking with other languages.
Optional define of HDF5 to enable write postprocess files using the HDF5 library

*Old changes*

2008/09/30
	* Updated cfortran to release 4.4

2008/02/25
	* changing to release 1.8
	* M build/Jamroot, source/gidpost.cpp source/gidpostInt.h: removed	C++ dependencies.
2008/01/29
	* source/gidpost.cpp:

	 - GiD_OpenPostMeshFile should not write a header. Thanks to	Janosch Stascheit <janosch.stascheit@rub.de>

	- invalid access in GiD_BeginMeshGroup: should ensure that a valid	mesh file object is available so etMeshFile() should be called
	before writting. Thanks to Janosch Stascheit	<janosch.stascheit@rub.de> for providing the fix.

2008/01/18
	* source/gidpost.cpp: GiD_BeginOnMeshGroup should be "OnGroup	\"%s\"" instead of "Group \"%s\"". Thanks to Janosch Stascheit
	<janosch.stascheit@rub.de>
2007/06/13
	* source/gidpost.cpp: GiD_WritePlainDefMatrix can be written in	binary files.

2007/06/11
	* GiD_Sphere and GiD_Circle
2007/06/11
	* source/gidpost.cc: Write2DMatrix is now available on binary	format, Z component is set to 0.

2007/01/23:
	* source/gidpostfor.c: bug reported by a user:	GiD_Begin3DMatResult must be ncomp = SetupComponents(6, not 4	Factoring fortran declaration in gidpostfor.h
	Defining symbols for both gcc and icc.

2006/06/25
	* doc/*:
	* sources/gidpost.cpp,gidpost.h: included pyramid element.

2006/06/25
	* tag 1.7: rel_1_7
	* examples/testpost.c: fixed a bug, AsciiZipped must has the mesh
	file in ascii format.

2006/05/25
	* gidpost.h: dllexport.

2006/05/18 
	* all: GiD_BeginMeshColor, GiD_BeginMeshGroup, GiD_EndMeshGroup,
	GiD_BeginOnMeshGroup, GiD_EndOnMeshGroup

2005/10/24
	* doc:
	* gidpost.h : documented GiD_ResultDescriptionDim
2005/10/21
	* gidpost* : added support to new type specification for result
	groups. The type could be Type:dim, where Type is one valid result
	type and dim is a number giving a valid dimension. Updated fortran
	interface for recently added functions.

2005/10/07
	* all: new version 1.60. "Result Groups" can now be written. The
	macro USE_CONST can be used to enable (if defined) or desable (if
	not defined) the use of const char* arguments.

2005/09/22
	* all: ha sido un error cambiar de 'char *' a 'const char *'

2005/09/22
	* gidpost.cpp: flag_begin_values = 0; must be set in
	GiD_BeginResultGroup. 
2005/09/21
	* gidpost.h:
	* gidpost.cpp: values_location_{min,max} not needed, CBufferValue
	* gidpostInt.cpp: now controls the types of results written within 
	* gidpostInt.h: a ResultGroup.
2005/09/16
	* gidpost.cpp: in GiD_WriteVectorModule, when writing if we are
	inside a ResultGroup block the module value is ignored.
	* gidpostInt.cpp: in CBufferValues::WriteValues, buffer overflow
	is checked after checking for change in location.

2005/09/15
	* gidpost.cpp:
	* gidpostInt.cpp:
	* gidpostInt.h: fixed a bug when writing results in a "Result
	Group", vectors can be 3 or 4 component long.

	* examples/testpost.c : bug when passing location id in result
	group sample code. 
2005/06/27
	* source/gidpost.{h,cpp} including Prism element
	* doc/gidpost.{html,subst} including Prism element documentation.
	* all: change to version 1.52
2005/05/09
	* all: constification, version change from 1.5 to 1.51

2005/01/04
	* source/gidpostInt.cpp: fixed a bug when writing 3D vectors.

2005/01/07
	* source/gidpost.cpp added const char and a filter to remove
	double quotes from names which can cause problems within gid. 

2003/07/29
	* doc/gidpost.html   : commented the GiD_ResultUnit
	* source/gidpost.h   : interface because it is not 
	* source/gidpostfor.c: supported yet inside GiD.
2003/07/28
	* doc/gidpost.html: updated documentation for release 1.5
	* binary/gidpost.lib: binary release for windows
2003/07/15
	* gidpost.h:
	* gidpost.cpp:
	* gidpostfor.c: new function 'GiD_WriteCoordinates2D' (to write
	coordinates in 2D but only in ASCII format, the function can be
	used in binary but the library provide the 'z' cordinate a zero.

2003/07/15
	* gidpost.h:
	* gidpost.cpp: removed 'char * UnitName' argument , new functions:
	 GiD_BeginResultHeader, GiD_ResultRange, GiD_ResultComponents,
	 GiD_ResultUnit, GiD_BeginResultGroup, GiD_ResultDescription,
	 GiD_ResultValues, GiD_FlushPostFile. Validation in debug mode.

	* gidpostInt.h:
	* gidpostInt.cpp: new member
	CPostFile::WriteValues( int id, int n, double * buffer), new class
	CBufferValues to write the values in a result group, validation in
	debug mode.


	* gidpostfor.c: updated fortran interface for the new
	functionality.

	* testpost.c:
	* testpostfor.f: added test code for the new functionality.