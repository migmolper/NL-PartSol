
void * AllocateArray(int, /* Dimension of the array */
		     int); /* Number of bytes for each element */

void ** AllocateMatrix(int, /* Number of rows */
		       int, /* Number of columns */
		       int); /* Number of bytes for each element */

void ** AllocateTableOfPointers(int, /* Dimension of the table */
				int);  /* Number of bytes for each element */

double CalculeDeterminant(double *, /* Input matrix linearized */
			  int); /* Shape of the matrix */
