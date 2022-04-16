
#ifndef _MACROS_H_
#define _MACROS_H_

/*
  Color text
*/
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


/*
  Constant macros
*/
#define MAXW 100
#define MAXC 1000
#define NumberDimensions 2
#define TOL_InOut 10E-23
#define TOL_NR 10E-6
#define TOL_zero 10E-23
#define PI__MatrixLib__ 3.14159265358979323846

/*
  Math macros
*/
static float sqr_arg;
#define SQR(a) ((sqr_arg=(a)) == 0.0 ? 0.0 : sqr_arg*sqr_arg)
static double dsqr_arg;
#define DSQR(a) ((dsqr_arg=(a)) == 0.0 ? 0.0 : dsqr_arg*dsqr_arg)
static double dmax_arg1, dmax_arg2;
#define DMAX(a,b) (dmax_arg1=(a),dmax_arg2=(b),(dmax_arg1) > (dmax_arg2) ? \
		   (dmax_arg1) : (dmax_arg2))
static double dmin_arg1, dmin_arg2;
#define DMIN(a,b) (dmin_arg1=(a),dmin_arg2=(b),(dmin_arg1) < (dmin_arg2) ? \
		   (dmin_arg1) : (dmin_arg2))
static float max_arg1, max_arg2;
#define FMAX(a,b) (max_arg1=(a),max_arg2=(b),(max_arg1) > (max_arg2) ? \
		   (max_arg1) : (max_arg2))
static float min_arg1, min_arg2;
#define FMIN(a,b) (min_arg1=(a),min_arg2=(b),(min_arg1) < (min_arg2) ? \
		   (min_arg1) : (min_arg2))
static long lmax_arg1, lmax_arg2;
#define LMAX(a,b) (lmax_arg1=(a),lmax_arg2=(b),(lmax_arg1) > (lmax_arg2) ? \
		   (lmax_arg1) : (lmax_arg2))
static long lmin_arg1, lmin_arg2;
#define LMIN(a,b) (lmin_arg1=(a),lmin_arg2=(b),(lmin_arg1) < (lmin_arg2) ? \
		   (lmin_arg1) : (lmin_arg2))
static int imax_arg1, imax_arg2;
#define IMAX(a,b) (imax_arg1=(a),imax_arg2=(b),(imax_arg1) > (imax_arg2) ?	\
		   (imax_arg1) : (imax_arg2))
static int imin_arg1, imin_arg2;
#define IMIN(a,b) (imin_arg1=(a),imin_arg2=(b),(imin_arg1) < (imin_arg2) ?	\
		   (imin_arg1) : (imin_arg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))s

#endif