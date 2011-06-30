/* This file contains routines (perhaps in modified form) by third parties.
 * Use and distribution of these works are determined by their respective copyrights. */

#ifndef THIRDPARTY_H_
#define THIRDPARTY_H_


/*============================================================================*/

/* these are the routines in thirdparty.cpp */

bool linfit(
	long n,
	double x[], /* x[n] */
	double y[], /* y[n] */
	double &a,
	double &siga,
	double &b,
	double &sigb
);

/** number of predefined factorials, so largest factorial
 * that can be represented as double is (NPRE_FACTORIAL-1)! */
static const int NPRE_FACTORIAL = 171;

/** factorial: compute n! by lookup in table of predefined factorials */
double factorial(long n);

/** lfactorial: compute log10(n!), this sroutine cahes its results for efficiency */
double lfactorial(long n);

complex<double> cdgamma(complex<double> x);

double bessel_j0(double x);
double bessel_y0(double x);
double bessel_j1(double x);
double bessel_y1(double x);
double bessel_jn(int n, double x);
double bessel_yn(int n, double x);

double bessel_k0(double x);
double bessel_k0_scaled(double x);
double bessel_k1(double x);
double bessel_k1_scaled(double x);
double bessel_i0(double x);
double bessel_i0_scaled(double x);
double bessel_i1(double x);
double bessel_i1_scaled(double x);

double ellpk(double x);

/** expn, returns exponential integral,
 \param n is order, 1 for first integral integral
 \param x is argument, must be positive
 */
double expn(int n, double x);

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void);

/* These real versions are due to Isaku Wada, 2002/01/09 added */
/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);

/*============================================================================*/

/* these are the routines in thirdparty_interpolate.cpp */

void spline_cubic_set( long n, const double t[], const double y[], double ypp[],
		       int ibcbeg, double ybcbeg, int ibcend, double ybcend );
void spline_cubic_val( long n, const double t[], double tval, const double y[], const double ypp[],
		       double *yval, double *ypval, double *yppval );

/** spline: set of routines to do spline interpolation
 * call spline first to set coefficients up,
 * then call splint to interpolate
 * spldrv gives dy/dx(x) rather than y(x)
 * splint_safe and spldrv_safe check whether x is within bounds
\param x[] - array of x values given
\param y[] - array of y values given
\param n - number of points in these arrays
\param yp1 some sort of boundary condition, set to 2e31
\param yp2 some sort of boundary condition, set to 2e31
\param y2a[] - array n long to save coefficients
*/
inline void spline(const double x[], 
		   const double y[], 
		   long int n, 
		   double yp1, 
		   double ypn, 
		   double y2a[])
{
	int ibcbeg = ( yp1 > 0.99999e30 ) ? 2 : 1;
	double ybcbeg = ( yp1 > 0.99999e30 ) ? 0. : yp1;
	int ibcend = ( ypn > 0.99999e30 ) ? 2 : 1;
	double ybcend = ( ypn > 0.99999e30 ) ? 0. : ypn;
	spline_cubic_set( n, x, y, y2a, ibcbeg, ybcbeg, ibcend, ybcend );
	return;
}

/** splint: do spline interpolation
\param x[] - array of x values given
\param y[] - array of y valus
\param y2a[] - array of save values found above
\param n - number of points in these arrays
\param x value where x is desired
\param y2a[] - array n long to save coefficients
*/
inline void splint(const double xa[], 
		   const double ya[], 
		   const double y2a[], 
		   long int n, 
		   double x, 
		   double *y)
{
	spline_cubic_val( n, xa, x, ya, y2a, y, NULL, NULL );
	return;
}

inline void spldrv(const double xa[], 
		   const double ya[], 
		   const double y2a[], 
		   long int n, 
		   double x, 
		   double *y)
{
	spline_cubic_val( n, xa, x, ya, y2a, NULL, y, NULL );
	return;
}

/* wrapper routine for splint that checks whether x-value is within bounds
 * if the x-value is out of bounds, a flag will be raised and the function
 * will be evaluated at the nearest boundary */
/* >>chng 03 jan 15, added splint_safe, PvH */
inline void splint_safe(const double xa[], 
			const double ya[], 
			const double y2a[], 
			long int n, 
			double x, 
			double *y,
			bool *lgOutOfBounds)
{
	double xsafe;

	const double lo_bound = MIN2(xa[0],xa[n-1]);
	const double hi_bound = MAX2(xa[0],xa[n-1]);
	const double SAFETY = MAX2(hi_bound-lo_bound,1.)*10.*DBL_EPSILON;

	DEBUG_ENTRY( "splint_safe()" );

	if( x < lo_bound-SAFETY )
	{
		xsafe = lo_bound;
		*lgOutOfBounds = true;
	}
	else if( x > hi_bound+SAFETY )
	{
		xsafe = hi_bound;
		*lgOutOfBounds = true;
	}
	else
	{
		xsafe = x;
		*lgOutOfBounds = false;
	}

	splint(xa,ya,y2a,n,xsafe,y);
	return;
}

/* wrapper routine for spldrv that checks whether x-value is within bounds
 * if the x-value is out of bounds, a flag will be raised and the function
 * will be evaluated at the nearest boundary */
/* >>chng 03 jan 15, added spldrv_safe, PvH */
inline void spldrv_safe(const double xa[],
			const double ya[],
			const double y2a[],
			long int n,
			double x,
			double *y,
			bool *lgOutOfBounds)
{
	double xsafe;

	const double lo_bound = MIN2(xa[0],xa[n-1]);
	const double hi_bound = MAX2(xa[0],xa[n-1]);
	const double SAFETY = MAX2(fabs(lo_bound),fabs(hi_bound))*10.*DBL_EPSILON;

	DEBUG_ENTRY( "spldrv_safe()" );

	if( x < lo_bound-SAFETY )
	{
		xsafe = lo_bound;
		*lgOutOfBounds = true;
	}
	else if( x > hi_bound+SAFETY )
	{
		xsafe = hi_bound;
		*lgOutOfBounds = true;
	}
	else
	{
		xsafe = x;
		*lgOutOfBounds = false;
	}

	spldrv(xa,ya,y2a,n,xsafe,y);
	return;
}

/** lagrange: do lagrange interpolation of order n on x[], y[]
 * use with caution, especialy for high order n!
 * using spline interpolation above is preferred
\param x[]
\param y[]
\param n
\param xval
*/
double lagrange(const double x[], /* x[n] */
		const double y[], /* y[n] */
		long n,
		double xval);

/** do linear interpolation on x[], y[];
 * it is assumed that x[] is strictly monotonically increasing
\param x[]
\param y[]
\param n
\param xval
*/
double linint(const double x[], /* x[n] */
	      const double y[], /* y[n] */
	      long n,
	      double xval);

/** find index ilo such that x[ilo] <= xval < x[ilo+1] using bisection
 * this version implicitly assumes that x is monotonically increasing */
template<class T>
inline long hunt_bisect(const T x[], /* x[n] */
			long n,
			T xval)
{
	/* do bisection hunt */
	long ilo = 0, ihi = n-1;
	while( ihi-ilo > 1 )
	{
		long imid = (ilo+ihi)/2;
		if( xval < x[imid] )
			ihi = imid;
		else
			ilo = imid;
	}
	return ilo;
}

/** find index ilo such that x[ilo] <= xval < x[ilo+1] using bisection
 * this version implicitly assumes that x is monotonically decreasing */
template<class T>
inline long hunt_bisect_reverse(const T x[], /* x[n] */
				long n,
				T xval)
{
	/* do bisection hunt */
	long ilo = 0, ihi = n-1;
	while( ihi-ilo > 1 )
	{
		long imid = (ilo+ihi)/2;
		if( xval <= x[imid] )
			ilo = imid;
		else
			ihi = imid;
	}
	return ilo;
}

/*============================================================================*/

/* these are the routines in thirdparty_lapack.cpp */

/* there are wrappers for lapack linear algebra routines.
 * there are two versions of the lapack routines - a fortran
 * version that I converted to C with forc to use if nothing else is available
 * (included in the Cloudy distribution),
 * and an option to link into an external lapack library that may be optimized
 * for your machine.  By default the tralated C routines will be used.
 * To use your machine's lapack library instead, define the macro
 * LAPACK and link into your library.  This is usually done with a command line
 * switch "-DLAPACK" on the compiler command, and the linker option "-llapack"
 */

/**getrf_wrapper return value is zero for success, non-zero is error condition 
\param M
\param N
\param *A
\param lda
\param *ipiv
\param *info
*/
void getrf_wrapper(long M, long N, double *A, long lda, int32 *ipiv, int32 *info);

/**getrs_wrapper return value is zero for success, non-zero is error condition 
\param trans
\param N
\param nrhs
\param *A
\param lda
\param *ipiv
\param *B
\param ldb
\param *info
*/
void getrs_wrapper(char trans, long N, long nrhs, double *A, long lda, int32 *ipiv, double *B, long ldb, int32 *info);

void humlik(int n, realnum x[], realnum y, realnum k[]);

/** the size of an MD5 checksum in characters */
static const unsigned int NMD5 = 32;

/** calculate the MD5 sum of a file */
string MD5file(const char* fnam, access_scheme scheme=AS_DATA_ONLY);
/** non-standard MD5 algorithm that skips eol characters and comments lines */
string MD5datafile(const char* fnam, access_scheme scheme=AS_DATA_ONLY);
/** calculate the MD5 sum of a string */
string MD5string(const string& str);

#endif /* THIRDPARTY_H_ */
