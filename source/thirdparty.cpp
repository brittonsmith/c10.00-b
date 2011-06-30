/* This file contains routines (perhaps in modified form) written by third parties.
 * Use and distribution of these works are determined by their respective copyrights. */
#include "cddefines.h"
#include "thirdparty.h"
#include "physconst.h"

inline double polevl(double x, const double coef[], int N);
inline double p1evl(double x, const double coef[], int N);
inline double chbevl(double, const double[], int);


/* the routine linfit was derived from the program slopes.f */

/********************************************************************/
/*                                                                  */
/*                      PROGRAM SLOPES                              */
/*                                                                  */
/*    PROGRAM TO COMPUTE THE THEORETICAL REGRESSION SLOPES          */
/*    AND UNCERTAINTIES (VIA DELTA METHOD), AND UNCERTAINTIES       */
/*    VIA BOOTSTRAP AND BIVARIATE NORMAL SIMULATION FOR A           */
/*    (X(I),Y(I)) DATA SET WITH UNKNOWN POPULATION DISTRIBUTION.    */
/*                                                                  */
/*       WRITTEN BY ERIC FEIGELSON, PENN STATE.  JUNE 1991          */
/*                                                                  */
/*                                                                  */
/*  A full description of these methods can be found in:            */
/*    Isobe, T., Feigelson, E. D., Akritas, M. and Babu, G. J.,     */
/*       Linear Regression in Astronomy I, Astrophys. J. 364, 104   */
/*       (1990)                                                     */
/*    Babu, G. J. and Feigelson, E. D., Analytical and Monte Carlo  */
/*       Comparisons of Six Different Linear Least Squares Fits,    */
/*       Communications in Statistics, Simulation & Computation,    */
/*       21, 533 (1992)                                             */
/*    Feigelson, E. D. and Babu, G. J., Linear Regression in        */
/*       Astronomy II, Astrophys. J. 397, 55 (1992).                */
/*                                                                  */
/********************************************************************/

/* this used to be sixlin, but only the first fit has been retained */

/********************************************************************/
/************************* routine linfit ***************************/
/********************************************************************/

bool linfit(
	long n,
	double x[], /* x[n] */
	double y[], /* y[n] */
	double &a,
	double &siga,
	double &b,
	double &sigb
)
{

/*
 *                       linear regression
 *     written by t. isobe, g. j. babu and e. d. feigelson
 *               center for space research, m.i.t.
 *                             and
 *              the pennsylvania state university
 *
 *                   rev. 1.0,   september 1990
 *
 *       this subroutine provides linear regression coefficients
 *    computed by six different methods described in isobe,
 *    feigelson, akritas, and babu 1990, astrophysical journal
 *    and babu and feigelson 1990, subm. to technometrics.
 *    the methods are ols(y/x), ols(x/y), ols bisector, orthogonal,
 *    reduced major axis, and mean-ols regressions.
 *
 *    [Modified and translated to C/C++ by Peter van Hoof, Royal
 *     Observatory of Belgium; only the first method has been retained]
 *     
 *
 *    input
 *         x[i] : independent variable
 *         y[i] : dependent variable
 *            n : number of data points
 *
 *    output
 *         a : intercept coefficients
 *         b : slope coefficients
 *      siga : standard deviations of intercepts
 *      sigb : standard deviations of slopes
 *
 *    error returns
 *         returns true when division by zero will
 *         occur (i.e. when sxy is zero)
 */

	DEBUG_ENTRY( "linfit()" );

	/* initializations */
	a = 0.0;
	siga = 0.0;
	b = 0.0;
	sigb = 0.0;

	/* compute averages and sums */
	double s1 = 0.0;
	double s2 = 0.0;
	for( long i=0; i < n; i++ )
	{
		s1 += x[i];
		s2 += y[i];
	}
	double rn = (double)n;
	double xavg = s1/rn;
	double yavg = s2/rn;
	double sxx = 0.0;
	double sxy = 0.0;
	for( long i=0; i < n; i++ )
	{
		x[i] -= xavg;
		y[i] -= yavg;
		sxx += pow2(x[i]);
		sxy += x[i]*y[i];
	}

	if( sxy == 0.0 )
	{
		return true;
	}

	/* compute the slope coefficient */
	b = sxy/sxx;

	/* compute intercept coefficient */
	a = yavg - b*xavg;

	/* prepare for computation of variances */
	double sum1 = 0.0;
	for( long i=0; i < n; i++ )
		sum1 += pow2(x[i]*(y[i] - b*x[i]));

	/* compute variance of the slope coefficient */
	sigb = sum1/pow2(sxx);

	/* compute variance of the intercept coefficient */
	for( long i=0; i < n; i++ )
		siga += pow2((y[i] - b*x[i])*(1.0 - rn*xavg*x[i]/sxx));

	/* convert variances to standard deviations */
	sigb = sqrt(sigb);
	siga = sqrt(siga)/rn;

	/* return data arrays to their original form */
	for( long i=0; i < n; i++ )
	{
		x[i] += xavg;
		y[i] += yavg;
	}
	return false;
}

/************************************************************************
 * This marks the end of the block of code from Isobe, Babu & Feigelson *
 ************************************************************************/


/* the routines factorial and lfactorial came originally from hydro_bauman.cpp
 * and were written by Robert Paul Bauman. lfactorial was modified by Peter van Hoof */

/************************************************************************************************/
/* pre-calculated factorials                                                                    */
/************************************************************************************************/

static const double pre_factorial[NPRE_FACTORIAL] =
{
	1.00000000000000000000e+00,
	1.00000000000000000000e+00,
	2.00000000000000000000e+00,
	6.00000000000000000000e+00,
	2.40000000000000000000e+01,
	1.20000000000000000000e+02,
	7.20000000000000000000e+02,
	5.04000000000000000000e+03,
	4.03200000000000000000e+04,
	3.62880000000000000000e+05,
	3.62880000000000000000e+06,
	3.99168000000000000000e+07,
	4.79001600000000000000e+08,
	6.22702080000000000000e+09,
	8.71782912000000000000e+10,
	1.30767436800000000000e+12,
	2.09227898880000000000e+13,
	3.55687428096000000000e+14,
	6.40237370572800000000e+15,
	1.21645100408832000000e+17,
	2.43290200817664000000e+18,
	5.10909421717094400000e+19,
	1.12400072777760768000e+21,
	2.58520167388849766400e+22,
	6.20448401733239439360e+23,
	1.55112100433309859840e+25,
	4.03291461126605635592e+26,
	1.08888694504183521614e+28,
	3.04888344611713860511e+29,
	8.84176199373970195470e+30,
	2.65252859812191058647e+32,
	8.22283865417792281807e+33,
	2.63130836933693530178e+35,
	8.68331761881188649615e+36,
	2.95232799039604140861e+38,
	1.03331479663861449300e+40,
	3.71993326789901217463e+41,
	1.37637530912263450457e+43,
	5.23022617466601111726e+44,
	2.03978820811974433568e+46,
	8.15915283247897734264e+47,
	3.34525266131638071044e+49,
	1.40500611775287989839e+51,
	6.04152630633738356321e+52,
	2.65827157478844876773e+54,
	1.19622220865480194551e+56,
	5.50262215981208894940e+57,
	2.58623241511168180614e+59,
	1.24139155925360726691e+61,
	6.08281864034267560801e+62,
	3.04140932017133780398e+64,
	1.55111875328738227999e+66,
	8.06581751709438785585e+67,
	4.27488328406002556374e+69,
	2.30843697339241380441e+71,
	1.26964033536582759243e+73,
	7.10998587804863451749e+74,
	4.05269195048772167487e+76,
	2.35056133128287857145e+78,
	1.38683118545689835713e+80,
	8.32098711274139014271e+81,
	5.07580213877224798711e+83,
	3.14699732603879375200e+85,
	1.98260831540444006372e+87,
	1.26886932185884164078e+89,
	8.24765059208247066512e+90,
	5.44344939077443063905e+92,
	3.64711109181886852801e+94,
	2.48003554243683059915e+96,
	1.71122452428141311337e+98,
	1.19785716699698917933e+100,
	8.50478588567862317347e+101,
	6.12344583768860868500e+103,
	4.47011546151268434004e+105,
	3.30788544151938641157e+107,
	2.48091408113953980872e+109,
	1.88549470166605025458e+111,
	1.45183092028285869606e+113,
	1.13242811782062978295e+115,
	8.94618213078297528506e+116,
	7.15694570462638022794e+118,
	5.79712602074736798470e+120,
	4.75364333701284174746e+122,
	3.94552396972065865030e+124,
	3.31424013456535326627e+126,
	2.81710411438055027626e+128,
	2.42270953836727323750e+130,
	2.10775729837952771662e+132,
	1.85482642257398439069e+134,
	1.65079551609084610774e+136,
	1.48571596448176149700e+138,
	1.35200152767840296226e+140,
	1.24384140546413072522e+142,
	1.15677250708164157442e+144,
	1.08736615665674307994e+146,
	1.03299784882390592592e+148,
	9.91677934870949688836e+149,
	9.61927596824821198159e+151,
	9.42689044888324774164e+153,
	9.33262154439441526381e+155,
	9.33262154439441526388e+157,
	9.42594775983835941673e+159,
	9.61446671503512660515e+161,
	9.90290071648618040340e+163,
	1.02990167451456276198e+166,
	1.08139675824029090008e+168,
	1.14628056373470835406e+170,
	1.22652020319613793888e+172,
	1.32464181945182897396e+174,
	1.44385958320249358163e+176,
	1.58824554152274293982e+178,
	1.76295255109024466316e+180,
	1.97450685722107402277e+182,
	2.23119274865981364576e+184,
	2.54355973347218755612e+186,
	2.92509369349301568964e+188,
	3.39310868445189820004e+190,
	3.96993716080872089396e+192,
	4.68452584975429065488e+194,
	5.57458576120760587943e+196,
	6.68950291344912705515e+198,
	8.09429852527344373681e+200,
	9.87504420083360135884e+202,
	1.21463043670253296712e+205,
	1.50614174151114087918e+207,
	1.88267717688892609901e+209,
	2.37217324288004688470e+211,
	3.01266001845765954361e+213,
	3.85620482362580421582e+215,
	4.97450422247728743840e+217,
	6.46685548922047366972e+219,
	8.47158069087882050755e+221,
	1.11824865119600430699e+224,
	1.48727070609068572828e+226,
	1.99294274616151887582e+228,
	2.69047270731805048244e+230,
	3.65904288195254865604e+232,
	5.01288874827499165889e+234,
	6.91778647261948848943e+236,
	9.61572319694108900019e+238,
	1.34620124757175246000e+241,
	1.89814375907617096864e+243,
	2.69536413788816277557e+245,
	3.85437071718007276916e+247,
	5.55029383273930478744e+249,
	8.04792605747199194159e+251,
	1.17499720439091082343e+254,
	1.72724589045463891049e+256,
	2.55632391787286558753e+258,
	3.80892263763056972532e+260,
	5.71338395644585458806e+262,
	8.62720977423324042775e+264,
	1.31133588568345254503e+267,
	2.00634390509568239384e+269,
	3.08976961384735088657e+271,
	4.78914290146339387432e+273,
	7.47106292628289444390e+275,
	1.17295687942641442768e+278,
	1.85327186949373479574e+280,
	2.94670227249503832518e+282,
	4.71472363599206132029e+284,
	7.59070505394721872577e+286,
	1.22969421873944943358e+289,
	2.00440157654530257674e+291,
	3.28721858553429622598e+293,
	5.42391066613158877297e+295,
	9.00369170577843736335e+297,
	1.50361651486499903974e+300,
	2.52607574497319838672e+302,
	4.26906800900470527345e+304,
	7.25741561530799896496e+306
};

double factorial( long n )
{
	DEBUG_ENTRY( "factorial()" );

	if( n < 0 || n >= NPRE_FACTORIAL )
	{
		fprintf( ioQQQ, "factorial: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return pre_factorial[n];
}

/* NB - this implementation is not thread-safe! */

class t_lfact : public Singleton<t_lfact>
{
	friend class Singleton<t_lfact>;
protected:
	t_lfact()
	{
		p_lf.reserve( 512 );
		p_lf.push_back( 0. ); /* log10( 0! ) */
		p_lf.push_back( 0. ); /* log10( 1! ) */
	}
private:
	vector<double> p_lf;
public:
	double get_lfact( unsigned long n )
	{
		if( n < p_lf.size() )
		{
			return p_lf[n];
		}
		else
		{
			for( unsigned long i=(unsigned long)p_lf.size(); i <= n; i++ )
				p_lf.push_back( p_lf[i-1] + log10(static_cast<double>(i)) );
			return p_lf[n];
		}
	}
};

double lfactorial( long n )
{
	/******************************/
	/*                 n          */
	/*                ---         */
	/*    log( n! ) =  >  log(j)  */
	/*                ---         */
	/*                j=1         */
	/******************************/

	DEBUG_ENTRY( "lfactorial()" );

	if( n < 0 )
	{
		fprintf( ioQQQ, "lfactorial: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return t_lfact::Inst().get_lfact( static_cast<unsigned long>(n) );
}

/*******************************************************************
 * This marks the end of the block of code from Robert Paul Bauman *
 *******************************************************************/


/* complex Gamma function in double precision */
/* this routine is a slightly modified version of the one found 
 * at http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html	
 * The following copyright applies: 
 *   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
 *   You may use, copy, modify this code for any purpose and 
 *   without fee. You may distribute this ORIGINAL package.	*/
complex<double> cdgamma(complex<double> x)
{
	double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

	DEBUG_ENTRY( "cdgamma()" );

	xr = x.real();
	xi = x.imag();
	if(xr < 0)
	{
		wr = 1. - xr;
		wi = -xi;
	}
	else
	{
		wr = xr;
		wi = xi;
	}
	ur = wr + 6.00009857740312429;
	vr = ur * (wr + 4.99999857982434025) - wi * wi;
	vi = wi * (wr + 4.99999857982434025) + ur * wi;
	yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 
		 0.293729529320536228;
	yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
	ur = vr * (wr + 4.00000003016801681) - vi * wi;
	ui = vi * (wr + 4.00000003016801681) + vr * wi;
	vr = ur * (wr + 2.99999999944915534) - ui * wi;
	vi = ui * (wr + 2.99999999944915534) + ur * wi;
	yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
	yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
	ur = vr * (wr + 2.00000000000603851) - vi * wi;
	ui = vi * (wr + 2.00000000000603851) + vr * wi;
	vr = ur * (wr + 0.999999999999975753) - ui * wi;
	vi = ui * (wr + 0.999999999999975753) + ur * wi;
	yr += ur * 10.5400280458730808 + vr;
	yi += ui * 10.5400280458730808 + vi;
	ur = vr * wr - vi * wi;
	ui = vi * wr + vr * wi;
	t = ur * ur + ui * ui;
	vr = yr * ur + yi * ui + t * 0.0327673720261526849;
	vi = yi * ur - yr * ui;
	yr = wr + 7.31790632447016203;
	ur = log(yr * yr + wi * wi) * 0.5 - 1;
	ui = atan2(wi, yr);
	yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
	yi = ui * (wr - 0.5) + ur * wi;
	ur = yr * cos(yi);
	ui = yr * sin(yi);
	yr = ur * vr - ui * vi;
	yi = ui * vr + ur * vi;
	if(xr < 0)
	{
		wr = xr * 3.14159265358979324;
		wi = exp(xi * 3.14159265358979324);
		vi = 1 / wi;
		ur = (vi + wi) * sin(wr);
		ui = (vi - wi) * cos(wr);
		vr = ur * yr + ui * yi;
		vi = ui * yr - ur * yi;
		ur = 6.2831853071795862 / (vr * vr + vi * vi);
		yr = ur * vr;
		yi = ur * vi;
	}
	return complex<double>( yr, yi );
}

/*************************************************************
 * This marks the end of the block of code from Takuya OOURA *
 *************************************************************/

/*====================================================================
 *
 * Below are routines from the Cephes library.
 *
 * This is the copyright statement included with the library:
 *
 *   Some software in this archive may be from the book _Methods and
 * Programs for Mathematical Functions_ (Prentice-Hall, 1989) or
 * from the Cephes Mathematical Library, a commercial product. In
 * either event, it is copyrighted by the author.  What you see here
 * may be used freely but it comes with no support or guarantee.
 *
 *   The two known misprints in the book are repaired here in the
 * source listings for the gamma function and the incomplete beta
 * integral.
 *
 *   Stephen L. Moshier
 *   moshier@world.std.com
 *
 *====================================================================*/

/*							j0.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */
/*							y0.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        9400       7.0e-17     7.9e-18
 *    IEEE      0, 30       30000       1.3e-15     1.6e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

static const double b0_PP[7] = {
	7.96936729297347051624e-4,
	8.28352392107440799803e-2,
	1.23953371646414299388e0,
	5.44725003058768775090e0,
	8.74716500199817011941e0,
	5.30324038235394892183e0,
	9.99999999999999997821e-1,
};

static const double b0_PQ[7] = {
	9.24408810558863637013e-4,
	8.56288474354474431428e-2,
	1.25352743901058953537e0,
	5.47097740330417105182e0,
	8.76190883237069594232e0,
	5.30605288235394617618e0,
	1.00000000000000000218e0,
};

static const double b0_QP[8] = {
	-1.13663838898469149931e-2,
	-1.28252718670509318512e0,
	-1.95539544257735972385e1,
	-9.32060152123768231369e1,
	-1.77681167980488050595e2,
	-1.47077505154951170175e2,
	-5.14105326766599330220e1,
	-6.05014350600728481186e0,
};

static const double b0_QQ[7] = {
	/* 1.00000000000000000000e0,*/
	6.43178256118178023184e1,
	8.56430025976980587198e2,
	3.88240183605401609683e3,
	7.24046774195652478189e3,
	5.93072701187316984827e3,
	2.06209331660327847417e3,
	2.42005740240291393179e2,
};

static const double b0_YP[8] = {
	1.55924367855235737965e4,
	-1.46639295903971606143e7,
	5.43526477051876500413e9,
	-9.82136065717911466409e11,
	8.75906394395366999549e13,
	-3.46628303384729719441e15,
	4.42733268572569800351e16,
	-1.84950800436986690637e16,
};

static const double b0_YQ[7] = {
	/* 1.00000000000000000000e0,*/
	1.04128353664259848412e3,
	6.26107330137134956842e5,
	2.68919633393814121987e8,
	8.64002487103935000337e10,
	2.02979612750105546709e13,
	3.17157752842975028269e15,
	2.50596256172653059228e17,
};

/*  5.783185962946784521175995758455807035071 */
static const double DR1 = 5.78318596294678452118e0;
/* 30.47126234366208639907816317502275584842 */
static const double DR2 = 3.04712623436620863991e1;

static double b0_RP[4] = {
	-4.79443220978201773821e9,
	1.95617491946556577543e12,
	-2.49248344360967716204e14,
	9.70862251047306323952e15,
};

static double b0_RQ[8] = {
	/* 1.00000000000000000000e0,*/
	4.99563147152651017219e2,
	1.73785401676374683123e5,
	4.84409658339962045305e7,
	1.11855537045356834862e10,
	2.11277520115489217587e12,
	3.10518229857422583814e14,
	3.18121955943204943306e16,
	1.71086294081043136091e18,
};

static const double TWOOPI = 2./PI;
static const double SQ2OPI = sqrt(2./PI);
static const double PIO4 = PI/4.;

double bessel_j0(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_j0()" );

	if( x < 0 )
		x = -x;

	if( x <= 5.0 )
	{
		z = x * x;
		if( x < 1.0e-5 )
			return 1.0 - z/4.0;

		p = (z - DR1) * (z - DR2);
		p = p * polevl( z, b0_RP, 3)/p1evl( z, b0_RQ, 8 );
		return p;
	}

	w = 5.0/x;
	q = 25.0/(x*x);
	p = polevl( q, b0_PP, 6)/polevl( q, b0_PQ, 6 );
	q = polevl( q, b0_QP, 7)/p1evl( q, b0_QQ, 7 );
	xn = x - PIO4;
	p = p * cos(xn) - w * q * sin(xn);
	return p * SQ2OPI / sqrt(x);
}

/*							y0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  y0(x)  -  2 * log(x) * j0(x) / PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EULER ) / PI
 * = 0.073804295108687225.
 */

double bessel_y0(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_y0()" );

	if( x <= 5.0 )
	{
		if( x <= 0.0 )
		{
			fprintf( ioQQQ, "bessel_y0: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		z = x * x;
		w = polevl( z, b0_YP, 7 ) / p1evl( z, b0_YQ, 7 );
		w += TWOOPI * log(x) * bessel_j0(x);
		return w;
	}

	w = 5.0/x;
	z = 25.0 / (x * x);
	p = polevl( z, b0_PP, 6)/polevl( z, b0_PQ, 6 );
	q = polevl( z, b0_QP, 7)/p1evl( z, b0_QQ, 7 );
	xn = x - PIO4;
	p = p * sin(xn) + w * q * cos(xn);
	return p * SQ2OPI / sqrt(x);
}

/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       8.6e-17     1.3e-17
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

static const double b1_RP[4] = {
	-8.99971225705559398224e8,
	4.52228297998194034323e11,
	-7.27494245221818276015e13,
	3.68295732863852883286e15,
};

static const double b1_RQ[8] = {
	/* 1.00000000000000000000E0,*/
	6.20836478118054335476e2,
	2.56987256757748830383e5,
	8.35146791431949253037e7,
	2.21511595479792499675e10,
	4.74914122079991414898e12,
	7.84369607876235854894e14,
	8.95222336184627338078e16,
	5.32278620332680085395e18,
};

static const double b1_PP[7] = {
	7.62125616208173112003e-4,
	7.31397056940917570436e-2,
	1.12719608129684925192e0,
	5.11207951146807644818e0,
	8.42404590141772420927e0,
	5.21451598682361504063e0,
	1.00000000000000000254e0,
};

static const double b1_PQ[7] = {
	5.71323128072548699714e-4,
	6.88455908754495404082e-2,
	1.10514232634061696926e0,
	5.07386386128601488557e0,
	8.39985554327604159757e0,
	5.20982848682361821619e0,
	9.99999999999999997461e-1,
};

static const double b1_QP[8] = {
	5.10862594750176621635e-2,
	4.98213872951233449420e0,
	7.58238284132545283818e1,
	3.66779609360150777800e2,
	7.10856304998926107277e2,
	5.97489612400613639965e2,
	2.11688757100572135698e2,
	2.52070205858023719784e1,
};

static const double b1_QQ[7] = {
	/* 1.00000000000000000000e0,*/
	7.42373277035675149943e1,
	1.05644886038262816351e3,
	4.98641058337653607651e3,
	9.56231892404756170795e3,
	7.99704160447350683650e3,
	2.82619278517639096600e3,
	3.36093607810698293419e2,
};

static const double b1_YP[6] = {
	1.26320474790178026440e9,
	-6.47355876379160291031e11,
	1.14509511541823727583e14,
	-8.12770255501325109621e15,
	2.02439475713594898196e17,
	-7.78877196265950026825e17,
};

static const double b1_YQ[8] = {
	/* 1.00000000000000000000E0,*/
	5.94301592346128195359E2,
	2.35564092943068577943E5,
	7.34811944459721705660E7,
	1.87601316108706159478E10,
	3.88231277496238566008E12,
	6.20557727146953693363E14,
	6.87141087355300489866E16,
	3.97270608116560655612E18,
};

static const double Z1 = 1.46819706421238932572E1;
static const double Z2 = 4.92184563216946036703E1;

static const double THPIO4 = 3.*PI/4.;

double bessel_j1(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_j1()" );

	w = x;
	if( x < 0 )
		w = -x;

	if( w <= 5.0 )
	{
		z = x * x;	
		w = polevl( z, b1_RP, 3 ) / p1evl( z, b1_RQ, 8 );
		w = w * x * (z - Z1) * (z - Z2);
		return w;
	}

	w = 5.0/x;
	z = w * w;
	p = polevl( z, b1_PP, 6)/polevl( z, b1_PQ, 6 );
	q = polevl( z, b1_QP, 7)/p1evl( z, b1_QQ, 7 );
	xn = x - THPIO4;
	p = p * cos(xn) - w * q * sin(xn);
	return p * SQ2OPI / sqrt(x);
}

double bessel_y1(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_y1()" );

	if( x <= 5.0 )
	{
		if( x <= 0.0 )
		{
			fprintf( ioQQQ, "bessel_y1: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		z = x * x;
		w = x * (polevl( z, b1_YP, 5 ) / p1evl( z, b1_YQ, 8 ));
		w += TWOOPI * ( bessel_j1(x) * log(x)  -  1.0/x );
		return w;
	}

	w = 5.0/x;
	z = w * w;
	p = polevl( z, b1_PP, 6 )/polevl( z, b1_PQ, 6 );
	q = polevl( z, b1_QP, 7 )/p1evl( z, b1_QQ, 7 );
	xn = x - THPIO4;
	p = p * sin(xn) + w * q * cos(xn);
	return p * SQ2OPI / sqrt(x);
}

/*							jn.c
 *
 *	Bessel function of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, jn();
 *
 * y = jn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The ratio of jn(x) to j0(x) is computed by backward
 * recurrence.  First the ratio jn/jn-1 is found by a
 * continued fraction expansion.  Then the recurrence
 * relating successive orders is applied until j0 or j1 is
 * reached.
 *
 * If n = 0 or 1 the routine for j0 or j1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   range      # trials      peak         rms
 *    DEC       0, 30        5500       6.9e-17     9.3e-18
 *    IEEE      0, 30        5000       4.4e-16     7.9e-17
 *
 *
 * Not suitable for large n or x. Use jv() instead.
 *
 */

/*							jn.c
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

double bessel_jn(int n, double x)
{
	double pkm2, pkm1, pk, xk, r, ans;
	int k, sign;

	DEBUG_ENTRY( "bessel_jn()" );

	if( n < 0 )
	{
		n = -n;
		if( (n & 1) == 0 )	/* -1**n */
			sign = 1;
		else
			sign = -1;
	}
	else
		sign = 1;

	if( x < 0.0 )
	{
		if( n & 1 )
			sign = -sign;
		x = -x;
	}

	if( x < DBL_EPSILON )
	{
		return sign * powi(x/2.,n)/factorial(n);
	}

	if( n == 0 )
	{
		return sign * bessel_j0(x);
	}
	if( n == 1 )
	{
		return sign * bessel_j1(x);
	}
	// avoid cancellation error for very small x
	if( n == 2 && x > 0.1 )
	{
		return sign * (2.0 * bessel_j1(x) / x  -  bessel_j0(x));
	}

	/* continued fraction */
	k = 52;

	pk = 2 * (n + k);
	ans = pk;
	xk = x * x;

	do
	{
		pk -= 2.0;
		ans = pk - (xk/ans);
	}
	while( --k > 0 );
	ans = x/ans;

	/* backward recurrence */
	pk = 1.0;
	pkm1 = 1.0/ans;
	k = n-1;
	r = 2 * k;

	do
	{
		pkm2 = (pkm1 * r  -  pk * x) / x;
		pk = pkm1;
		pkm1 = pkm2;
		r -= 2.0;
	}
	while( --k > 0 );

	if( fabs(pk) > fabs(pkm1) )
		ans = bessel_j1(x)/pk;
	else
		ans = bessel_j0(x)/pkm1;
	return sign * ans;
}

/*							yn.c
 *
 *	Bessel function of second kind of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, yn();
 * int n;
 *
 * y = yn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The function is evaluated by forward recurrence on
 * n, starting with values computed by the routines
 * y0() and y1().
 *
 * If n = 0 or 1 the routine for y0 or y1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *
 *                      Absolute error, except relative
 *                      when y > 1:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        2200       2.9e-16     5.3e-17
 *    IEEE      0, 30       30000       3.4e-15     4.3e-16
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * yn singularity   x = 0              MAXNUM
 * yn overflow                         MAXNUM
 *
 * Spot checked against tables for x, n between 0 and 100.
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

double bessel_yn(int n, double x)
{
	double an, anm1, anm2, r;
	int k, sign;

	DEBUG_ENTRY( "bessel_yn()" );

	if( n < 0 )
	{
		n = -n;
		if( (n & 1) == 0 )	/* -1**n */
			sign = 1;
		else
			sign = -1;
	}
	else
		sign = 1;

	if( n == 0 )
	{
		return sign * bessel_y0(x);
	}
	if( n == 1 )
	{
		return sign * bessel_y1(x);
	}

	/* test for overflow */
	if( x <= 0.0 )
	{
		fprintf( ioQQQ, "bessel_yn: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* forward recurrence on n */
	anm2 = bessel_y0(x);
	anm1 = bessel_y1(x);
	k = 1;
	r = 2 * k;
	do
	{
		an = r * anm1 / x  -  anm2;
		anm2 = anm1;
		anm1 = an;
		r += 2.0;
		++k;
	}
	while( k < n );
	return sign * an;
}

/*							k0.c
 *
 *	Modified Bessel function, third kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0();
 *
 * y = k0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order zero of the argument.
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at 2000 random points between 0 and 8.  Peak absolute
 * error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3100       1.3e-16     2.1e-17
 *    IEEE      0, 30       30000       1.2e-15     1.6e-16
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 *  K0 domain          x <= 0          MAXNUM
 *
 */
/*							k0e()
 *
 *	Modified Bessel function, third kind, order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0e();
 *
 * y = k0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order zero of the argument.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.4e-15     1.4e-16
 * See k0().
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
 * in the interval [0,2].  The odd order coefficients are all
 * zero; only the even order coefficients are listed.
 * 
 * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EULER.
 */

static const double k0_A[] =
{
	1.37446543561352307156e-16,
	4.25981614279661018399e-14,
	1.03496952576338420167e-11,
	1.90451637722020886025e-9,
	2.53479107902614945675e-7,
	2.28621210311945178607e-5,
	1.26461541144692592338e-3,
	3.59799365153615016266e-2,
	3.44289899924628486886e-1,
	-5.35327393233902768720e-1
};

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 * 
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */

static const double k0_B[] = {
	5.30043377268626276149e-18,
	-1.64758043015242134646e-17,
	5.21039150503902756861e-17,
	-1.67823109680541210385e-16,
	5.51205597852431940784e-16,
	-1.84859337734377901440e-15,
	6.34007647740507060557e-15,
	-2.22751332699166985548e-14,
	8.03289077536357521100e-14,
	-2.98009692317273043925e-13,
	1.14034058820847496303e-12,
	-4.51459788337394416547e-12,
	1.85594911495471785253e-11,
	-7.95748924447710747776e-11,
	3.57739728140030116597e-10,
	-1.69753450938905987466e-9,
	8.57403401741422608519e-9,
	-4.66048989768794782956e-8,
	2.76681363944501510342e-7,
	-1.83175552271911948767e-6,
	1.39498137188764993662e-5,
	-1.28495495816278026384e-4,
	1.56988388573005337491e-3,
	-3.14481013119645005427e-2,
	2.44030308206595545468e0
};

double bessel_k0(double x)
{
	double y, z;

	DEBUG_ENTRY( "bessel_k0()" );

	if( x <= 0.0 )
	{
		fprintf( ioQQQ, "bessel_k0: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x <= 2.0 )
	{
		y = x * x - 2.0;
		y = chbevl( y, k0_A, 10 ) - log( 0.5 * x ) * bessel_i0(x);
		return y;
	}
	z = 8.0/x - 2.0;
	y = exp(-x) * chbevl( z, k0_B, 25 ) / sqrt(x);
	return y;
}

double bessel_k0_scaled(double x)
{
	double y;

	DEBUG_ENTRY( "bessel_k0_scaled()" );

	if( x <= 0.0 )
	{
		fprintf( ioQQQ, "bessel_k0_scaled: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x <= 2.0 )
	{
		y = x * x - 2.0;
		y = chbevl( y, k0_A, 10 ) - log( 0.5 * x ) * bessel_i0(x);
		return y * exp(x);
	}
	return chbevl( 8.0/x - 2.0, k0_B, 25 ) / sqrt(x);
}

/*							k1.c
 *
 *	Modified Bessel function, third kind, order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1();
 *
 * y = k1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the modified Bessel function of the third kind
 * of order one of the argument.
 *
 * The range is partitioned into the two intervals [0,2] and
 * (2, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3300       8.9e-17     2.2e-17
 *    IEEE      0, 30       30000       1.2e-15     1.6e-16
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * k1 domain          x <= 0          MAXNUM
 *
 */
/*							k1e.c
 *
 *	Modified Bessel function, third kind, order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1e();
 *
 * y = k1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order one of the argument:
 *
 *      k1e(x) = exp(x) * k1(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       7.8e-16     1.2e-16
 * See k1().
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
 * in the interval [0,2].
 * 
 * lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
 */

static const double k1_A[] =
{
	-7.02386347938628759343e-18,
	-2.42744985051936593393e-15,
	-6.66690169419932900609e-13,
	-1.41148839263352776110e-10,
	-2.21338763073472585583e-8,
	-2.43340614156596823496e-6,
	-1.73028895751305206302e-4,
	-6.97572385963986435018e-3,
	-1.22611180822657148235e-1,
	-3.53155960776544875667e-1,
	1.52530022733894777053e0
};

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */

static const double k1_B[] =
{
	-5.75674448366501715755e-18,
	1.79405087314755922667e-17,
	-5.68946255844285935196e-17,
	1.83809354436663880070e-16,
	-6.05704724837331885336e-16,
	2.03870316562433424052e-15,
	-7.01983709041831346144e-15,
	2.47715442448130437068e-14,
	-8.97670518232499435011e-14,
	3.34841966607842919884e-13,
	-1.28917396095102890680e-12,
	5.13963967348173025100e-12,
	-2.12996783842756842877e-11,
	9.21831518760500529508e-11,
	-4.19035475934189648750e-10,
	2.01504975519703286596e-9,
	-1.03457624656780970260e-8,
	5.74108412545004946722e-8,
	-3.50196060308781257119e-7,
	2.40648494783721712015e-6,
	-1.93619797416608296024e-5,
	1.95215518471351631108e-4,
	-2.85781685962277938680e-3,
	1.03923736576817238437e-1,
	2.72062619048444266945e0
};

double bessel_k1(double x)
{
	double y, z;

	DEBUG_ENTRY( "bessel_k1()" );

	z = 0.5 * x;
	if( z <= 0.0 )
	{
		fprintf( ioQQQ, "bessel_k1: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x <= 2.0 )
	{
		y = x * x - 2.0;
		y =  log(z) * bessel_i1(x)  +  chbevl( y, k1_A, 11 ) / x;
		return y;
	}
	return exp(-x) * chbevl( 8.0/x - 2.0, k1_B, 25 ) / sqrt(x);
}

double bessel_k1_scaled(double x)
{
	double y;

	DEBUG_ENTRY( "bessel_k1_scaled()" );

	if( x <= 0.0 )
	{
		fprintf( ioQQQ, "bessel_k1_scaled: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x <= 2.0 )
	{
		y = x * x - 2.0;
		y =  log( 0.5 * x ) * bessel_i1(x)  +  chbevl( y, k1_A, 11 ) / x;
		return y * exp(x);
	}
	return chbevl( 8.0/x - 2.0, k1_B, 25 ) / sqrt(x);
}

/*							i0.c
 *
 *	Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0();
 *
 * y = i0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         6000       8.2e-17     1.9e-17
 *    IEEE      0,30        30000       5.8e-16     1.4e-16
 *
 */
/*							i0e.c
 *
 *	Modified Bessel function of order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0e();
 *
 * y = i0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order zero of the argument.
 *
 * The function is defined as i0e(x) = exp(-|x|) j0( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       5.4e-16     1.2e-16
 * See i0().
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */

static const double i0_A[] =
{
	-4.41534164647933937950e-18,
	3.33079451882223809783e-17,
	-2.43127984654795469359e-16,
	1.71539128555513303061e-15,
	-1.16853328779934516808e-14,
	7.67618549860493561688e-14,
	-4.85644678311192946090e-13,
	2.95505266312963983461e-12,
	-1.72682629144155570723e-11,
	9.67580903537323691224e-11,
	-5.18979560163526290666e-10,
	2.65982372468238665035e-9,
	-1.30002500998624804212e-8,
	6.04699502254191894932e-8,
	-2.67079385394061173391e-7,
	1.11738753912010371815e-6,
	-4.41673835845875056359e-6,
	1.64484480707288970893e-5,
	-5.75419501008210370398e-5,
	1.88502885095841655729e-4,
	-5.76375574538582365885e-4,
	1.63947561694133579842e-3,
	-4.32430999505057594430e-3,
	1.05464603945949983183e-2,
	-2.37374148058994688156e-2,
	4.93052842396707084878e-2,
	-9.49010970480476444210e-2,
	1.71620901522208775349e-1,
	-3.04682672343198398683e-1,
	6.76795274409476084995e-1
};

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

static const double i0_B[] =
{
	-7.23318048787475395456e-18,
	-4.83050448594418207126e-18,
	4.46562142029675999901e-17,
	3.46122286769746109310e-17,
	-2.82762398051658348494e-16,
	-3.42548561967721913462e-16,
	1.77256013305652638360e-15,
	3.81168066935262242075e-15,
	-9.55484669882830764870e-15,
	-4.15056934728722208663e-14,
	1.54008621752140982691e-14,
	3.85277838274214270114e-13,
	7.18012445138366623367e-13,
	-1.79417853150680611778e-12,
	-1.32158118404477131188e-11,
	-3.14991652796324136454e-11,
	1.18891471078464383424e-11,
	4.94060238822496958910e-10,
	3.39623202570838634515e-9,
	2.26666899049817806459e-8,
	2.04891858946906374183e-7,
	2.89137052083475648297e-6,
	6.88975834691682398426e-5,
	3.36911647825569408990e-3,
	8.04490411014108831608e-1
};

double bessel_i0(double x)
{
	double y;

	DEBUG_ENTRY( "bessel_i0()" );

	if( x < 0 )
		x = -x;

	if( x <= 8.0 )
	{
		y = (x/2.0) - 2.0;
		return exp(x) * chbevl( y, i0_A, 30 );
	}
	return exp(x) * chbevl( 32.0/x - 2.0, i0_B, 25 ) / sqrt(x);
}

double bessel_i0_scaled(double x)
{
	double y;

	DEBUG_ENTRY( "bessel_i0_scaled()" );

	if( x < 0 )
		x = -x;

	if( x <= 8.0 )
	{
		y = (x/2.0) - 2.0;
		return chbevl( y, i0_A, 30 );
	}
	return chbevl( 32.0/x - 2.0, i0_B, 25 ) / sqrt(x);
}

/*							i1.c
 *
 *	Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1();
 *
 * y = i1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3400       1.2e-16     2.3e-17
 *    IEEE      0, 30       30000       1.9e-15     2.1e-16
 *
 *
 */
/*							i1e.c
 *
 *	Modified Bessel function of order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1e();
 *
 * y = i1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order one of the argument.
 *
 * The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       2.0e-15     2.0e-16
 * See i1().
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */

static double i1_A[] =
{
	2.77791411276104639959e-18,
	-2.11142121435816608115e-17,
	1.55363195773620046921e-16,
	-1.10559694773538630805e-15,
	7.60068429473540693410e-15,
	-5.04218550472791168711e-14,
	3.22379336594557470981e-13,
	-1.98397439776494371520e-12,
	1.17361862988909016308e-11,
	-6.66348972350202774223e-11,
	3.62559028155211703701e-10,
	-1.88724975172282928790e-9,
	9.38153738649577178388e-9,
	-4.44505912879632808065e-8,
	2.00329475355213526229e-7,
	-8.56872026469545474066e-7,
	3.47025130813767847674e-6,
	-1.32731636560394358279e-5,
	4.78156510755005422638e-5,
	-1.61760815825896745588e-4,
	5.12285956168575772895e-4,
	-1.51357245063125314899e-3,
	4.15642294431288815669e-3,
	-1.05640848946261981558e-2,
	2.47264490306265168283e-2,
	-5.29459812080949914269e-2,
	1.02643658689847095384e-1,
	-1.76416518357834055153e-1,
	2.52587186443633654823e-1
};

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */

static double i1_B[] =
{
	7.51729631084210481353e-18,
	4.41434832307170791151e-18,
	-4.65030536848935832153e-17,
	-3.20952592199342395980e-17,
	2.96262899764595013876e-16,
	3.30820231092092828324e-16,
	-1.88035477551078244854e-15,
	-3.81440307243700780478e-15,
	1.04202769841288027642e-14,
	4.27244001671195135429e-14,
	-2.10154184277266431302e-14,
	-4.08355111109219731823e-13,
	-7.19855177624590851209e-13,
	2.03562854414708950722e-12,
	1.41258074366137813316e-11,
	3.25260358301548823856e-11,
	-1.89749581235054123450e-11,
	-5.58974346219658380687e-10,
	-3.83538038596423702205e-9,
	-2.63146884688951950684e-8,
	-2.51223623787020892529e-7,
	-3.88256480887769039346e-6,
	-1.10588938762623716291e-4,
	-9.76109749136146840777e-3,
	7.78576235018280120474e-1
};

double bessel_i1(double x)
{ 
	double y, z;

	DEBUG_ENTRY( "bessel_i1()" );

	z = fabs(x);
	if( z <= 8.0 )
	{
		y = (z/2.0) - 2.0;
		z = chbevl( y, i1_A, 29 ) * z * exp(z);
	}
	else
	{
		z = exp(z) * chbevl( 32.0/z - 2.0, i1_B, 25 ) / sqrt(z);
	}
	if( x < 0.0 )
		z = -z;
	return z;
}

double bessel_i1_scaled(double x)
{ 
	double y, z;

	DEBUG_ENTRY( "bessel_i1_scaled()" );

	z = fabs(x);
	if( z <= 8.0 )
	{
		y = (z/2.0) - 2.0;
		z = chbevl( y, i1_A, 29 ) * z;
	}
	else
	{
		z = chbevl( 32.0/z - 2.0, i1_B, 25 ) / sqrt(z);
	}
	if( x < 0.0 )
		z = -z;
	return z;
}

/*							ellpk.c
 *
 *	Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m1, y, ellpk();
 *
 * y = ellpk( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0,1        16000       3.5e-17     1.1e-17
 *    IEEE       0,1        30000       2.5e-16     6.8e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpk domain       x<0, x>1           0.0
 *
 */

/*
Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

static const double elk_P[] =
{
	1.37982864606273237150e-4,
	2.28025724005875567385e-3,
	7.97404013220415179367e-3,
	9.85821379021226008714e-3,
	6.87489687449949877925e-3,
	6.18901033637687613229e-3,
	8.79078273952743772254e-3,
	1.49380448916805252718e-2,
	3.08851465246711995998e-2,
	9.65735902811690126535e-2,
	1.38629436111989062502e0
};

static const double elk_Q[] =
{
	2.94078955048598507511e-5,
	9.14184723865917226571e-4,
	5.94058303753167793257e-3,
	1.54850516649762399335e-2,
	2.39089602715924892727e-2,
	3.01204715227604046988e-2,
	3.73774314173823228969e-2,
	4.88280347570998239232e-2,
	7.03124996963957469739e-2,
	1.24999999999870820058e-1,
	4.99999999999999999821e-1
};

static const double C1 = 1.3862943611198906188e0; /* log(4) */

double ellpk(double x)
{
	DEBUG_ENTRY( "ellpk()" );

	if( x < 0.0 || x > 1.0 )
	{
		fprintf( ioQQQ, "ellpk: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x > DBL_EPSILON )
	{
		return polevl(x,elk_P,10) - log(x) * polevl(x,elk_Q,10);
	}
	else
	{
		if( x == 0.0 )
		{
			fprintf( ioQQQ, "ellpk: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			return C1 - 0.5 * log(x);
		}
	}
}

/*							expn.c
 *
 *		Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, expn();
 *
 * y = expn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                  inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                 1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        5000       2.0e-16     4.6e-17
 *    IEEE      0, 30       10000       1.7e-15     3.6e-16
 *
 */

/* Cephes Math Library Release 2.8:  June, 2000
   Copyright 1985, 2000 by Stephen L. Moshier */

static const double MAXLOG = log(DBL_MAX);
static const double BIG = 1.44115188075855872E+17; /* 2^57 */

/*expn exponential intergral for any n */
double expn(int n, double x)
{
	double ans, r, t, yk, xk;
	double pk, pkm1, pkm2, qk, qkm1, qkm2;
	double psi, z;
	int i, k;

	DEBUG_ENTRY( "expn()" );

	if( n < 0 || x < 0. )
	{
		fprintf( ioQQQ, "expn: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x > MAXLOG )
	{
		return 0.0;
	}

	if( x == 0.0 )
	{
		if( n < 2 )
		{
			fprintf( ioQQQ, "expn: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			return 1.0/((double)n-1.0);
		}
	}

	if( n == 0 )
	{
		return exp(-x)/x;
	}

	/*		Expansion for large n		*/
	if( n > 5000 )
	{
		xk = x + n;
		yk = 1.0 / (xk * xk);
		t = n;
		ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
		ans = yk * (ans + t * (t  -  2.0 * x));
		ans = yk * (ans + t);
		ans = (ans + 1.0) * exp( -x ) / xk;
		return ans;
	}

	if( x <= 1.0 )
	{
		/*		Power series expansion		*/
		psi = -EULER - log(x);
		for( i=1; i < n; i++ )
			psi = psi + 1.0/i;

		z = -x;
		xk = 0.0;
		yk = 1.0;
		pk = 1.0 - n;
		if( n == 1 )
			ans = 0.0;
		else
			ans = 1.0/pk;
		do
		{
			xk += 1.0;
			yk *= z/xk;
			pk += 1.0;
			if( pk != 0.0 )
			{
				ans += yk/pk;
			}
			if( ans != 0.0 )
				t = fabs(yk/ans);
			else
				t = 1.0;
		}
		while( t > DBL_EPSILON );
		ans = powi(z,n-1)*psi/factorial(n-1) - ans;
		return ans;
	}
	else
	{
		/*		continued fraction		*/
		k = 1;
		pkm2 = 1.0;
		qkm2 = x;
		pkm1 = 1.0;
		qkm1 = x + n;
		ans = pkm1/qkm1;

		do
		{
			k += 1;
			if( is_odd(k) )
			{
				yk = 1.0;
				xk = static_cast<double>(n + (k-1)/2);
			}
			else
			{
				yk = x;
				xk = static_cast<double>(k/2);
			}
			pk = pkm1 * yk  +  pkm2 * xk;
			qk = qkm1 * yk  +  qkm2 * xk;
			if( qk != 0 )
			{
				r = pk/qk;
				t = fabs( (ans - r)/r );
				ans = r;
			}
			else
				t = 1.0;
			pkm2 = pkm1;
			pkm1 = pk;
			qkm2 = qkm1;
			qkm1 = qk;
			if( fabs(pk) > BIG )
			{
				pkm2 /= BIG;
				pkm1 /= BIG;
				qkm2 /= BIG;
				qkm1 /= BIG;
			}
		}
		while( t > DBL_EPSILON );

		ans *= exp( -x );
		return ans;
	}
}

/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

inline double polevl(double x, const double coef[], int N)
{
	double ans;
	int i;
	const double *p = coef;

	ans = *p++;
	i = N;

	do
		ans = ans * x  +  *p++;
	while( --i );

	return ans;
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

inline double p1evl(double x, const double coef[], int N)
{
	double ans;
	const double *p = coef;
	int i;

	ans = x + *p++;
	i = N-1;

	do
		ans = ans * x  + *p++;
	while( --i );

	return ans;
}

/*							chbevl.c
 *
 *	Evaluate Chebyshev series
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N], chebevl();
 *
 * y = chbevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the series
 *
 *        N-1
 *         - '
 *  y  =   >   coef[i] T (x/2)
 *         -            i
 *        i=0
 *
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero
 * order term is last in the array.  Note N is the number of
 * coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must
 * have been transformed to x -> 2(2x - b - a)/(b-a) before
 * entering the routine.  This maps x from (a, b) to (-1, 1),
 * over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in
 * which (a, b) is mapped to (1/b, 1/a), the transformation
 * required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
 * this becomes x -> 4a/x - 1.
 *
 *
 *
 * SPEED:
 *
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

inline double chbevl(double x, const double array[], int n)
{
	double b0, b1, b2;
	const double *p = array;
	int i;

	b0 = *p++;
	b1 = 0.0;
	i = n - 1;

	do
	{
		b2 = b1;
		b1 = b0;
		b0 = x * b1  -  b2  + *p++;
	}
	while( --i );

	return 0.5*(b0-b2);
}

/*******************************************************************
 * This marks the end of the block of code from the Cephes library *
 *******************************************************************/



/* >>refer Mersenne Twister Matsumoto, M. & Nishimura, T. 1998, ACM Transactions on Modeling
 * >>refercon	and Computer Simulation (TOMACS), 8, 1998 */

/********************************************************************
 * This copyright notice must accompany the following block of code *
 *******************************************************************/

/* 
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Before using, initialize the state by using init_genrand(seed) 
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* Period parameters */  
static const int N = 624;
static const int M = 397;
static const unsigned long MATRIX_A = 0x9908b0dfUL;   /* constant vector a */
static const unsigned long UMASK = 0x80000000UL; /* most significant w-r bits */
static const unsigned long LMASK = 0x7fffffffUL; /* least significant r bits */
inline unsigned long MIXBITS(unsigned long u, unsigned long v)
{
	return (u & UMASK) | (v & LMASK);
}
inline unsigned long TWIST(unsigned long u, unsigned long v)
{
	return (MIXBITS(u,v) >> 1) ^ (v&1UL ? MATRIX_A : 0UL);
}

static unsigned long state[N]; /* the array for the state vector  */
static int nleft = 1;
static int initf = 0;
static unsigned long *nexxt;

/* initializes state[N] with a seed */
void init_genrand(unsigned long s)
{
	int j;
	state[0]= s & 0xffffffffUL;
	for (j=1; j<N; j++) {
		state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array state[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		state[j] &= 0xffffffffUL;  /* for >32 bit machines */
	}
	nleft = 1; initf = 1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
	int i, j, k;
	init_genrand(19650218UL);
	i=1; j=0;
	k = (N>key_length ? N : key_length);
	for (; k; k--) {
		state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
			+ init_key[j] + j; /* non linear */
		state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++; j++;
		if (i>=N) { state[0] = state[N-1]; i=1; }
		if (j>=key_length) j=0;
	}
	for (k=N-1; k; k--) {
		state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
			- i; /* non linear */
		state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++;
		if (i>=N) { state[0] = state[N-1]; i=1; }
	}

	state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
	nleft = 1; initf = 1;
}

static void next_state()
{
	unsigned long *p=state;
	int j;

	/* if init_genrand() has not been called, */
	/* a default initial seed is used         */
	if (initf==0) init_genrand(5489UL);

	nleft = N;
	nexxt = state;
    
	for (j=N-M+1; --j; p++) 
		*p = p[M] ^ TWIST(p[0], p[1]);

	for (j=M; --j; p++) 
		*p = p[M-N] ^ TWIST(p[0], p[1]);

	*p = p[M-N] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32()
{
	unsigned long y;

	if (--nleft == 0) next_state();
	y = *nexxt++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31()
{
	unsigned long y;

	if (--nleft == 0) next_state();
	y = *nexxt++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1()
{
	unsigned long y;

	if (--nleft == 0) next_state();
	y = *nexxt++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return (double)y * (1.0/4294967295.0); 
	/* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2()
{
	unsigned long y;

	if (--nleft == 0) next_state();
	y = *nexxt++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return (double)y * (1.0/4294967296.0); 
	/* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3()
{
	unsigned long y;

	if (--nleft == 0) next_state();
	y = *nexxt++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return ((double)y + 0.5) * (1.0/4294967296.0); 
	/* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53() 
{ 
	unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
	return (a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

/* These real versions are due to Isaku Wada, 2002/01/09 added */

/************************************************************************
 * This marks the end of the block of code from Matsumoto and Nishimura *
 ************************************************************************/

/*
  Calculate the Voigt profile aka Faddeeva function with relative
  error less than 10^(-R).  
  
  R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the
  user subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
  
  Translated from a Fortran version by R.J. Wells,
  see http://www.atm.ox.ac.uk/user/wells/voigt.html
*/
void humlik(int n, realnum x[], realnum y, realnum k[])
{
	/* n	IN  Number of points 
		x  IN  Input x array 
		y  IN  Input y value >=0.0
		k  OUT (Voigt) array */
	/* Initialized data */

	const realnum c[6] = { 1.0117281f,-.75197147f,.012557727f,.010022008f,
								  -2.4206814e-4f,5.0084806e-7f };
	const realnum s[6] = { 1.393237,.23115241f,-.15535147f,.0062183662f,
								9.1908299e-5f,-6.2752596e-7f };
	const realnum t[6] = { .31424038f,.94778839f,1.5976826f,2.2795071f,
								3.020637f,3.8897249f };

	const realnum R0 = 146.7f, R1 = 14.67f;

	/* Local variables */
	realnum d;
	int i, j;
	realnum a0, d0, e0, d2, e2, h0, e4, h2, h4, h6, p0, 
		p2, p4, p6, p8, z0, z2, z4, z6, z8, mf[6], pf[6], 
		mq[6], yf, pq[6], xm[6], ym[6], xp[6], xq, yq, yp[6];
	bool rg1, rg2, rg3;
	realnum abx, ypy0, xlim0, xlim1, xlim2, xlim3, xlim4, ypy0q, yrrtpi;

	rg1 = rg2 = rg3 = true;
	// Initialization to quiet warnings
	z0 = z2 = z4 = z6 = z8 = 0.;
	p0 = p2 = p4 = p6 = p8 = 0.;
	h0 = h2 = h4 = h6 = 0.;
	a0 = d0 = d2 = e0 = e2 = e4 = 0.;

	yq = y * y;
	yrrtpi = y * .56418958; // 1/SQRT(pi)
	/* Region boundaries */
	xlim0 = R0 - y;
	xlim1 = R1 - y;
	xlim3 = y * 3.097f - .45f;
	xlim2 = 6.8f - y;
	xlim4 = y * 18.1f + 1.65f;
	if (y <= 1e-6f) 
	{ /* avoid W4 algorithm */
		xlim1 = xlim0;
		xlim2 = xlim0;
	}

	for (i = 0; i < n; ++i) 
	{
		abx = fabs(x[i]);
		xq = abx * abx;
		if (abx > xlim0) 
		{ /* Region 0 algorithm */
			k[i] = yrrtpi / (xq + yq);
		} 
		else if (abx > xlim1) 
		{	/* Humlicek W4 Region 1 */
			if (rg1) 
			{
				/* First point in Region 1 */
				rg1 = false;
				a0 = yq + .5f;
				d0 = a0 * a0;
				d2 = yq + yq - 1.f;
			}
			d = .56418958f / (d0 + xq * (d2 + xq));
			k[i] = d * y * (a0 + xq);
		} 
		else if (abx > xlim2) 
		{ /* Humlicek W4 Region 2 */
			if (rg2) 
			{ /* First point in Region 2 */
				rg2 = false;
				h0 = yq * (yq * (yq * (yq + 6.f) + 10.5f) + 4.5f) + .5625f;
				h2 = yq * (yq * (yq * 4.f + 6.f) + 9.f) - 4.5f;
				h4 = 10.5f - yq * (6.f - yq * 6.f);
				h6 = yq * 4. - 6.;
				e0 = yq * (yq * (yq + 5.5f) + 8.25f) + 1.875f;
				e2 = yq * (yq * 3.f + 1.f) + 5.25f;
				e4 = h6 * .75f;
			}
			d = .56418958f / (h0 + xq * (h2 + xq * (h4 + xq * (h6 + xq))));
			k[i] = d * y * (e0 + xq * (e2 + xq * (e4 + xq)));
		} 
		else if (abx < xlim3) 
		{ /* Humlicek W4 Region 3 */
			if (rg3) 
			{
				/* First point in Region 3 */
				rg3 = false;
				z0 = y * (y * (y * (y * (y * (y * 
														(y * (y * (y * (y + 13.3988f) 
																	  + 88.26741f) 
																+ 369.1989f) + 1074.409f) + 
														2256.981f) 
												 + 3447.629f) + 3764.966f) + 2802.87f) + 
							 1280.829f) + 272.1014f;
				z2 = y * (y * (y * (y * (y * 
												 (y * (y * (y * 5.f + 53.59518f)
														 + 266.2987f) + 793.4273f) 
												 + 1549.675f) + 
										  2037.31f) + 1758.336f) + 902.3066f) + 211.678f;
				z4 = y * (y * (y * (y * (y * (y * 10.f + 80.39278f) + 
												 269.2916f) + 479.2576f) + 497.3014f) 
							 + 308.1852f) + 78.86585f;
				z6 = y * (y * (y * (y * 10.f + 53.59518f) + 92.75679f) + 
							 55.02933f) + 22.03523f;
				z8 = y * (y * 5.f + 13.3988f) + 1.49646f;
				p0 = y * (y * (y * (y * (y * (y * 
														(y * (y * (y * 
																	  .3183291f + 4.264678f) 
																+ 27.93941f) + 115.3772f) + 
														328.2151f) + 662.8097f) 
										  + 946.897f) + 919.4955f) + 
							 549.3954f) + 153.5168f;
				p2 = y * (y * (y * (y * (y * (y * (y * 1.2733163f + 
															  12.79458f) + 56.81652f) 
												 + 139.4665f) + 189.773f) + 
									124.5975f) - 1.322256f) - 34.16955f;
				p4 = y * (y * (y * (y * (y * 1.9099744f + 12.79568f) + 
										  29.81482f) + 24.01655f)
							 + 10.46332f) + 2.584042f;
				p6 = y * (y * (y * 1.273316f + 4.266322f) + .9377051f) - 
					.07272979f;
				p8 = y * .3183291f + 5.480304e-4f;
			}
			d = 1.7724538f / (z0 + xq * (z2 + xq * (z4 + xq * 
																 (z6 + xq * ( z8 + xq)))));
			k[i] = d * (p0 + xq * (p2 + xq * (p4 + xq * (p6 + xq * p8))));
		} 
		else 
		{ /* Humlicek CPF12 algorithm */
			ypy0 = y + 1.5f;
			ypy0q = ypy0 * ypy0;
			k[i] = 0.f;
			for (j = 0; j <= 5; ++j) 
			{
				d = x[i] - t[j];
				mq[j] = d * d;
				mf[j] = 1.f / (mq[j] + ypy0q);
				xm[j] = mf[j] * d;
				ym[j] = mf[j] * ypy0;
				d = x[i] + t[j];
				pq[j] = d * d;
				pf[j] = 1. / (pq[j] + ypy0q);
				xp[j] = pf[j] * d;
				yp[j] = pf[j] * ypy0;
			}
			if (abx <= xlim4) 
			{ /* Humlicek CPF12 Region I */
				for (j = 0; j <= 5; ++j) 
				{
					k[i] = k[i] + c[j] * (ym[j] + yp[j]) - s[j] * (xm[j] - xp[j]);
				}
			} 
			else 
			{ /* Humlicek CPF12 Region II */
				yf = y + 3.f;
				for (j = 0; j <= 5; ++j) 
				{
					k[i] = k[i] + (c[j] * (mq[j] * mf[j] - ym[j] * 1.5f)
										+ s[j] * yf * xm[j]) / (mq[j] + 2.25f) + 
						(c[j] * (pq[j] * pf[j] - yp[j] * 1.5f) - s[j] * yf * xp[j]) 
						/ (pq[j] + 2.25f);
				}
				k[i] = y * k[i] + exp(-xq);
			}
		}
	}
}

/******************************************************************
 * This marks the end of the block of code for the Voigt function *
 ******************************************************************/

STATIC uint32 MD5swap( uint32 word );
STATIC void MD5_Transform (uint32 *digest, const uint32 *in);

//
// The routines MD5file(), MD5textfile(), MD5string() and MD5swap() were written by Peter van Hoof
//
// this version returns the md5sum of a file and is identical to the well known md5sum algorithm
string MD5file(const char* fnam, access_scheme scheme)
{
	DEBUG_ENTRY( "MD5file()" );

	fstream ioFile;
	open_data( ioFile, fnam, mode_r, scheme );

	char c;
	string content;
	while( ioFile.get(c) )
		content += c;

	return MD5string( content );
}


// this version returns the md5sum of a text file. It filters out the eol characters,
// which makes it incompatible with the md5sum algorithm, but it is OS agnostic...
// comment lines that start with the hash symbol are also skipped
string MD5datafile(const char* fnam, access_scheme scheme)
{
	DEBUG_ENTRY( "MD5datafile()" );

	fstream ioFile;
	open_data( ioFile, fnam, mode_r, scheme );

	string line, content;
	while( getline( ioFile, line ) )
		if( line[0] != '#' )
			content += line;

	return MD5string( content );
}


// calculate the md5sum of an arbitrary string
string MD5string(const string& str)
{
	DEBUG_ENTRY( "MD5string()" );

	uint32 state[4];

	state[0] = 0x67452301L;
	state[1] = 0xefcdab89L;
	state[2] = 0x98badcfeL;
	state[3] = 0x10325476L;

	string lstr = str;

	// pad the string following RFC 1321 3.1 Step 1.
	uint32 bytes = str.length()%64;
	uint32 padlen = ( bytes >= 56 ) ? 64 + 56 - bytes : 56 - bytes;
	lstr += '\x80';
	for( uint32 i=1; i < padlen; ++i )
		lstr += '\0';

	ASSERT( lstr.length()%64 == 56 );

	uint32 i;
	union {
		uint32 in[16];
		unsigned char chr[64];
	} u;

	for( i=0; i < lstr.length()/64; ++i )
	{
		for( uint32 j=0; j < 64; ++j )
		{
			if( cpu.little_endian() )
				u.chr[j] = lstr[i*64+j];
			else
			{
				uint32 jr = j%4;
				uint32 j4 = j-jr;
				u.chr[j4+3-jr] = lstr[i*64+j];
			}
		}
		MD5_Transform( state, u.in );
	}
	for( uint32 j=0; j < 56; ++j )
	{
		if( cpu.little_endian() )
			u.chr[j] = lstr[i*64+j];
		else
		{
			uint32 jr = j%4;
			uint32 j4 = j-jr;
			u.chr[j4+3-jr] = lstr[i*64+j];
		}
	}
	// append the length of the string in _bits_ following RFC 1321 3.1 Step 2.
	u.in[14] = (str.length()<<3)&0xffffffff;
	u.in[15] = (str.length()>>29)&0xffffffff;

	MD5_Transform( state, u.in );

	ostringstream hash;
	for( uint32 i=0; i < 4; ++i )
		hash << hex << setfill('0') << setw(8) << MD5swap(state[i]);

	return hash.str();
}

STATIC uint32 MD5swap( uint32 word )
{
	DEBUG_ENTRY( "MD5swap()" );

	union {
		uint32 word;
		unsigned char byte[4];
	} ui, uo;

	ui.word = word;
	uo.byte[0] = ui.byte[3];
	uo.byte[1] = ui.byte[2];
	uo.byte[2] = ui.byte[1];
	uo.byte[3] = ui.byte[0];

	return uo.word;
}

//
// The following implementation of the MD5 algorithm was taken from the
// Crypto++ library (http://www.cryptopp.com/) and is subject to the
// following license:
//
//
// Compilation Copyright (c) 1995-2010 by Wei Dai.  All rights reserved.
// This copyright applies only to this software distribution package 
// as a compilation, and does not imply a copyright on any particular 
// file in the package.
// 
// All individual files in this compilation are placed in the public domain by
// Wei Dai and other contributors.
// 
// I would like to thank the following authors for placing their works into
// the public domain:
// 
// Joan Daemen - 3way.cpp
// Leonard Janke - cast.cpp, seal.cpp
// Steve Reid - cast.cpp
// Phil Karn - des.cpp
// Andrew M. Kuchling - md2.cpp, md4.cpp
// Colin Plumb - md5.cpp
// Seal Woods - rc6.cpp
// Chris Morgan - rijndael.cpp
// Paulo Baretto - rijndael.cpp, skipjack.cpp, square.cpp
// Richard De Moliner - safer.cpp
// Matthew Skala - twofish.cpp
// Kevin Springle - camellia.cpp, shacal2.cpp, ttmac.cpp, whrlpool.cpp, ripemd.cpp
// 
// Permission to use, copy, modify, and distribute this compilation for
// any purpose, including commercial applications, is hereby granted
// without fee, subject to the following restrictions:
// 
// 1. Any copy or modification of this compilation in any form, except
// in object code form as part of an application software, must include
// the above copyright notice and this license.
// 
// 2. Users of this software agree that any modification or extension
// they provide to Wei Dai will be considered public domain and not
// copyrighted unless it includes an explicit copyright notice.
// 
// 3. Wei Dai makes no warranty or representation that the operation of the
// software in this compilation will be error-free, and Wei Dai is under no
// obligation to provide any services, by way of maintenance, update, or
// otherwise.  THE SOFTWARE AND ANY DOCUMENTATION ARE PROVIDED "AS IS"
// WITHOUT EXPRESS OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE. IN NO EVENT WILL WEI DAI OR ANY OTHER CONTRIBUTOR BE LIABLE FOR
// DIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
// 
// 4. Users will not use Wei Dai or any other contributor's name in any 
// publicity or advertising, without prior written consent in each case.
// 
// 5. Export of this software from the United States may require a
// specific license from the United States Government.  It is the
// responsibility of any person or organization contemplating export
// to obtain such a license before exporting.
// 
// 6. Certain parts of this software may be protected by patents.  It
// is the users' responsibility to obtain the appropriate
// licenses before using those parts.
// 
// If this compilation is used in object code form in an application
// software, acknowledgement of the author is not required but would be
// appreciated. The contribution of any useful modifications or extensions
// to Wei Dai is not required but would also be appreciated.
// 

// md5.cpp - modified by Wei Dai from Colin Plumb's public domain md5.c
// any modifications are placed in the public domain

inline uint32 rotlFixed(uint32 x, unsigned int y)
{
	return uint32((x<<y) | (x>>(32-y)));
}

STATIC void MD5_Transform (uint32 *digest, const uint32 *in)
{
	DEBUG_ENTRY( "MD5_Transform()" );

#define F1(x, y, z) (z ^ (x & (y ^ z)))
#define F2(x, y, z) F1(z, x, y)
#define F3(x, y, z) (x ^ y ^ z)
#define F4(x, y, z) (y ^ (x | ~z))

#define MD5STEP(f, w, x, y, z, data, s) \
	w = rotlFixed(w + f(x, y, z) + data, s) + x

	uint32 a, b, c, d;

	a=digest[0];
	b=digest[1];
	c=digest[2];
	d=digest[3];

	MD5STEP(F1, a, b, c, d, in[0] + 0xd76aa478, 7);
	MD5STEP(F1, d, a, b, c, in[1] + 0xe8c7b756, 12);
	MD5STEP(F1, c, d, a, b, in[2] + 0x242070db, 17);
	MD5STEP(F1, b, c, d, a, in[3] + 0xc1bdceee, 22);
	MD5STEP(F1, a, b, c, d, in[4] + 0xf57c0faf, 7);
	MD5STEP(F1, d, a, b, c, in[5] + 0x4787c62a, 12);
	MD5STEP(F1, c, d, a, b, in[6] + 0xa8304613, 17);
	MD5STEP(F1, b, c, d, a, in[7] + 0xfd469501, 22);
	MD5STEP(F1, a, b, c, d, in[8] + 0x698098d8, 7);
	MD5STEP(F1, d, a, b, c, in[9] + 0x8b44f7af, 12);
	MD5STEP(F1, c, d, a, b, in[10] + 0xffff5bb1, 17);
	MD5STEP(F1, b, c, d, a, in[11] + 0x895cd7be, 22);
	MD5STEP(F1, a, b, c, d, in[12] + 0x6b901122, 7);
	MD5STEP(F1, d, a, b, c, in[13] + 0xfd987193, 12);
	MD5STEP(F1, c, d, a, b, in[14] + 0xa679438e, 17);
	MD5STEP(F1, b, c, d, a, in[15] + 0x49b40821, 22);

	MD5STEP(F2, a, b, c, d, in[1] + 0xf61e2562, 5);
	MD5STEP(F2, d, a, b, c, in[6] + 0xc040b340, 9);
	MD5STEP(F2, c, d, a, b, in[11] + 0x265e5a51, 14);
	MD5STEP(F2, b, c, d, a, in[0] + 0xe9b6c7aa, 20);
	MD5STEP(F2, a, b, c, d, in[5] + 0xd62f105d, 5);
	MD5STEP(F2, d, a, b, c, in[10] + 0x02441453, 9);
	MD5STEP(F2, c, d, a, b, in[15] + 0xd8a1e681, 14);
	MD5STEP(F2, b, c, d, a, in[4] + 0xe7d3fbc8, 20);
	MD5STEP(F2, a, b, c, d, in[9] + 0x21e1cde6, 5);
	MD5STEP(F2, d, a, b, c, in[14] + 0xc33707d6, 9);
	MD5STEP(F2, c, d, a, b, in[3] + 0xf4d50d87, 14);
	MD5STEP(F2, b, c, d, a, in[8] + 0x455a14ed, 20);
	MD5STEP(F2, a, b, c, d, in[13] + 0xa9e3e905, 5);
	MD5STEP(F2, d, a, b, c, in[2] + 0xfcefa3f8, 9);
	MD5STEP(F2, c, d, a, b, in[7] + 0x676f02d9, 14);
	MD5STEP(F2, b, c, d, a, in[12] + 0x8d2a4c8a, 20);

	MD5STEP(F3, a, b, c, d, in[5] + 0xfffa3942, 4);
	MD5STEP(F3, d, a, b, c, in[8] + 0x8771f681, 11);
	MD5STEP(F3, c, d, a, b, in[11] + 0x6d9d6122, 16);
	MD5STEP(F3, b, c, d, a, in[14] + 0xfde5380c, 23);
	MD5STEP(F3, a, b, c, d, in[1] + 0xa4beea44, 4);
	MD5STEP(F3, d, a, b, c, in[4] + 0x4bdecfa9, 11);
	MD5STEP(F3, c, d, a, b, in[7] + 0xf6bb4b60, 16);
	MD5STEP(F3, b, c, d, a, in[10] + 0xbebfbc70, 23);
	MD5STEP(F3, a, b, c, d, in[13] + 0x289b7ec6, 4);
	MD5STEP(F3, d, a, b, c, in[0] + 0xeaa127fa, 11);
	MD5STEP(F3, c, d, a, b, in[3] + 0xd4ef3085, 16);
	MD5STEP(F3, b, c, d, a, in[6] + 0x04881d05, 23);
	MD5STEP(F3, a, b, c, d, in[9] + 0xd9d4d039, 4);
	MD5STEP(F3, d, a, b, c, in[12] + 0xe6db99e5, 11);
	MD5STEP(F3, c, d, a, b, in[15] + 0x1fa27cf8, 16);
	MD5STEP(F3, b, c, d, a, in[2] + 0xc4ac5665, 23);

	MD5STEP(F4, a, b, c, d, in[0] + 0xf4292244, 6);
	MD5STEP(F4, d, a, b, c, in[7] + 0x432aff97, 10);
	MD5STEP(F4, c, d, a, b, in[14] + 0xab9423a7, 15);
	MD5STEP(F4, b, c, d, a, in[5] + 0xfc93a039, 21);
	MD5STEP(F4, a, b, c, d, in[12] + 0x655b59c3, 6);
	MD5STEP(F4, d, a, b, c, in[3] + 0x8f0ccc92, 10);
	MD5STEP(F4, c, d, a, b, in[10] + 0xffeff47d, 15);
	MD5STEP(F4, b, c, d, a, in[1] + 0x85845dd1, 21);
	MD5STEP(F4, a, b, c, d, in[8] + 0x6fa87e4f, 6);
	MD5STEP(F4, d, a, b, c, in[15] + 0xfe2ce6e0, 10);
	MD5STEP(F4, c, d, a, b, in[6] + 0xa3014314, 15);
	MD5STEP(F4, b, c, d, a, in[13] + 0x4e0811a1, 21);
	MD5STEP(F4, a, b, c, d, in[4] + 0xf7537e82, 6);
	MD5STEP(F4, d, a, b, c, in[11] + 0xbd3af235, 10);
	MD5STEP(F4, c, d, a, b, in[2] + 0x2ad7d2bb, 15);
	MD5STEP(F4, b, c, d, a, in[9] + 0xeb86d391, 21);

	digest[0] += a;
	digest[1] += b;
	digest[2] += c;
	digest[3] += d;
}

/***************************************************************
 * This marks the end of the block of code from Wei Dai et al. *
 ***************************************************************/
