/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*phfit derive photoionization cross sections for first 30 elements */
#include "cddefines.h"
#include "physconst.h"
#include "atmdat.h"
#include "iso.h"
#include "taulines.h"

/** constructor: read in all the ADfA data files */
t_ADfA::t_ADfA()
{
	DEBUG_ENTRY( "t_ADfA()" );

	/* this option, use the new atmdat_rad_rec recombination rates */
	version = PHFIT_UNDEF;

	double help[9];
	const long VERSION_MAGIC = 20061204L;

	static const char chFile[] = "phfit.dat";

	FILE *io = open_data( chFile, "r" );

	bool lgErr = false;
	long i=-1, j=-1, k=-1, n;

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( n=0; n < 7; n++ )
		lgErr = lgErr || ( fscanf( io, "%ld", &L[n] ) != 1 );
	for( n=0; n < 30; n++ )
		lgErr = lgErr || ( fscanf( io, "%ld", &NINN[n] ) != 1 );
	for( n=0; n < 30; n++ )
		lgErr = lgErr || ( fscanf( io, "%ld", &NTOT[n] ) != 1 );
	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld %ld", &i, &j, &k ) != 3 );
		if( i == -1 && j == -1 && k == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le", &help[0], &help[1],
					   &help[2], &help[3], &help[4], &help[5] ) != 6 );
		for( int l=0; l < 6; ++l )
			PH1[i][j][k][l] = (realnum)help[l];
	}
	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le", &help[0], &help[1],
					   &help[2], &help[3], &help[4], &help[5], &help[6] ) != 7 );
		for( int l=0; l < 7; ++l )
			PH2[i][j][l]  = (realnum)help[l];
	}
	fclose( io );

	ASSERT( !lgErr );

	static const char chFile2[] = "hpfit.dat";

	io = open_data( chFile2, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile2, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < NHYDRO_MAX_LEVEL; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le", &help[0], &help[1],
					   &help[2], &help[3], &help[4] ) != 5 );
		for( int l=0; l < 5; ++l )
			PHH[i][l] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile3[] = "rec_lines.dat";

	io = open_data( chFile3, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile3, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < 110; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le", &help[0], &help[1], &help[2],
					   &help[3], &help[4], &help[5], &help[6], &help[7] ) != 8 );
		for( int l=0; l < 8; ++l )
			P[l][i] = (realnum)help[l];
	}


	for( i=0; i < 405; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le %le", &help[0], &help[1], &help[2],
					   &help[3], &help[4], &help[5], &help[6], &help[7], &help[8] ) != 9 );
		for( int l=0; l < 9; ++l )
			ST[l][i] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile4[] = "rad_rec.dat";

	io = open_data( chFile4, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile4, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le", &help[0], &help[1] ) != 2 );
		for( int l=0; l < 2; ++l )
			rrec[i][j][l] = (realnum)help[l];
	}
	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le", &help[0], &help[1],
					   &help[2], &help[3] ) != 4 );
		for( int l=0; l < 4; ++l )
			rnew[i][j][l] = (realnum)help[l];
	}
	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
		if( i == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le", &help[0], &help[1], &help[2] ) != 3 );
		for( int l=0; l < 3; ++l )
			fe[i][l] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile5[] = "h_rad_rec.dat";

	io = open_data( chFile5, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile5, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < NHYDRO_MAX_LEVEL; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le %le", &help[0], &help[1], &help[2],
					   &help[3], &help[4], &help[5], &help[6], &help[7], &help[8] ) != 9 );
		for( int l=0; l < 9; ++l )
			HRF[i][l] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile6[] = "h_phot_cs.dat";

	io = open_data( chFile6, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile6, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < NHYDRO_MAX_LEVEL; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le", &help[0] ) != 1 );
		STH[i] = (realnum)help[0];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile7[] = "coll_ion.dat";

	io = open_data( chFile7, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile7, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le", &CF[i][j][0], &CF[i][j][1],
					   &CF[i][j][2], &CF[i][j][3], &CF[i][j][4] ) != 5 );
	}

	fclose( io );

	ASSERT( !lgErr );

	/*refer	HI	cs	Anderson, H., Ballance, C.P., Badnell, N.R., 
	 *refercon	& Summers, H.P  2000, J Phys B, 33, 1255 */
	static const char chFile8[] = "h_coll_str.dat";

	io = open_data( chFile8, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile8, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le", &HCS[i-2][j-1][0], &HCS[i-2][j-1][1],
					   &HCS[i-2][j-1][2], &HCS[i-2][j-1][3], &HCS[i-2][j-1][4], &HCS[i-2][j-1][5],
					   &HCS[i-2][j-1][6], &HCS[i-2][j-1][7] ) != 8 );
	}

	fclose( io );

	ASSERT( !lgErr );
}

double t_ADfA::phfit(long int nz, 
		     long int ne,
		     long int is, 
		     double e)
{
	long int nint, 
	  nout;
	double a, 
	  b, 
	  crs, 
	  einn, 
	  p1, 
	  q, 
	  x, 
	  y, 
	  z;

	DEBUG_ENTRY( "phfit()" );

	/*** Version 3. October 8, 1996.
	 *** Written by D. A. Verner, verner@pa.uky.edu
	 *** Inner-shell ionization energies of some low-ionized species are slightly
	 *** improved to fit smoothly the experimental inner-shell ionization energies 
	 *** of neutral atoms.
	 ******************************************************************************
	 *** This subroutine calculates partial photoionization cross sections
	 *** for all ionization stages of all atoms from H to Zn (Z=30) by use of
	 *** the following fit parameters:
	 *** Outer shells of the Opacity Project (OP) elements:
	 *** >>refer	all	photocs	Verner, Ferland, Korista, Yakovlev, 1996, ApJ, in press.
	 *** Inner shells of all elements, and outer shells of the non-OP elements:
	 ***  Verner and Yakovlev, 1995, A&AS, 109, 125
	 *** Input parameters:  nz - atomic number from 1 to 30 (integer) 
	 ***          ne - number of electrons from 1 to iz (integer)
	 ***          is - shell number (integer)
	 ***          e - photon energy, eV 
	 ***          version - enum, PHFIT96 (default): calculates 
	 ***                 new cross sections, PHFIT95: calculates
	 ***                 only old Hartree-Slater cross sections
	 *** Output parameter:  photoionization cross section, Mb
	 *** Shell numbers:
	 *** 1 - 1s, 2 - 2s, 3 - 2p, 4 - 3s, 5 - 3p, 6 - 3d, 7 - 4s. 
	 *** If a species in the ground state has no electrons on the given shell,
	 *** the subroutine returns 0.
	 ****************************************************************************** */

	crs = 0.0;
	if( nz < 1 || nz > 30 )
	{ 
		return(crs);
	}

	if( ne < 1 || ne > nz )
	{ 
		return(crs);
	}

	nout = NTOT[ne-1];
	if( nz == ne && nz > 18 )
		nout = 7;
	if( nz == (ne + 1) && ((((nz == 20 || nz == 21) || nz == 22) || 
	  nz == 25) || nz == 26) )
		nout = 7;
	if( is > nout )
	{ 
		return(crs);
	}

	if( (is == 6 && (nz == 20 || nz == 19)) && ne >= 19 )
	{ 
		return(crs);
	}

	ASSERT( is >= 1 && is <= 7 );

	if( e < PH1[is-1][ne-1][nz-1][0] )
	{ 
		return(crs);
	}

	nint = NINN[ne-1];
	if( ((nz == 15 || nz == 17) || nz == 19) || (nz > 20 && nz != 26) )
	{
		einn = 0.0;
	}
	else
	{
		if( ne < 3 )
		{
			einn = 1.0e30;
		}
		else
		{
			einn = PH1[nint-1][ne-1][nz-1][0];
		}
	}

	if( (is <= nint || e >= einn) || version == PHFIT95 )
	{
		p1 = -PH1[is-1][ne-1][nz-1][4];
		y = e/PH1[is-1][ne-1][nz-1][1];
		q = -0.5*p1 - L[is-1] - 5.5;
		a = PH1[is-1][ne-1][nz-1][2]*(POW2(y - 1.0) + 
			POW2(PH1[is-1][ne-1][nz-1][5]));
		b = sqrt(y/PH1[is-1][ne-1][nz-1][3]) + 1.0;
		crs = a*pow(y,q)*pow(b,p1);
	}
	else
	{
		if( (is < nout && is > nint) && e < einn )
		{ 
			return(crs);
		}
		p1 = -PH2[ne-1][nz-1][3];
		q = -0.5*p1 - 5.5;
		x = e/PH2[ne-1][nz-1][0] - PH2[ne-1][nz-1][5];
		z = sqrt(x*x+POW2(PH2[ne-1][nz-1][6]));
		a = PH2[ne-1][nz-1][1]*(POW2(x - 1.0) + 
			POW2(PH2[ne-1][nz-1][4]));
		b = 1.0 + sqrt(z/PH2[ne-1][nz-1][2]);
		crs = a*pow(z,q)*pow(b,p1);
	}
	return(crs);
}

double t_ADfA::hpfit(long int iz,
		     long int n,
		     double e)
{
	long int l, 
	  m;
	double cs,
	  eth, 
	  ex, 
	  q, 
	  x;

	DEBUG_ENTRY( "hpfit()" );

	/*state specific photoionization cross sections for model hydrogen atom
	 * Version 1, September 23, 1997
	 ******************************************************************************
	 *** This subroutine calculates state-specific photoionization cross sections
	 *** for hydrogen and hydrogen-like ions.
	 *** Input parameters:  iz - atomic number 
	 ***          n  - shell number, from 0 to 400:
	 ***                                    0 - 1s
	 ***                                    1 - 2s
	 ***                                    2 - 2p
	 ***                                    3 - 3 
	 ***                                    ......
	 ***          e  - photon energy, eV
	 *** return value - cross section, cm^(-2)     
	 *******************************************************************************/

	cs = 0.0;

	ASSERT( iz > 0 && e>0. );

	if( n >= NHYDRO_MAX_LEVEL )
	{ 
		fprintf( ioQQQ, " hpfit called with too large n, =%li\n" , n );
		cdEXIT(EXIT_FAILURE);
	}

	l = 0;
	if( n == 2 )
	{
		l = 1;
	}
	q = 3.5 + l - 0.5*PHH[n][1];

	if( n == 0 )
	{
		m = 1;
	}
	else
	{
		if( n == 1 )
		{
			m = 2;
		}
		else
		{
			m = n;
		}
	}

	eth = ph1(0,0,iz-1,0)/POW2((double)m);
	ex = MAX2(1. , e/eth );

	/* Don't just force to be at least one...make sure e/eth is close to one or greater.	*/
	ASSERT( e/eth > 0.95 );

	if( ex < 1.0 )
	{ 
		return(0.);
	}

	x = ex/PHH[n][0];
	cs = (PHH[n][4]*pow(1.0 + ((double)PHH[n][2]/x),(double)PHH[n][3])/
	  pow(x,q)/pow(1.0 + sqrt(x),(double)PHH[n][1])*8.79737e-17/
	  POW2((double)iz));
	return(cs);
}

void t_ADfA::rec_lines(double t, 
		       realnum r[][471])
{
	long int i, 
	  j, 
	  ipj1, 
	  ipj2;

	double a, 
	  a1, 
	  dr[4][405], 
	  p1, 
	  p2, 
	  p3, 
	  p4, 
	  p5, 
	  p6, 
	  rr[4][110], 
	  te, 
	  x, 
	  z;

	static long jd[6]={143,145,157,360,376,379};

	static long ip[38]={7,9,12,13,14,16,18,19,20,21,22,44,45,49,50,
	  52,53,54,55,56,57,58,59,60,66,67,78,83,84,87,88,95,96,97,100,
	  101,103,104};

	static long id[38]={7,3,10,27,23,49,51,64,38,47,60,103,101,112,
	  120,114,143,145,157,152,169,183,200,163,225,223,237,232,235,
	  249,247,300,276,278,376,360,379,384};

	DEBUG_ENTRY( "rec_lines()" );

	/*effective recombination coefficients for lines of C, N, O, by D. Verner
	 * Version 2, April 30, 1997
	 ******************************************************************************
	 *** This subroutine calculates effective recombination coefficients
	 *** for 110 permitted recombination lines of C, N, O (Pequignot, Petitjean,
	 *** & Boisson, 1991, A&A, 251, 680) and 405 permitted dielectronic 
	 *** recombination lines (Nussbaumer & Storey, 1984, A&AS, 56, 293)
	 *** Input parameter:   t  - temperature, K
	 *** Output parameters: r(i,j), i=1,471
	 ***          r(i,1) - atomic number
	 ***          r(i,2) - number of electrons
	 ***          r(i,3) - wavelength, angstrom
	 ***          r(i,4) - rate coefficient, cm^3 s^(-1)
	 ****************************************************************************** */

	for( i=0; i < 110; i++ )
	{
		rr[0][i] = P[0][i];
		rr[1][i] = P[1][i];
		rr[2][i] = P[2][i];
		z = P[0][i] - P[1][i] + 1.0;
		te = 1.0e-04*t/z/z;
		p1 = P[3][i];
		p2 = P[4][i];
		p3 = P[5][i];
		p4 = P[6][i];
		if( te < 0.004 )
		{
			a1 = p1*pow(0.004,p2)/(1.0 + p3*pow(0.004,p4));
			a = a1/sqrt(te/0.004);
		}
		else
		{
			if( te > 2.0 )
			{
				a1 = p1*pow(2.0,p2)/(1.0 + p3*pow(2.0,p4));
				a = a1/pow(te/2.0,1.5);
			}
			else
			{
				a = p1*pow(te,p2)/(1.0 + p3*pow(te,p4));
			}
		}
		rr[3][i] = 1.0e-13*z*a*P[7][i];
	}

	for( i=0; i < 405; i++ )
	{
		dr[0][i] = ST[0][i];
		dr[1][i] = ST[1][i];
		dr[2][i] = ST[2][i];
		te = 1.0e-04*t;
		p1 = ST[3][i];
		p2 = ST[4][i];
		p3 = ST[5][i];
		p4 = ST[6][i];
		p5 = ST[7][i];
		p6 = ST[8][i];
		if( te < p6 )
		{
			x = p5*(1.0/te - 1.0/p6);
			if( x > 80.0 )
			{
				a = 0.0;
			}
			else
			{
				a1 = (p1/p6 + p2 + p3*p6 + p4*p6*p6)/pow(p6,1.5)/exp(p5/
				  p6);
				a = a1/exp(x);
			}
		}
		else
		{
			if( te > 6.0 )
			{
				a1 = (p1/6.0 + p2 + p3*6.0 + p4*36.0)/pow(6.0,1.5)/
				  exp(p5/6.0);
				a = a1/pow(te/6.0,1.5);
			}
			else
			{
				a = (p1/te + p2 + p3*te + p4*te*te)/pow(te,1.5)/exp(p5/
				  te);
			}
		}
		dr[3][i] = 1.0e-12*a;
	}

	for( i=0; i < 6; i++ )
	{
		ipj1 = jd[i];
		ipj2 = ipj1 + 1;
		dr[3][ipj1-1] += dr[3][ipj2-1];
		dr[0][ipj2-1] = 0.0;
	}

	for( i=0; i < 38; i++ )
	{
		ipj1 = ip[i];
		ipj2 = id[i];
		rr[3][ipj1-1] += dr[3][ipj2-1];
		dr[0][ipj2-1] = 0.0;
	}

	for( i=0; i < 110; i++ )
	{
		r[0][i] = (realnum)rr[0][i];
		r[1][i] = (realnum)rr[1][i];
		r[2][i] = (realnum)rr[2][i];
		r[3][i] = (realnum)rr[3][i];
	}

	j = 110;
	for( i=0; i < 405; i++ )
	{
		if( dr[0][i] > 1.0 )
		{
			j += 1;
			r[0][j-1] = (realnum)dr[0][i];
			r[1][j-1] = (realnum)dr[1][i];
			r[2][j-1] = (realnum)dr[2][i];
			r[3][j-1] = (realnum)dr[3][i];
		}
	}
	return;
}

double t_ADfA::rad_rec(long int iz, 
		       long int in, 
		       double t)
{
	/*
	 *** Version 4. June 29, 1999.
	 *** Written by D. A. Verner, verner@pa.uky.edu 
	 ******************************************************************************
	 *** This subroutine calculates rates of radiative recombination for all ions
	 *** of all elements from H through Zn by use of the following fits:
	 *** H-like, He-like, Li-like, Na-like - 
	 *** >>refer	all	reccoef	Verner & Ferland, 1996, ApJS, 103, 467
	 *** Other ions of C, N, O, Ne - Pequignot et al. 1991, A&A, 251, 680,
	 ***    refitted by Verner & Ferland formula to ensure correct asymptotes
	 *** Fe XVII-XXIII - 
	 *** >>refer	Fe17-23	recom	Arnaud & Raymond, 1992, ApJ, 398, 394
	 *** Fe I-XV - refitted by Verner & Ferland formula to ensure correct asymptotes
	 *** Other ions of Mg, Si, S, Ar, Ca, Fe, Ni - 
	 ***                      -
	 *** >>refer	all	recom	Shull & Van Steenberg, 1982, ApJS, 48, 95
	 *** Other ions of Na, Al - 
	 *** >>refer	Na, Al	recom	Landini & Monsignori Fossi, 1990, A&AS, 82, 229
	 *** Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV,
	 *** Mn I-V, Co I)        - 
	 *** >>refer	many	recom	Landini & Monsignori Fossi, 1991, A&AS, 91, 183
	 *** All other species    - interpolations of the power-law fits
	 *** Input parameters:  iz - atomic number 
	 ***                    in - number of electrons from 1 to iz 
	 ***                    t  - temperature, K
	 *** return result:  - rate coefficient, cm^3 s^(-1)
	 ******************************************************************************
	 */
	double tt;
	double rate;

	DEBUG_ENTRY( "rad_rec()" );

	rate = 0.0;
	if( iz < 1 || iz > 30 )
	{
		fprintf( ioQQQ, " rad_rec called with insane atomic number, =%4ld\n", 
		  iz );
		cdEXIT(EXIT_FAILURE);
	}
	if( in < 1 || in > iz )
	{
		fprintf( ioQQQ, " rad_rec called with insane number elec =%4ld\n", 
		  in );
		cdEXIT(EXIT_FAILURE);
	}
	if( (((in <= 3 || in == 11) || (iz > 5 && iz < 9)) || iz == 10) || 
	  (iz == 26 && in > 11) )
	{
		tt = sqrt(t/rnew[in-1][iz-1][2]);
		rate = 
		  rnew[in-1][iz-1][0]/(tt*pow(tt + 1.0,1.0 - rnew[in-1][iz-1][1])*
		  pow(1.0 + sqrt(t/rnew[in-1][iz-1][3]),1.0 + rnew[in-1][iz-1][1]));
	}
	else
	{
		tt = t*1.0e-04;
		if( iz == 26 && in <= 13 )
		{
			rate = fe[in-1][0]/pow(tt,fe[in-1][1] + 
			  fe[in-1][2]*log10(tt));
		}
		else
		{
			rate = rrec[in-1][iz-1][0]/pow(tt,(double)rrec[in-1][iz-1][1]);
		}
	}

	return rate;
}

double t_ADfA::H_rad_rec(long int iz,
			 long int n,
			 double t)
{
	/*
	 * Version 4, October 9, 1997
	 ******************************************************************************
	 *** This subroutine calculates state-specific recombination rates 
	 *** for hydrogen and hydrogen-like ions.
	 *** Input parameters:  iz - atomic number 
	 ***          n  - shell number, from 0 to 400:
	 ***                                    0 - 1s
	 ***                                    1 - 2s
	 ***                                    2 - 2p
	 ***                                    3 - 3 
	 ***                                    ......
	 ***          t  - temperature, K
	 *** Output parameter:  r  - rate coefficient, cm^3 s^(-1)
	 *** If n is negative, the subroutine returns the total recombination 
	 *** rate coefficient
	 ******************************************************************************
	 */
	double rate,
	  TeScaled, 
	  x, 
	  x1, 
	  x2;

	DEBUG_ENTRY( "H_rad_rec()" );

	rate = 0.0;

	/* iz is charge, must be 1 or greater */
	ASSERT( iz > 0 );

	/* n is level number, must be less than dim or hydro vectors */
	ASSERT( n < NHYDRO_MAX_LEVEL );

	TeScaled = t/POW2((double)iz);

	if( n < 0 )
	{
		x1 = sqrt(TeScaled/3.148);
		x2 = sqrt(TeScaled/7.036e05);
		rate = 7.982e-11/x1/pow(1.0 + x1,0.252)/pow(1.0 + x2,1.748);
	}
	else
	{
		x = log10(TeScaled);
		rate = (HRF[n][0] + HRF[n][2]*x + HRF[n][4]*
		  x*x + HRF[n][6]*powi(x,3) + HRF[n][8]*powi(x,4))/
		  (1.0 + HRF[n][1]*x + HRF[n][3]*x*x + HRF[n][5]*
		  powi(x,3) + HRF[n][7]*powi(x,4));
		rate = pow(10.0,rate)/TeScaled;
	}
	rate *= iz;

	return rate;
}

/*coll_ion D Verner's routine to compute collisional ionization rate coefficients,
 * returns collisional ionization rate coefficient cm^3 s^-1*/
double t_ADfA::coll_ion(
	/* atomic number, 1 for hydrogen */
	long int iz, 
	/* stage of ionization, 1 for atom */
	long int in, 
	/* temperature */
	double t)
{
	double rate, te, u;

	DEBUG_ENTRY( "coll_ion()" );
	/*D Verner's routine to compute collisional ionization rate coefficients
	 * Version 3, April 21, 1997
	 * Cu (Z=29) and Zn (Z=30) are added (fits from Ni, correct thresholds).
	 ******************************************************************************
	 *** This subroutine calculates rates of direct collisional ionization 
	 *** for all ionization stages of all elements from H to Ni (Z=28)
	 *** by use of the fits from
	 *>>refer	all	collion	Voronov, G. S., 1997, At. Data Nucl. Data Tables, 65, 1
	 *** Input parameters:  iz - atomic number on pphysical scale, H is 1
	 ***          in - number of electrons from 1 to iz 
	 ***          t  - temperature, K
	 *** Output parameter:  c  - rate coefficient, cm^3 s^(-1)
	 ****************************************************************************** */
	rate = 0.0;

	if( iz < 1 || iz > 30 )
	{ 
		/* return zero rate is atomic number outside range of code */
		return( 0. );
	}

	if( in < 1 || in > iz )
	{ 
		/* return zero rate is ion stage is impossible */
		return( 0. );
	}

	te = t*EVRYD/TE1RYD;
	u = CF[in-1][iz-1][0]/te;
	if( u > 80.0 )
	{ 
		return( 0. );
	}

	rate = (CF[in-1][iz-1][2]*(1.0 + CF[in-1][iz-1][1]*
	  sqrt(u))/(CF[in-1][iz-1][3] + u)*pow(u,CF[in-1][iz-1][4])*
	  exp(-u));

	return(rate);
}

realnum t_ADfA::h_coll_str( long ipLo, long ipHi, long ipTe )
{
	realnum rate;

	DEBUG_ENTRY( "h_coll_str()" );

	ASSERT( ipLo < ipHi );

#	if !defined NDEBUG
	long ipISO = ipH_LIKE;
	long nelem = ipHYDROGEN;
#	endif
	ASSERT( N_(ipLo) < N_(ipHi) );
	ASSERT( N_(ipHi) <= 5 );

	rate = (realnum)HCS[ipHi-1][ipLo][ipTe];

	return rate;
}
