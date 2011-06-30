/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonSilic determine ionization balance of Silicon */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonSilic(void)
{
	const int NDIM = ipSILICON+1;

	static const double dicoef[2][NDIM] = {
		{1.10e-3,5.87e-3,5.03e-3,5.43e-3,8.86e-3,1.68e-2,2.49e-2,3.13e-2,
		 4.25e-2,6.18e-2,1.38e-2,.327,.189,0.},
		{0.,.753,.188,.450,0.,1.80,1.88,2.01,1.22,.303,1.42,.306,.286,0.}
	};
	static const double  dite[2][NDIM] = {
		{7.70e4,9.63e4,8.75e4,1.05e6,1.14e6,4.85e5,4.15e5,3.66e5,3.63e5,
		 3.88e5,2.51e5,1.88e7,1.99e7,0.},
		{0.,6.46e4,4.71e4,7.98e5,0.,1.03e6,1.91e6,2.11e6,2.14e6,1.12e6,
		 3.93e6,3.60e6,4.14e6,0.}
	};
	static const double ditcrt[NDIM] = {1.1e4,1.1e4,1.1e4,1.7e5,9.5e4,8.0e4,
	  7.4e4,6.8e4,6.6e4,6.5e4,4.5e4,3.7e6,6.3e6,1e20};
	static const double aa[NDIM] = {-0.0219,3.2163,0.1203,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.4364,-12.0571,-2.6900,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.0684,16.2118,19.1943,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {-0.0032,-0.5886,-0.1479,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.1342,0.5613,0.1118,0.1,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.};

	bool lgPrtDebug=false;

	DEBUG_ENTRY( "IonSilic()" );

	/* silicon, atomic number 14 */
	if( !dense.lgElmtOn[ipSILICON] )
		return;

	ion_zero(ipSILICON);

	ion_photo(ipSILICON,false);

	/* find collisional ionization rates */
	ion_collis(ipSILICON);

	/* get recombination coefficients */
	/*lint -e64 type mismatch */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipSILICON);
	/*lint +e64 type mismatch */

	/* solve for ionization balance */
	/*if( nzone>700 ) lgPrtDebug = true;*/
	ion_solver(ipSILICON,lgPrtDebug);

	return;
}
