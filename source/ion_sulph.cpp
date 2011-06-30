/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonSulph compute ionization balance for sulphur */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonSulph(void)
{
	const int NDIM = ipSULPHUR+1;

	static const double dicoef[2][NDIM] = {
		{1.62e-3,1.09e-2,3.35e-2,3.14e-2,1.27e-2,1.47e-2,1.34e-2,
		 2.38e-2,3.19e-2,7.13e-2,.08,.0796,.0134,.402,.241,0.},
		{0.,1.20e-2,6.59e-2,6.89e-2,.187,.129,1.04,1.12,1.40,1.00,
		 .555,1.63,.304,.298,.281,0.}
	};
	static const double dite[2][NDIM] = {
		{1.25e5,1.92e5,1.89e5,1.68e5,1.38e5,1.80e6,6.90e5,5.84e5,
		 5.17e5,6.66e5,6.00e5,5.09e5,2.91e5,2.41e7,2.54e7,0.},
		{0.,1.80e4,1.59e5,8.04e4,1.71e5,1.75e6,2.15e6,2.59e6,2.91e6,
		 2.32e6,2.41e6,6.37e6,1.04e6,4.67e6,5.30e6,0.}
	};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.};
	static const double ditcrt[NDIM] = {2.2e4,1.2e4,1.4e4,1.5e4,1.4e4,2.9e5,
	  1.3e5,1.1e5,9.0e4,9.0e4,9.0e4,8.3e4,6.0e4,5.0e6,9.0e6,1e20};

	bool lgPrtDebug=false;

	DEBUG_ENTRY( "IonSulph()" );

	/* sulphur, atomic number 16 */
	if( !dense.lgElmtOn[ipSULPHUR] )
	{
		return;
	}

	ion_zero(ipSULPHUR);

	ion_photo(ipSULPHUR,lgPrtDebug);

	/* find collisional ionization rates */
	ion_collis(ipSULPHUR);

	/* get recombination coefficients */
	/*lint -e64 type mismatch */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipSULPHUR);
	/*lint +e64 type mismatch */

	/* solve for ionization balance */
	/*if( nzone>700 ) lgPrtDebug = true;*/
	ion_solver(ipSULPHUR,lgPrtDebug);

	return;
}
