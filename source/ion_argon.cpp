/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonArgon compute ionization balance of argon */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonArgon(void)
{
	const int NDIM = ipARGON+1;

	static const double dicoef[2][NDIM] = {
		{1.00e-3,1.10e-2,3.40e-2,6.85e-2,9.00e-2,6.35e-2,2.60e-2,1.70e-2,
		 2.10e-2,3.50e-2,4.30e-2,7.13e-2,9.60e-2,8.50e-2,1.70e-2,.476,.297,0.},
		{.005,.045,.057,.087,.0769,.140,.120,.1,1.92,1.66,1.67,1.40,1.31,
		 1.02,.245,.294,.277,0.}
	};
	static const double dite[2][NDIM] = {
		{3.20e5,2.90e5,2.39e5,2.56e5,2.50e5,2.10e5,1.80e5,2.70e6,8.30e5,
		 6.95e5,6.05e5,6.68e5,6.50e5,5.30e5,3.55e5,3.01e7,3.13e7,0.},
		{3.10e5,5.50e5,6.00e5,3.81e5,3.30e5,2.15e5,2.15e5,3.30e6,3.50e6,
		 3.60e6,3.80e6,2.90e6,3.60e6,2.80e6,1.10e6,6.05e6,6.54e6,0.}
	};
	static const double ditcrt[NDIM] = {2.5e4,3.0e4,2.5e4,2.5e4,1.8e4,1.8e4,2.2e4,
	  5.0e5,1.6e5,1.5e5,1.5e5,1.3e5,1.3e5,1.1e5,7.6e4,6.5e6,1.4e7,1e20};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonArgon()" );

	/* argon nelem=18
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */

	/* rec from +15, 16, 17 from Arnauld et al 85 */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipARGON] )
	{
		return;
	}

	ion_zero(ipARGON);

	ion_photo(ipARGON,false);

	/* find collisional ionization rates */
	ion_collis(ipARGON);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipARGON);

	/* solve for ionization balance */
	ion_solver(ipARGON,false);
	return;
}
