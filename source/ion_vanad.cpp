/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonVanad do ionization balance for vanadium */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonVanad(void)
{
	const int NDIM = ipVANADIUM+1;

	static const double dicoef[2][NDIM] = {
		{1.58e-3,8.38e-3,1.54e-2,3.75e-2,0.117,0.254,0.291,0.150,0.140,0.100,
		 0.200,0.240,0.260,0.190,0.120,0.350,0.066,0.10,0.13,0.23,0.14,0.11,0.},
		{.456,.323,.310,.411,.359,.0975,.229,4.20,3.30,5.30,1.50,0.700,.600,
		 .5,1.,0.,7.8,6.3,5.5,3.6,4.9,1.6,0.}
	};
	static const double dite[2][NDIM] = {
		{6.00e4,1.94e5,3.31e5,4.32e5,6.28e5,7.50e5,7.73e5,2.62e5,2.50e5,2.57e5,2.84e5,
		 8.69e5,4.21e5,4.57e5,2.85e5,8.18e5,1.51e6,1.30e6,1.19e6,1.09e6,9.62e5,7.23e5,0.},
		{8.97e4,1.71e5,2.73e5,3.49e5,5.29e5,4.69e5,6.54e5,1.32e6,1.33e6,1.41e6,1.52e6,
		 1.51e6,1.82e6,1.84e6,2.31e6,0.,9.98e6,9.98e6,1.00e7,1.10e7,8.34e6,1.01e7,0.}
	};
	static const double ditcrt[NDIM] = {6e3,2e4,4e4,5e4,7e4,8e4,8e4,3e4,
	  3e4,3e4,3e4,9e4,4e4,5e4,3e4,9e5,2e5,2e5,2e5,2e5,1e5,7e4,1e20};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonVanad()" );

	/* chromium nelem=23
	 * based on iron
	 *
	 * Fe rates from Woods et al. Ap.J. 249, 399.
	 * rec from +23, 24 25 from Arnauld et al 85 */

	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipVANADIUM] )
	{
		return;
	}

	ion_zero(ipVANADIUM);

	ion_photo(ipVANADIUM,false);

	/* find collisional ionization rates */
	ion_collis(ipVANADIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipVANADIUM);

	/* solve for ionization balance */
	ion_solver(ipVANADIUM,false);
	return;
}
