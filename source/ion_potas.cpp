/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonPotas compute ionization equilibrium for Potassium */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonPotas(void)
{
	const int NDIM = ipPOTASSIUM+1;

	static const double dicoef[2][NDIM] = {
		{3.28e-4,.0584,.122,.132,.133,.126,.139,.0955,.0402,.0419,.0257,
		 .0445,.0548,.0713,.0903,.110,.0205,.549,0.},
		{.0907,.110,.0174,.132,.114,.162,.0878,.263,.0627,.0616,2.77,2.23,
		 2.00,1.82,.424,.243,.185,.292,0.}
	};
	static const double dite[2][NDIM] = {
		{3.46e4,3.84e5,4.08e5,3.82e5,3.53e5,3.19e5,3.22e5,2.47e5,2.29e5,
		 3.73e6,9.26e5,7.96e5,6.90e5,6.70e5,4.72e5,5.67e5,4.21e5,3.65e7,0.},
		{1.64e4,2.45e5,4.27e5,6.92e5,8.78e5,7.43e5,6.99e5,4.43e5,2.81e5,
		 5.84e6,4.89e6,4.62e6,4.52e6,3.32e6,1.37e6,4.41e6,2.27e6,7.25e6,0.}
	};
	static const double ditcrt[NDIM] = {6e3,2e4,4e4,5e4,7e4,8e4,8e4,3e4,
	  3e4,3e4,3e4,9e4,4e4,5e4,3e4,9e5,2e5,2e5,1e20};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonPotas()" );

	/* potassium nelem=19
	 * rec data from calcium
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */

	/* rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipPOTASSIUM] )
	{
		return;
	}

	ion_zero(ipPOTASSIUM);

	ion_photo(ipPOTASSIUM,false);

	/* find collisional ionization rates */
	ion_collis(ipPOTASSIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipPOTASSIUM);

	/* >>chng 03 nov 02, rm pl and rec components, were old */
#	if 0
	/* correct low temp rad rec coef */
	if( phycon.te < 1e3 )
	{
		ionbal.RateRecomTot[ipPOTASSIUM][1] += ((5.49e-10*(pow(phycon.te,-0.647f)) - 
		  rec[1]*(pow((double)phycon.te,pl[1])))*dense.eden);
	}
#	endif

	/* solve for ionization balance */
	ion_solver(ipPOTASSIUM,false);
	return;
}
