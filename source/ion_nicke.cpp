/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonNicke ionization balance for nickel */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonNicke(void)
{
	const int NDIM = ipNICKEL+1;

	static const double dicoef[2][NDIM] = {
		{1.41e-3,5.2e-3,.0138,.023,.0419,.0683,.122,.300,.15,.697,.709,.644,.525,
		 .446,.363,.302,.102,.27,.0467,.0835,.0996,.199,.24,.115,.0316,.803,.575,0.},
		{.469,.357,.281,.128,.0417,.0558,.0346,0.,1.9,.277,.135,.134,.192,.332,.337,
		 .121,.0514,.183,7.56,4.55,4.87,2.19,1.15,1.23,.132,.289,.286,0.}
	};
	static const double dite[2][NDIM] = {
		{9.82e4,2.01e5,3.05e5,4.20e5,5.56e5,6.72e5,7.93e5,9.00e5,1.00e6,7.81e5,7.64e5,
		 7.44e5,6.65e5,5.97e5,5.24e5,4.96e5,4.46e5,8.46e6,1.36e6,1.23e6,1.06e6,1.25e6,
		 1.23e6,3.32e5,6.45e5,6.65e7,6.81e7,0.},
		{1.01e5,1.91e5,2.32e5,3.18e5,4.55e5,5.51e5,5.28e5,0.,5.50e5,8.87e5,1.80e6,1.25e6,
		 1.89e6,8.84e5,1.29e6,6.24e5,1.59e6,8.01e6,9.32e6,9.45e6,9.45e6,8.01e6,7.57e6,
		 2.64e6,1.93e6,1.19e7,9.08e6,0.}
	};
	static const double ditcrt[NDIM] = {6e3,2e4,4e4,5e4,7e4,8e4,8e4,3e4,3e4,3e4,3e4,
	  9e4,4e4,5e4,3e4,9e5,2e5,2e5,2e5,2e5,1e5,7e4,4e4,6e6,6e6,1e20,1e20,1e20};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonNicke()" );

	/* nickel, nelem=28
	 *
	 * all rates from Shull and Van Steenberg apj sup 48, 95. */

	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipNICKEL] )
	{
		return;
	}

	ion_zero(ipNICKEL);

	ion_photo(ipNICKEL,false);

	/* find collisional ionization rates */
	ion_collis(ipNICKEL);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipNICKEL);

	/* solve for ionization balance */
	ion_solver(ipNICKEL,false);
	return;
}
