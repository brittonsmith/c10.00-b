/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonCobal do ionization balance for cobalt */
#include "cddefines.h"
#include "dense.h"
#include "ionbal.h"

void IonCobal(void)
{
	const int NDIM = ipCOBALT+1;

	static const double dicoef[2][NDIM] = {
		{5.2e-3,1.38e-3,2.3e-2,4.19e-2,6.83e-2,0.122,0.30,0.15,0.697,0.709,0.644,0.525,0.446,
		 0.363,0.302,0.102,0.270,0.0467,0.0835,0.0996,0.199,0.240,0.115,0.0316,0.803,0.575,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double dite[2][NDIM] = {
		{2.01e5,3.05e5,4.20e5,5.56e5,6.72e5,7.93e5,9.00e5,1.00e6,7.81e5,7.64e5,7.44e5,6.65e5,
		 5.97e5,5.24e5,4.96e5,4.46e5,8.49e6,1.36e6,1.23e6,1.06e6,1.25e6,1.23e6,3.32e5,6.45e5,
		 6.65e7,6.81e7,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double ditcrt[NDIM] = {2e4,4e4,5e4,7e4,8e4,8e4,3e4,3e4,
	  3e4,3e4,9e4,4e4,5e4,3e4,9e5,2e5,2e5,2e5,2e5,1e5,7e4,4e4,6e6,
	  6e6,1e20,1e20,1e20};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonCobal()" );

	/* cobalt nelem=27
	 *
	 * rates for ni from Shull and van Steenberg, Ap.J. Sup 48, 95. */

	/* DATA FYIELD/1*.34,2*.35,.36,2*.37,.38,.39,.40,.41,.42,.43,.44,
	 *  4.45,.46,2*.47,2*.48,2*.49, .11,.75,5*0./
	 * above fluorescent yields quoted in Krolik+Kallman Ap.J.(Let) 320, L5.
	 * they are correct for iron, not ni
	 *
	 * GRDEFF is fraction of recombinations to ground state, used for
	 * DATA GRDEFF/27*0.2/
	 *
	 * all rates from Shull and Van Steenberg apj sup 48, 95. for Ni
	 * assumed to have same cs per stage of ionization */

	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipCOBALT] )
	{
		return;
	}

	ion_zero(ipCOBALT);

	ion_photo(ipCOBALT,false);

	/* find collisional ionization rates */
	ion_collis(ipCOBALT);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipCOBALT);

	/* solve for ionization balance */
	ion_solver(ipCOBALT,false);
	return;
}
