/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonChlor ionization balance for chlorine */
#include "cddefines.h"
#include "dense.h"
#include "atmdat.h"
#include "ionbal.h"

void IonChlor(void)
{
	const int NDIM = ipCHLORINE+1;

	static const double dicoef[2][NDIM] = {
		{5.5e-5,1.0e-2,1.1e-2,1.0e-2,5.0e-2,3.2e-2,3.4e-2,1.6e-2,2.4e-2,4.0e-2,
		 4.0e-2,3.8e-2,6.8e-2,2.6e-2,4.6e-1,11.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double dite[2][NDIM] = {
		{1.3e5,1.4e5,1.4e5,1.4e5,2.0e5,1.8e5,2.3e5,7.2e5,6.4e5,6.0e5,4.9e5,
		 4.6e5,5.3e5,3.2e5,2.6e7,2.8e7,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double ditcrt[NDIM] = {3.0e4,2.5e4,2.5e4,1.8e4,1.8e4,2.2e4,
	  5.0e5,1.6e5,1.5e5,1.5e5,1.3e5,1.3e5,1.1e5,7.6e4,6.5e6,1.4e7,1e20};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	double ct_save;

	DEBUG_ENTRY( "IonChlor()" );

	/* chlorine, element number 17
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95.
	 * these are for argon, but no better can be done today */

	/* fits to Aldrovandi and Pequignot Rev Bras Fisica 4, 491 */

	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipCHLORINE] )
	{
		return;
	}

	ion_zero(ipCHLORINE);

	ion_photo(ipCHLORINE,false);

	/* find collisional ionization rates */
	ion_collis(ipCHLORINE);

	/* save the CT recom rate for H only */
	ct_save = atmdat.HCharExcRecTo[ipCHLORINE][0];

	/* >>chng 05 mar 23 add process Cl+ + H2 -> Cl0 + several molecules 
	 * rate is rate for first step with H2 from UMIST
	atmdat.HCharExcRecTo[ipCHLORINE][0] += 
		1e-9 * hmi.H2_total / SDIV(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop); */

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipCHLORINE);

	/* solve for ionization balance */
	ion_solver(ipCHLORINE,false);

	/* reset the rate */
	atmdat.HCharExcRecTo[ipCHLORINE][0] = ct_save;
	return;
}
