/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonBoron ionization balance for boron */
#include "cddefines.h"
#include "dense.h"
#include "trace.h"
#include "ionbal.h"

void IonBoron(void)
{
	const int NDIM = ipBORON+1;

	static const double dicoef[2][NDIM] = {
		{2.54e-3,6.15e-3,1.62e-3,4.78e-2,0.}, {4.42e-2,5.88e-2,0.343,0.362,0.}
	};
	static const double dite[2][NDIM] = {
		{1.57e5,1.41e5,8.19e4,3.44e6,0.}, {3.74e5,1.41e5,1.59e5,5.87e5,0.}
	};
	static const double ditcrt[NDIM] = {1.2e4,1.2e4,1.1e4,4.4e5,1e20};
	static const double aa[NDIM] = {1.8267,2.3196,0.,0.,0.};
	static const double bb[NDIM] = {4.1012,10.7328,0.,0.,0.};
	static const double cc[NDIM] = {4.8443,6.8830,0.,0.,0.};
	static const double dd[NDIM] = {.2261,-0.1824,0.,0.,0.};
	static const double ff[NDIM] = {0.5960,0.4101,0.1,0.1,0.};

	DEBUG_ENTRY( "IonBoron()" );

	/* boron nelem=5
	 * data are for carbon
	 *
	 * real CollidRate(nelem,2)
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */
	/* DATA GRDEFF/0.10,0.10,0.10,0.053,0.10/
	 * GRDEFF is fraction of recombinations to ground state, used for
	 * outward diffuse fields
	 *
	 * rec from +3, +4 from Arnaud et al Ast Ap Sup 60 425. (1985)
	 * rec from fully ionized uses Seaton '79 in ionrat */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipBORON] )
	{
		return;
	}

	/* zero out ionization balance arrays */
	ion_zero(ipBORON);

	ion_photo(ipBORON,false);

	/* find collisional ionization rates */
	ion_collis(ipBORON);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipBORON);

	/* solve for ionization balance */
	ion_solver(ipBORON,false);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     Boroni returns; frac=" );
		for( int i=0; i < ipBORON+2; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipBORON][i]/
			  dense.gas_phase[ipBORON] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
