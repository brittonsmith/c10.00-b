/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonBeryl ionization balance for beryllium */
#include "cddefines.h"
#include "dense.h"
#include "trace.h"
#include "ionbal.h"

void IonBeryl(void)
{
	const int NDIM = ipBERYLLIUM+1;

	static const double dicoef[2][NDIM] = {
		{2.54e-3,6.15e-3,1.62e-3,0.}, {4.42e-2,5.88e-2,0.343,0.}
	};
	static const double dite[2][NDIM] = {
		{1.57e5,1.41e5,8.19e4,0.}, {3.74e5,1.41e5,1.59e5,0.}
	};
	static const double ditcrt[NDIM] = {1.2e4,1.2e4,1.1e4,1e20};
	static const double aa[NDIM] = {2.3196,0.,0.,0.};
	static const double bb[NDIM] = {10.7328,0.,0.,0.};
	static const double cc[NDIM] = {6.8830,0.,0.,0.};
	static const double dd[NDIM] = {-0.1824,0.,0.,0.};
	static const double ff[NDIM] = {0.4101,0.1,0.1,0.};

	DEBUG_ENTRY( "IonBeryl()" );

	/* boron nelem=4
	 * data are for carbon
	 *
	 * real CollidRate(nelem,2)
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */
	/* DATA GRDEFF/0.10,0.10,0.10,0.053/
	 * GRDEFF is fraction of recombinations to ground state, used for
	 * outward diffuse fields
	 *
	 * rec from +3, +4 from Arnaud et al Ast Ap Sup 60 425. (1985)
	 * rec from fully ionized uses Seaton '79 in ionrat */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipBERYLLIUM] )
	{
		return;
	}

	/* zero out ionization balance arrays */
	ion_zero(ipBERYLLIUM);

	ion_photo(ipBERYLLIUM,false);

	/* find collisional ionization rates */
	ion_collis(ipBERYLLIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipBERYLLIUM);

	/* solve for ionization balance */
	ion_solver(ipBERYLLIUM,false);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     Beryli returns; frac=" );
		for( int i=0; i < ipBERYLLIUM+2; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipBERYLLIUM][i]/
			  dense.gas_phase[ipBERYLLIUM] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
