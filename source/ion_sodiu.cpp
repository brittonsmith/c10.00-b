/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonSodiu ionization balance for sodium */
#include "cddefines.h"
#include "trace.h"
#include "dense.h"
#include "ionbal.h"

void IonSodiu(void)
{
	const int NDIM = ipSODIUM+1;

	static const double dicoef[2][NDIM] = {
		{1.0e-3,2.6e-3,6.0e-3,1.1e-2,8.0e-3,1.0e-2,3.2e-2,1.2e-2,2.2e-1,1.8e-1,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double dite[2][NDIM] = {
		{3.6e5,3.8e5,3.4e5,3.0e5,2.7e5,2.8e5,3.2e5,1.9e5,1.2e7,1.3e7,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double ditcrt[NDIM] = {0.,1.3e5,4.4e5,7.2e4,2.8e5,4.0e5,
	  1.3e6,1.6e6,1.9e6,3.4e4,0.};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonSodiu()" );

	/* sodium nelem=11
	 *
	 * rates taken from Mg Shull and van Steenberg, Ap.J. Sup 48, 95. */

	/* from Aldrovandi and Pequignot Revista Brasileira de Fisica, 4, 491, */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipSODIUM] )
	{
		return;
	}

	ion_zero(ipSODIUM);

	ion_photo(ipSODIUM,false);

	/* find collisional ionization rates */
	ion_collis(ipSODIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipSODIUM);

	/* solve for ionization balance */
	ion_solver(ipSODIUM,false);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     IonSodiu returns; frac=" );
		for( int i=0; i < 10; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipSODIUM][i]/
			  dense.gas_phase[ipSODIUM] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
