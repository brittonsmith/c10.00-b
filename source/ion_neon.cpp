/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonNeon ionization balance for neon */
#include "cddefines.h"
#include "dense.h"
#include "trace.h"
#include "ionbal.h"

void IonNeon(void)
{
	const int NDIM = ipNEON+1;

	static const double dicoef[2][NDIM] = {
		{9.77e-4,2.65e-3,3.69e-3,1.12e-2,2.44e-2,3.02e-2,6.10e-3,2.52e-1,1.44e-1,0.},
		{.073,.242,1.01,.391,2.52,.445,.254,.304,.296,0.}
	};
	static const double dite[2][NDIM] = {
		{3.11e5,2.84e5,2.24e5,2.7e5,3.09e5,2.83e5,1.68e5,1.4e7,1.5e7,0.},
		{2.06e5,3.07e5,2.94e5,5.50e5,9.91e5,1.73e6,6.13e5,1.80e6,2.24e6,0.}
	};
	static const double ditcrt[NDIM] = {3.0e4,3.3e4,3.3e4,3.5e4,3.6e4,3.6e4,2.9e4,1.5e6,3.8e6,1e20};
	static const double aa[NDIM] = {0.,0.0129,3.6781,-0.0254,-0.0141,19.9280,5.4751,0.,0.,0.};
	static const double bb[NDIM] = {0.,-0.1779,14.1481,5.5365,33.8479,235.0536,203.9751,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.9353,17.1175,17.0727,43.1608,152.5096,86.9016,0.,0.,0.};
	static const double dd[NDIM] = {0.,-0.0682,-0.5017,-0.7225,-1.6072,9.1413,-7.4568,0.,0.,0.};
	static const double ff[NDIM] = {0.1,0.4516,0.2313,0.1702,0.1942,0.1282,2.5145,0.,0.,0.};

	DEBUG_ENTRY( "IonNeon()" );

	/* neon, atomic number 10 */
	if( !dense.lgElmtOn[ipNEON] )
	{
		return;
	}

	ion_zero(ipNEON);

	ion_photo(ipNEON,false);

	/* find collisional ionization rates */
	ion_collis(ipNEON);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipNEON);

	/* solve for ionization balance */
	ion_solver(ipNEON,false);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     IonNeon   returns; frac=" );
		for( int i=0; i < 10; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipNEON][i]/
			  dense.gas_phase[ipNEON] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
