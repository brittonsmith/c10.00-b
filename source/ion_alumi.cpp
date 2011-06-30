/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonAlumi ionization balance for aluminum */
#include "cddefines.h"
#include "trace.h"
#include "dense.h"
#include "ionbal.h"

void IonAlumi(void)
{
	const int NDIM = ipALUMINIUM+1;

	static const double dicoef[2][NDIM] = {
		{5.5e-3,6.0e-3,7.5e-3,5.6e-3,1.2e-2,1.9e-2,1.6e-2,1.7e-2,4.4e-2,1.7e-2,3.0e-1,3.4e-0,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double dite[2][NDIM] = {
		{9.0e4,7.5e4,8.5e5,5.0e5,4.4e5,3.8e5,3.4e5,3.4e5,3.9e5,2.3e5,1.6e7,1.7e7,1e20},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};
	static const double ditcrt[NDIM] = {8.0e3,7.0e3,1.2e5,7.0e4,6.5e4,5.5e4,
	  5.5e4,5.5e4,5.8e4,4.0e4,2.9e6,5.5e6,1e20};
	static const double aa[NDIM] = {0.0219,0.7086,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {-0.4528,-3.1083,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {2.5427,7.0422,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {-0.1678,0.5998,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {0.2276,0.4194,0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonAlumi()" );

	/* aluminium nelem=13*/

	/* rates from Shull and van Steenberg, Ap.J. Sup 48, 95.
	 * use silicon since not in shull + van steenberg */

	/* rec coef from Aldrovandi and Pequignot Revista Bras de Fisica 4, 491. */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */


	/* from Nussbaumer and Storey, Mg, Al, Si: n.b., al1 uses table 2(b) */

	if( !dense.lgElmtOn[ipALUMINIUM] )
	{
		return;
	}

	ion_zero(ipALUMINIUM);

	ion_photo(ipALUMINIUM,false);

	/* find collisional ionization rates */
	ion_collis(ipALUMINIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipALUMINIUM);

	/* solve for ionization balance */
	ion_solver(ipALUMINIUM,false);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     IonAlumi returns; frac=" );
		for( int i=0; i < 10; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipALUMINIUM][i]/
			  dense.gas_phase[ipALUMINIUM] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
