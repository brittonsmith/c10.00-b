/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonMagne ionization balance for magnesium */
#include "cddefines.h"
#include "atoms.h"
#include "trace.h"
#include "iso.h"
#include "dense.h"
#include "mole.h"
#include "opacity.h"
#include "thermal.h"
#include "gammas.h"
#include "ionbal.h"

void IonMagne(void)
{
	const int NDIM = ipMAGNESIUM+1;

	static const double dicoef[2][NDIM] = {
		{4.49e-4,1.95e-3,5.12e-3,7.74e-3,1.17e-2,3.69e-2,3.63e-2,4.15e-2,8.86e-3,.252,.144,0.},
		{.021,.074,.323,.636,.807,.351,.548,.233,.318,.315,.291,0.}
	};
	static const double dite[2][NDIM] = {
		{5.01e4,6.06e5,4.69e5,3.74e5,3.28e5,4.80e5,3.88e5,3.39e5,2.11e5,1.40e7,1.50e7,0.},
		{2.81e4,1.44e6,7.55e5,7.88e5,1.02e6,9.73e5,7.38e5,3.82e5,1.54e6,2.64e6,3.09e6,0.}
	};
	static const double ditcrt[NDIM] = {4.0e3,7.4e4,6.6e4,5.5e4,4.4e4,4.5e4,
	  4.5e4,5.0e5,3.4e4,2.4e6,4.0e6,1e20};
	static const double aa[NDIM] = {1.2044,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {-4.6836,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {7.6620,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {-0.5930,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double ff[NDIM] = {1.6260,0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "IonMagne()" );

	/* magnesium nelem=12
	 *
	 * rates from Shull and van Steenberg, Ap.J. Sup 48, 95. */

	/* rec from +9, +10, +11 from Arnauld et al '85 */
	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	/* mg i from nuss+storey, who say that =>Mgii is very small */

	if( !dense.lgElmtOn[ipMAGNESIUM] )
	{
		atoms.xMg2Max = 0.;
		return;
	}

	ion_zero(ipMAGNESIUM);

	ion_photo(ipMAGNESIUM,false);
	/* debugging printout for shell photo rates - 0 for atom, last true, also print details */
	/*GammaPrtRate( ioQQQ , 0 , ipMAGNESIUM , true );*/

	/* find collisional ionization rates */
	ion_collis(ipMAGNESIUM);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipMAGNESIUM);

	// photoionization of Mg+ from the excited upper level of Mg II 2798
	if( dense.IonLow[ipMAGNESIUM] <= 1 )
	{
		t_phoHeat dummy;
		/* photoionization from excited upper state of 2798 */
		realnum rmg2l = (realnum)GammaK(opac.ipmgex,
			iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][0],opac.ipOpMgEx,1., &dummy );
		ionbal.PhotoRate_Shell[ipMAGNESIUM][1][3][0] += rmg2l*atoms.popmg2;

		if( nzone <= 1 )
		{
			atoms.xMg2Max = 0.;
		}
		else if( ionbal.PhotoRate_Shell[ipMAGNESIUM][1][3][0] > 1e-30 )
		{
			/* remember max relative photoionization rate for possible comment */
			atoms.xMg2Max = (realnum)(MAX2(atoms.xMg2Max,rmg2l*atoms.popmg2/
			  ionbal.PhotoRate_Shell[ipMAGNESIUM][1][3][0]));
		}
	}
	else
	{
		atoms.xMg2Max = 0.;
	}

	/* charge transfer with the atom */
	if( dense.IonLow[ipMAGNESIUM] <= 0 )
	{
		/* >>chng 06 Feb 28 -- NPA.  Add in charge transfer ionization of Fe with S+, Si+, and C+ */
		/* only include this if molecular network is enabled - otherwise no feedback onto
		 * Si+, S+, and C+ soln */
		if( !co.lgNoCOMole )
		{
			/* Use sink rate from last completed state of network, including _all_ present & future Mg sinks */
			ionbal.PhotoRate_Shell[ipMAGNESIUM][0][3][0] +=	
				CO_findrk("S+,Mg=>S,Mg+")*dense.xIonDense[ipSULPHUR][1] +
				CO_findrk("Si+,Mg=>Si,Mg+")*dense.xIonDense[ipSILICON][1] + 
				CO_findrk("C+,Mg=>C,Mg+")*dense.xIonDense[ipCARBON][1];
			/* CO_sink_rate("Mg"); */
		}
	}

	/* solve for ionization balance */
	ion_solver(ipMAGNESIUM,false);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     IonMagne returns; frac=" );
		for( int i=0; i < 10; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipMAGNESIUM][i]/
			  dense.gas_phase[ipMAGNESIUM] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
