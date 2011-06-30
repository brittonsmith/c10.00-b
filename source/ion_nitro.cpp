/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonNitro ionization balance for nitrogen */
#include "cddefines.h"
#include "dense.h"
#include "thermal.h"
#include "atoms.h"
#include "opacity.h"
#include "trace.h"
#include "conv.h"
#include "gammas.h"
#include "ionbal.h"

void IonNitro(void)
{
	const int NDIM = ipNITROGEN+1;

	static const double dicoef[2][NDIM] = {
		{2.98e-3,7.41e-3,1.13e-2,2.62e-3,7.5e-2,4.61e-2,0.},
		{0.,.0764,.164,.243,.35,.309,0.}
	};
	static const double dite[2][NDIM] = {
		{2.2e5,2.01e5,1.72e5,1.02e5,4.75e6,5.44e6,0.},
		{0.,7.37e4,2.25e5,1.25e5,8.35e5,1.14e6,0.}
	};
	static const double ditcrt[NDIM] = {1.8e4,1.8e4,2.4e4,1.5e4,6.8e5,1.0e6,1e20};
	static const double aa[NDIM] = {0.0,0.0320,-0.8806,0.4134,0.,0.,0.};
	static const double bb[NDIM] = {0.6310,-0.6624,11.2406,-4.6319,0.,0.,0.};
	static const double cc[NDIM] = {0.1990,4.3191,30.7066,25.9172,0.,0.,0.};
	static const double dd[NDIM] = {-0.0197,0.0003,-1.1721,-2.2290,0.,0.,0.};
	static const double ff[NDIM] = {0.4398,0.5946,0.6127,0.2360,0.1,0.,0.};

	bool lgDebug;

	DEBUG_ENTRY( "IonNitro()" );

	/* nitrogen, atomic number 7 */
	if( !dense.lgElmtOn[ipNITROGEN] )
	{
		return;
	}

	ion_zero(ipNITROGEN);

	ion_photo(ipNITROGEN,false);

	/* find collisional ionization rates */
	ion_collis(ipNITROGEN);

	/* get recombination coefficients */
	/*lint -e64 type mismatch */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipNITROGEN);
	/*lint +e64 type mismatch */

	/* photoionization from 2D of NI, is atomic nitrogen present? */
	if( dense.xIonDense[ipNITROGEN][0] > 0. )
	{
		t_phoHeat photoHeat;
		// photo rate, population atoms.p2nit evaluated in cooling
		atoms.d5200r = (realnum)GammaK(opac.in1[0],opac.in1[1],opac.in1[2],1.,&photoHeat);
		/* valence shell photoionization, followed by heating; [0][6] => atomic nitrogen */
		ionbal.PhotoRate_Shell[ipNITROGEN][0][2][0] = ionbal.PhotoRate_Shell[ipNITROGEN][0][2][0]*
		  (1. - atoms.p2nit) + atoms.p2nit*atoms.d5200r;
		ionbal.PhotoRate_Shell[ipNITROGEN][0][2][1] = ionbal.PhotoRate_Shell[ipNITROGEN][0][2][1]*
		  (1. - atoms.p2nit) + photoHeat.HeatNet*atoms.p2nit;
	}
	else
	{
		atoms.p2nit = 0.;
		atoms.d5200r = 0.;
	}

	/* solve for ionization balance */
#	if 0
	broken();/* rm following true */
	if(nzone == 1 )
		lgDebug = true;
	else
		lgDebug = false;
#	endif
	lgDebug = false;
	ion_solver(ipNITROGEN,lgDebug);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     IonNitro retun; frac=" );
		for( int i=0; i < 8; i++ )
		{
			fprintf( ioQQQ, "%10.3e", dense.xIonDense[ipNITROGEN][i]/
			  dense.gas_phase[ipNITROGEN] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
