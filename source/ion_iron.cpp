/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonIron ionization balance for iron */
#include "cddefines.h"
#include "dense.h"
#include "fe.h"
#include "thermal.h"
#include "gammas.h"
#include "opacity.h"
#include "trace.h"
#include "grainvar.h"
#include "ionbal.h"
#include "mole.h"

/*IonIron ionization balance for iron */
void IonIron(void)
{
	const int NDIM = ipIRON+1;

	static const double dicoef[2][NDIM] = {
		{1.58e-3,8.38e-3,1.54e-2,3.75e-2,0.117,0.254,0.291,0.150,0.140,0.100,0.200,0.240,
		 0.260,0.190,0.120,0.350,0.066,0.10,0.13,0.23,0.14,0.11,0.041,0.747,0.519,0.},
		{.456,.323,.310,.411,.359,.0975,.229,4.20,3.30,5.30,1.50,0.700,.600,.5,1.,0.,7.8,
		 6.3,5.5,3.6,4.9,1.6,4.2,.284,.279,0.}
	};
	static const double dite[2][NDIM] = {
		{6.00e4,1.94e5,3.31e5,4.32e5,6.28e5,7.50e5,7.73e5,2.62e5,2.50e5,2.57e5,2.84e5,
		 8.69e5,4.21e5,4.57e5,2.85e5,8.18e5,1.51e6,1.30e6,1.19e6,1.09e6,9.62e5,7.23e5,
		 4.23e5,5.84e7,6.01e7,0.},
		{8.97e4,1.71e5,2.73e5,3.49e5,5.29e5,4.69e5,6.54e5,1.32e6,1.33e6,1.41e6,1.52e6,
		 1.51e6,1.82e6,1.84e6,2.31e6,0.,9.98e6,9.98e6,1.00e7,1.10e7,8.34e6,1.01e7,1.07e7,
		 1.17e7,9.97e6,0.}
	};
	static const double ditcrt[NDIM] = {1e11,1e11,1e11,1e11,1e11,1e11,1e11,
	  1e11,1e11,1e11,1e11,1e11,1e11,1e11,1e11,1e11,1e11,1e11,1e11,
	  1e11,1e11,1e11,1e11,1e11,1e11,1e11};
	static const double aa[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double bb[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double cc[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double dd[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	/* >>chng 03 aug 15, zero out this array.  having non-zero values in
	 * this array had the effect of using 0 for the guess to the low-T dr
	 * in makerecomb */
	static const double ff[NDIM] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	static const double fyield[NDIM+1] = {.34,.34,.35,.35,.36,.37,.37,.38,.39,.40,
	  .41,.42,.43,.44,.45,.46,.47,.47,.48,.48,.49,.49,.11,.75,0.,0.,0.};

	long int i, limit, limit2;

	DEBUG_ENTRY( "IonIron()" );

	/* iron nelem=26
	 *
	 * Fe rates from Woods et al. Ap.J. 249, 399.
	 * rec from +23, 24 25 from Arnauld et al 85 */

	/* Pequignot and Aldrovandi Ast Ap 161, 169. */

	if( !dense.lgElmtOn[ipIRON] )
	{
		fe.fekcld = 0.;
		fe.fekhot = 0.;
		fe.fegrain = 0.;
		return;
	}

	ion_zero(ipIRON);

	ion_photo(ipIRON,false);
	/* debugging printout for shell photo rates */
	/*GammaPrtRate( ioQQQ , 0 , ipIRON , true);*/
	/*GammaPrtShells( ipIRON , 13 );
	GammaPrtShells( ipIRON , 12 );
	GammaPrtShells( ipIRON , 0 );*/
	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			long int iplow , iphi , ipop , ns , ion;
			double rate;
			ns = 6;
			ion = 0;
			/* show what some of the ionization sources are */
			iplow = opac.ipElement[ipIRON][ion][ns][0];
			iphi = opac.ipElement[ipIRON][ion][ns][1];
			ipop = opac.ipElement[ipIRON][ion][ns][2];
			rate = ionbal.PhotoRate_Shell[ipIRON][ion][ns][0];
			GammaPrt( iplow , iphi , ipop , ioQQQ, rate , rate*0.1 );
		}
	}

	/* find collisional ionization rates */
	ion_collis(ipIRON);

	/* get recombination coefficients */
	ion_recomb(false,(const double*)dicoef,(const double*)dite,ditcrt,aa,bb,cc,dd,ff,ipIRON);

	/* >>chng 06 Feb 28 -- NPA.  Add in charge transfer ionization of Mg with S+, Si+, and C+ */
	/* only include this if molecular network is enabled - otherwise no feedback onto
	 * Si+, S+, and C+ soln */
	if( !co.lgNoCOMole )
	{
		ionbal.PhotoRate_Shell[ipIRON][0][6][0] +=
			CO_findrk("S+,Fe=>S,Fe+")*dense.xIonDense[ipSULPHUR][1] +
			CO_findrk("Si+,Fe=>Si,Fe+")*dense.xIonDense[ipSILICON][1] + 
			CO_findrk("C+,Fe=>C,Fe+")*dense.xIonDense[ipCARBON][1];
		/* CO_sink_rate("Fe");*/

	}

	/* solve for ionization balance */
	ion_solver(ipIRON,false);

	/* now find total Auger yield of K-alphas
	 * "cold" iron has M-shell electrons, up to Fe 18 */
	fe.fekcld = 0.;
	limit = MIN2(18,dense.IonHigh[ipIRON]);

	for( i=dense.IonLow[ipIRON]; i < limit; i++ )
	{
		ASSERT( i < NDIM + 1 );
		fe.fekcld += 
			(realnum)(ionbal.PhotoRate_Shell[ipIRON][i][0][0]*dense.xIonDense[ipIRON][i]*
		  fyield[i]);
	}

	/* same sum for hot iron */
	fe.fekhot = 0.;
	limit = MAX2(18,dense.IonLow[ipIRON]);

	limit2 = MIN2(ipIRON+1,dense.IonHigh[ipIRON]);
	ASSERT( limit2 <= LIMELM + 1 );

	for( i=limit; i < limit2; i++ )
	{
		ASSERT( i < NDIM + 1 );
		fe.fekhot += 
			(realnum)(ionbal.PhotoRate_Shell[ipIRON][i][0][0]*dense.xIonDense[ipIRON][i]*
		  fyield[i]);
	}

	/* Fe Ka from grains - Fe in grains assumed to be atomic
	 * gv.elmSumAbund[ipIRON] is number density of iron added over all grain species */
	i = 0;
	/* fyield is 0.34 for atomic fe */
	fe.fegrain = ( gv.lgWD01 ) ? 0.f : (realnum)(ionbal.PhotoRate_Shell[ipIRON][i][0][0]*fyield[i]*
				 gv.elmSumAbund[ipIRON]);

	if( trace.lgTrace && trace.lgHeavyBug )
	{
		/* print densities of each stage of ionization */
		fprintf( ioQQQ, "     Fe" );
		for( i=0; i < 15; i++ )
		{
			fprintf( ioQQQ, "\t%.1e", dense.xIonDense[ipIRON][i]/dense.gas_phase[ipIRON] );
		}
		fprintf( ioQQQ, "\n" );
	}

	if( trace.lgTrace && trace.lgFeBug )
	{
		fprintf( ioQQQ, " IonIron-Abund:" );
		for( i=0; i < 27; i++ )
		{
			fprintf( ioQQQ, "%10.2e", dense.xIonDense[ipIRON][i] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " IonIron - Ka(Auger)%10.2e\n", fe.fekcld + 
		  fe.fekhot );
	}
	return;
}
