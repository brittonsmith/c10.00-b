/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_stark compute stark broadening escape probabilities using Puetter formalism */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "dense.h"
#include "hydrogenic.h"
#include "phycon.h"
#include "rt.h"

void RT_stark(void)
{
	long int ipLo, 
	  ipHi,
	  nelem,
	  ipISO;

	double aa , ah, 
	  stark, 
	  strkla;

	DEBUG_ENTRY( "RT_stark()" );

	/* only evaluate one time per zone */
	static long int nZoneEval=-1;
	if( nzone==nZoneEval )
		return;
	nZoneEval = nzone;

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( nelem >= 2 && !dense.lgElmtOn[nelem] )
				continue;

			if( !rt.lgStarkON || dense.eden < 1e8 )
			{
				for( ipHi=0; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					for( ipLo=0; ipLo < iso.numLevels_max[ipISO][nelem]; ipLo++ )
					{
						iso.pestrk[ipISO][nelem][ipHi][ipLo] = 0.;
						iso.pestrk[ipISO][nelem][ipLo][ipHi] = 0.;
					}
				}
				continue;
			}

			/* evaluate Stark escape probability from 
			 * >>ref Puetter Ap.J. 251, 446. */

			/* coefficients for Stark broadening escape probability
			 * to be Puetters AH, equation 9b, needs factor of (Z^-4.5 * (nu*nl)^3 * xl) */
			ah = 6.9e-6*1000./1e12/(phycon.sqrte*phycon.te10*phycon.te10*
			  phycon.te03*phycon.te01*phycon.te01)*dense.eden;

			/* include Z factor */
			ah *= pow( (double)(nelem+1), -4.5 );

			/* coefficient for all lines except Ly alpha */
			/* equation 10b, except missing tau^-0.6 */
			stark = 0.264*pow(ah,0.4);

			/* coefficient for Ly alpha */
			/* first few factors resemble equation 13c...what about the rest? */
			strkla = 0.538*ah*4.*9.875*(phycon.sqrte/phycon.te10/phycon.te03);

			/* Lyman lines always have outer optical depths */
			/*ASSERT( Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].TauIn> 0. );*/
			/* >>chng 02 mar 31, put in max, crashed on some first iteration 
			 * with negative optical depths,
			 * NB did not understand why neg optical depths started */
			aa = (realnum)SDIV(Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->TauIn);
			aa = pow( aa ,-0.75);
			iso.pestrk[ipISO][nelem][0][iso.nLyaLevel[ipISO]] = strkla/2.*MAX2(1.,aa);

			/**\todo	2	- Stark is disabled for now since Lya escape causes density dependent
			 * feedback on the radiative transfer.  Would need to redo the escape
			 * probs every time the electron density is updated - see blr89.in for an 
			 * example */
			iso.pestrk[ipISO][nelem][0][iso.nLyaLevel[ipISO]] = MIN2(.01,iso.pestrk[ipISO][nelem][0][iso.nLyaLevel[ipISO]]);
			iso.pestrk[ipISO][nelem][0][iso.nLyaLevel[ipISO]] = 0.;
			iso.pestrk[ipISO][nelem][iso.nLyaLevel[ipISO]][0] = 
				iso.pestrk[ipISO][nelem][0][iso.nLyaLevel[ipISO]]*Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Emis->Aul;


			/* >>chng 06 aug 28, from numLevels_max to _local. */
			for( ipHi=3; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
			{
				if( Transitions[ipISO][nelem][ipHi][ipH1s].ipCont <= 0 )
					continue;

				iso.pestrk[ipISO][nelem][0][ipHi] = stark*iso.strkar[ipISO][nelem][0][ipHi]/2.*pow(MAX2(1.,
				  Transitions[ipISO][nelem][ipHi][ipH1s].Emis->TauIn),-0.75);

				iso.pestrk[ipISO][nelem][0][ipHi] = MIN2(.01,iso.pestrk[ipISO][nelem][0][ipHi]);
				iso.pestrk[ipISO][nelem][ipHi][0] = Transitions[ipISO][nelem][ipHi][ipH1s].Emis->Aul*
				  iso.pestrk[ipISO][nelem][0][ipHi];
			}

			/* zero out rates above iso.numLevels_local */
			for( ipHi=iso.numLevels_local[ipISO][nelem]; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				iso.pestrk[ipISO][nelem][0][ipHi] = 0.;
				iso.pestrk[ipISO][nelem][ipHi][0] = 0.;
			}

			/* all other lines */
			for( ipLo=ipH2s; ipLo < (iso.numLevels_local[ipISO][nelem] - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
				{
					if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont <= 0 )
						continue;

					aa = stark*iso.strkar[ipISO][nelem][ipLo][ipHi]*
					  pow(MAX2(1.,Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn),-0.75);
					iso.pestrk[ipISO][nelem][ipLo][ipHi] = 
						(realnum)MIN2(.01,aa);

					iso.pestrk[ipISO][nelem][ipHi][ipLo] = Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul*
					  iso.pestrk[ipISO][nelem][ipLo][ipHi];
				}
			}

			/* zero out rates above iso.numLevels_local */
			for( ipLo=(iso.numLevels_local[ipISO][nelem] - 1); ipLo<(iso.numLevels_max[ipISO][nelem] - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					iso.pestrk[ipISO][nelem][ipLo][ipHi] = 0.;
					iso.pestrk[ipISO][nelem][ipHi][ipLo] = 0.;
				}
			}
		}
	}

	return;
}
