/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_ionize_recombine find state specific creation and destruction rates for iso sequences */
/*ChargeTransferUpdate update rate of ct ionization and recombination for H atoms */
#include "cddefines.h"
#include "ionbal.h"
#include "conv.h"
#include "atmdat.h"
#include "dense.h"
#include "hmi.h"
#include "mole.h"
#include "secondaries.h"
#include "iso.h"
#include "trace.h"
/*lint -e778 constant expression evaluates to zero - in HMRATE macro */
/*iso_charge_transfer_update update rate of ct ionization and recombination for H atoms */
void iso_charge_transfer_update(void)
{
	long int ion;
	long int nelem;

	DEBUG_ENTRY( "iso_charge_transfer_update()" );

	/* find total rate for charge transfer ionization of hydrogen,
	 * recombination for other element is ionization of hydrogen */
	atmdat.HCharExcIonTotal = 0.;
	for( nelem=ipHELIUM; nelem<LIMELM; ++nelem)
	{
		/* ion on the C scale, 0 is atom, so goes up to nelem+1,
		 * for helium nelem=1, ions go from 1 to 2 */
		for( ion=1; ion<=nelem+1; ++ion )
		{
			double one;
			/* we intentionally skip CT with O+ since this is in hmole */
			/* >>chng 05 aug 17, add logic on !ionbal.lgHO_ct_chem
			 * this is flag that is true when H O charge transfer is done in
			 * chemistry, not here.  set false with set HO charge transfer ionization
			 * command, so we do it here - introduced by NA for the Leiden hacks */
			if( (nelem == ipOXYGEN) && (ion == 1) && ionbal.lgHO_ct_chem ) 
				continue;
			/* charge transfer ionization of H, recombination for other species */
			one = atmdat.HCharExcRecTo[nelem][ion-1]*dense.xIonDense[nelem][ion];  //units s-1
			/* this, when multiplied by dense.xIonDense[ipHYDROGEN][0], is the charge transfer
			 * ionization RATE of neutral hydrogen, units cm-3 s-1 */
			atmdat.HCharExcIonTotal += one;		 
		}
	}

	/* >>chng 01 may 07,  add this in */
	/* charge transfer recombination of hydrogen,
	 * which is ionization of the heavy element */
	atmdat.HCharExcRecTotal = 0.;
	for( nelem=ipHELIUM; nelem<LIMELM; ++nelem)
	{
		/* this is ion on the abundances scale, 1 is atom, so goes up to nelem+1,
		 * for helium nelem=1, ion must go up to 3 */
		for( ion=0; ion<=nelem; ++ion )
		{
			/* option to skip Oo => O+ */
			/* >>chng 05 aug 17, add logic on !ionbal.lgHO_ct_chem
			 * this is flag that is true when H O charge transfer is done in
			 * chemistry, not here.  set false with set HO charge transfer ionization
			 * command, so we do it here - introduced by NA for the Leiden hacks */
			if( (nelem == ipOXYGEN) && (ion == 0) && ionbal.lgHO_ct_chem ) 
				continue;
			/* charge transfer ionization of H, recombination for other species */
			/* this, when multiplied by dense.xIonDense[ipHYDROGEN][1], is the charge transfer
			 * recombination RATE of neutral hydrogen, units cm-3 s-1 */
			atmdat.HCharExcRecTotal += atmdat.HCharExcIonOf[nelem][ion]*dense.xIonDense[nelem][ion];
		}
	}

	/* >>logic checked 02 may 02, think it's right */
	/* add on charge transfer ionization of helium,
	 * recombination for other element is ionization of helium */
	atmdat.HeCharExcIonTotal = 0.;
	/* loop up from Li, the next heavier element */
	for( nelem=ipLITHIUM; nelem<LIMELM; ++nelem)
	{
		double hold_one = 0.;
		/* check that array bounds not exceeded */
		/* ion on the C scale, 0 is atom, so goes up to nelem+1,
		 * for helium nelem=1, ions go from 1 to 2 */
		for( ion=1; ion<=nelem+1; ++ion )
		{
			/* charge transfer ionization of He, recombination for other species */
			hold_one = atmdat.HeCharExcRecTo[nelem][ion-1]*dense.xIonDense[nelem][ion];
			atmdat.HeCharExcIonTotal += hold_one;
		}
	}

	/* >>chng 04 jul 02,
	 * add on charge transfer reactions of He-H */
	atmdat.HeCharExcIonTotal += atmdat.HCharExcIonOf[ipHELIUM][0]*dense.xIonDense[ipHYDROGEN][1];

	/* charge transfer recombination of helium,
	 * which is ionization of the heavy element */
	atmdat.HeCharExcRecTotal = 0.;
	for( nelem=ipLITHIUM; nelem<LIMELM; ++nelem)
	{
		/* this is ion on the physics scale, 1 is atom, so goes up to nelem+1,
		 * for helium nelem=1, ion must go up to 3 */
		for( ion=0; ion<=nelem; ++ion )
		{
			/* charge transfer recombination of He, ionization for other species */
			atmdat.HeCharExcRecTotal += 
				atmdat.HeCharExcIonOf[nelem][ion]*dense.xIonDense[nelem][ion];
		}
	}
	/* >>chng 04 jul 02
	 * Add on charge transfer reactions of He+ +H0 -> He0 + H+ */
	atmdat.HeCharExcRecTotal += atmdat.HCharExcRecTo[ipHELIUM][0]*dense.xIonDense[ipHYDROGEN][0];

	return;
}

/*iso_ionize_recombine find state specific creation and destruction rates for iso sequences */
void iso_ionize_recombine(
	/* iso sequence, hydrogen or helium for now */
	long ipISO , 
	/* the chemical element, 0 for hydrogen */
	long int nelem )
{
	long int level;
	double Recom3Body;

	DEBUG_ENTRY( "iso_ionize_recombine()" );

	ASSERT( ipISO >= 0 && ipISO < NISO );
	ASSERT( nelem >= 0 && nelem < LIMELM );

	/* find recombination and ionization elements, will use to get simple estimate
	 * of ionization ratio below */
	/* >>chng 06 jul 20, level should go to numLevels_local instead of numLevels_max */

	fixit(); /* must apply error to iso.gamnc */

	for( level=ipH1s; level< iso.numLevels_local[ipISO][nelem]; ++level)
	{
		/* all process moving level to continuum, units s-1 */
		iso.RateLevel2Cont[ipISO][nelem][level] = iso.gamnc[ipISO][nelem][level] + 
		  iso.ColIoniz[ipISO][nelem][level]* dense.EdenHCorr + 
		  secondaries.csupra[nelem][nelem-ipISO]*iso.lgColl_ionize[ipISO];

		/* all processes from continuum to level n, units s-1 */
		iso.RateCont2Level[ipISO][nelem][level] = (
			/* radiative recombination */
			iso.RadRecomb[ipISO][nelem][level][ipRecRad]*
			iso.RadRecomb[ipISO][nelem][level][ipRecNetEsc] + 

			/* dielectronic recombination */
			iso.DielecRecomb[ipISO][nelem][level] +

			/* induced recombination */
			iso.RecomInducRate[ipISO][nelem][level]*iso.PopLTE[ipISO][nelem][level] + 

			/* collisional or three body recombination */
			/* PopLTE(level,nelem) is only LTE pop when mult by n_e n_H */
			iso.ColIoniz[ipISO][nelem][level]*dense.EdenHCorr*iso.PopLTE[ipISO][nelem][level]
			) * dense.eden;

			if( iso.lgRandErrGen[ipISO] )
			{
				iso.RateCont2Level[ipISO][nelem][level] *=
					iso.ErrorFactor[ipISO][nelem][ iso.numLevels_max[ipISO][nelem] ][level][IPRAD];
			}
	}

	/* now charge transfer - all into/from ground, two cases, H and not H */
	if( nelem==ipHYDROGEN )
	{
		/* charge transfer, hydrogen onto everything else */
		/* charge exchange ionization, these come out of every level */
		for( level=0; level< iso.numLevels_local[ipISO][nelem]; ++level)
		{
			iso.RateLevel2Cont[ipISO][nelem][level] += atmdat.HCharExcIonTotal;
		}
		/* charge exchange recombination */
		iso.RateCont2Level[ipISO][nelem][0] += atmdat.HCharExcRecTotal;
	}
	else if( nelem==ipHELIUM && ipISO==ipHE_LIKE )
	{
		/* add CO here, only for he itself, 
		 * destruction of CO but recombination for he */

		/* The following terms destroy He+ and thereby need to be included 
		 * here.  Also included are the contributions
		 * to the recombination due to hmole_step.*/

		fixit(); // should these chemistry rates be accounted for in mole.source and mole.sink?
		double COsource = CO_source_rate("He+");
		for( level=0; level< iso.numLevels_local[ipISO][nelem]; ++level)
		{
			iso.RateLevel2Cont[ipISO][nelem][level] += COsource;
			fixit(); //need reverse rates for hmi rates just below
		}
		iso.RateCont2Level[ipISO][nelem][0] += CO_sink_rate("He+") +
			hmi.H2_total*hmi.rheph2hpheh + hmi.heph2heh2p*hmi.H2_total;

		/* this is ioniz of He due to ct with all other gas constituents */
		for( level=0; level< iso.numLevels_local[ipISO][nelem]; ++level)
		{
			iso.RateLevel2Cont[ipISO][nelem][level] += atmdat.HeCharExcIonTotal;
		}
		/* this is recom of He due to ct with all other gas constituents */
		iso.RateCont2Level[ipISO][nelem][0] += atmdat.HeCharExcRecTotal;
	}
	else
	{
		for( level=0; level< iso.numLevels_local[ipISO][nelem]; ++level)
		{
			iso.RateLevel2Cont[ipISO][nelem][level] +=
				atmdat.HeCharExcIonOf[nelem][nelem-ipISO]*dense.xIonDense[ipHELIUM][1] +
				atmdat.HCharExcIonOf[nelem][nelem-ipISO]*dense.xIonDense[ipHYDROGEN][1];
		}

		iso.RateCont2Level[ipISO][nelem][0] +=
			atmdat.HeCharExcRecTo[nelem][nelem-ipISO]*dense.xIonDense[ipHELIUM][0] + 
			atmdat.HCharExcRecTo[nelem][nelem-ipISO]*dense.xIonDense[ipHYDROGEN][0];
	}


	/* now create sums of recombination and ionization rates for ISO species */
	ionbal.RateRecomTot[nelem][nelem-ipISO] = 0.;
	ionbal.RR_rate_coef_used[nelem][nelem-ipISO] = 0.;
	Recom3Body = 0.;
	/* >>chng 06 jul 20, level should go to numLevels_local instead of numLevels_max */
	for( level=0; level< iso.numLevels_local[ipISO][nelem]; ++level)
	{

		/* units of ionbal.RateRecomTot are s-1,
		 * equivalent ionization term is done after level populations are known */
		ionbal.RateRecomTot[nelem][nelem-ipISO] += iso.RateCont2Level[ipISO][nelem][level];

		/* just the radiative recombination rate coef, cm3 s-1 */
		ionbal.RR_rate_coef_used[nelem][nelem-ipISO] += iso.RadRecomb[ipISO][nelem][level][ipRecRad]*
			iso.RadRecomb[ipISO][nelem][level][ipRecNetEsc];

		/* >>chng 05 jul 11, from > to >=, some very opt thick sims did block escape to zero */
		ASSERT( ionbal.RR_rate_coef_used[nelem][nelem-ipISO]>= 0. );

		/* this is three-body recombination rate coef by itself - 
		 * need factor of eden to become rate */
		Recom3Body += iso.ColIoniz[ipISO][nelem][level]*dense.EdenHCorr*iso.PopLTE[ipISO][nelem][level];
	}

	/* fraction of total recombinations due to three body - when most are due to this
	 * small changes in temperature can produce large changes in recombination coefficient,
	 * and so in ionization */
	iso.RecomCollisFrac[ipISO][nelem] = Recom3Body* dense.eden / ionbal.RateRecomTot[nelem][nelem-ipISO];

	/* very first pass through here rate RateIoniz not yet evaluated */
	if( conv.nTotalIoniz==0 )
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = iso.RateLevel2Cont[ipISO][nelem][0];

	/* get simple estimate of atom to ion ratio */
	if( ionbal.RateRecomTot[nelem][nelem-ipISO] > 0. )
	{
		iso.xIonSimple[ipISO][nelem] = ionbal.RateIonizTot(nelem,nelem-ipISO)/ionbal.RateRecomTot[nelem][nelem-ipISO];
	}
	else
	{
		iso.xIonSimple[ipISO][nelem] = 0.;
	}

	if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
	{
		fprintf( ioQQQ, "     iso_ionize_recombine iso=%2ld Z=%2ld Level2Cont[0] %10.2e RateRecomTot %10.2e xIonSimple %10.2e\n", 
			ipISO, nelem, iso.RateLevel2Cont[ipISO][nelem][0], ionbal.RateRecomTot[nelem][nelem-ipISO], iso.xIonSimple[ipISO][nelem] );
	}

	return;
}
/*lint +e778 constant expression evaluates to zero - in HMRATE macro */
