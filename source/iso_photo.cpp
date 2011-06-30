/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_photo do photoionization rates for element nelem on the ipISO isoelectronic sequence */
#include "cddefines.h"
#include "hydrogenic.h"
#include "rfield.h"
#include "opacity.h"
#include "trace.h"
#include "ionbal.h"
#include "thermal.h"
#include "gammas.h"
#include "iso.h"
#include "taulines.h"

void iso_photo(
	/* iso sequence, hydrogen or helium for now */
	long ipISO , 
	/* the chemical element, 0 for hydrogen */
	long int nelem)
{
	long int limit ,
		n;

	t_phoHeat photoHeat;

	DEBUG_ENTRY( "iso_photo()" );

	/* check that we were called with valid charge */
	ASSERT( nelem >= 0 && nelem < LIMELM );
	ASSERT( ipISO < NISO );

	/* do photoionization rates */
	/* induced recombination; FINDUC is integral of
	 * pho rate times EXP(-hn/kt) for induc rec
	 * CIND is this times hnu-hnu0 to get ind rec cooling
	 * ionbal.lgPhotoIoniz_On is 1, set to 0 with 'no photoionization' command
	 * ipSecIon points to 7.353 Ryd, lowest energy where secondary ioniz
	 * of hydrogen is possible */

	/* photoionization of ground, this is different from excited states because 
	 * there will eventually be more than one shell, (when li-like done)
	 * and because upper limit is high-energy limit to code, so secondaries are 
	 * included */
	iso.gamnc[ipISO][nelem][0] = GammaBn(iso.ipIsoLevNIonCon[ipISO][nelem][0],
		rfield.nflux,
		iso.ipOpac[ipISO][nelem][0],
		iso.xIsoLevNIonRyd[ipISO][nelem][0],
		&iso.RecomInducRate[ipISO][nelem][0],
		&iso.RecomInducCool_Coef[ipISO][nelem][0],
		&photoHeat)*
		ionbal.lgPhotoIoniz_On;

	/* heating due to photo of ground */
	iso.PhotoHeat[ipISO][nelem][0] = photoHeat.HeatNet*ionbal.lgPhotoIoniz_On;

	/* save these rates into ground photo rate vector */
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][0] = iso.gamnc[ipISO][nelem][ipH1s];
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][1] = photoHeat.HeatLowEnr*ionbal.lgPhotoIoniz_On;
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][2] = photoHeat.HeatHiEnr*ionbal.lgPhotoIoniz_On;

	/* CompRecoilIonRate is direct photioniz rate due to 
	 * bound compton scattering of very hard x-rays+Compton scat */
	/* now add in compton recoil, to ground state, save heating as high energy heat */
	ASSERT( ionbal.CompRecoilIonRate[nelem][nelem-ipISO]>=0. &&
		ionbal.CompRecoilHeatRate[nelem][nelem-ipISO]>= 0. );
	iso.gamnc[ipISO][nelem][0] += ionbal.CompRecoilIonRate[nelem][nelem-ipISO];
	iso.PhotoHeat[ipISO][nelem][0] += ionbal.CompRecoilHeatRate[nelem][nelem-ipISO];

	/* now add bound compton scattering to ground term photoionization rate */
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][0] += ionbal.CompRecoilIonRate[nelem][nelem-ipISO];
	/* add heat to high energy heating term */
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][2] += ionbal.CompRecoilHeatRate[nelem][nelem-ipISO];

	/* option to print ground state photoionization rates */
	if( trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO]) )
	{
		GammaPrt(iso.ipIsoLevNIonCon[ipISO][nelem][0],
			rfield.nflux,
			iso.ipOpac[ipISO][nelem][0],
			ioQQQ,
			iso.gamnc[ipISO][nelem][0],
			iso.gamnc[ipISO][nelem][0]*0.05);
	}

	/* for excited states upper limit to integration is ground threshold */
	limit = iso.ipIsoLevNIonCon[ipISO][nelem][0]-1;
	/* >>chng 06 aug 17, to numLevels_local instead of _max. */
	for( n=1; n < iso.numLevels_local[ipISO][nelem]; n++ )
	{
		/* continuously update rates for n <=3, but only update
		 * rates for higher levels when redoing static opacities */
		if( StatesElemNEW[nelem][nelem-ipISO][n].n>4 && !opac.lgRedoStatic )
			break;
		/** \todo	2	- hydro.lgHInducImp should depend on iso and nelem,
		 * even better - just call one gamnc and within that code
		 * check to see whether induced is important by looking
		 * at occnum near threshold */
		if( hydro.lgHInducImp )
		{
			iso.gamnc[ipISO][nelem][n] = 
				GammaBn(
				iso.ipIsoLevNIonCon[ipISO][nelem][n],
				limit,
				iso.ipOpac[ipISO][nelem][n],
				iso.xIsoLevNIonRyd[ipISO][nelem][n],
				&iso.RecomInducRate[ipISO][nelem][n],
				&iso.RecomInducCool_Coef[ipISO][nelem][n],
				&photoHeat)*
				ionbal.lgPhotoIoniz_On;
		}
		else
		{
			iso.gamnc[ipISO][nelem][n] = 
				GammaK(iso.ipIsoLevNIonCon[ipISO][nelem][n],
				limit,
				iso.ipOpac[ipISO][nelem][n],1.,
				&photoHeat)*
				ionbal.lgPhotoIoniz_On;

			/* these are zero */
			iso.RecomInducRate[ipISO][nelem][n] = 0.;
			iso.RecomInducCool_Coef[ipISO][nelem][n] = 0.;
		}
		iso.PhotoHeat[ipISO][nelem][n] = photoHeat.HeatNet*ionbal.lgPhotoIoniz_On;

		ASSERT( iso.gamnc[ipISO][nelem][n]>= 0. );
		ASSERT( iso.PhotoHeat[ipISO][nelem][n]>= 0. );
		/* this loop only has excited states */
	}

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			if( nelem==ipHYDROGEN )
			{
				fprintf(ioQQQ," buggbugg hphotodebugg%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
					nzone,
				  iso.gamnc[ipISO][nelem][0],
				  iso.gamnc[ipISO][nelem][1],
				  iso.gamnc[ipISO][nelem][3],
				  iso.gamnc[ipISO][nelem][4],
				  iso.gamnc[ipISO][nelem][5],
				  iso.gamnc[ipISO][nelem][6]);
			}
		}
	}

	/* >>chng 02 jan 19, kill excited state photoionization with case b no photo */
	/* option for case b conditions, kill all excited state photoionizations */
	if( opac.lgCaseB_no_photo )
	{
		for( n=1; n < iso.numLevels_max[ipISO][nelem]; n++ )
		{
			iso.gamnc[ipISO][nelem][n] = 0.;
			iso.RecomInducRate[ipISO][nelem][n] = 0.;
			iso.RecomInducCool_Coef[ipISO][nelem][n] = 0.;
		}
	}
	{
		/* this block turns off induced recom for some element */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && ipISO==1 && nelem==5)
		{
			/* this debugging block is normally not active */
			for( n=0; n < iso.numLevels_max[ipISO][nelem]; n++ )
			{
				iso.RecomInducRate[ipISO][nelem][n] = 0.;
			}
		}
	}

	if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
	{
		fprintf( ioQQQ, "     iso_photo, ipISO%2ld nelem%2ld low, hi=",ipISO,nelem);
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso.gamnc[ipISO][nelem][ipH1s]));
		ASSERT(nelem>=ipISO);
		fprintf( ioQQQ,PrintEfmt("%9.2e", ionbal.CompRecoilIonRate[nelem][nelem-ipISO]));
		fprintf( ioQQQ, " total=");
		fprintf( ioQQQ,PrintEfmt("%9.2e",iso.gamnc[ipISO][nelem][ipH1s] ));
		fprintf( ioQQQ, "\n");
	}
	return;
}
