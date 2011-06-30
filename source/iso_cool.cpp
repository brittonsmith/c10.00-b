/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_cool compute net cooling due to species in iso-sequences */
#include "cddefines.h"
#include "physconst.h"
#include "taulines.h"
#include "hydrogenic.h"
#include "elementnames.h"
#include "phycon.h"
#include "dense.h"
#include "thermal.h"
#include "cooling.h"
#include "iso.h"

/* HP cc cannot compile this routine with any optimization */ 
#if defined(__HP_aCC)
#	pragma OPT_LEVEL 1
#endif

// set to true to enable debug print of contributors to collisional ionization cooling
const bool lgPrintIonizCooling = false;

void iso_cool(
	   /* iso sequence, 0 for hydrogenic */
	   long int ipISO , 
		/* nelem is charge -1, so 0 for H itself */
		long int nelem)
{
	long int ipHi, 
	  ipbig,
	  ipLo, 
	  n;
	double RecCoolExtra,
	  biggest = 0.,
	  dCdT_all, 
	  edenIonAbund, 
	  CollisIonizatCoolingTotal, 
	  dCollisIonizatCoolingTotalDT, 
	  HeatExcited,
	  heat_max,
	  CollisIonizatCooling, 
	  CollisIonizatCoolingDT, 
	  hlone,
	  thin,
	  ThinCoolingCaseB, 
	  ThinCoolingSum;

	valarray<double> CollisIonizatCoolingArr, 
	  CollisIonizatCoolingDTArr,
	  SavePhotoHeat,
	  SaveInducCool,
	  SaveRadRecCool;

	long int nlo_heat_max  , nhi_heat_max;

	/* place to put string labels for iso lines */
	char chLabel[NCOLNT_LAB_LEN+1];

	DEBUG_ENTRY( "iso_cool()" );

	/* validate the incoming data */
	ASSERT( nelem >= ipISO );
	ASSERT( ipISO < NISO );
	ASSERT( nelem < LIMELM );
	/* local number of levels may be less than malloced number if continuum
	 * has been lowered due to high density */
	ASSERT( iso.numLevels_local[ipISO][nelem] <= iso.numLevels_max[ipISO][nelem] );

	if( dense.xIonDense[nelem][nelem-ipISO]<=0. ||
		!dense.lgElmtOn[nelem] )
	{
		/* all global variables must be zeroed here or set below */
		iso.coll_ion[ipISO][nelem] = 0.;
		iso.cLya_cool[ipISO][nelem] = 0.;
		iso.cLyrest_cool[ipISO][nelem] = 0.;
		iso.cBal_cool[ipISO][nelem] = 0.;
		iso.cRest_cool[ipISO][nelem] = 0.;
		iso.xLineTotCool[ipISO][nelem] = 0.;
		iso.RadRecCool[ipISO][nelem] = 0.;
		iso.FreeBnd_net_Cool_Rate[ipISO][nelem] = 0.;
		iso.dLTot[ipISO][nelem] = 0.;
		iso.RecomInducCool_Rate[ipISO][nelem] = 0.;
		return;
	}

	/* create some space, these go to numLevels_local instead of _max
	 * since continuum may have been lowered by density */
	if( lgPrintIonizCooling )
	{
		CollisIonizatCoolingArr.resize( iso.numLevels_local[ipISO][nelem] );
		CollisIonizatCoolingDTArr.resize( iso.numLevels_local[ipISO][nelem] );
	}
	SavePhotoHeat.resize( iso.numLevels_local[ipISO][nelem] );
	SaveInducCool.resize( iso.numLevels_local[ipISO][nelem] );
	SaveRadRecCool.resize( iso.numLevels_local[ipISO][nelem] );

	/***********************************************************************
	 *                                                                     *
	 * collisional ionization cooling, less three-body recombination  heat *
	 *                                                                     *
	 ***********************************************************************/

	/* will be net collisional ionization cooling, units erg/cm^3/s */
	CollisIonizatCoolingTotal = 0.;
	dCollisIonizatCoolingTotalDT = 0.;

	/* collisional ionization cooling minus three body heating 
	 * depending on how topoff is done, highest level can have large population
	 * and its coupling to continuum can be large, at various times code 
	 * had to ignore effects of very highest level, but starting in mid
	 * 20006 all levels have been included 
	 * 2008 apr 18, do not include highest - when only 2 collapsed levels
	  * are used several density 13 BLR sims have serious convergence
	  * problems */
	for( n=0; n < iso.numLevels_local[ipISO][nelem]-1; ++n )
	{
		/* total collisional ionization cooling less three body heating */
		CollisIonizatCooling = 
			EN1RYD*iso.xIsoLevNIonRyd[ipISO][nelem][n]*iso.ColIoniz[ipISO][nelem][n]*dense.EdenHCorr*
		  	(StatesElemNEW[nelem][nelem-ipISO][n].Pop -iso.PopLTE[ipISO][nelem][n]*dense.eden*
			dense.xIonDense[nelem][nelem+1-ipISO]);
		CollisIonizatCoolingTotal += CollisIonizatCooling;

		/* the derivative of the cooling */
		CollisIonizatCoolingDT = CollisIonizatCooling*
			(iso.xIsoLevNIonRyd[ipISO][nelem][n]*TE1RYD/POW2(phycon.te)- thermal.halfte);


		dCollisIonizatCoolingTotalDT += CollisIonizatCoolingDT;
		// save values for debug printout
		if( lgPrintIonizCooling )
		{
			CollisIonizatCoolingArr[n] = CollisIonizatCooling;
			CollisIonizatCoolingDTArr[n] = CollisIonizatCoolingDT;
		}
	}

	/* save net collisional ionization cooling less H-3 body heating
	 * for inclusion in printout */
	iso.coll_ion[ipISO][nelem] = CollisIonizatCoolingTotal;

	/* add this derivative to total */
	thermal.dCooldT += dCollisIonizatCoolingTotalDT;

	/* create a meaningful label */
	sprintf(chLabel , "IScion%2s%2s" , elementnames.chElementSym[ipISO] , 
		elementnames.chElementSym[nelem] );
	/* dump the coolant onto the cooling stack */
	CoolAdd(chLabel,(realnum)nelem,MAX2(0.,iso.coll_ion[ipISO][nelem]));

	/* heating[0][3] is all three body heating, opposite of collisional 
	 * ionization cooling,
	 * would be unusual for this to be non-zero since would require excited
	 * state departure coefficients to be greater than unity */
	thermal.heating[0][3] += MAX2(0.,-iso.coll_ion[ipISO][nelem]);

	/* debug block printing individual contributors to collisional ionization cooling */
	if( lgPrintIonizCooling && nelem==1 && ipISO==1 )
	{
		fprintf(ioQQQ,"DEBUG coll ioniz cool contributors:");
		for( n=0; n < iso.numLevels_local[ipISO][nelem]; n++ )
		{
			if( CollisIonizatCoolingArr[n] / SDIV( CollisIonizatCoolingTotal ) > 0.1 )
				fprintf(ioQQQ," %2li %.1e",
					n,
					CollisIonizatCoolingArr[n]/ SDIV( CollisIonizatCoolingTotal ) );
		}
		fprintf(ioQQQ,"\n");
		fprintf(ioQQQ,"DEBUG coll ioniz derivcontributors:");
		for( n=0; n < iso.numLevels_local[ipISO][nelem]; n++ )
		{
			if( CollisIonizatCoolingDTArr[n] / SDIV( dCollisIonizatCoolingTotalDT ) > 0.1 )
				fprintf(ioQQQ," %2li %.1e",
					n,
					CollisIonizatCoolingDTArr[n]/ SDIV( dCollisIonizatCoolingTotalDT ) );
		}
		fprintf(ioQQQ,"\n");
	}

	/***********************************************************************
	 *                                                                     *
	 * hydrogen recombination free-bound free bound cooling                *
	 *                                                                     *
	 ***********************************************************************/

	/* this is the product of the ion abundance times the electron density */
	edenIonAbund = dense.eden*dense.xIonDense[nelem][nelem+1-ipISO];

	/* now do case b sum to compare with exact value below */
	iso.RadRecCool[ipISO][nelem] = 0.;
	ThinCoolingSum = 0.;

	if( ipISO == ipH_LIKE )
	{
		/* do ground with special approximate fits to Ferland et al. '92 */
		thin = HydroRecCool(
			/* n is the prin quantum number on the physical scale */
			1 , 
			/* nelem is the charge on the C scale, 0 is hydrogen */
			nelem);
	}
	else
	{
		/* this is the cooling before correction for optical depths */
		 thin = iso.RadRecomb[ipISO][nelem][0][ipRecRad]*
			/* arg is the scaled temperature, T * n^2 / Z^2, 
			 * n is principal quantum number, Z is charge, 1 for H */
			HCoolRatio( 
			phycon.te * POW2( (double)StatesElemNEW[nelem][nelem-ipISO][0].n / (double)(nelem+1-ipISO) ))*
			/* convert results to energy per unit vol */
			phycon.te * BOLTZMANN;
	}
	/* the cooling, corrected for optical depth */
	SaveRadRecCool[0] = iso.RadRecomb[ipISO][nelem][0][ipRecNetEsc] * thin;
	/* this is now total free-bound cooling */
	iso.RadRecCool[ipISO][nelem] += SaveRadRecCool[0] * edenIonAbund;

	/* radiative recombination cooling for all excited states */
	for( n=1; n < iso.numLevels_local[ipISO][nelem]; n++ )
	{
		/* this is the cooling before correction for optical depths */
		 thin = iso.RadRecomb[ipISO][nelem][n][ipRecRad]*
			/* arg is the scaled temperature, T * n^2 / Z^2, 
			 * n is principal quantum number, Z is charge, 1 for H */
			HCoolRatio( 
			phycon.te * POW2( (double)StatesElemNEW[nelem][nelem-ipISO][n].n / (double)(nelem+1-ipISO) ))*
			/* convert results to energy per unit vol */
			phycon.te * BOLTZMANN;

		/* the cooling, corrected for optical depth */
		SaveRadRecCool[n] = iso.RadRecomb[ipISO][nelem][n][ipRecNetEsc] * thin;
		/* this is now total free-bound cooling */
		iso.RadRecCool[ipISO][nelem] += SaveRadRecCool[n] * edenIonAbund;

		/* keep track of case b sum for topoff below */
		ThinCoolingSum += thin;
	}
	{
		/* debug block for state specific recombination cooling */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC  )
		{
			if( nelem==ipISO )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( n=0; n < (iso.numLevels_local[ipISO][nelem] - 1); n++ )
				{
					fprintf(ioQQQ,"\t%.2f",SaveRadRecCool[n]/ThinCoolingSum);
				}
				fprintf(ioQQQ,"\n");
			}
		}
	}

	/* Case b sum of optically thin radiative recombination cooling. 
	 * add any remainder to the sum from above - high precision is needed 
	 * to get STE result to converge close to equilibrium - only done for
	 * H-like ions where exact result is known */
	if( ipISO == ipH_LIKE )
	{
		/* these expressions are only valid for hydrogenic sequence */
		if( nelem == 0 )
		{
			/*expansion for hydrogen itself */
			ThinCoolingCaseB = (-25.859117 + 
			0.16229407*phycon.telogn[0] + 
			0.34912863*phycon.telogn[1] - 
			0.10615964*phycon.telogn[2])/(1. + 
			0.050866793*phycon.telogn[0] - 
			0.014118924*phycon.telogn[1] + 
			0.0044980897*phycon.telogn[2] + 
			6.0969594e-5*phycon.telogn[3]);
		}
		else
		{
			/* same expansion but for hydrogen ions */
			ThinCoolingCaseB = (-25.859117 + 
			0.16229407*(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) + 
			0.34912863*POW2(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) - 
			0.10615964*powi( (phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]),3) )/(1. + 
			0.050866793*(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) - 
			0.014118924*POW2(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) + 
			0.0044980897*powi( (phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]),3) + 
			6.0969594e-5*powi( (phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]),4) );
		}

		/* now convert to linear cooling coefficient */
		ThinCoolingCaseB = POW3(1.+nelem-ipISO)*pow(10.,ThinCoolingCaseB)/(phycon.te/POW2(1.+nelem-ipISO) );

		/* this is the error, expect positive since do not include infinite number of levels */
		RecCoolExtra = ThinCoolingCaseB - ThinCoolingSum;
	}
	else
	{
		ThinCoolingCaseB = 0.;
		RecCoolExtra = 0.;
	}

	/* don't let the extra be negative - should be positive if H-like, negative for
	 * he-like only due to real difference in recombination coefficients */
	RecCoolExtra = MAX2(0., RecCoolExtra );

	/* add error onto total - this is significant for approach to STE */
	iso.RadRecCool[ipISO][nelem] += RecCoolExtra* edenIonAbund *iso.RadRecomb[ipISO][nelem][iso.numLevels_local[ipISO][nelem]-1][ipRecNetEsc];

	/***********************************************************************
	 *                                                                     
	 * heating  by photoionization of                                      
	 * excited states of all species                            
	 *                                                                     
	 ***********************************************************************/

	/* photoionization of excited levels */
	HeatExcited = 0.;
	ipbig = -1000;
	for( n=1; n < (iso.numLevels_local[ipISO][nelem] - 1); ++n )
	{
		ASSERT( iso.PhotoHeat[ipISO][nelem][n] >= 0. );
		ASSERT( StatesElemNEW[nelem][nelem-ipISO][n].Pop >= 0. );

		SavePhotoHeat[n] = StatesElemNEW[nelem][nelem-ipISO][n].Pop*
			iso.PhotoHeat[ipISO][nelem][n];
		HeatExcited += SavePhotoHeat[n];
		if( SavePhotoHeat[n] > biggest )
		{
			biggest = SavePhotoHeat[n];
			ipbig = (int)n;
		}
	}
	{
		/* debug block for heating due to photo of each n */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC  && ipISO==0 && nelem==0  && nzone > 700)
		{
			/* this was not done above */
			SavePhotoHeat[ipH1s] = StatesElemNEW[nelem][nelem-ipISO][ipH1s].Pop*
				iso.PhotoHeat[ipISO][nelem][ipH1s];
			fprintf(ioQQQ,"ipISO = %li nelem=%li ipbig=%li biggest=%g HeatExcited=%.2e ctot=%.2e\n",
				ipISO , nelem,
				ipbig , 
				biggest,
				HeatExcited ,
				thermal.ctot);
			/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
			for(n=ipH1s; n< (iso.numLevels_local[ipISO][nelem] - 1); ++n )
			{
				fprintf(ioQQQ,"DEBUG phot heat%2li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
					n,
					SavePhotoHeat[n]/HeatExcited,
					dense.xIonDense[nelem][nelem+1-ipISO],
					StatesElemNEW[nelem][nelem-ipISO][n].Pop,
					iso.PhotoHeat[ipISO][nelem][n],
					iso.gamnc[ipISO][nelem][n]);
			}
		}
	}

	/* FreeBnd_net_Cool_Rate is net cooling due to recombination 
	 * RadRecCool is total radiative recombination cooling sum to all levels,
	 * with n>=2 photoionization heating subtracted */
	iso.FreeBnd_net_Cool_Rate[ipISO][nelem] = iso.RadRecCool[ipISO][nelem] - HeatExcited;
	/*fprintf(ioQQQ,"DEBUG Hn2\t%.3e\t%.3e\n",
		-iso.RadRecCool[ipISO][nelem]/SDIV(thermal.htot),
		HeatExcited/SDIV(thermal.htot));*/

	/* heating[0][1] is all excited state photoionization heating from ALL 
	 * species, this is set to zero in CoolEvaluate before loop where this 
	 * is called. */
	thermal.heating[0][1] += MAX2(0.,-iso.FreeBnd_net_Cool_Rate[ipISO][nelem]);

	/* net free-bound cooling minus free-free heating */
	/* create a meaningful label */
	sprintf(chLabel , "ISrcol%2s%2s" , elementnames.chElementSym[ipISO]  , 
		elementnames.chElementSym[nelem]);
	CoolAdd(chLabel, (realnum)nelem, MAX2(0.,iso.FreeBnd_net_Cool_Rate[ipISO][nelem]));

	/* if rec coef goes as T^-.8, but times T, so propto t^.2 */
	thermal.dCooldT += 0.2*iso.FreeBnd_net_Cool_Rate[ipISO][nelem]*phycon.teinv;

	/***********************************************************************
	 *                                                                     *
	 * induced recombination cooling                                       *
	 *                                                                     *
	 ***********************************************************************/

	iso.RecomInducCool_Rate[ipISO][nelem] = 0.;
	/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
	for( n=0; n < (iso.numLevels_local[ipISO][nelem] - 1); ++n )
	{
		/* >>chng 02 jan 22, removed cinduc, replace with RecomInducCool */
		SaveInducCool[n] = iso.RecomInducCool_Coef[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n]*edenIonAbund;
		iso.RecomInducCool_Rate[ipISO][nelem] += SaveInducCool[n];
		thermal.dCooldT += SaveInducCool[n]*
			(iso.xIsoLevNIonRyd[ipISO][nelem][n]/phycon.te_ryd - 1.5)*phycon.teinv;
	}

	{
		/* print rec cool, induced rec cool, photo heat */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && ipISO==0 && nelem==5 /**/ )
		{
			fprintf(ioQQQ," ipISO=%li nelem=%li ctot = %.2e\n",
				ipISO,
				nelem,
				thermal.ctot);
			fprintf(ioQQQ,"sum\t%.2e\t%.2e\t%.2e\n", 
				HeatExcited,
				iso.RadRecCool[ipISO][nelem],
				iso.RecomInducCool_Rate[ipISO][nelem]);
			fprintf(ioQQQ,"sum\tp ht\tr cl\ti cl\n");

			/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
			for(n=0; n< (iso.numLevels_local[ipISO][nelem] - 1); ++n )
			{
				fprintf(ioQQQ,"%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e  \n", 
					n,
					SavePhotoHeat[n],
					SaveRadRecCool[n],
					SaveInducCool[n] ,
					iso.RecomInducCool_Coef[ipISO][nelem][n],
					iso.PopLTE[ipISO][nelem][n],
					iso.RecomInducRate[ipISO][nelem][n]);
			}
			fprintf(ioQQQ," \n");
		}
	}
	/* create a meaningful label - induced rec cooling */
	sprintf(chLabel , "ISicol%2s%2s" , elementnames.chElementSym[ipISO]  , 
		elementnames.chElementSym[nelem]);
	/* induced rec cooling */
	CoolAdd(chLabel,(realnum)nelem,iso.RecomInducCool_Rate[ipISO][nelem]);

	/* find total collisional energy exchange due to bound-bound */
	iso.xLineTotCool[ipISO][nelem] = 0.;
	dCdT_all = 0.;
	heat_max = 0.;
	nlo_heat_max = -1;
	nhi_heat_max = -1;

	/* loop does not include highest levels - their population may be
	 * affected by topoff */
	for( ipLo=0; ipLo < iso.numLevels_local[ipISO][nelem]-2; ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < iso.numLevels_local[ipISO][nelem]-1; ipHi++ )
		{
			/* collisional energy exchange between ipHi and ipLo - net cool */
			hlone = 
			  Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL*
			  (StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop*
			  iso.Boltzmann[ipISO][nelem][ipHi][ipLo]*
			  StatesElemNEW[nelem][nelem-ipISO][ipHi].g/StatesElemNEW[nelem][nelem-ipISO][ipLo].g - 
			  StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop)*dense.EdenHCorr*
			  Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg;

			if( hlone > 0. )
				Transitions[ipISO][nelem][ipHi][ipLo].Coll.cool = hlone;
			else
				Transitions[ipISO][nelem][ipHi][ipLo].Coll.heat = -1.*hlone;

			iso.xLineTotCool[ipISO][nelem] += hlone;

			/* next get derivative */
			if( hlone > 0. )
			{
				/* usual case, this line was a net coolant - for derivative 
				 * taking the exponential from ground gave better
				 * representation of effects */
				dCdT_all += hlone*
				  (Transitions[ipISO][nelem][ipHi][ipH1s].EnergyK*thermal.tsq1 - thermal.halfte);
			}
			else
			{
				/* this line heats the gas, remember which one it was,
				 * the following three vars never are used, but could be for
				 * debugging */
				if( hlone < heat_max )
				{
					heat_max = hlone;
					nlo_heat_max = ipLo;
					nhi_heat_max = ipHi;
				} 
				dCdT_all -= hlone*thermal.halfte;
			}
		}
	}
	{
		/* debug block announcing which line was most important */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			if( nelem==ipISO )
				fprintf(ioQQQ,"%li %li %.2e\n", nlo_heat_max, nhi_heat_max, heat_max );
		}
	}

	/* Lyman line collisional heating/cooling */
	/* Lya itself */
	iso.cLya_cool[ipISO][nelem] = 0.;
	/* lines higher than Lya */
	iso.cLyrest_cool[ipISO][nelem] = 0.;

	for( ipHi=1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
	{
		hlone = Transitions[ipISO][nelem][ipHi][ipH1s].Coll.ColUL*
		  (StatesElemNEW[nelem][nelem-ipISO][0].Pop*iso.Boltzmann[ipISO][nelem][ipHi][0]*
		  StatesElemNEW[nelem][nelem-ipISO][ipHi].g/StatesElemNEW[nelem][nelem-ipISO][0].g - 
		  StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop)* dense.EdenHCorr*
		  Transitions[ipISO][nelem][ipHi][0].EnergyErg;

		if( ipHi == iso.nLyaLevel[ipISO] )
		{
			iso.cLya_cool[ipISO][nelem] = hlone;

		}
		else
		{
			/* sum energy in higher lyman lines */
			iso.cLyrest_cool[ipISO][nelem] += hlone;
		}
	}

	/* Balmer line collisional heating/cooling and derivative 
	 * only used for print out, incorrect if not H */
	iso.cBal_cool[ipISO][nelem] = 0.;
	for( ipHi=3; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
	{
		/* single line cooling */
		ipLo = ipH2s;
		hlone = 
		  Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL*(
		  StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop*
		  iso.Boltzmann[ipISO][nelem][ipHi][ipLo]*
		  StatesElemNEW[nelem][nelem-ipISO][ipHi].g/StatesElemNEW[nelem][nelem-ipISO][ipLo].g - 
		  StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop)*dense.EdenHCorr*
		  Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg;

		ipLo = ipH2p;
		hlone += 
		  Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL*(
		  StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop*
		  iso.Boltzmann[ipISO][nelem][ipHi][ipLo]*
		  StatesElemNEW[nelem][nelem-ipISO][ipHi].g/StatesElemNEW[nelem][nelem-ipISO][ipLo].g - 
		  StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop)*dense.EdenHCorr*
		  Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg;

		/* this is only used to add to line array */
		iso.cBal_cool[ipISO][nelem] += hlone;
	}

	/* all hydrogen lines from Paschen on up, but not Balmer 
	 * only used for printout, incorrect if not H */
	iso.cRest_cool[ipISO][nelem] = 0.;
	for( ipLo=3; ipLo < iso.numLevels_local[ipISO][nelem]-1; ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
		{
			hlone = 
			  Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL*(
			  StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop*
			  iso.Boltzmann[ipISO][nelem][ipHi][ipLo]*
			  StatesElemNEW[nelem][nelem-ipISO][ipHi].g/StatesElemNEW[nelem][nelem-ipISO][ipLo].g - 
			  StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop)*dense.EdenHCorr*
			  Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg;

			iso.cRest_cool[ipISO][nelem] += hlone;
		}
	}

	/* add total line heating or cooling into stacks, derivatives */
	/* line energy exchange can be either heating or coolant
	 * must add this to total heating or cooling, as appropriate */
	/* create a meaningful label */
	sprintf(chLabel , "ISclin%2s%2s" , elementnames.chElementSym[ipISO] , 
		elementnames.chElementSym[nelem]);
	if( iso.xLineTotCool[ipISO][nelem] > 0. )
	{
		/* species is a net coolant label gives iso sequence, "wavelength" gives element */
		CoolAdd(chLabel,(realnum)nelem,iso.xLineTotCool[ipISO][nelem]);
		thermal.dCooldT += dCdT_all;
		iso.dLTot[ipISO][nelem] = 0.;
	}
	else
	{
		/* species is a net heat source, thermal.heating[0][23]was set to 0 in CoolEvaluate*/
		thermal.heating[0][23] -= iso.xLineTotCool[ipISO][nelem];
		CoolAdd(chLabel,(realnum)nelem,0.);
		iso.dLTot[ipISO][nelem] = -dCdT_all;
	}

	{
		/* debug print for understanding contributors to HLineTotCool */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			if( nelem == 0 )
			{
				fprintf(ioQQQ,"%.2e la %.2f restly %.2f barest %.2f hrest %.2f\n",
					iso.xLineTotCool[ipISO][nelem] ,
					iso.cLya_cool[ipISO][nelem]/iso.xLineTotCool[ipISO][nelem] ,
					iso.cLyrest_cool[ipISO][nelem]/iso.xLineTotCool[ipISO][nelem] ,
					iso.cBal_cool[ipISO][nelem]/iso.xLineTotCool[ipISO][nelem] ,
					iso.cRest_cool[ipISO][nelem]/iso.xLineTotCool[ipISO][nelem] );
			}
		}
	}
	{
		/* this is an option to print out each rate affecting each level in strict ste,
		 * this is a check that the rates are indeed in detailed balance */
		enum {DEBUG_LOC=false};
		enum {LTEPOP=true};
		/* special debug print for gas near STE */
		if( DEBUG_LOC  && (nelem==1 || nelem==0) )
		{
			/* LTEPOP is option to assume STE for rates */
			if( LTEPOP )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for(n=ipH1s; n<iso.numLevels_local[ipISO][nelem]-1; ++n )
				{
					fprintf(ioQQQ,"%li\t%li\t%g\t%g\t%g\t%g\tT=\t%g\t%g\t%g\t%g\n", nelem,n,
						iso.gamnc[ipISO][nelem][n] *iso.PopLTE[ipISO][nelem][n], 
						/* induced recom, intergral times hlte */
						/*iso.RadRecomb[ipISO][nelem][n][ipRecRad]+iso.rinduc[ipISO][nelem][n]  ,*/
						/* >>chng 02 jan 22, remove rinduc, replace with RecomInducRate */
						iso.RadRecomb[ipISO][nelem][n][ipRecRad]+ 
							iso.RecomInducRate[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n]  ,
						/* induced rec */
						iso.RecomInducRate[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n]  ,
						/* spontaneous recombination */
						iso.RadRecomb[ipISO][nelem][n][ipRecRad] ,
						/* heating, followed by two processes that must balance it */
						iso.PhotoHeat[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n], 
						iso.RecomInducCool_Coef[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n]+SaveRadRecCool[n] ,
						/* induced rec cooling, integral times hlte */
						iso.RecomInducCool_Coef[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n] ,
						/* free-bound cooling per unit vol */
						SaveRadRecCool[n] );
				}
			}
			else
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for(n=ipH1s; n<iso.numLevels_local[ipISO][nelem]-1; ++n )
				{
					fprintf(ioQQQ,"%li\t%li\t%g\t%g\t%g\t%g\tT=\t%g\t%g\t%g\t%g\n", nelem,n,
						iso.gamnc[ipISO][nelem][n]*StatesElemNEW[nelem][nelem-ipISO][n].Pop, 
						/* induced recom, intergral times hlte */
						iso.RadRecomb[ipISO][nelem][n][ipRecRad]*edenIonAbund+
							iso.RecomInducRate[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n] *edenIonAbund ,
						iso.RadRecomb[ipISO][nelem][n][ipRecRad]*edenIonAbund ,
						iso.RecomInducRate[ipISO][nelem][n]*iso.PopLTE[ipISO][nelem][n] *edenIonAbund ,
						/* heating, followed by two processes that must balance it */
						SavePhotoHeat[n], 
						SaveInducCool[n]+SaveRadRecCool[n]*edenIonAbund ,
						/* induced rec cooling, integral times hlte */
						SaveInducCool[n] ,
						/* free-bound cooling per unit vol */
						SaveRadRecCool[n]*edenIonAbund );
				}
			}
		}
	}
	return;
}
#if defined(__HP_aCC)
#pragma OPTIMIZE OFF
#pragma OPTIMIZE ON
#endif
