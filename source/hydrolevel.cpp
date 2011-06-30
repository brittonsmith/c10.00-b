/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroLevel calls iso_level to solve for ionization balance 
 * and level populations of model hydrogen atom */
/*PrtHydroTrace1 print trace info for hydrogen-like species */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "dense.h"
#include "secondaries.h"
#include "trace.h"
#include "phycon.h"
#include "ionbal.h"
#include "hydrogenic.h"

/*PrtHydroTrace1a print trace info for hydrogen-like species */
STATIC void PrtHydroTrace1a(long nelem )
{
	double colfrc, 
	  phtfrc, 
	  secfrc;

	DEBUG_ENTRY( "PrtHydroTrace1a()" );

	/* identify how atom is ionized for full trace */
	if( iso.xIonSimple[ipH_LIKE][nelem] > 0. )
	{
		/* fraction of ionization due to photoionization */
		phtfrc = iso.gamnc[ipH_LIKE][nelem][ipH1s]/((dense.eden*(iso.RadRec_effec[ipH_LIKE][nelem] + 
			ionbal.CotaRate[nelem]) )*
			iso.xIonSimple[ipH_LIKE][nelem]);

		/* fraction of ionization due to collisional ionization */
		colfrc = (iso.ColIoniz[ipH_LIKE][nelem][ipH1s]*dense.EdenHCorr )/
			((dense.eden*(iso.RadRec_effec[ipH_LIKE][nelem] + 
			ionbal.CotaRate[0]) )*
			iso.xIonSimple[ipH_LIKE][nelem]);

		/* fraction of ionization due to secondary ionization */
		secfrc = secondaries.csupra[nelem][nelem]/((dense.eden*(iso.RadRec_effec[ipH_LIKE][nelem] + 
			ionbal.CotaRate[0]) )*
			iso.xIonSimple[ipH_LIKE][nelem]);
	}
	else
	{
		phtfrc = 0.;
		colfrc = 0.;
		secfrc = 0.;
	}

	fprintf( ioQQQ, "     HydroLevel Z=%2ld called, simple II/I=",nelem);
	PrintE93( ioQQQ, iso.xIonSimple[ipH_LIKE][nelem]);
	fprintf( ioQQQ," PhotFrc:");
	PrintE93( ioQQQ,phtfrc);
	fprintf(ioQQQ," ColFrc:");
	PrintE93( ioQQQ,colfrc);
	fprintf( ioQQQ," SecFrc");
	PrintE93( ioQQQ, secfrc);
	fprintf( ioQQQ,"  Te:");
	PrintE93( ioQQQ,phycon.te);
	fprintf( ioQQQ," eden:");
	PrintE93( ioQQQ,dense.eden);
	fprintf( ioQQQ,"\n"); 
	return;
}

/*PrtHydroTrace1 print trace info for hydrogen-like species */
STATIC void PrtHydroTrace1(long nelem )
{
	long int ipHi , ipLo , i;

	DEBUG_ENTRY( "PrtHydroTrace1()" );

	fprintf( ioQQQ, 
		"       HydroLevel%3ld finds arrays, with optical depths defined? %li induced 2ph=%12.3e\n", 
		nelem, iteration, Transitions[ipH_LIKE][nelem][ipH2s][ipH1s].Emis->pump );
	/* 06 aug 28, from numLevels_max to _local. */
	for( ipHi=ipH2s; ipHi < iso.numLevels_local[ipH_LIKE][nelem]; ipHi++ )
	{
		fprintf( ioQQQ, "up:%2ld", ipHi );
		fprintf( ioQQQ, "lo" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ, "%9ld", ipLo );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " A*esc" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Aul*
				Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pesc ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " A*ees" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Aul*
				Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pelec_esc ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " tauin" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->TauIn ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " t tot" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->TauTot ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " Esc  " );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pesc ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " Eesc " );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pelec_esc ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " Dest " );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pdest) );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " A*dst" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Aul*
				Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pdest ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " StrkE" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso.pestrk[ipH_LIKE][nelem][ipLo][ipHi] ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " B(ul)" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->pump*
				StatesElemNEW[nelem][nelem-ipH_LIKE][ipLo].g/StatesElemNEW[nelem][nelem-ipH_LIKE][ipHi].g ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " tcont" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->TauCon ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " C(ul)" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", 
				Transitions[ipH_LIKE][nelem][ipHi][ipLo].Coll.ColUL*dense.eden ));
		}
		fprintf( ioQQQ, "\n" );

		if( ipHi == 2 )
		{
			fprintf( ioQQQ, "    FeIIo");
			fprintf( ioQQQ,PrintEfmt("%9.2e", 
				hydro.dstfe2lya* Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Aul ));
			fprintf( ioQQQ, "\n");
		}
	}

	fprintf( ioQQQ, "         " );
	/* 06 aug 28, from numLevels_max to _local. */
	for( i=1; i < iso.numLevels_local[ipH_LIKE][nelem]; i++ )
	{
		fprintf( ioQQQ, "%9ld", i );
	}
	fprintf( ioQQQ, "\n" );
	return;
}

/*HydroLevel calls iso_level to solve for ionization balance 
 * and level populations of model hydrogen atom */
void HydroLevel(long int nelem)
{
	long int i; 
	int ipISO = ipH_LIKE;

	DEBUG_ENTRY( "HydroLevel()" );

	/* check that we were called with valid charge */
	ASSERT( nelem >= 0);
	ASSERT( nelem < LIMELM );

	/* option to print some rates */
	if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO]) && (nelem == trace.ipIsoTrace[ipISO]) )
		PrtHydroTrace1(nelem);

	if( trace.lgHBug && trace.lgTrace )
		PrtHydroTrace1a(nelem);

	/* this is main trace h-like printout */
	if( (trace.lgIsoTraceFull[ipISO] && trace.lgTrace) && (nelem == trace.ipIsoTrace[ipISO]) )
	{
		fprintf( ioQQQ, "       HLEV HGAMNC" );
		PrintE93( ioQQQ, iso.gamnc[ipISO][nelem][ipH1s] );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH2s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso.gamnc[ipISO][nelem][i] ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV TOTCAP" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=1; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso.RateCont2Level[ipISO][nelem][i]/dense.eden ));
		}
		fprintf( ioQQQ," tot");
		fprintf( ioQQQ,PrintEfmt("%10.2e", ionbal.RateRecomTot[nelem][nelem-ipISO]/dense.eden ) );
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV IND Rc" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso.RecomInducRate[ipISO][nelem][i]*iso.PopLTE[ipISO][nelem][i] ));
		}
		fprintf( ioQQQ, "\n" );

		/* incuded recombination rate coefficients */
		fprintf( ioQQQ, "       IND Rc LTE " );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e",
				iso.gamnc[ipISO][nelem][i]*iso.PopLTE[ipISO][nelem][i] ));
		}
		fprintf( ioQQQ, "\n" );

		/* LTE level populations */
		fprintf( ioQQQ, "       HLEV   HLTE" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso.PopLTE[ipISO][nelem][i] ));
		}
		fprintf( ioQQQ, "\n" );

		/* fraction of total ionization due to collisions from given level */
		fprintf( ioQQQ, "       HLEVfr cion" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", 
				iso.ColIoniz[ipISO][nelem][i]*
			  StatesElemNEW[nelem][nelem-ipISO][i].Pop*dense.EdenHCorr/MAX2(1e-30,iso.RateLevel2Cont[ipISO][nelem][i]) ) );
		}
		fprintf( ioQQQ, "\n" );

		/* fraction of total ionization due to photoionization from given level */
		if( ionbal.RateRecomTot[nelem][nelem]> 0. )
		{
			fprintf( ioQQQ, "       HLEVfrPhIon" );
			/* 06 aug 28, from numLevels_max to _local. */
			for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e", 
					iso.gamnc[ipISO][nelem][i]*StatesElemNEW[nelem][nelem-ipISO][i].Pop/MAX2(1e-30,iso.RateLevel2Cont[ipISO][nelem][i]) ) );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, "       HLEV     HN" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", StatesElemNEW[nelem][nelem-ipISO][i].Pop ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV   b(n)" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso.numLevels_local[ipISO][nelem]; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso.DepartCoef[ipISO][nelem][i] ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV   X12tot");
		fprintf(ioQQQ,PrintEfmt("%9.2e", secondaries.x12tot ));
		fprintf( ioQQQ," Grn dest:");
		fprintf(ioQQQ,PrintEfmt("%9.2e",
		  ionbal.RateIoniz[nelem][nelem][nelem+1] ));
		fprintf(ioQQQ, "\n"); 
	}


	/* >>chng 05 mar 24,
	 * renormalize the populations and emission of H atom to agree with chemistry 
	 * these were not being kept parallel with chemistry, and caused large changes in O+
	 * abundance when finally done */
	fixit();  /* this probably does not need to be called here. */
	if( nelem == ipHYDROGEN )
		HydroRenorm();

	if( trace.lgTrace )
	{
		/* iso.RecomTotal[nelem] is gross rec coef, computed here while filling in matrix
		 * elements, all physical processes included. 
		 * RadRec_effec is total effective radiative only */
		fprintf( ioQQQ, "       HydroLevel Z:%2ld return %s te=",
			nelem, 
			iso.chTypeAtomUsed[ipISO][nelem] );
		PrintE93( ioQQQ,phycon.te);
		fprintf( ioQQQ," density=%.4e", dense.xIonDense[nelem][nelem-ipISO] );

		fprintf( ioQQQ," simple=%.4e",iso.xIonSimple[ipISO][nelem]);

		fprintf( ioQQQ," b1=%.2e",iso.DepartCoef[ipISO][nelem][ipH1s]);

		fprintf( ioQQQ," ion rate=%.4e",ionbal.RateIonizTot(nelem,nelem-ipISO) );

		fprintf( ioQQQ," TotRec=%.4e",ionbal.RateRecomTot[nelem][nelem-ipISO]);

		fprintf( ioQQQ," RadRec=%.4e",iso.RadRec_effec[ipISO][nelem]);
		fprintf( ioQQQ, "\n");
	}
	return;
}
