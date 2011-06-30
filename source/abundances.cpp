/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*AbundancesPrt print all abundances, both gas phase and grains */
/*AbundancesSet sets initial abundances after parameters are entered by reading input */
/*AbundancesTable interpolate on table of points to do 'element table' command, */
/*PrtElem print chemical composition at start of calculation */
#include "cddefines.h"
#include "physconst.h"
#include "phycon.h"
#include "called.h"
#include "stopcalc.h"
#include "thermal.h"
#include "trace.h"
#include "elementnames.h"
#include "dense.h"
#include "radius.h"
#include "grainvar.h"
#include "abund.h"

/*PrtElem print chemical composition at start of calculation */
STATIC void PrtElem(
  /* the job to do, the options are "init", "fill", "flus" */
  const char *chJob, 
  /* label for the element */
  const char *chLabl, 
  /* its abundance */
  double abund_prt);

/*AbundancesPrt print all abundances, both gas phase and grains */
void AbundancesPrt( void )
{
	long int i;
	double GrainNumRelHydrSilicate ,
		GrainNumRelHydrCarbonaceous ,
		GrainNumRelHydr_PAH,
		GrainMassRelHydrSilicate,
		GrainMassRelHydrCarbonaceous,
		GrainMassRelHydr_PAH;

	DEBUG_ENTRY( "AbundancesPrt()" );

	/* this is main loop to print abundances of each element */
	if( called.lgTalk )
	{
		PrtElem("initG","  ",0.);/* initialize print routine for gas*/
		for( i=0; i < LIMELM; i++ )
		{
			if( dense.lgElmtOn[i] )
			{
				/* fill in print buffer with abundances */
				PrtElem("fill",(char*)elementnames.chElementSym[i],
				  abund.solar[i]);
			}
		}

		/* flush the print buffer */
		PrtElem("flus","  ",0.);
		/* final carriage return */
		fprintf( ioQQQ, " \n" );

		/* now grains if present */
		if( gv.lgDustOn() )
		{
			/* we will first print the total abundances of each element locked up in grains */
			/* initialize print routine for dust*/
			PrtElem("initD","  ",0.);
			for( i=0; i < LIMELM; i++ )
			{
				if( gv.elmSumAbund[i]>SMALLFLOAT )
				{
					/* fill in print buffer with abundances */
					PrtElem("fill",(char*)elementnames.chElementSym[i],
						gv.elmSumAbund[i]/dense.gas_phase[ipHYDROGEN]);
				}
			}
			/* flush the print buffer */
			PrtElem("flus","  ",0.);
			/* final carriage return */
			fprintf( ioQQQ, " \n" );

			/* this is used to store grain number density per hydrogen */
			GrainNumRelHydrSilicate = 0.;
			GrainNumRelHydrCarbonaceous = 0;
			GrainNumRelHydr_PAH = 0.;
			GrainMassRelHydrSilicate = 0.;
			GrainMassRelHydrCarbonaceous = 0;
			GrainMassRelHydr_PAH = 0.;

			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{

				/* number density of grains per hydrogen, the ratio
				 * gv.bin[nd]->IntVol/gv.bin[nd]->AvVol is the number of grain particles 
				 * per H at standard grain abundance*/
				realnum DensityNumberPerHydrogen = 
					(gv.bin[nd]->IntVol/gv.bin[nd]->AvVol)*gv.bin[nd]->dstAbund / 
					gv.bin[nd]->GrnDpth;
				/* mass of grains per hydrogen */
				realnum DensityMassPerHydrogen = 
					gv.bin[nd]->IntVol*gv.bin[nd]->dustp[0]*gv.bin[nd]->dstAbund/
					(realnum)ATOMIC_MASS_UNIT / gv.bin[nd]->GrnDpth;

				/* >>chng 06 mar 05, fix expression for calculating grain number density, PvH */
				if( gv.bin[nd]->matType == MAT_CAR || gv.bin[nd]->matType == MAT_CAR2 )
				{
					/* carbonaceous grains */
					GrainNumRelHydrCarbonaceous += DensityNumberPerHydrogen;
					GrainMassRelHydrCarbonaceous += DensityMassPerHydrogen;
				}
				else if( gv.bin[nd]->matType == MAT_SIL || gv.bin[nd]->matType == MAT_SIL2 )
				{
					/* silicate grains */
					GrainNumRelHydrSilicate += DensityNumberPerHydrogen;
					GrainMassRelHydrSilicate += DensityMassPerHydrogen;
				}
				else if( gv.bin[nd]->matType == MAT_PAH || gv.bin[nd]->matType == MAT_PAH2 )
				{
					/* PAHs - full abundance - remove possible factor accounting for
					 * variation of abundances with physical conditions - this will 
					 * be the PAH abundance with scale factor of unity */
					GrainNumRelHydr_PAH += DensityNumberPerHydrogen;
					GrainMassRelHydr_PAH += DensityMassPerHydrogen;
				}
				else
					TotalInsanity();
			}

			/* now print total number of grains of each type */
			fprintf(ioQQQ,"              Number of grains per hydrogen (scale=1)                         Mass of grains per hydrogen (scale=1)\n");
			fprintf(ioQQQ,"        Carbonaceous: %.3f  Silicate: %.3f  PAH: %.3f         Carbonaceous: %.3f  Silicate: %.3f  PAH: %.3f\n\n" ,
				log10( MAX2( 1e-30, GrainNumRelHydrCarbonaceous ) ) ,
				log10( MAX2( 1e-30, GrainNumRelHydrSilicate ) ) ,
				log10( MAX2( 1e-30, GrainNumRelHydr_PAH ) ) ,
				log10( MAX2( 1e-30, GrainMassRelHydrCarbonaceous ) ) ,
				log10( MAX2( 1e-30, GrainMassRelHydrSilicate ) ) ,
				log10( MAX2( 1e-30, GrainMassRelHydr_PAH ) )  );
		}
	}
	return;
}

/*AbundancesSet print all abundances, both gas phase and grains */
void AbundancesSet(void)
{
	long int i, 
	  nelem;
	double fac;
	static bool lgFirstCall=true;
	static bool lgElOnOff[LIMELM];

	DEBUG_ENTRY( "AbundancesSet()" );

	/* if this is the first call to this routine in this core load, 
	 * save the state of the lgElmOn array, so that it is possible
	 * to turn off elements in later models, but not turn on an
	 * element that was initially turned off.  This is necessary since
	 * the Create... routines that create space for elements will
	 * not be revisited in later models.  You can turn off an initially
	 * enabled element, but not turn a disabled one on.  */

	if( lgFirstCall )
	{
		/* first call - save the initial state of the lgElmtOn vector */
		for( i=0; i<LIMELM; ++i )
		{
			lgElOnOff[i] = dense.lgElmtOn[i];
		}
	}
	lgFirstCall = false;

	/* make sure that initially false elements remain off, while letting 
	 * enabled elements be turned off */
	for( i=ipHYDROGEN; i<LIMELM; ++i )
	{
		dense.lgElmtOn[i] = lgElOnOff[i] && dense.lgElmtOn[i];
	}

	/* rescale so that abundances are H=1 */
	for( i=ipHELIUM; i < LIMELM; i++ )
	{
		abund.solar[i] /= abund.solar[0];
	}
	abund.solar[ipHYDROGEN] = 1.;

	/* set current abundances to "solar" times metals scale factor
	 * and grain depletion factor */
	abund.solar[ipHELIUM] *= abund.depset[1]*abund.ScaleElement[1];

	/* option for density or abundance variations, this flag is true by default,
	 * set in zero, but set false if variations are enabled AND these
	 * are not density variations, but rather abundances */
	if( dense.lgDenFlucOn )
	{
		/* usual case - either density fluctuations or none at all */
		fac = 1.;
	}
	else
	{
		/* abundance fluctuations enabled, set initial value */
		fac = dense.cfirst*cos(dense.flcPhase) + dense.csecnd;
	}

	for( i=ipLITHIUM; i < LIMELM; i++ )
	{
		abund.solar[i] *= (realnum)(abund.ScaleMetals*abund.depset[i]*
		  abund.ScaleElement[i]*fac);
	}

	/* now fix abundance of any element with element table set */
	if( abund.lgAbTaON )
	{
		for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
		{
			if( abund.lgAbunTabl[nelem] )
			{
				abund.solar[nelem] = (realnum)(AbundancesTable(radius.Radius,
				  radius.depth,nelem+1));
			}
		}
	}

	/* dense.gas_phase[nelem] contains total abundance of element */
	/* the density of hydrogen itself has already been set at this point -
	 * it is set when commands parsed, most likely by the hden command -
	 * set all heavier elements */
	for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
	{
		/* this implements the element off limit xxx command, where
		 * xxx is the limit to the smallest n(A)/n(H) that will remain on */
		if( abund.solar[nelem] < dense.AbundanceLimit )
			dense.lgElmtOn[nelem] = false;

		if( dense.lgElmtOn[nelem] )
		{
			dense.gas_phase[nelem] = abund.solar[nelem]*dense.gas_phase[ipHYDROGEN];
			if( dense.gas_phase[nelem] <= 0. )
			{
				fprintf( ioQQQ, " Abundances must be greater than zero.  "
					"Check entered abundance for element%3ld  = %2.2s\n", 
				nelem, elementnames.chElementSym[nelem] );
				cdEXIT(EXIT_FAILURE);
			}
			else if( dense.gas_phase[nelem] < SMALLFLOAT )
			{
				fprintf(ioQQQ," Abundance for %s is %.2e, less than lower "
					"limit of %.3e, so turning element off.\n",
					elementnames.chElementSym[nelem],
					dense.gas_phase[nelem],
					SMALLFLOAT );
				dense.lgElmtOn[nelem] = false;
			}
			else if( dense.gas_phase[nelem] > MAX_DENSITY )
			{
				fprintf(ioQQQ," Abundance for %s is %.2e.  This version of Cloudy does not "
					"permit densities greater than %e cm-3.\n",
					elementnames.chElementSym[nelem],
					dense.gas_phase[nelem],
					MAX_DENSITY );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			/* >>chng 04 apr 20, set to zero if element is off */
			dense.gas_phase[nelem] = 0.;
		}
	}

	/* if stop temp set below default then we are going into cold and possibly 
	 * molecular gas - check some parameters in this case */
	if( called.lgTalk && (StopCalc.TempLoStopZone < phycon.TEMP_STOP_DEFAULT || 
		/* thermal.ConstTemp def is zero, set pos when used */
		(thermal.ConstTemp > 0. && thermal.ConstTemp < phycon.TEMP_STOP_DEFAULT ) ) )
	{

		/* print warning if temperature set below default but C > O */
		if( dense.lgElmtOn[ipOXYGEN] && dense.gas_phase[ipCARBON]/SDIV( dense.gas_phase[ipOXYGEN]) >= 1. )
		{
			fprintf( ioQQQ, "\n >>> \n"
							" >>> The simulation is going into possibly molecular gas but the carbon/oxygen abundance ratio is greater than unity.\n" );
			fprintf( ioQQQ, " >>> Standard interstellar chemistry networks are designed for environments with C/O < 1.\n" );
			fprintf( ioQQQ, " >>> The chemistry network may (or may not) collapse deep in molecular regions where CO is fully formed.\n" );
			fprintf( ioQQQ, " >>> \n\n\n\n\n" );
		}
	}

	if( trace.lgTrace )
	{
		realnum sumx , sumy , sumz = 0.;

		sumx = dense.gas_phase[ipHYDROGEN]*dense.AtomicWeight[ipHYDROGEN];
		sumy = dense.gas_phase[ipHELIUM]*dense.AtomicWeight[ipHELIUM];

		fprintf( ioQQQ, "\n AbundancesSet sets following densities (cm^-3); \n" );
		for( i=0; i<3; i++ )
		{
			for( nelem=i*10; nelem < i*10+10; nelem++ )
			{
				fprintf( ioQQQ, " %2.2s", elementnames.chElementSym[nelem] );
				PrintE82( ioQQQ, dense.gas_phase[nelem] );
				if( nelem>ipHELIUM )
					sumz += dense.gas_phase[nelem]*dense.AtomicWeight[nelem];
			}
			fprintf( ioQQQ, " \n" );
		}
		fprintf( ioQQQ, "\n AbundancesSet sets following abundances rel to H; \n" );
		for( i=0; i<3; i++ )
		{
			for( nelem=i*10; nelem < i*10+10; nelem++ )
			{
				fprintf( ioQQQ, " %2.2s", elementnames.chElementSym[nelem] );
				PrintE82( ioQQQ, dense.gas_phase[nelem]/dense.gas_phase[ipHYDROGEN] );
			}
			fprintf( ioQQQ, " \n" );
		}
		fprintf( ioQQQ, " \n" );
		fprintf(ioQQQ," Gas-phase mass fractions, X:%.3e Y:%.3e Z:%.3e\n\n",
			sumx/SDIV(sumx+sumy+sumz) , 
			sumy/SDIV(sumx+sumy+sumz) , 
			sumz/SDIV(sumx+sumy+sumz) );
	}
	return;
}

/* this is number of elements across one line */
#define	NELEM1LINE	9

/*PrtElem print chemical composition at start of calculation */
STATIC void PrtElem(
  /* the job to do, the options are "init", "fill", "flus" */
  const char *chJob, 
  /* label for the element */
  const char *chLabl, 
  /* its abundance */
  double abund_prt)
{
	static char chAllLabels[NELEM1LINE][14];/* buffer where elements will be stored*/
	long int i, 
	  noffset;
	static long int nelem;  /* counter for number of elements read in*/

	DEBUG_ENTRY( "PrtElem()" );

	if( strcmp(chJob,"initG") == 0 )
	{
		/* gas phase abundances */
		nelem = 0;
		fprintf( ioQQQ, 
			"                                                  Gas Phase Chemical Composition\n" );
	}
	else if( strcmp(chJob,"initD") == 0 )
	{
		/* abundances in grains */
		nelem = 0;
		fprintf( ioQQQ, 
			"                                                    Grain Chemical Composition\n" );
	}

	else if( strcmp(chJob,"fill") == 0 )
	{
		/* print log of abundance to avoid exponential output */
		abund_prt = log10( abund_prt );
		/* stuff in labels and abundances */
		sprintf( chAllLabels[nelem], "  %2.2s:%8.4f", chLabl, abund_prt );
		if( nelem == NELEM1LINE-1 )
		{
			/* we hit as many as it will hold - print it out and reset*/
			fprintf( ioQQQ, "      " );
			for( i=0; i < NELEM1LINE; i++ )
			{
				fprintf( ioQQQ, "%13.13s", chAllLabels[i] );
			}
			fprintf( ioQQQ, "\n" );
			/* reset counter to zero */
			nelem = 0;
		}
		else
		{
			/* just increment */
			++nelem;
		}
	}

#	if 0
	/* Do this if you want to know about PAH number abundance */
	else if( strcmp(chJob,"fillp") == 0 )
	{
		/* print log of abundance to avoid exponential output */
		abund_prt = log10( abund_prt );

		/* stuff in labels and abundances */
		sprintf( chAllLabels[nelem], "  %2.2s:%8.4f", chLabl, abund_prt );
		if( nelem == NELEM1LINE-1 )
		{
			/* we hit as many as it will hold - print it out and reset*/
			fprintf( ioQQQ, "      " );
			for( i=0; i < NELEM1LINE; i++ )
			{
				fprintf( ioQQQ, "%13.13s", chAllLabels[i] );
			}
			fprintf( ioQQQ, "\n" );
			/* reset counter to zero */
			nelem = 0;
		}
		else
		{
			/* just increment */
			++nelem;
		}
	}
#	endif

	else if( strcmp(chJob,"flus") == 0 )
	{
		/* flush the stack */
		i = NELEM1LINE - (nelem - 2);
		noffset = i/2-1;
		/* make format pretty */
		fprintf( ioQQQ, "      " );

		for(i=0; i < noffset; i++)
		{
			/* skip out this many fields */
			fprintf( ioQQQ, "             " );
		}

		/* if nelem is even we need to space out another 8 */
		if( !(nelem%2) && nelem > 0)
			fprintf( ioQQQ,"        ");

		for( i=0; i < nelem; i++ )
		{
			fprintf( ioQQQ, "%13.13s", chAllLabels[i] );
		}

		fprintf( ioQQQ, "\n" );
	}
	else
	{
		fprintf( ioQQQ, " PrtElem does not understand job=%4.4s\n", 
		  chJob );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}


/*AbundancesTable interpolate on table of points to do 'element table' command, */
double AbundancesTable(double r0, 
  double depth, 
  long int iel)
{
	bool lgHit;
	long int j;
	double frac, 
	  tababun_v, 
	  x;

	DEBUG_ENTRY( "AbundancesTable()" );
	/* interpolate on table of points to do 'element table' command, based 
	 * on code by K Volk, each line is log radius and abundance. */

	/* interpolate on radius or depth? */
	if( abund.lgAbTaDepth[iel-1] )
	{
		/* depth key appeared = we want depth */
		x = log10(depth);
	}
	else
	{
		/* use radius */
		x = log10(r0);
	}

	/* this will be reset below, but is here as a safety check */
	tababun_v = -DBL_MAX;

	if( x < abund.AbTabRad[0][iel-1] || x >= abund.AbTabRad[abund.nAbunTabl-1][iel-1] )
	{
		fprintf( ioQQQ, " requested radius outside range of AbundancesTable\n" );
		fprintf( ioQQQ, " radius was%10.2e min, max=%10.2e%10.2e\n", 
		  x, abund.AbTabRad[0][iel-1], abund.AbTabRad[abund.nAbunTabl-1][iel-1] );
		cdEXIT(EXIT_FAILURE);
	}

	else
	{
		lgHit = false;
		j = 1;

		while( !lgHit && j <= abund.nAbunTabl - 1 )
		{
			if( abund.AbTabRad[j-1][iel-1] <= (realnum)x && 
				abund.AbTabRad[j][iel-1] > (realnum)x )
			{
				frac = (x - abund.AbTabRad[j-1][iel-1])/(abund.AbTabRad[j][iel-1] - 
				  abund.AbTabRad[j-1][iel-1]);
				tababun_v = abund.AbTabFac[j-1][iel-1] + frac*
				  (abund.AbTabFac[j][iel-1] - abund.AbTabFac[j-1][iel-1]);
				lgHit = true;
			}
			++j;
		}

		if( !lgHit )
		{
			fprintf( ioQQQ, " radius outran dlaw table scale, requested=%6.2f largest=%6.2f\n", 
			  x, abund.AbTabRad[abund.nAbunTabl-1][iel-1] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* got it, now return value, not log of density */
	tababun_v = pow(10.,tababun_v);
	return( tababun_v );
}

#ifdef _MSC_VER
#	pragma warning( disable : 4305 )/* disable const double to float warning in MS VS - 
									   very large number of warns result */
#endif
/*AbundancesZero set initial abundances for different mixes */
void AbundancesZero(void)
{
	long int i;

	DEBUG_ENTRY( "AbundancesZero()" );

	// option to turn off an element by default, many have no abundances
	// and assume 1e-30 - better to turn the element off and same time
	for( int nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		abund.lgElmONapn[nelem] = true;
		abund.lgElmONahii[nelem] = true;
		abund.lgElmONaism[nelem] = true;
		abund.lgElmONaCrab[nelem] = true;
	}

	// Crab Nebula abundances from constant density model M1 of Table 3
	//>>refer	Crab	abundances	Pequignot, D. & Dennefeld, M. 1983, A&A, 120, 249
	abund.aCrab[ipHYDROGEN] = 1.0;
	abund.aCrab[ipHELIUM] = 7000e-4;
	abund.aCrab[ipLITHIUM] = 2.04e-9;
	abund.lgElmONaCrab[ipLITHIUM] = false;
	abund.aCrab[ipBERYLLIUM] = 2.63e-11;
	abund.lgElmONaCrab[ipBERYLLIUM] = false;
	abund.aCrab[ipBORON] = 7.59e-10;
	abund.lgElmONaCrab[ipBORON] = false;
	abund.aCrab[ipCARBON] = 70e-4;
	abund.aCrab[ipNITROGEN] = 1.23e-4;
	abund.aCrab[ipOXYGEN] = 17e-4;
	abund.aCrab[ipFLUORINE] = 3.02e-8;
	abund.lgElmONaCrab[ipFLUORINE] = false;
	abund.aCrab[ipNEON] = 5.33e-4;
	abund.aCrab[ipSODIUM] = 2.06e-6;
	abund.lgElmONaCrab[ipSODIUM] = false;
	abund.aCrab[ipMAGNESIUM] = 0.80e-4;
	abund.aCrab[ipALUMINIUM] = 2.95e-6;
	abund.lgElmONaCrab[ipALUMINIUM] = false;
	abund.aCrab[ipSILICON] = 3.55e-5;
	abund.lgElmONaCrab[ipSILICON] = false;
	abund.aCrab[ipPHOSPHORUS] = 3.73e-7;
	abund.lgElmONaCrab[ipPHOSPHORUS] = false;
	abund.aCrab[ipSULPHUR] = 0.35e-4;
	abund.aCrab[ipCHLORINE] = 1.88e-7;
	abund.lgElmONaCrab[ipCHLORINE] = false;
	abund.aCrab[ipARGON] = 3.98e-6;
	abund.lgElmONaCrab[ipARGON] = false;
	abund.aCrab[ipPOTASSIUM] = 1.35e-7;
	abund.lgElmONaCrab[ipPOTASSIUM] = false;
	abund.aCrab[ipCALCIUM] = 2.29e-6;
	abund.lgElmONaCrab[ipCALCIUM] = false;
	abund.aCrab[ipSCANDIUM] = 1.58e-9;
	abund.lgElmONaCrab[ipSCANDIUM] = false;
	abund.aCrab[ipTITANIUM] = 1.10e-7;
	abund.lgElmONaCrab[ipTITANIUM] = false;
	abund.aCrab[ipVANADIUM] = 1.05e-8;
	abund.lgElmONaCrab[ipVANADIUM] = false;
	abund.aCrab[ipCHROMIUM] = 4.84e-7;
	abund.lgElmONaCrab[ipCHROMIUM] = false;
	abund.aCrab[ipMANGANESE] = 3.42e-7;
	abund.lgElmONaCrab[ipMANGANESE] = false;
	abund.aCrab[ipIRON] = 3.24e-5;
	abund.aCrab[ipCOBALT] = 8.32e-8;
	abund.lgElmONaCrab[ipCOBALT] = false;
	abund.aCrab[ipNICKEL] = 1.76e-6;
	abund.lgElmONaCrab[ipNICKEL] = false;
	abund.aCrab[ipCOPPER] = 1.87e-8;
	abund.lgElmONaCrab[ipCOPPER] = false;
	abund.aCrab[ipZINC] = 4.52e-8;
	abund.lgElmONaCrab[ipZINC] = false;

	/* solar abundances 
	 * >>refer	solar	abund	Grevesse, N., & Sauval, A.J., 2001, Space Science Review, 85, 161-174 */
	/* >>chng 02 aug 20, update to these values */
	abund.SolarSave[ipHYDROGEN] = 1.0f;
	abund.SolarSave[ipHELIUM] = 0.100f;
	abund.SolarSave[ipLITHIUM] = 2.04e-9f;
	abund.SolarSave[ipBERYLLIUM] = 2.63e-11f;
	abund.SolarSave[ipBORON] = 6.17E-10f;
	/* >>chng 02 jul 30, from 3.55 to 2.45, */
	/* >>refer	C	abund	Allende Prieto, C., 
	 * >>refercon	Lambert, D.L., & Asplund, M., 2002, ApJ, 573, L137 */
	abund.SolarSave[ipCARBON] = 2.45e-4f;
	/* >>chng 02 jul 30, from 9.33 to 8.53, */
	/* >>refer	N	abund	Holweger, H., 2001, Joint SOHO/ACE workshop `Solar and Galactic 
	 * >>refercon	Composition'. Edited by Robert F. Wimmer-Schweingruber. Publisher: 
	 * >>refercon	American Institute of Physics Conference proceedings vol. 598 location: Bern, 
	 * >>refercon	Switzerland, March 6 - 9, 2001., p.23 */
	abund.SolarSave[ipNITROGEN] = 8.51e-5f;
	/* >>chng 02 jul 30, from 7.41 to 4.90, */
	/* >>refer	O	abund	Allende Prieto, C., 
	 * >>refercon	Lambert, D.L., & Asplund, M., 2001, ApJ, 556, L63 */
	abund.SolarSave[ipOXYGEN] = 4.90e-4f;
	abund.SolarSave[ipFLUORINE] = 3.02e-8f;
	/* >>chng 02 jul 30, from 1.17 to 1.00, */
	/* >>refer	Ne	abund	Holweger, H., 2001, Joint SOHO/ACE workshop `Solar and Galactic 
	 * >>refercon	Composition'. Edited by Robert F. Wimmer-Schweingruber. Publisher: 
	 * >>refercon	American Institute of Physics Conference proceedings vol. 598 location: Bern, 
	 * >>refercon	Switzerland, March 6 - 9, 2001., p.23 */
	abund.SolarSave[ipNEON] = 1.00e-4f;
	abund.SolarSave[ipSODIUM] = 2.14e-6f;
	/* >>chng 02 jul 30, from 3.80 to 3.45, */
	/* >>refer	Mg	abund	Holweger, H., 2001, Joint SOHO/ACE workshop `Solar and Galactic 
	 * >>refercon	Composition'. Edited by Robert F. Wimmer-Schweingruber. Publisher: 
	 * >>refercon	American Institute of Physics Conference proceedings vol. 598 location: Bern, 
	 * >>refercon	Switzerland, March 6 - 9, 2001., p.23 */
	abund.SolarSave[ipMAGNESIUM] = 3.47e-5f;
	abund.SolarSave[ipALUMINIUM] = 2.95e-6f;
	/* >>chng 02 jul 30, from 3.55 to 3.44, */
	/* >>refer	Si	abund	Holweger, H., 2001, Joint SOHO/ACE workshop `Solar and Galactic 
	 * >>refercon	Composition'. Edited by Robert F. Wimmer-Schweingruber. Publisher: 
	 * >>refercon	American Institute of Physics Conference proceedings vol. 598 location: Bern, 
	 * >>refercon	Switzerland, March 6 - 9, 2001., p.23 */
	abund.SolarSave[ipSILICON] = 3.47e-5f;
	abund.SolarSave[ipPHOSPHORUS] = 3.20e-7f;
	abund.SolarSave[ipSULPHUR] = 1.84e-5f;
	abund.SolarSave[ipCHLORINE] = 1.91e-7f;
	abund.SolarSave[ipARGON] = 2.51e-6f;
	abund.SolarSave[ipPOTASSIUM] = 1.32e-7f;
	abund.SolarSave[ipCALCIUM] = 2.29e-6f;
	abund.SolarSave[ipSCANDIUM] = 1.48e-9f;
	abund.SolarSave[ipTITANIUM] = 1.05e-7f;
	abund.SolarSave[ipVANADIUM] = 1.00e-8f;
	abund.SolarSave[ipCHROMIUM] = 4.68e-7f;
	abund.SolarSave[ipMANGANESE] = 2.88e-7f;
	/* >>chng 02 jul 30, from 3.24 to 2.81, */
	/* >>refer	Fe	abund	Holweger, H., 2001, Joint SOHO/ACE workshop `Solar and Galactic 
	 * >>refercon	Composition'. Edited by Robert F. Wimmer-Schweingruber. Publisher: 
	 * >>refercon	American Institute of Physics Conference proceedings vol. 598 location: Bern, 
	 * >>refercon	Switzerland, March 6 - 9, 2001., p.23 */
	abund.SolarSave[ipIRON] = 2.82e-5f;
	abund.SolarSave[ipCOBALT] = 8.32e-8f;
	abund.SolarSave[ipNICKEL] = 1.78e-6f;
	abund.SolarSave[ipCOPPER] = 1.62e-8f;
	abund.SolarSave[ipZINC] = 3.98e-8f;

	/* abundance set from pre-c96 */
	/* solar abundances Grevesse and Anders 1989, Grevesse and Noel 1993 */
	abund.OldSolar84[ipHYDROGEN] = 1.0;
	abund.OldSolar84[ipHELIUM] = 0.100;
	abund.OldSolar84[ipLITHIUM] = 2.04e-9;
	abund.OldSolar84[ipBERYLLIUM] = 2.63e-11;
	abund.OldSolar84[ipBORON] = 7.59e-10;
	abund.OldSolar84[ipCARBON] = 3.55e-4;
	abund.OldSolar84[ipNITROGEN] = 9.33e-5;
	abund.OldSolar84[ipOXYGEN] = 7.41e-4;
	abund.OldSolar84[ipFLUORINE] = 3.02e-8;
	abund.OldSolar84[ipNEON] = 1.17e-4;
	abund.OldSolar84[ipSODIUM] = 2.06e-6;
	abund.OldSolar84[ipMAGNESIUM] = 3.80e-5;
	abund.OldSolar84[ipALUMINIUM] = 2.95e-6;
	abund.OldSolar84[ipSILICON] = 3.55e-5;
	abund.OldSolar84[ipPHOSPHORUS] = 3.73e-7;
	abund.OldSolar84[ipSULPHUR] = 1.62e-5;
	abund.OldSolar84[ipCHLORINE] = 1.88e-7;
	abund.OldSolar84[ipARGON] = 3.98e-6;
	abund.OldSolar84[ipPOTASSIUM] = 1.35e-7;
	abund.OldSolar84[ipCALCIUM] = 2.29e-6;
	abund.OldSolar84[ipSCANDIUM] = 1.58e-9;
	abund.OldSolar84[ipTITANIUM] = 1.10e-7;
	abund.OldSolar84[ipVANADIUM] = 1.05e-8;
	abund.OldSolar84[ipCHROMIUM] = 4.84e-7;
	abund.OldSolar84[ipMANGANESE] = 3.42e-7;
	abund.OldSolar84[ipIRON] = 3.24e-5;
	abund.OldSolar84[ipCOBALT] = 8.32e-8;
	abund.OldSolar84[ipNICKEL] = 1.76e-6;
	abund.OldSolar84[ipCOPPER] = 1.87e-8;
	abund.OldSolar84[ipZINC] = 4.52e-8;

	/* Nova Cyg 75 abundances, C, O, NE UP 20, NIT UP 100, REST SOLAR AR */
	abund.anova[ipHYDROGEN] = 1.0;
	abund.anova[ipHELIUM] = 0.098;
	abund.anova[ipLITHIUM] = 2.04e-9;
	abund.anova[ipBERYLLIUM] = 2.6e-11;
	abund.anova[ipBORON] = 7.60e-9;
	abund.anova[ipCARBON] = 9.4e-4;
	abund.anova[ipNITROGEN] = 9.8e-3;
	abund.anova[ipOXYGEN] = 1.7e-2;
	abund.anova[ipFLUORINE] = 3.02e-8;
	abund.anova[ipNEON] = 2.03e-3;
	abund.anova[ipSODIUM] = 2.06e-6;
	abund.anova[ipMAGNESIUM] = 3.80e-5;
	abund.anova[ipALUMINIUM] = 2.95e-6;
	abund.anova[ipSILICON] = 3.55e-5;
	abund.anova[ipPHOSPHORUS] = 3.73e-7;
	abund.anova[ipSULPHUR] = 1.62e-5;
	abund.anova[ipCHLORINE] = 1.88e-7;
	abund.anova[ipARGON] = 3.63e-6;
	abund.anova[ipPOTASSIUM] = 1.35e-7;
	abund.anova[ipCALCIUM] = 2.29e-6;
	abund.anova[ipSCANDIUM] = 1.22e-9;
	abund.anova[ipTITANIUM] = 8.60e-8;
	abund.anova[ipVANADIUM] = 1.05e-8;
	abund.anova[ipCHROMIUM] = 4.84e-7;
	abund.anova[ipMANGANESE] = 3.42e-7;
	abund.anova[ipIRON] = 4.68e-5;
	abund.anova[ipCOBALT] = 2.24e-9;
	abund.anova[ipNICKEL] = 1.76e-6;
	abund.anova[ipCOPPER] = 1.87e-8;
	abund.anova[ipZINC] = 4.52e-8;

	/* primordial abundances */
	abund.aprim[ipHYDROGEN] = 1.0;
	abund.aprim[ipHELIUM] = 0.072;
	abund.aprim[ipLITHIUM] = 1e-10;
	abund.aprim[ipBERYLLIUM] = 1e-16;

	for( i=4; i < LIMELM; i++ )
	{
		abund.aprim[i] = 1e-25;
	}

	/* typical ISM abundances, mean of Table 3, Cowie+Songaila, Ann Rev '86
	 * also Table 5, Savage and Sembach, Ann Rev 1996 */
	abund.aism[ipHYDROGEN] = 1.;
	abund.aism[ipHELIUM] = 0.098;
	abund.aism[ipLITHIUM] = 5.4e-11;
	abund.aism[ipBERYLLIUM] = 1e-20;
	abund.lgElmONaism[ipBERYLLIUM] = false;
	abund.aism[ipBORON] = 8.9e-11;
	abund.aism[ipCARBON] = 2.51e-4;
	abund.aism[ipNITROGEN] = 7.94e-5;
	/* >>chng >>01 feb 19, from 5.01e-4 to 3.19e-4, value from */
	/* >>refer	O	abundance	Meyers, D.M., Jura, M., & Cardelli, J.A., 1998, ApJ, 493, 222-229 */
	/* they quote 3.19 +/- 0.14 e-4 */
	abund.aism[ipOXYGEN] = 3.19e-4;
	/* >>chng 10 jul 22 -- NPA.  F abundance slightly depleted (0.1 to -0.6 dex) relative to the solar
	 * F/H ratio of 3.6e-8 (reference below).  Use value of 3e-8 */
	/* >>refer	F	abundance	Snow, T. P., Destree, J. D., & Jensen, A. G. 2007, ApJ, 655, 285 */
	abund.aism[ipFLUORINE] = 2.0e-8;
	abund.aism[ipNEON] = 1.23e-4;
	abund.aism[ipSODIUM] = 3.16e-7;
	abund.aism[ipMAGNESIUM] = 1.26e-5;
	abund.aism[ipALUMINIUM] = 7.94e-8;
	abund.aism[ipSILICON] = 3.16e-6;
	abund.aism[ipPHOSPHORUS] = 1.6e-7;
	abund.aism[ipSULPHUR] = 3.24e-5;
	abund.aism[ipCHLORINE] = 1e-7;
	abund.aism[ipARGON] = 2.82e-6;
	abund.aism[ipPOTASSIUM] = 1.1e-8;
	abund.aism[ipCALCIUM] = 4.1e-10;
	abund.aism[ipSCANDIUM] = 1e-20;
	abund.lgElmONaism[ipSCANDIUM] = false;
	abund.aism[ipTITANIUM] = 5.8e-10;
	abund.aism[ipVANADIUM] = 1.0e-10;
	abund.aism[ipCHROMIUM] = 1.0e-8;
	abund.aism[ipMANGANESE] = 2.3e-8;
	abund.aism[ipIRON] = 6.31e-7;
	/* >>chng 10 jul 22 -- NPA.  Co/H ratio in cold and warm phases towards rho Oph are derived by reference
	 * below to be 6.4e-10 and 1.1e-8, respectively.  Average of two comes out to 5.8e-9 */
	/* >>refer	Co	abundance	Mullman, K. L. et al. 1998, ApJ, 500, 1064 */
	abund.aism[ipCOBALT] = 5.9e-9;
	abund.aism[ipNICKEL] = 1.82e-8;
	abund.aism[ipCOPPER] = 1.5e-9;
	abund.aism[ipZINC] = 2.0e-8;

	/* HII region abundances, Orion mean of Baldwin et al, Rubin et al,
	 * and DEO et al, all 1991 apj
	 * also Table 5, Savage and Sembach, Ann Rev 1996 for ism */
	abund.ahii[ipHYDROGEN] = 1.;
	abund.ahii[ipHELIUM] = 0.095;
	abund.ahii[ipLITHIUM] = 5.4e-11;
	abund.ahii[ipBERYLLIUM] = 1e-20;
	abund.lgElmONahii[ipBERYLLIUM] = false;
	abund.ahii[ipBORON] = 8.9e-11;
	abund.ahii[ipCARBON] = 3.e-4;
	abund.ahii[ipNITROGEN] = 7.0e-5;
	abund.ahii[ipOXYGEN] = 4.0e-4;
	abund.ahii[ipFLUORINE] = 1e-20;
	abund.lgElmONahii[ipFLUORINE] = false;
	abund.ahii[ipNEON] = 6e-5;
	abund.ahii[ipSODIUM] = 3e-7;
	abund.ahii[ipMAGNESIUM] = 3.e-6;
	abund.ahii[ipALUMINIUM] = 2.e-7;
	abund.ahii[ipSILICON] = 4.e-6;
	abund.ahii[ipPHOSPHORUS] = 1.6e-7;
	abund.ahii[ipSULPHUR] = 1.0e-5;
	abund.ahii[ipCHLORINE] = 1.e-7;
	abund.ahii[ipARGON] = 3.e-6;
	abund.ahii[ipPOTASSIUM] = 1.1e-8;
	abund.ahii[ipCALCIUM] = 2.e-8;
	abund.ahii[ipSCANDIUM] = 1e-20;
	abund.lgElmONahii[ipSCANDIUM] = false;
	abund.ahii[ipTITANIUM] = 5.8e-10;
	abund.ahii[ipVANADIUM] = 1.0e-10;
	abund.ahii[ipCHROMIUM] = 1.0e-8;
	abund.ahii[ipMANGANESE] = 2.3e-8;
	abund.ahii[ipIRON] = 3.0e-6;
	abund.ahii[ipCOBALT] = 1e-20;
	abund.lgElmONahii[ipCOBALT] = false;
	abund.ahii[ipNICKEL] = 1e-7;
	abund.ahii[ipCOPPER] = 1.5e-9;
	abund.ahii[ipZINC] = 2.0e-8;

	/* PN abund from  */
	/* >>refer	PN	abundances	Aller+Czyzak, ApJ Sup 51, 211 */
	abund.apn[ipHYDROGEN] = 1.;
	abund.apn[ipHELIUM] = 0.1;
	abund.apn[ipLITHIUM] = 1e-20;
	abund.lgElmONapn[ipLITHIUM] = false;
	abund.apn[ipBERYLLIUM] = 1e-20;
	abund.lgElmONapn[ipBERYLLIUM] = false;
	abund.apn[ipBORON] = 1e-20;
	abund.lgElmONapn[ipBORON] = false;
	abund.apn[ipCARBON] = 7.8e-4;
	abund.apn[ipNITROGEN] = 1.8e-4;
	abund.apn[ipOXYGEN] = 4.4e-4;
	abund.apn[ipFLUORINE] = 3e-7;
	abund.apn[ipNEON] = 1.1e-4;
	abund.apn[ipSODIUM] = 1.9e-6;
	abund.apn[ipMAGNESIUM] = 1.6e-6;
	abund.apn[ipALUMINIUM] = 2.7e-7;
	abund.apn[ipSILICON] = 1e-5;
	abund.apn[ipPHOSPHORUS] = 2e-7;
	abund.apn[ipSULPHUR] = 1e-5;
	abund.apn[ipCHLORINE] = 1.7e-7;
	abund.apn[ipARGON] = 2.7e-6;
	abund.apn[ipPOTASSIUM] = 1.2e-7;
	abund.apn[ipCALCIUM] = 1.2e-8;
	abund.apn[ipSCANDIUM] = 1e-20;
	abund.lgElmONapn[ipSCANDIUM] = false;
	abund.apn[ipTITANIUM] = 1e-20;
	abund.lgElmONapn[ipTITANIUM] = false;
	abund.apn[ipVANADIUM] = 1e-20;
	abund.lgElmONapn[ipVANADIUM] = false;
	abund.apn[ipCHROMIUM] = 1e-20;
	abund.lgElmONapn[ipCHROMIUM] = false;
	abund.apn[ipMANGANESE] = 1e-20;
	abund.lgElmONapn[ipMANGANESE] = false;
	abund.apn[ipIRON] = 5.0e-7;
	abund.apn[ipCOBALT] = 1e-20;
	abund.lgElmONapn[ipCOBALT] = false;
	abund.apn[ipNICKEL] = 1.8e-8;
	abund.apn[ipCOPPER] = 1e-20;
	abund.lgElmONapn[ipCOPPER] = false;
	abund.apn[ipZINC] = 1e-20;
	abund.lgElmONapn[ipZINC] = false;

	/* mix from Cameron 1982, in "Essays on Nuclear Astro" */
	abund.camern[ipHYDROGEN] = 1.;
	abund.camern[ipHELIUM] = .0677;
	abund.camern[ipLITHIUM] = 2.2e-9;
	abund.camern[ipBERYLLIUM] = 4.5e-11;
	abund.camern[ipBORON] = 3.4e-10;
	abund.camern[ipCARBON] = 4.22e-4;
	abund.camern[ipNITROGEN] = 8.72e-5;
	abund.camern[ipOXYGEN] = 6.93e-4;
	abund.camern[ipFLUORINE] = 2.9e-8;
	abund.camern[ipNEON] = 9.77e-5;
	abund.camern[ipSODIUM] = 2.25e-6;
	abund.camern[ipMAGNESIUM] = 3.98e-5;
	abund.camern[ipALUMINIUM] = 3.20e-6;
	abund.camern[ipSILICON] = 3.76e-5;
	abund.camern[ipPHOSPHORUS] = 2.4e-7;
	abund.camern[ipSULPHUR] = 1.88e-5;
	abund.camern[ipCHLORINE] = 1.78e-7;
	abund.camern[ipARGON] = 3.99e-6;
	abund.camern[ipPOTASSIUM] = 1.3e-7;
	abund.camern[ipCALCIUM] = 2.35e-6;
	abund.camern[ipSCANDIUM] = 1.16e-9;
	abund.camern[ipTITANIUM] = 9.0e-8;
	abund.camern[ipVANADIUM] = 9.5e-9;
	abund.camern[ipCHROMIUM] = 4.8e-7;
	abund.camern[ipMANGANESE] = 3.5e-7;
	abund.camern[ipIRON] = 3.38e-5;
	abund.camern[ipCOBALT] = 8.27e-8;
	abund.camern[ipNICKEL] = 1.80e-6;
	abund.camern[ipCOPPER] = 2.0e-8;
	abund.camern[ipZINC] = 4.7e-8;

	/* set logical flags saying whether to include element in AGN tables */
	/* first set all false, since most not included */
	for( i=0; i < LIMELM; i++ )
	{
		abund.lgAGN[i] = false;
	}
	abund.lgAGN[ipHYDROGEN] = true;
	abund.lgAGN[ipHELIUM] = true;
	abund.lgAGN[ipCARBON] = true;
	abund.lgAGN[ipNITROGEN] = true;
	abund.lgAGN[ipOXYGEN] = true;
	abund.lgAGN[ipNEON] = true;
	abund.lgAGN[ipMAGNESIUM] = true;
	abund.lgAGN[ipSILICON] = true;
	abund.lgAGN[ipSULPHUR] = true;
	abund.lgAGN[ipARGON] = true;
	abund.lgAGN[ipIRON] = true;
	return;
}

