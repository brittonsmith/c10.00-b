/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtContinuum print information about continuum if requested with PRINT CONTINUUM command,
 * called by PrtFinal */
#include "cddefines.h"
#include "rfield.h"
#include "iso.h"
#include "radius.h"
#include "opacity.h"
#include "prt.h"

static const int NFLUXSV = 360;
static const int NXBD = 9;

void PrtContinuum(void)
{
	long int i, 
	  i1, 
	  j, 
	  ninc, 
	  nline;
	realnum fluxsv[NFLUXSV], 
	  xbdsav[NXBD];

	/* energies for the x-ray bands */
	static double xraybd[NXBD + 1]={
		7.3498,
		36.8,
		73.5,
		110.3,
		147.1,
		183.8,
		220.6,
		367.6,
		551.5,
		735.3};

	DEBUG_ENTRY( "PrtContinuum()" );
	/*print information about continuum if requested with PRINT CONTINUUM command */
	/* this is number of ranges for adding x-ray bands */

	/* this stuff was always printed before version 84, and is now an option
	 * to turn on information with PRINT CONTINUUM command */

	/* the default - not turned on, just return */
	if( !prt.lgPrtCont )
	{ 
		return;
	}

	/* derive x-ray fluxes in various bands for later printout
	 *
	 * >>chng 97 jun 08, get rid of statement labels */
	if( xraybd[0] < rfield.anu[rfield.nflux-1] )
	{
		i1 = opac.ipCKshell - 10;
		/* get out to lower bound of first range */
		while( (double)rfield.anu[i1-1] < xraybd[0] && i1 < rfield.nflux )
		{
			i1 += 1;
		}
		/* now add up over that range */
		for( i=0; i < NXBD; i++ )
		{
			xbdsav[i] = 0.;
			while( (double)rfield.anu[i1-1] < xraybd[i+1] && i1 < rfield.nflux )
			{
				xbdsav[i] += rfield.flux[0][i1-1] + 
				  rfield.outlin[0][i1-1] +rfield.outlin_noplot[i1-1] + rfield.ConInterOut[i1-1];
				i1 += 1;
			}
			xbdsav[i] *= (realnum)radius.r1r0sq;
		}
	}
	else
	{
		/* continuum not defined at x-ray energies, but zero it for safety */
		for( i=0; i < NXBD; i++ )
		{
			xbdsav[i] = 0.;
		}
	}

	/* now print x-ray flux (photons s^-1, flux or lum) in various bands */
	if( xbdsav[0] > 0 )
	{
		fprintf( ioQQQ, "\n 0.1-0.5kev:" );
		for(i=0; i < NXBD; i++)
			fprintf( ioQQQ, "%8.2e", xbdsav[i] );
		fprintf( ioQQQ, " 0.5-1.0kev:\n" );
	}

	/* renorm fLUX array before printing it */
	if( iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s] - iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2] + 1 > NFLUXSV )
	{
		fprintf( ioQQQ, " PCONTN: not enough cells in flux_total_incident, need%4ld\n", 
		  iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s] - iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2] + 1 );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]-1; i < iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]; i++ )
	{
		if( rfield.flux_total_incident[0][i] > 0. )
		{
			fluxsv[i-iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]+1] = (realnum)((rfield.flux[0][i] +rfield.outlin[0][i] +  rfield.outlin_noplot[i] + 
			  rfield.ConInterOut[i] )*radius.r1r0sq/
			  rfield.flux_total_incident[0][i]);
		}
		else
		{
			fluxsv[i-iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]+1] = 0.;
		}
	}
	/* renormalize whole continuum */
	for( i=0; i < rfield.nflux; i++ )
	{
		rfield.flux[0][i] = (realnum)(((rfield.flux[0][i] + 
		  rfield.ConInterOut[i]/rfield.anu[i])/rfield.widflx[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i])*
		  radius.r1r0sq);
	}

	fprintf( ioQQQ, 
		"\n\n                                                        Normalised continuum\n" );

	for( i=iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]; i <= iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]; i += 3 )
	{
		fprintf( ioQQQ, "%7.3f%6.3f", rfield.anu[i-1], fluxsv[i-iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]] );
	}
	fprintf( ioQQQ, "\n" );

	fprintf( ioQQQ, 
		"\n                                                  Emergent continuum - phot/ryd/cm2/s (r in)\n" );

	nline = ((rfield.nflux - iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2] - 1)/4 + 7)/7;
	ninc = nline*4;
	for( j=0; j < nline; j++ )
	{
		i1 = j*4 + iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2];
		fprintf( ioQQQ, " " );

		for( i=i1; i<rfield.nflux; i = i + ninc)
		{
			fprintf( ioQQQ, "%6.3f%10.2e", rfield.anu[i], 
			  rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i] +rfield.ConInterOut[i] );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
