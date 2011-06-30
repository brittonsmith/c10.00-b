/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityAdd1Element enter total photo cross section for all subshells into opacity array */
#include "cddefines.h"
#include "iso.h"
#include "rfield.h"
#include "dense.h"
#include "heavy.h"
#include "opacity.h"
#include "taulines.h"

void OpacityAdd1Element(
		/* nelem is 0 for H, 1 for He, etc */
		long int nelem)
{
	long int ipHi, 
	  ipop, 
	  limit,
	  low, 
	  n, 
	  ion, 
	  nshell;
	char chStat;
	double abundance;

	DEBUG_ENTRY( "OpacityAdd1Element()" );

	/* this routine drives OpacityAdd1Subshell to put in total opacities for all shells*/

	/*begin sanity check */
	ASSERT( (nelem >=0 ) && (nelem < LIMELM) );

	/* first do simple two-level systems -
	 * this is number of species that are not treated on common iso-electronic series */
	limit = nelem + 1 - NISO;
	/* this can be called with hydrogen itself, in which case nelem is 0, and limit is
	 * -1 - do not do any of the simple ions */
	limit = MAX2( 0 , limit );

	/* do not include the ion stages that have complete atoms,
	 * currently H and He like iso sequences */
	for( ion=0; ion < limit; ion++ )
	{
		if( dense.xIonDense[nelem][ion] > 0. )
		{
			/*start with static opacities, then do volatile*/

			chStat = 's';
			/* number of bound electrons */
			for( nshell=0; nshell < Heavy.nsShells[nelem][ion]; nshell++ )
			{
				/* highest shell will be volatile*/
				if( nshell== Heavy.nsShells[nelem][ion]-1 )
					chStat = 'v';
				/* set lower and upper limits to this range */
				low = opac.ipElement[nelem][ion][nshell][0];
				ipHi = opac.ipElement[nelem][ion][nshell][1];
				ipop = opac.ipElement[nelem][ion][nshell][2];
				/* OpacityAdd1Subshell will not do anything if static opacities do not need to be reset*/
				OpacityAdd1Subshell(ipop,low,ipHi,dense.xIonDense[nelem][ion] , chStat );
			}
		}
	}

	/* now loop over all species done as large multi-level systems */
	/* >>chng 02 jan 17, add loop over H and He like */
	/* ion is on the c scale, =0 for HI, =1 for HeII */
	for( ion=limit; ion<nelem+1; ++ion )
	{
		/* ipISO is 0 for H-like, 1 for He-like */
		long int ipISO = nelem-ion;

		/* do multi level systems, but only if present 
		 * test for nelem+1 in case atom present but not ion, test is whether the
		 * abundance of the recombined species is present */
		/* >>chng 02 jan 17, sec dim had been nelem+1, change to ion+1 */
		/*if( dense.xIonDense[nelem][nelem] > 0. )*/
		if( dense.xIonDense[nelem][ion] > 0. )
		{
			ASSERT( ipISO < NISO );

			/* do ground first, then all excited states */
			n = 0;
			/* abundance of recombined species, which can be zero if no ion present */
			abundance = StatesElemNEW[nelem][nelem-ipISO][n].Pop;

			/* >>chng 02 may 06, add second test, had been just the chck on helium,
			 * with no option to use new soln */
			if( abundance == 0.  )
			{
				/* no ionized species, assume everything in ground */
				abundance = dense.xIonDense[nelem][ion];
			}

			/* >>chng 02 jan 17, to arbitrary iso sequence */
			/* use computed opacities and departure coef for level */
			OpacityAdd1SubshellInduc(
				iso.ipOpac[ipISO][nelem][n],
				iso.ipIsoLevNIonCon[ipISO][nelem][n],
				/* the upper limit to the integration, 
				* ground opacity goes up to the high-energy limit of code*/
				rfield.nflux,
				/* the abundance of the ion */
				abundance,
				/* departure coef, volatile opac, always reevaluate */
				iso.DepartCoef[ipISO][nelem][n] , 'v' );

			/* do excited levvels,
			 * this loop only if upper levels have finite population*/
			if( StatesElemNEW[nelem][nelem-ipISO][3].Pop > 0. )
			{
				char chType = 'v';
				/* always want to evaluate all opacities for n=3, 4, use static opacities for higher levels */
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max */
				for( long level =1; level < iso.numLevels_local[ipISO][nelem]; level++ )
				{
					/* above 4 is static */
					if( StatesElemNEW[nelem][nelem-ipISO][level].n >= 5 )
						chType = 's';

					/* include correction for stimulated emission */
					OpacityAdd1SubshellInduc(
						iso.ipOpac[ipISO][nelem][level],
						iso.ipIsoLevNIonCon[ipISO][nelem][level],
						/* the high energy bound of excited states is the 
						* edge of the Lyman continuum */
						iso.ipIsoLevNIonCon[ipISO][nelem][0],
						StatesElemNEW[nelem][nelem-ipISO][level].Pop,
						/* departure coef, volitile opacities */
						iso.DepartCoef[ipISO][nelem][level] , chType );
				}
			}
		}
	}
	return;
}
