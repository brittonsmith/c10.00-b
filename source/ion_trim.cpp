/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_trim raise or lower most extreme stages of ionization considered,
 * called by ConvBase - ion limits were originally set by  */
#include "cddefines.h"
#include "elementnames.h"
#include "radius.h"
#include "heavy.h"
#include "conv.h"
#include "rfield.h"
#include "phycon.h"
#include "mole.h"
#include "thermal.h"
#include "iso.h"
#include "dense.h"
#include "conv.h"
#include "struc.h"
#include "ionbal.h"
#include "taulines.h"

void ion_trim(
	/* nelem is on the C scale, 0 for H, 5 for C */
	long int nelem )
{

	/* this will remember that higher stages trimed up or down */
	bool lgHi_Down = false;
	bool lgHi_Up = false;
	bool lgHi_Up_enable;
	/* this will remember that lower stages trimmed up or own*/
	bool lgLo_Up = false;
	bool lgLo_Down = false;
	long int ion_lo_save = dense.IonLow[nelem],
		ion_hi_save = dense.IonHigh[nelem];
	long int ion;
	realnum trimhi , trimlo;

	/* this is debugging code that turns on print after certain number of calls */
	/*realnum xlimit = SMALLFLOAT *10.;*/
	/*static int ncall=0;
	if( nelem==5 )
		++ncall;*/

	DEBUG_ENTRY( "ion_trim()" );

	/* this routine should not be called early in the simulation, before
	 * things have settled down */
	ASSERT( conv.nTotalIoniz>2 );

	/*confirm all ionization stages are within their range of validity */
	ASSERT( nelem >= ipHELIUM && nelem < LIMELM );
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonHigh[nelem] <= nelem+1 );
	/* IonHigh can be equal to IonLow if both are zero - so no ionization,
	 * or is ionization has been set with "element ionization" command */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] ||
		(dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 )  || 
		dense.lgSetIoniz[nelem] );

	/* during search phase of mostly neutral matter the electron density
	 * can be vastly too large, and the ionization suppressed as a result.
	 * during search do not trim down or up as much */
	if( conv.lgSearch )
	{
		trimhi = (realnum)ionbal.trimhi * 1e-4f;
		trimlo = (realnum)ionbal.trimlo * 1e-4f;
	}
	else
	{
		trimhi = (realnum)ionbal.trimhi;
		trimlo = (realnum)ionbal.trimlo;
	}

	/* helium is special case since abundance is so high, and He+ CT with 
	 * CO is the dominant CO destruction process in molecular regions */
	if( nelem == ipHELIUM )
	{
		/* never want to trip up a lower stage of ionization */
		trimlo = SMALLFLOAT;

		/* if He+ is highest stage of ionization, probably want to keep it
		 * since important for CO chemistry in molecular clouds */
		if( dense.IonHigh[ipHELIUM] == 1 )
		{
			trimhi = MIN2( trimhi , 1e-20f );
		}
		else if( dense.IonHigh[ipHELIUM] == 2 )
		{
			if( conv.lgSearch )
			{
				/* during search phase we may be quite far from solution, in the
				 * case of a PDR sim the electron density may be vastly higher 
				 * than it will be with stable solution.  He++ can be very important
				 * for the chemistry and we want to consider it.  Make the
				 * threshold for ignoring He++ much higher than normal to prevent
				 * a premature removal of the ion */
				trimhi = MIN2( trimhi , 1e-17f ); 
			}
			else
			{
				/* similar smaller upper limit for ion*/
				trimhi = MIN2( trimhi , 1e-12f ); 
			}
		}
	}

	/* logic for PDRs, for elements included in chemistry, need stable solutions, 
	 * keep 3 ion stages in most cases - added by NA to do HII/PDR sims */
	if( mole.lgElem_in_chemistry[nelem] )
	{
		trimlo = SMALLFLOAT;
		if(dense.IonHigh[nelem] ==2)
		{
			trimhi = MIN2(trimhi, 1e-20f);
		}
	}

	/* raise or lower highest and lowest stages of ionization to
	 * consider only significant stages
	 * IonHigh[nelem]  stage of ionization  */

	/* this is a special block for initialization only - it checks on absurd abundances
	 * and can trim multiple stages of ionization at one time. */
	if( conv.lgSearch )
	{
		/* special - trim down higher stages if they have essentially zero abundance */
		while( 
			(dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] < dense.density_low_limit || 
			dense.xIonDense[nelem][dense.IonHigh[nelem]] < dense.density_low_limit ) && 
			/* >>chng 02 may 12, rm +1 since this had effect of not allowing fully atomic */
			dense.IonHigh[nelem] > dense.IonLow[nelem] )
		{
			/* dense.xIonDense[nelem][i] is density of i+1 th ionization stage (cm^-3)
			 * the -1 is correct for the heating, -1 since heating is caused by ionization of stage below it */
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
			thermal.heating[nelem][dense.IonHigh[nelem]-1] = 0.;
			if( dense.IonHigh[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonHigh[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso.numLevels_max[ipISO][nelem]; level++ )
					StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
			}

			/* decrement high stage limit */
			--dense.IonHigh[nelem];
			ASSERT( dense.IonHigh[nelem] >= 0);
			/* remember this happened */
			lgHi_Down = true;
		}

		/* special - trim up lower stages trim if they have essentially zero abundance */
		while( 
			(dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] < dense.density_low_limit || 
			dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit ) && 
			dense.IonLow[nelem] < dense.IonHigh[nelem] - 1 )
		{
			/* dense.xIonDense[nelem][i] is density of ith ionization stage (cm^-3)
			 * there is no-1 since we are removing the agent that heats */
			dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;
			/* >>chng 01 aug 04, remove -1 which clobbers thermal.heating when IonLow == 0 */
			/*thermal.heating[nelem][dense.IonLow[nelem]-1] = 0.;*/
			thermal.heating[nelem][dense.IonLow[nelem]] = 0.;
			if( dense.IonLow[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonLow[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso.numLevels_max[ipISO][nelem]; level++ )
					StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
			}

			/* increment low stage limit */
			++dense.IonLow[nelem];
			lgLo_Up = true;
		}
	}

	/* sanity check */
	/* IonHigh can be equal to IonLow if both are zero - no ionization*/
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] ||
		(dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 )  || 
		dense.lgSetIoniz[nelem] );

	/* this checks on error condition where the gas is stupendously highly ionized,
	 * the low limit is already one less than high, and we need to raise the low
	 * limit further */
	if( dense.IonHigh[nelem] > 1 &&
		(dense.IonLow[nelem]==dense.IonHigh[nelem]-1) &&
		(dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit) )
	{
#		if 0
		fprintf(ioQQQ,
			"PROBLEM DISASTER the density of ion %li of element %s is too small to be computed on this cpu.\n",
			dense.IonLow[nelem]+1,
			elementnames.chElementName[nelem]);
		fprintf(ioQQQ,
			"Turn off the element with the command ELEMENT XXX OFF or do not consider "
			"gas with low density, the current hydrogen density is %.2e cm-3.\n",
			dense.gas_phase[ipHYDROGEN]);
		fprintf(ioQQQ,
			"The outer radius of the cloud is %.2e cm - should a stopping "
			"radius be set?\n",
			radius.Radius );
		fprintf(ioQQQ,
			"The abort flag is being set - the calculation is stopping.\n");
#		endif
		//lgAbort = true;
		dense.xIonDense[nelem][dense.IonLow[nelem]] = (realnum)dense.density_low_limit;
		return;
	}

	/* trim down high stages that have too small an abundance */
	if( 
		dense.IonHigh[nelem] > dense.IonLow[nelem]  && 
		( (dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] <= 
		trimhi ) ||
		(dense.xIonDense[nelem][dense.IonHigh[nelem]] <= dense.density_low_limit ) )
		) 
	{
		/* >>chng 03 sep 30, the atom and its first ion are a special case
		 * since we want to compute even trivial ions in molecular clouds */
		if( dense.IonHigh[nelem]>1 ||
			(dense.IonHigh[nelem]==1&&dense.xIonDense[nelem][1]<100.*dense.density_low_limit) )
		{
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
			thermal.heating[nelem][dense.IonHigh[nelem]-1] = 0.;
			if( dense.IonHigh[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonHigh[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso.numLevels_max[ipISO][nelem]; level++ )
					StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
			}
			--dense.IonHigh[nelem];
			lgHi_Down = true;
		}
	}

	/* trim up highest stages - will this be done? */
	lgHi_Up_enable = true;
	/* trimming can be turned off */
	if( !ionbal.lgTrimhiOn )
		lgHi_Up_enable = false;
	/* high stage is already fully stripped */
	if( dense.IonHigh[nelem]>=nelem+1 )
		lgHi_Up_enable = false;
	/* we have previously trimmed this stage down */
	if( lgHi_Down )
		lgHi_Up_enable = false;
	/* we are in neutral gas */
	if( dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN]<0.9  )
		lgHi_Up_enable = false;

	if( lgHi_Up_enable )
	{
		realnum abundold = 0;
		if( nzone>2 ) 
			abundold = struc.xIonDense[nelem][dense.IonHigh[nelem]][nzone-3]/
			SDIV( struc.gas_phase[nelem][nzone-3]);
		realnum abundnew = dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem];
		/* only raise highest stage if ionization potential of next highest stage is within
		 * continuum array and the abundance of the highest stage is significant */
		if( (Heavy.Valence_IP_Ryd[nelem][dense.IonHigh[nelem]] < rfield.anu[rfield.nflux-1] ||
			Heavy.Valence_IP_Ryd[nelem][dense.IonHigh[nelem]] < phycon.te_ryd*10.) &&
			dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] > 1e-4f &&
			/* this checks that abundance of highest stage is increasing */
			abundnew > abundold*1.01 )
		{
			/*fprintf(ioQQQ,"uuppp %li %li \n", nelem, dense.IonHigh[nelem] );*/
			/* raise highest level of ionization */
			++dense.IonHigh[nelem];
			lgHi_Up = true;
			/* set abundance to almost zero so that sanity check elsewhere will be ok */
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = (realnum)(dense.density_low_limit*10.);
		}
	}

	/* sanity check 
	 * IonHigh can be equal to IonLow if both are zero - no ionization*/
	ASSERT( dense.IonLow[nelem] < dense.IonHigh[nelem] || 
		(dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 )  || 
		dense.lgSetIoniz[nelem] );

	/* lower lowest stage of ionization if we have significant abundance at current lowest */
	if( dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] > 1e-3f && 
		dense.IonLow[nelem] > 0 )
	{
		/* lower lowest level of ionization */
		--dense.IonLow[nelem];
		lgLo_Down = true;
		/* >>chng 01 aug 02, set this to zero so that sanity check elsewhere will be ok */
		dense.xIonDense[nelem][dense.IonLow[nelem]] = (realnum)dense.density_low_limit;
	}

	/* raise lowest stage of ionization, but only if we are near illuminated face of cloud*/
	else if( nzone < 10 &&
		(dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] <= (realnum)trimlo) && 
		(dense.IonLow[nelem] < (dense.IonHigh[nelem] - 2) ) )
	{
		/* raise lowest level of ionization */
		dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;
		/* no minus one in below since this is low bound, already bounds at atom */
		thermal.heating[nelem][dense.IonLow[nelem]] = 0.;
		if( dense.IonLow[nelem] == nelem+1-NISO )
		{
			long int ipISO = nelem - dense.IonLow[nelem];
			ASSERT( ipISO>=0 && ipISO<NISO );
			for( long level = 0; level < iso.numLevels_max[ipISO][nelem]; level++ )
				StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
		}
		++dense.IonLow[nelem];
		lgLo_Up = true;
	}
	/* test on zero */
	else if( (dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit) && 
		/*>>chng 06 may 24, from IonLow < IonHigh to <IonHigh-1 because
		 * old test would allow low to be raised too high, which then
		 * leads to insanity */
		(dense.IonLow[nelem] < dense.IonHigh[nelem]-1) )
	{
		while(dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit && 
			/* >>chng 07 feb 27 from < IonHigh to < IonHigh-1 to prevent raising
			 * low to high */
			(dense.IonLow[nelem] < dense.IonHigh[nelem]-1) )
		{
			/* raise lowest level of ionization */
			dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;
			/* no minus one in below since this is low bound, already bounds at atom */
			thermal.heating[nelem][dense.IonLow[nelem]] = 0.;
			if( dense.IonLow[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonLow[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso.numLevels_max[ipISO][nelem]; level++ )
					StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
			}
			++dense.IonLow[nelem];
			lgLo_Up = true;
		}
	}

	/* sanity check */
	/* IonHigh can be equal to IonLow if both are zero - no ionization*/
	ASSERT( dense.IonLow[nelem] < dense.IonHigh[nelem] || 
		(dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 )  || 
		dense.lgSetIoniz[nelem] );

	/* these are standard bounds checks that appear throughout this routine
	 * dense.xIonDense[IonLow] is first one with positive abundances
	 * 
	 * in case where lower ionization stage was just lowered the abundance
	 * was set to dense.density_low_limit so test must be < dense.density_low_limit */
	ASSERT( dense.IonLow[nelem] <= 1 ||
		dense.xIonDense[nelem][dense.IonLow[nelem]-1] == 0. );

	ASSERT( (dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 ) || lgLo_Up ||
		dense.xIonDense[nelem][dense.IonLow[nelem]] >= dense.density_low_limit ||
		dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] >= dense.density_low_limit ||
		/*>>chng 06 may 24, include this to cover case where cap is set by not allowing
		 * lower limit to come up to upper limit */
		(dense.IonLow[nelem]==dense.IonHigh[nelem]-1 ));

	{
		/* option to print out what has happened so far .... */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nelem == ipHELIUM )
		{
			if( lgHi_Down ||lgHi_Up ||lgLo_Up ||lgLo_Down )
			{
				fprintf(ioQQQ,"DEBUG TrimZone\t%li\t",nzone );
				if(  lgHi_Down )
				{
					fprintf(ioQQQ,"high dn %li to %li",
						ion_hi_save , 
						dense.IonHigh[nelem] );
				}
				if( lgHi_Up )
				{
					fprintf(ioQQQ,"high up %li to %li",
						ion_hi_save , 
						dense.IonHigh[nelem] );
				}
				if( lgLo_Up )
				{
					fprintf(ioQQQ,"low up %li to %li",
						ion_lo_save , 
						dense.IonLow[nelem] );
				}
				if( lgLo_Down )
				{
					fprintf(ioQQQ,"low dn %li to %li",
						ion_lo_save , 
						dense.IonLow[nelem] );
				}
				fprintf(ioQQQ,"\n" );
			}
		}
	}

	/* only assert that the stage above the highest has a population of
	 * zero if that stage exists */
	if(dense.IonHigh[nelem] < nelem+1 )
		ASSERT( dense.xIonDense[nelem][dense.IonHigh[nelem]+1] == 0. );

	/* option to print if any trimming occurred */
	if( lgHi_Down || lgHi_Up || lgLo_Up || lgLo_Down )
	{
		conv.lgIonStageTrimed = true;
		{
			/* option to print out what has happened so far .... */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && nelem==ipHELIUM )
			{
				fprintf(ioQQQ,"DEBUG ion_trim zone\t%.2f \t%li\t", fnzone, nelem);
				if( lgHi_Down  )
					fprintf(ioQQQ,"\tHi_Down");
				if( lgHi_Up )
					fprintf(ioQQQ,"\tHi_Up");
				if( lgLo_Up )
					fprintf(ioQQQ,"\tLo_Up");
				if( lgLo_Down )
					fprintf(ioQQQ,"\tLo_Down");
				fprintf(ioQQQ,"\n");
			}
		}
	}

	/* asserts that that densities of ion stages that are now turned 
	 * off are zero */
	for( ion=0; ion<dense.IonLow[nelem]; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}
	for( ion=dense.IonHigh[nelem]+1; ion<nelem+1; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}

	for( ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
	{
		/* in case where ionization stage was changed the
		 * density is zero, we set it to SMALLFLOAT since a density
		 * of zero is not possible */
		dense.xIonDense[nelem][ion] = MAX2(dense.xIonDense[nelem][ion], SMALLFLOAT);
	}
	return;
}
