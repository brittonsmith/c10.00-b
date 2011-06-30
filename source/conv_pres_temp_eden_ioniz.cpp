/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvPresTempEdenIoniz solve for current pressure, calls PressureChange, ConvTempEdenIonize,
 * called by cloudy */
/*ConvFail handle conergece failure */
#include "cddefines.h"
#include "phycon.h"
#include "rt.h"
#include "dense.h"
#include "pressure.h"
#include "trace.h"
#include "conv.h"
#include "thermal.h"
#include "grainvar.h"
#include "grains.h"

/* the limit to the number of loops */
/* >>chng 02 jun 13, from 40 to 50 */
static const int LOOPMAX = 100;

/*ConvPresTempEdenIoniz solve for current pressure, calls PressureChange, ConvTempEdenIoniz,
 * called by cloudy 
 * returns 0 if ok, 1 if disaster */
int ConvPresTempEdenIoniz(void)
{
	long int loop,
		LoopMax=LOOPMAX;
	double 
		hden_old ,
		hden_chng_old ,
		hden_chng,
		/*pres_old ,*/
		/*pres_chng_old ,*/
		/*pres_chng,*/
		/*old_slope,*/
		/*rel_slope,*/
		dP_chng_factor;
	bool lgPresOscil;
	long int nloop_pres_oscil;
	double TemperatureInitial;

	DEBUG_ENTRY( "ConvPresTempEdenIoniz()" );

	/* this will count number of times we call ConvBase in this zone,
	 * counter is incremented there 
	 * zero indicates first pass through solvers on this zone */
	conv.nPres2Ioniz = 0;
	conv.lgFirstSweepThisZone = true;
	conv.lgLastSweepThisZone = false;
	loop = 0;

	/* this will be the limit, which we will increase if no oscillations occur */
	LoopMax = LOOPMAX;
	/* set the initial temperature to the current value, so we will know
	 * if we are trying to jump over a thermal front */
	TemperatureInitial = phycon.te;

	/* this will be flag to check for pressure oscillations */
	lgPresOscil = false;
	/* this is loop where it happened */
	nloop_pres_oscil = 0;
	/* should still be true at end */
	conv.lgConvPops = true;

	/* we will use these to check whether hden oscillating - would need to decrease step size */
	hden_old = dense.gas_phase[ipHYDROGEN];
	hden_chng = 0.;
	/*pres_old = pressure.PresTotlCurr;*/
	/*pres_chng = 0.;*/

	/* this is dP_chng_factor 
	 * cut in half when pressure error changes sign */
	dP_chng_factor = 1.;
	/*rel_slope = 0.;*/

	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			" ConvPresTempEdenIoniz1 entered, will call ConvIoniz to initialize\n");
	}

	/* converge the ionization first, so that we know where we are, and have
	 * a valid foundation to begin the search */
	/* the true electron density dense.EdenTrue is set in eden_sum called by ConvBase */

	/* chng 02 dec 11 rjrw -- ConvIoniz => ConvTempEdenIoniz() here for consistency with inside loop */
	/* ConvIoniz; */
	if( ConvTempEdenIoniz() )
	{
		return 1;
	}

	/* this evaluates current pressure, and returns whether or not 
	 * it is within tolerance of correct pressure */
	conv.lgConvPres = false; /* lgConvPres(); */

	/* convergence trace at this level */
	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			" ConvPresTempEdenIoniz1 ConvIoniz found following converged: Pres;%c, Temp;%c, Eden;%c, Ion:%c, Pops:%c\n", 
			TorF(conv.lgConvPres), 
			TorF(conv.lgConvTemp),
			TorF(conv.lgConvEden),
			TorF(conv.lgConvIoniz),
			TorF(conv.lgConvPops));
	}

	/* >>chng 01 apr 01, add test for at least 2 loops to get better pressure convergence */
	/* >>chng 01 oct 31, add test for number of times converged, for constant
	 * pressure will get two valid solutions */
	/*while( (loop < LoopMax) && !(conv.lgConvPres  &&nPresConverged > 1) &&  !lgAbort )*/
	/* >>chng 02 dec 12, do not demand two constant pressure being valid - should not be
	 * necessary if first one really is valid, and this caused unneeded second evaluation
	 * in the constant density cases */

	/* trace convergence print at this level */
	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			"\n ConvPresTempEdenIoniz1 entering main pressure loop.\n");
	}
	while( (loop < LoopMax) && !conv.lgConvPres &&  !lgAbort )
	{
		/* there can be a pressure or density oscillation early in the search - if not persistent
		 * ok to clear flag */
		/* >>chng 01 aug 24, if large change in temperature allow lots more loops */
		if( fabs( TemperatureInitial - phycon.te )/phycon.te > 0.3 )
			LoopMax = 2*LOOPMAX;
		/* if start of calculation and we are solving for set pressure,
		 * allow a lot more iterations */
		if( nzone <= 1 && pressure.lgPressureInitialSpecified )
			LoopMax = 10*LOOPMAX;

		/* change current densities of all constituents if necessary, 
		 * PressureChange evaluates lgPresOK, true if pressure is now ok
		 * sets CurrentPressure and CorrectPressure */
		hden_old = dense.gas_phase[ipHYDROGEN];
		/*pres_old = pressure.PresTotlCurr;*/

		/* this will evaluate current pressure, update the densities, 
		 * determine the wind velocity, and set conv.lgConvPres,
		 * return value is true if density was changed, false if no changes were necessary 
		 * if density changes then we must redo the temperature and ionization 
		 * PressureChange contains the logic that determines how to change the density to get
		 * the right pressure */
		if( PressureChange( dP_chng_factor ) ) 
		{
			/* heating cooling balance while doing ionization,
			 * this is where the heavy lifting is done, this calls PresTotCurrent,
			 * which sets pressure.PresTotlCurr */
			if( ConvTempEdenIoniz() )
			{
				return 1;
			}
		}

		/* if product of these two is negative then hden is oscillating */
		hden_chng_old = hden_chng;
		/*pres_chng_old = pres_chng;*/
		hden_chng = dense.gas_phase[ipHYDROGEN] - hden_old;
		if( fabs(hden_chng) < SMALLFLOAT )
			hden_chng = sign( (double)SMALLFLOAT, hden_chng );
		/*pres_chng = pressure.PresTotlCurr - pres_old;*/
		/*old_slope = rel_slope;*/
		/*rel_slope = (pres_chng/pressure.PresTotlCurr) / (hden_chng/dense.gas_phase[ipHYDROGEN]);*/

		{
			/*@-redef@*/
			enum{DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC && nzone > 150 && iteration > 1 )
			{
				fprintf(ioQQQ,"%li\t%.2e\t%.2e\t%.2e\n", 
					nzone,
					pressure.PresTotlCurr, 
					pressure.PresTotlCorrect,
					(pressure.PresTotlCorrect - pressure.PresTotlCurr)*100./pressure.PresTotlCorrect
					);
			}
		}

		/* check whether pressure is oscillating */
		/* >>chng 02 may 31, add check on sign on hden changes */
		/* >>chng 02 jun 05, add loop > 1 so that don't trigger off old vals that have
		 * not stabilized yet */
		/* >>chng 04 dec 11, rm test of pressure changing, only check whether hden is changing
		 * very small changes in hden resulted in changes in pressure with some noise due to 
		 * finite precision in te solver - that makes very small oscillations that are not
		 * actually oscillations, only noise */
		if( ( /*( pres_chng*pres_chng_old < 0. )||*/ ( hden_chng*hden_chng_old < 0. ) ) && 
			loop > 1 )
		{
			/*fprintf(ioQQQ,"DEBUG\t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,pres_chng,pres_chng_old ,hden_chng,hden_chng_old);*/
			/* the sign of the change in pressure has changed, so things
			 * are oscillating.  This would be a problem */
			lgPresOscil = true;
			nloop_pres_oscil = loop;
			/* >>chng 04 dec 09, go to this factor becoming smaller every time oscillation occurs */
			dP_chng_factor = MAX2(0.1 , dP_chng_factor * 0.75 );
			/* dP_chng_factor is how pressure changes with density - pass this to
			 * changing routine if it is stable */
			/* >>chng 04 dec 09, take average of old and new dP_chng_factor to avoid pressure
			 * failures in orion_hii_pdr_pp.in 
			if( loop > 4 && old_slope*rel_slope > 0. )
				dP_chng_factor = (1.-FRACNEW)*dP_chng_factor + FRACNEW*rel_slope;*/

			/*fprintf(ioQQQ,"oscilll %li %.2e %.2e %.2e %.2e dP_chng_factor %.2e\n", 
				loop ,
				pres_chng, 
				pres_chng_old,
				hden_chng , 
				hden_chng_old ,
				rel_slope);*/
		}

		/* convergence trace at this level */
		if( trace.nTrConvg>=1  )
		{
			fprintf( ioQQQ, 
				" ConvPresTempEdenIoniz1 %.2f l:%3li nH:%.4e ne:%.4e PCurnt:%.4e PCorct:%.4e err:%6.3f%% dP/dn:%.2e Te:%.4e Osc:%c\n", 
			  fnzone,
			  loop, 
			  dense.gas_phase[ipHYDROGEN], 
			  dense.eden,
			  pressure.PresTotlCurr, 
			  pressure.PresTotlCorrect, 
			  /* this is percentage error */
			  100.*(pressure.PresTotlCurr - pressure.PresTotlCorrect )/pressure.PresTotlCorrect,
			  dP_chng_factor ,
			  phycon.te,
			  TorF(lgPresOscil)  );
		}

		/* increment loop counter */
		++loop;

		/* there can be a pressure or density oscillation early in the search - if not persistent
		 * ok to clear flag 
		 * >>chng 04 sep 22, add this logic */
		if( loop - nloop_pres_oscil > 4 )
			lgPresOscil = false;

		/* if we hit limit of loop, but no oscillations have happened, then we are
		 * making progress, and can keep going */
		if( loop == LoopMax && !lgPresOscil )
		{
			LoopMax = MIN2( 100 , LoopMax*2 );
		}
	}

	/* keep track of minimum and maximum temperature */
	thermal.thist = max((realnum)phycon.te,thermal.thist);
	thermal.tlowst = min((realnum)phycon.te,thermal.tlowst);

	/* >>chng 04 jan 31, now that all of the physics is converged, determine grain drift velocity */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
		GrainDrift();

	/* >>chng 01 mar 14, all ConvFail one more time, no matter how
	 * many failures occurred below.  Had been series of if, so multiple
	 * calls per failure possible. */
	/* >>chng 04 au 07, only announce pres fail here,
	 * we did not converge the pressure */
	if( !conv.lgConvIoniz )
		ConvFail("ioni","");
	else if( !conv.lgConvEden )
		ConvFail("eden","");
	else if( !conv.lgConvTemp )
		ConvFail("temp","");
	else if( !conv.lgConvPres )
		ConvFail("pres","");

	/* this is only a sanity check that the summed continua accurately reflect
	 * all of the individual components.  Only include this when NDEBUG is not set,
	 * we are in not debug compile */
#	if !defined(NDEBUG)
	RT_OTS_ChkSum(0);
#	endif

	return 0;
}
