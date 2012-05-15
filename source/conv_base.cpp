/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvBase main routine to drive ionization solution for all species, find total opacity
 * called by ConvIoniz */
/*lgConverg check whether ionization of element nelem has converged */
#include "cddefines.h"
#include "dynamics.h"
#include "trace.h"
#include "elementnames.h"
#include "save.h"
#include "phycon.h"
#include "secondaries.h"
#include "stopcalc.h"
#include "grainvar.h"
#include "highen.h"
#include "dense.h"
#include "hmi.h"
#include "rfield.h"
#include "pressure.h"
#include "taulines.h"
#include "rt.h"
#include "grains.h"
#include "atmdat.h"
#include "ionbal.h"
#include "opacity.h"
#include "cooling.h"
#include "thermal.h"
#include "mole.h"
#include "iso.h"
#include "conv.h"

/*lgIonizConverg check whether ionization of element nelem has converged, called by ionize,
 * returns true if element is converged, false if not */
STATIC bool lgIonizConverg(
	/* atomic number on the C scale, 0 for H, 25 for Fe */
	long int nelem ,
	/* this is allowed error as a fractional change.  Most are 0.15 */
	double delta ,
	/* option to print abundances */
	bool lgPrint )
{
	bool lgConverg_v;
	long int ion;
	double Abund, 
	  bigchange ,
	  change ,
	  one;
	double abundold=0. , abundnew=0.;

	/* this is fractions [ion stage][nelem], ion stage = 0 for atom, nelem=0 for H*/
	static realnum OldFracs[LIMELM+1][LIMELM+1];

	if( !dense.lgElmtOn[nelem] )
		return true;

	DEBUG_ENTRY( "lgIonizConverg()" );

	/* arguments are atomic number, ionization stage
	 * OldFracs[nelem][0] is abundance of nelem (cm^-3) */

	/* this function returns true if ionization of element
	 * with atomic number nelem has not changed by more than delta,*/

	/* check whether abundances exist yet, only do this for after first zone */
	/*if( nzone > 0 )*/
	/* >>chng 03 sep 02, check on changed ionization after first call through,
	 * to insure converged constant eden / temperature models */
	if( conv.nPres2Ioniz )
	{
		/* >>chng 04 aug 31, this had been static, caused last hits on large changes
		 * in ionization */
		bigchange = 0.;
		change = 0.;
		Abund = dense.gas_phase[nelem];

		/* loop over all ionization stages, loop over all ions, not just active ones,
		 * since this also sets old values, and so will zero out non-existant ones */
		for( ion=0; ion <= (nelem+1); ++ion )
		{
			/*lint -e727 OlsdFracs not initialized */
			if( OldFracs[nelem][ion]/Abund > 1e-4 && 
				dense.xIonDense[nelem][ion]/Abund > 1e-4 )
			/*lint +e727 OlsdFracs not initialized */
			{
				/* check change in old to new value */
				one = fabs(dense.xIonDense[nelem][ion]-OldFracs[nelem][ion])/
					OldFracs[nelem][ion];
				change = MAX2(change, one );
				/* remember abundances for largest change */
				if( change>bigchange )
				{
					bigchange = change;
					abundold = OldFracs[nelem][ion]/Abund;
					abundnew = dense.xIonDense[nelem][ion]/Abund;
				}
			}
			/* now set new value */
			OldFracs[nelem][ion] = dense.xIonDense[nelem][ion];
		}

		if( change < delta )
		{
			lgConverg_v = true;
		}
		else
		{
			lgConverg_v = false;
			conv.BadConvIoniz[0] = abundold;
			conv.BadConvIoniz[1] = abundnew;
			ASSERT( abundold>0. && abundnew>0. );
		}
	}
	else
	{
		for( ion=0; ion <= (nelem+1); ++ion )
		{
			OldFracs[nelem][ion] = dense.xIonDense[nelem][ion];
		}

		lgConverg_v = true;
	}

	/* option to print abundances */
	if( lgPrint )
	{
		fprintf(ioQQQ," element %li converged? %c ",nelem, TorF(lgConverg_v));
		for( ion=0; ion<(nelem+1); ++ion )
		{
			fprintf(ioQQQ,"\t%.2e", dense.xIonDense[nelem][ion]/dense.gas_phase[nelem]);
		}
		fprintf(ioQQQ,"\n");
	}
	return( lgConverg_v );
}

/*ConvBase main routine to drive ionization solution for all species, find total opacity
 * called by ConvIoniz
 * return 0 if ok, 1 if abort */
int ConvBase( 
	/* this is zero the first time ConvBase is called by convIoniz,
	 * counts number of call thereafter */
	 long loopi )
{
	double O_HIonRate_New , O_HIonRate_Old;
	double HeatOld,
		EdenTrue_old,
		EdenFromMolecOld,
		EdenFromGrainsOld,
		HeatingOld ,
		CoolingOld;
	realnum hmole_save[N_H_MOLEC];
	static double SecondOld;
	static long int nzoneOTS=-1;
#	define LOOP_ION_LIMIT	10
	long loop_ion;
	static double SumOTS=0. , OldSumOTS[2]={0.,0.};
	double save_iso_grnd[NISO][LIMELM];

	DEBUG_ENTRY( "ConvBase()" );

#if	0
	if( loopi == 0 )
	{
		OldSumOTS[0] = 0.;
		OldSumOTS[1] = 0.;
		SumOTS = 0.;
		nzoneOTS = -1;
	}
#endif

	/* this is set to phycon.te in tfidle, is used to insure that all temp
	 * vars are properly updated when conv_ionizeopacitydo is called 
	 * NB must be same type as phycon.te */
	ASSERT( fp_equal( phycon.te, thermal.te_update ) );

	/* this allows zone number to be printed with slight increment as zone converged
	 * conv.nPres2Ioniz is incremented at the bottom of this routine */
	fnzone = (double)nzone + (double)conv.nPres2Ioniz/100.;

	/* reevaluate pressure */
	/* this sets values of pressure.PresTotlCurr, also calls tfidle,
	 * and sets the total energy content of gas, which may be important for acvection */
	PresTotCurrent();

	/* >>chng 04 sep 15, move EdenTrue_old up here, and will redo at bottom
	 * to find change 
	 * find and save current true electron density - if this changes by more than the
	 * tolerance then ionization solution is not converged */
	/* >>chng 04 jul 27, move eden_sum here from after this loop, so that change in EdenTrue
	 * can be monitored */
	/* update EdenTrue, eden itself is actually changed in ConvEdenIoniz */
	/* must not call eden_sum on very first time since for classic PDR total
	 * ionization may still be zero on first call */
	if( conv.nTotalIoniz )
	{
		if( eden_sum() )
		{
			/* non-zero return indicates abort condition */
			++conv.nTotalIoniz;
			return 1;
		}
	}

	/* the following two were evaluated in eden_sum 
	 * will confirm that these are converged */
	EdenTrue_old = dense.EdenTrue;
	EdenFromMolecOld = co.comole_eden;
	EdenFromGrainsOld = gv.TotalEden;
	HeatingOld = thermal.htot;
	CoolingOld = thermal.ctot;

	/* remember current ground state population - will check if converged */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem<LIMELM;++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* save the ground state population */
				save_iso_grnd[ipISO][nelem] = StatesElemNEW[nelem][nelem-ipISO][0].Pop;
			}
		}
	}

	for( long i=0; i < mole.num_comole_calc; i++ )
	{
		COmole[i]->comole_save = COmole[i]->hevmol;
	}

	for( long i=0; i < N_H_MOLEC; i++ )
	{
		hmole_save[i] = hmi.Hmolec[i];
	}

	if( loopi==0 )
	{
		/* these will be used to look for oscillating ots rates */
		OldSumOTS[0] = 0.;
		OldSumOTS[1] = 0.;
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"   ConvBase called. %.2f Te:%.3e  HI:%.3e  HII:%.3e  H2:%.3e  Ne:%.3e  htot:%.3e  CSUP:%.2e Conv?%c\n", 
		  fnzone,
		  phycon.te, 
		  dense.xIonDense[ipHYDROGEN][0], 
		  dense.xIonDense[ipHYDROGEN][1], 
		  hmi.H2_total,
		  dense.eden, 
		  thermal.htot, 
		  secondaries.csupra[ipHYDROGEN][0] ,
		  TorF(conv.lgConvIoniz) );
	}
	/* want this flag to be true when we exit, various problems will set falst */
	conv.lgConvIoniz = true;

	/* this routine is in heatsum.c, and zeros out the heating array */
	HeatZero();

	/* if this is very first call, say not converged, since no meaningful way to
	 * check on changes in quantities.  this counter is false only on very first
	 * call, when equal to zero */
	if( !conv.nTotalIoniz )
	{
		conv.lgConvIoniz = false;
		strcpy( conv.chConvIoniz, "first call" );
	}

	/* this will be flag to check whether any ionization stages
	 * were trimmed down */
	conv.lgIonStageTrimed = false;

	/* must redo photoionization rates for all species on second and later tries */
	/* always reevaluate the rates when . . . */
	/* this flag says not to redo opac and photo rates, and following test will set
	 * true if one of several tests not done*/
	opac.lgRedoStatic = false;
	if( 
		/* opac.lgOpacStatic (usually true), is set false with no static opacity command */
		!opac.lgOpacStatic || 
		/* we are in search mode */
		conv.lgSearch || 
		/* this is the first call to this zone */
		conv.nPres2Ioniz == 0 )
	{
		/* we need to redo ALL photoionization rates */
		opac.lgRedoStatic = true;
	}
		
	/* deduce dielectronic suppression factors */
	atmdat_DielSupres();

	/* calculate radiative and dielectronic recombination rate coefficients */
	ion_recom_calculate();

	/* this adjusts the lowest and highest stages of ionization we will consider,
	 * only safe to call when lgRedoStatic is true since this could lower the 
	 * lowest stage of ionization, which needs all its photo rates */

	/* conv.nTotalIoniz is only 0 (false) on the very first call to here,
	 * when the ionization distribution is not yet done */
	if( conv.nTotalIoniz )
	{
		bool lgIonizTrimCalled = false;
		static long int nZoneCalled = 0;

		fixit(); // nZoneCalled should be reinitialized for each grid point?

		/* ionization trimming only used for He and heavier, not H */
		/* only do this one time per zone since up and down cycle can occur */
		/* >>chng 05 jan 15, increasing temperature above default first conditions, also
		 * no trim during search - this fixed major logical error when sim is
		 * totally mechanically heated to coronal temperatures -
		 * problem discovered by Ronnie Hoogerwerf */
		/* do not keep changing the trim after the first call within
		 * this zone once we are deep into layer - doing so introduced very 
		 * small level of noise as some stages
		 * went up and down - O+2 - O+3 was good example, when small H+3 after He+ i-front 
		 * limit to one increase per element per zone */
		if( conv.nTotalIoniz>2 &&
			/* only call one time per zone except during search phase,
			 * when only call after 20 times only if temperature has changed */
			( conv.lgSearch || nZoneCalled!=nzone) )
		{
			lgIonizTrimCalled = true;
			for( long nelem=ipHELIUM; nelem<LIMELM; ++nelem )
			{
				if( dense.lgElmtOn[nelem] )
				{
					/* ion_trim will set conv.lgIonStageTrimed true is any ion has its
					 * lowest stage of ionization dropped or trimmed */
					ion_trim(nelem);
				}
			}
			nZoneCalled = nzone;
		}

		/* following block only set of asserts */
#		if !defined(NDEBUG)
		/* check that proper abundances are either positive or zero */
		for( long nelem=ipHELIUM; nelem<LIMELM; ++nelem)
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long ion=0; ion<dense.IonLow[nelem]; ++ion )
				{
					ASSERT( dense.xIonDense[nelem][ion] == 0. );
				}
				/*if( nelem==5 ) fprintf(ioQQQ,"carbbb\t%li\n", dense.IonHigh[nelem]);*/
				for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
				{
					/* >>chng 02 feb 06, had been > o., chng to > SMALLFLOAT to
					* trip over VERY small floats that failed on alphas, but not 386
					* 
					* in case where lower ionization stage was just lowered or
					* trimmed down the abundance
					* was set to SMALLFLOAT so test must be < SMALLFLOAT */
					/* >>chng 02 feb 19, add check for search phase.  During this search
					* models with extreme ionization (all neutral or all ionized) can
					* have extreme but non-zero abundances far from the ionization peak for
					* element with lots of electrons.  These will go away once the model
					* becomes stable */
					/* >>chng 03 dec 01, add check on whether ion trim was called 
					 * conserve.in threw assert when iontrim not called and abund grew small */
					ASSERT( conv.lgSearch || !lgIonizTrimCalled ||
						/* this can happen if all C is in the form of CO */
						(ion==0 && dense.IonHigh[nelem]==0 ) ||
						dense.xIonDense[nelem][ion] >= SMALLFLOAT || 
						dense.xIonDense[nelem][ion]/dense.gas_phase[nelem] >= SMALLFLOAT );
				}
				for( long ion=dense.IonHigh[nelem]+1; ion<nelem+1; ++ion )
				{
					ASSERT( ion >= 0 );
					ASSERT( dense.xIonDense[nelem][ion] == 0. );
				}
			}
		}
#		endif
	}

	/* now check if anything trimmed down */
	if( conv.lgIonStageTrimed )
	{
		/* something was trimmed down, so say that ionization not yet stable */
		/* say that ionization has not converged, secondaries changed by too much */
		conv.lgConvIoniz = false;
		strcpy( conv.chConvIoniz, "IonTrimmed" );
	}

	/* reevaluate advective terms if turned on */
	if( dynamics.lgAdvection )
		DynaIonize();

	/* >>chng 04 feb 15, add loop over ionization until converged.  Non-convergence
	 * in this case means ionization ladder pop changed, probably because of way
	 * that Auger ejection is treated - loop exits when conditions tested at end are met */
	loop_ion = 0;
	do
	{
		if( trace.nTrConvg>=4 )
		{
			/* trace ionization  */
			fprintf( ioQQQ, 
				"    ConvBase4 ionization driver loop_ion %li converged? %c reason not converged %s\n" ,
				loop_ion , 
				TorF( conv.lgConvIoniz) ,
				conv.chConvIoniz );
		}

		/* set the convergence flag to true, 
		 * if anything changes in ConvBase, it will be set false */
		if( loop_ion )
			conv.lgConvIoniz = true;

		/* charge transfer evaluation needs to be here so that same rate 
		 * coefficient used for H ion and other recombination */
		/* fill in master charge transfer array, and zero out some arrays that track effects,
		 * O_HIonRateis rate oxygen ionizes hydrogen - will need to converge this */
		ChargTranEval( &O_HIonRate_Old );

		/* Update densities of species controlled outside the chemical network */
 		CO_update_species_cache();  

 		hmole_reactions();

 		co.nitro_dissoc_rate = 0;

 		/* Update the rate constants -- generic version */
 		CO_update_rks();

		/* find grain temperature, charge, and gas heating rate */
		/* >>chng 04 apr 15, moved here from after call to HeLike(), PvH */
		GrainDrive();

		ASSERT( dense.gas_phase[ipHYDROGEN] < MAX_DENSITY );
		ASSERT( dense.xNucleiTotal < MAX_DENSITY );

		/* update rates and solve level pops for all iso species */
		iso_drive();

		for( long nelem=ipHYDROGEN; nelem<LIMELM; nelem++ )
			ion_wrapper( nelem );

		/* evaluate Compton heating, bound E Compton, cosmic rays */
		highen();

		/* find corrections for three-body rec - collisional ionization */
		atmdat_3body();

		/* deduce dielectronic suppression factors */
		atmdat_DielSupres();

		/*>>chng 04 may 09, add option to abort here, inspired by H2 pop failures
		 * which cascaded downstream if we did not abort */
		/* return if one of above solvers has declared abort condition */
		if( lgAbort )
		{
			++conv.nTotalIoniz;
			return 1;
		}

		/* only reevaluate UTA ionization on first pass through on this zone */
		if( conv.lgFirstSweepThisZone )
		{
			/* inner shell ionization */
			for( long nelem=0; nelem< LIMELM; ++nelem )
			{
				for( long ion=0; ion<nelem+1; ++ion )
				{
					ionbal.UTA_ionize_rate[nelem][ion] = 0.;
					ionbal.UTA_heat_rate[nelem][ion] = 0.;
				}
			}
			/* inner shell ionization by UTA lines */
			/* this flag turned off with no UTA command */
			if( ionbal.lgInnerShellLine_on )
			{
				for( long i=0; i < nUTA; ++i )
				{
					/* rateone is inverse lifetime of level against autoionization */
					double rateone = UTALines[i].Emis->pump * UTALines[i].Emis->AutoIonizFrac;
					ionbal.UTA_ionize_rate[UTALines[i].Hi->nelem-1][UTALines[i].Hi->IonStg-1] += rateone;
					/* heating rate in erg atom-1 s-1 */
					ionbal.UTA_heat_rate[UTALines[i].Hi->nelem-1][UTALines[i].Hi->IonStg-1] += rateone*UTALines[i].Coll.heat;
					{
						/* DEBUG OTS rates - turn this on to debug line, continuum or both rates */
						enum {DEBUG_LOC=false};
						if( DEBUG_LOC /*&& UTALines[i].nelem==ipIRON+1 && (UTALines[i].IonStg==15||UTALines[i].IonStg==14)*/ )
						{
							fprintf(ioQQQ,"DEBUG UTA %3i %3i %.3f %.2e %.2e %.2e\n",
								UTALines[i].Hi->nelem , UTALines[i].Hi->IonStg , UTALines[i].WLAng ,
								rateone, UTALines[i].Coll.heat, 
								UTALines[i].Coll.heat*dense.xIonDense[UTALines[i].Hi->nelem-1][UTALines[i].Hi->IonStg-1] );
						}
					}
				}
			}
		}
	
		/* drive CO chemistry */
		CO_drive();

		/* all elements have now had ionization reevaluated - in some cases we may have upset
		 * the ionization that was forced with an "element ionization" command - here we will 
		 * re-impose that set ionization */
		/* >>chng 04 oct 13, add this logic */
		for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
		{
			if( dense.lgSetIoniz[nelem] )
			{
				dense.IonLow[nelem] = 0;
				dense.IonHigh[nelem] = nelem + 1;
				while( dense.SetIoniz[nelem][dense.IonLow[nelem]] == 0. )
					++dense.IonLow[nelem];
				while( dense.SetIoniz[nelem][dense.IonHigh[nelem]] == 0. )
					--dense.IonHigh[nelem];
			}
		}

		/* redo Oxygen ct ion rate of H to see if it changed */
		ChargTranEval( &O_HIonRate_New );

		/* lgIonizConverg is a function to check whether ionization has converged
		 * check whether ionization changed by more than relative error
		 * given by second number */
		/* >>chng 04 feb 14, loop over all elements rather than just a few */
		for( long nelem=ipHELIUM; nelem<LIMELM; ++nelem )
		{
			if( !lgIonizConverg(nelem, 0.05  ,false ) )
			{
				conv.lgConvIoniz = false;
				sprintf( conv.chConvIoniz , "%2s ion" , elementnames.chElementSym[nelem] );
			}
		}

		if( fabs(EdenTrue_old - dense.EdenTrue) > conv.EdenErrorAllowed/2.*fabs(dense.EdenTrue) )
		{
			conv.lgConvIoniz = false;
			sprintf( conv.chConvIoniz , "EdTr cng" );
			conv.BadConvIoniz[0] = EdenTrue_old;
			conv.BadConvIoniz[1] = dense.EdenTrue;
		}

		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( long nelem=ipISO; nelem<LIMELM; ++nelem )
			{
				if( conv.lgConvIoniz && 
					dense.lgElmtOn[nelem] && 
					dense.xIonDense[nelem][nelem-ipISO]/dense.eden > conv.EdenErrorAllowed/10. &&
					dense.xIonDense[nelem][nelem-ipISO] > SMALLFLOAT )
				{
					/* sim blr_hizqso has largest errors in following */
					double abund = 0.;
					for( long n=0; n<iso.numLevels_local[ipISO][nelem]; ++n )
					{
						abund += StatesElemNEW[nelem][nelem-ipISO][n].Pop;
					}
					double rel_err = fabs(abund-dense.xIonDense[nelem][nelem-ipISO])/
						SDIV(dense.xIonDense[nelem][nelem-ipISO]);

					fixit();// the following test should pass cleanly with much smaller tolerance 
					fixit();// also need to remove the nelem==1

					const realnum err_tol = 0.3f;
					if( nelem==1 && rel_err > err_tol ) 
					{
						if( 0 && loop_ion > 0 )
						{	
							fprintf(ioQQQ,"PROBLEM Inconsistent iso loop_ion %2ld nzone %3ld ipISO %1ld nelem %2ld %e %e %e %s\n",
								loop_ion,nzone,ipISO,nelem,
								abund,
								dense.xIonDense[nelem][nelem-ipISO],
								rel_err,
								iso.chTypeAtomUsed[ipISO][nelem] );
						}
						conv.lgConvIoniz = false;
						sprintf( conv.chConvIoniz , "Iso!=ion" );
						conv.BadConvIoniz[0] = abund;
						conv.BadConvIoniz[1] = dense.xIonDense[nelem][nelem-ipISO];
					}
				}
			}
		}

		++loop_ion;
	}
	/* loop is not converged, less than loop limit, and we are reevaluating */
	while( !conv.lgConvIoniz && loop_ion < LOOP_ION_LIMIT && rfield.lgIonizReevaluate);

	/* >>chng 05 oct 29, move CT heating here from heat_sum since this sometimes contributes
	 * cooling rather than heat and so needs to be sorted out before either heating or cooling 
	 * are derived first find net heating - */
	thermal.char_tran_heat = ChargTranSumHeat();
	/* net is cooling if negative */
	thermal.char_tran_cool = MAX2(0. , -thermal.char_tran_heat );
	thermal.char_tran_heat = MAX2(0. ,  thermal.char_tran_heat );

	/* get total cooling, thermal.ctot = does not occur since passes as pointer.  This can add heat.
	 * it calls coolSum at end to sum up the total cooling */
	CoolEvaluate(&thermal.ctot );

	/* get total heating rate - first save old quantities to check how much it changes */
	HeatOld = thermal.htot;

	/* HeatSum will update ElecFrac, 
	 * secondary electron ionization and excitation efficiencies,
	 * and sum up current secondary rates - remember old values before we enter */
	SecondOld = secondaries.csupra[ipHYDROGEN][0];

	/* update the total heating - it was all set to zero in HeatZero at top of this routine,
	 * occurs before secondaries bit below, since this updates electron fracs */
	HeatSum();

	/* test whether we can just set the rate to the new one, or whether we should
	 * take average and do it again.  secondaries.sec2total was set in hydrogenic, and
	 * is ratio of secondary to total hydrogen destruction rates */
	/* >>chng 02 nov 20, add test on size of csupra - primal had very close to underflow */
	if( (secondaries.csupra[ipHYDROGEN][0]>SMALLFLOAT) && (secondaries.sec2total>0.001) && 
		fabs( 1. - SecondOld/SDIV(secondaries.csupra[ipHYDROGEN][0]) ) > 0.1 && 
		SecondOld > 0. &&  secondaries.csupra[ipHYDROGEN][0] > 0.)
	{
		/* say that ionization has not converged, secondaries changed by too much */
		conv.lgConvIoniz = false;
		strcpy( conv.chConvIoniz, "SecIonRate" );
		conv.BadConvIoniz[0] = SecondOld;
		conv.BadConvIoniz[1] = secondaries.csupra[ipHYDROGEN][0];
	}

#	if 0
	static realnum hminus_old=0.;
	/* >>chng 04 apr 15, add this convergence test */
	if( conv.nTotalIoniz )
	{
		if( fabs( hminus_old-hmi.Hmolec[ipMHm])/SDIV(hmi.Hmolec[ipMHm]) >
			conv.EdenErrorAllowed/4. )
			conv.lgConvIoniz = false;
			strcpy( conv.chConvIoniz, "Big H- chn" );
			conv.BadConvIoniz[0] = hminus_old;
			conv.BadConvIoniz[1] = hmi.Hmolec[ipMHm];
	}
	hminus_old = hmi.Hmolec[ipMHm];
#	endif

	if( HeatOld > 0. && !thermal.lgTemperatureConstant )
	{
		/* check if heating has converged - tolerance is final match */
		if( fabs(1.-thermal.htot/HeatOld) > conv.HeatCoolRelErrorAllowed*.5 )
		{
			conv.lgConvIoniz = false;
			strcpy( conv.chConvIoniz, "Big d Heat" );
			conv.BadConvIoniz[0] = HeatOld;
			conv.BadConvIoniz[1] = thermal.htot;
		}
	}

	/* abort flag may have already been set - if so bail */
	if( lgAbort )
	{

		return 1;
	}

	/* evaluate current opacities 
	 * rfield.lgOpacityReevaluate normally true, 
	 * set false with no opacity reevaluate command, an option to only 
	 * evaluate opacity one time per zone */
	if( !conv.nPres2Ioniz || rfield.lgOpacityReevaluate  || conv.lgIonStageTrimed )
		OpacityAddTotal();

	/* >>chng 02 jun 11, call even first time that this routine is called -
	 * this seemed to help convergence */

	/* do OTS rates for all lines and all continua since
	* we now have ionization balance of all species.  Note that this is not
	* entirely self-consistent, since destruction probabilities here are not the same as
	* the ones used in the model atoms.  Problems??  if near convergence
	* then should be nearly identical */
	if( !conv.nPres2Ioniz || rfield.lgOpacityReevaluate || rfield.lgIonizReevaluate ||
		conv.lgIonStageTrimed || conv.lgSearch || nzone!=nzoneOTS )
	{
		RT_OTS();
		nzoneOTS = nzone;

		/* remember old ots rates */
		OldSumOTS[0] = OldSumOTS[1];
		OldSumOTS[1] = SumOTS;
		/*fprintf(ioQQQ," calling RT_OTS zone %.2f SumOTS is %.2e\n",fnzone,SumOTS);*/

		/* now update several components of the continuum, this only happens after
		 * we have gone through entire solution for this zone at least one time. 
		 * there can be wild ots oscillation on first call */
		/* the rel change of 0.2 was chosen by running hizqso - smaller increased
		 * itrzn but larger did not further decrease it. */
		RT_OTS_Update(&SumOTS);
		/*fprintf(ioQQQ,"RT_OTS_Updateee\t%.3f\t%.2e\t%.2e\n", fnzone,SumOTS , OldSumOTS[1] );*/
	}
	else
		SumOTS = 0.;

	/* now check whether the ots rates changed */
	if( SumOTS> 0. )
	{
		/* the ots rate must be converged to the error in the electron density */
		/* how fine should this be converged?? originally had above, 10%, but take
		 * smaller ratio?? */
		if( fabs(1.-OldSumOTS[1]/SumOTS) > conv.EdenErrorAllowed )
		{
			/* this branch, ionization not converged due to large change in ots rates.
			 * check whether ots rates are oscillating, if loopi > 1 so we have enough info*/
			if( loopi > 1 )
			{
				/* here we have three values, are they changing sign? */
				if( (OldSumOTS[0]-OldSumOTS[1]) * ( OldSumOTS[1] - SumOTS ) < 0. )
				{
					/* ots rates are oscillating, if so then use small fraction of new destruction rates 
					 * when updating Lyman line destruction probabilities in HydroPesc*/
					conv.lgOscilOTS = true;
				}
			}

			conv.lgConvIoniz = false;
			strcpy( conv.chConvIoniz, "OTSRatChng" );
			conv.BadConvIoniz[0] = OldSumOTS[1];
			conv.BadConvIoniz[1] = SumOTS;
		}

		/* produce info on the ots fields if either "trace ots" or 
		 * "trace convergence xxx ots " was entered */
		if( ( trace.lgTrace && trace.lgOTSBug ) ||
			( trace.nTrConvg && trace.lgOTSBug ) )
		{
			RT_OTS_PrtRate(SumOTS*0.05 , 'b' );
		}
		/*fprintf(ioQQQ,"DEBUG opac\t%.2f\t%.3e\t%.3e\n",fnzone,
			dense.xIonDense[ipNICKEL][0] ,
			dense.xIonDense[ipZINC][0] );*/
		{
			/* DEBUG OTS rates - turn this on to debug line, continuum or both rates */
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC && (nzone>110)/**/ )
			{
#				if 0
#				include "lines_service.h"
				DumpLine( &Transitions[ipH_LIKE][ipHYDROGEN][2][0] );
#				endif
				/* last par 'l' for line, 'c' for continua, 'b' for both,
				 * the numbers printed are:
				 * cell i, [i], so 1 less than ipoint
				 * anu[i], 
				 * otslin[i], 
				 * opacity_abs[i],
				 * otslin[i]*opacity_abs[i],
				 * rfield.chLineLabel[i] ,
				 * rfield.line_count[i] */
			}
		}
	}
	{
		/* DEBUG OTS rates - turn this on to debug line, continuum or both rates */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && (nzone>200)/**/ )
		{
			fprintf(ioQQQ,"debug otsss\t%li\t%.3e\t%.3e\t%.3e\n",
				nzone,
				Transitions[0][1][15][3].Emis->ots,
				TauLines[26].Emis->ots,
				opac.opacity_abs[2069]);
		}
	}

	/* option to print OTS continuum with TRACE OTS */
	if( trace.lgTrace && trace.lgOTSBug )
	{
		/* find ots rates here, so we only print fraction of it,
		 * SumOTS is both line and continuum contributing to ots, and is multiplied by opacity */
		/* number is weakest rate to print */
		RT_OTS_PrtRate(SumOTS*0.05 , 'b' );
	}

	// all RT routines called 
	RT_line_all(  );

	/* >>chng 01 mar 16, evaluate pressure here since changing and other values needed */
	/* reevaluate pressure */
	/* this sets values of pressure.PresTotlCurr, also calls tfidle */
	PresTotCurrent();

	/* this checks that all molecules and ions add up to gas phase abundance
	 * but only if we have a converged solution */
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		/* check that species add up if element is turned on, and
		 * ionization is not set */
		if( dense.lgElmtOn[nelem] && !dense.lgSetIoniz[nelem])
		{
			/* this sum of over the molecules did not include the atom and first
			 * ion that is calculated in the molecular solver */
			double sum = dense.xMolecules[nelem];
			for( long ion=0; ion<nelem+2; ++ion )
			{
				sum += dense.xIonDense[nelem][ion];
			}
			ASSERT( sum > SMALLFLOAT );
			/* check that sum agrees with total gas phase abundance */
			/** \todo	2	this assert is not passed if error made much smaller.  This
			 * error should be related to a check on convergence of the molecular
			 * networks and their sum rules, with a criteria used here and there */
			/*if( fabs(sum-dense.gas_phase[nelem])/sum > 1e-2 )
			{
				fprintf(ioQQQ,"DEBUG non conservation element %li sum %.3e gas phase %.3e rel error %.3f\n",
					nelem, sum , dense.gas_phase[nelem],(sum-dense.gas_phase[nelem])/sum);
			}*/
			/* >>chng 04 jul 23, needed to increase to 2e-2 from 1e-2 to get conserve.in to work,
			 * not clear what caused increase in error */
			if( fabs(sum-dense.gas_phase[nelem])/SDIV(sum) > 2e-2 )
			{
				/*fprintf(ioQQQ,"DEBUG non-conv nelem\t%li\tsum\t%e\ttot gas\t%e\trel err\t%.3e\n",
					nelem,
					sum,
					dense.gas_phase[nelem],
					(sum-dense.gas_phase[nelem])/sum);*/
				conv.lgConvIoniz = false;
				strcpy( conv.chConvIoniz, "CO con4" );
				sprintf( conv.chConvIoniz, "COnelem%2li",nelem );
				conv.BadConvIoniz[0] = sum;
				conv.BadConvIoniz[1] = dense.gas_phase[nelem];
			}
		}
	}

	/* update some counters that keep track of how many times this routine
	 * has been called */

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"   ConvBase return. fnzone %.2f nPres2Ioniz %li Te:%.3e  HI:%.3e  HII:%.3e  H2:%.3e  Ne:%.3e  htot:%.3e  CSUP:%.2e Conv?%c reason:%s\n", 
			fnzone,
			conv.nPres2Ioniz,
			phycon.te, 
			dense.xIonDense[ipHYDROGEN][0], 
			dense.xIonDense[ipHYDROGEN][1], 
			hmi.H2_total,
			dense.eden, 
			thermal.htot, 
			secondaries.csupra[ipHYDROGEN][0] ,
			TorF(conv.lgConvIoniz) ,
			conv.chConvIoniz);
	}

	/* this counts number of times we are called by ConvPresTempEdenIoniz, 
	 * number of calls in this zone so first call is zero
	 * reset to zero each time ConvPresTempEdenIoniz is called */
	++conv.nPres2Ioniz;

	/* this is abort option set with SET PRESIONIZ command,
	 * test on nzone since init can take many iterations
	 * this is seldom used except in special cases */
	if( conv.nPres2Ioniz > conv.limPres2Ioniz && nzone > 0)
	{
		fprintf(ioQQQ,"PROBLEM ConvBase sets lgAbort since nPres2Ioniz exceeds limPres2Ioniz.\n");
		fprintf(ioQQQ,"Their values are %li and %li \n",conv.nPres2Ioniz , conv.limPres2Ioniz);
		lgAbort = true;

		return 1;
	}

	/* various checks on the convergence of the current solution */
	if( eden_sum() )
	{
		/* non-zero return indicates abort condition */
		return 1;
	}

	ASSERT( dense.EdenTrue != 0. );
	fixit();  
	// There is no guarantee at this point that dense.EdenTrue is non-zero. 
	// Some models experience FPE here in newmole.  
	// Compelling reason to divide by dense.EdenTrue here rather than multiply on other side?
 
	/* is electron density converged? */
	/** \todo	0	PvH prefers test against err/10 */
	if( fabs(EdenTrue_old - dense.EdenTrue)/dense.EdenTrue > conv.EdenErrorAllowed/2. )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "eden chng" );
		conv.BadConvIoniz[0] = EdenTrue_old;
		conv.BadConvIoniz[1] = dense.EdenTrue;
	}

	/* check on molecular electron den */
	if( fabs(EdenFromMolecOld - co.comole_eden)/dense.EdenTrue > conv.EdenErrorAllowed/2. )
	{
		conv.lgConvIoniz = false;
		sprintf( conv.chConvIoniz, "edn chnCO" );
		conv.BadConvIoniz[0] = EdenFromMolecOld;
		conv.BadConvIoniz[1] = dense.EdenTrue;
	}

	if( gv.lgGrainElectrons )
	{
		/* check on grain electron den */
		if( fabs(EdenFromGrainsOld - gv.TotalEden)/dense.EdenTrue > conv.EdenErrorAllowed/2. )
		{
			conv.lgConvIoniz = false;
			sprintf( conv.chConvIoniz, "edn grn e" );
			conv.BadConvIoniz[0] = EdenFromGrainsOld;
			conv.BadConvIoniz[1] = gv.TotalEden;
		}

		/* check on sum of grain and molecular electron den - often two large numbers that cancel */
		if( fabs( (EdenFromMolecOld-EdenFromGrainsOld) - (co.comole_eden-gv.TotalEden) )/dense.EdenTrue > 
		    conv.EdenErrorAllowed/4. )
		{
			conv.lgConvIoniz = false;
			sprintf( conv.chConvIoniz, "edn mole-grn" );
			conv.BadConvIoniz[0] = (EdenFromMolecOld-EdenFromGrainsOld);
			conv.BadConvIoniz[1] = (co.comole_eden-gv.TotalEden);
		}
	}

	/* check on heating and cooling if vary temp model 
	 * >>chng 08 jul 01, over the code's entire history it had tested whether
	 * this is a constant temperature simulation and did not do this test if
	 * the thermal solution was not done.  There are some cases where we do
	 * want to specify the temperature and then find the heating or cooling -
	 * this is done in calculations of cooling curves for instance.  With this
	 * change the heating/cooling are converged even in a constant temperature
	 * sim.  this does make CT sims run more slowly but with greater accuracy
	 * if heating or cooling is reported */
	if( !thermal.lgTemperatureConstant )
	{
		if( fabs(HeatingOld - thermal.htot)/thermal.htot > conv.HeatCoolRelErrorAllowed/2. )
		{
			conv.lgConvIoniz = false;
			sprintf( conv.chConvIoniz, "heat chg" );
			conv.BadConvIoniz[0] = HeatingOld;
			conv.BadConvIoniz[1] = thermal.htot;
		}

		if( fabs(CoolingOld - thermal.ctot)/thermal.ctot > conv.HeatCoolRelErrorAllowed/2. )
		{
			conv.lgConvIoniz = false;
			sprintf( conv.chConvIoniz, "cool chg" );
			conv.BadConvIoniz[0] = CoolingOld;
			conv.BadConvIoniz[1] = thermal.ctot;
		}
	}

	/* check whether molecular abundances are stable - these were evaluated in CO_drive */
	for( long i=0; i < mole.num_comole_calc; ++i )
	{
		if( dense.lgElmtOn[COmole[i]->nelem_hevmol] && COmole[i]->hevmol>SMALLFLOAT )
		{
			if( fabs(COmole[i]->hevmol-COmole[i]->comole_save)/dense.gas_phase[COmole[i]->nelem_hevmol]-1. >
				conv.EdenErrorAllowed/2.)
			{
				conv.lgConvIoniz = false;
				sprintf( conv.chConvIoniz, "ch %s",COmole[i]->label );
				conv.BadConvIoniz[0] = COmole[i]->comole_save/dense.gas_phase[COmole[i]->nelem_hevmol];
				conv.BadConvIoniz[1] = COmole[i]->hevmol/dense.gas_phase[COmole[i]->nelem_hevmol];
			}
		}
	}

	/* check whether H molecular abundances are stable */
	for( long i=0; i < N_H_MOLEC; ++i )
	{
		if( fabs(hmi.Hmolec[i]-hmole_save[i])/dense.gas_phase[ipHYDROGEN]-1. >
			conv.EdenErrorAllowed/2.)
		{
			conv.lgConvIoniz = false;
			sprintf( conv.chConvIoniz, "ch %s",hmi.chLab[i] );
			conv.BadConvIoniz[0] = hmole_save[i]/dense.gas_phase[ipHYDROGEN];
			conv.BadConvIoniz[1] = hmi.Hmolec[i]/dense.gas_phase[ipHYDROGEN];
		}
	}

	/* >>chng 05 mar 26, add this convergence test - important for molecular or advective
	 * sims since iso ion solver must sync up with chemistry */
	/* remember current ground state population - will check if converged */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem<LIMELM;++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* only do this check for "significant" levels of ionization */
				/*lint -e644 var possibly not init */
				if( dense.xIonDense[nelem][nelem-ipISO]/dense.gas_phase[nelem] > 1e-5 )
				{
					if( fabs(StatesElemNEW[nelem][nelem-ipISO][0].Pop-save_iso_grnd[ipISO][nelem])/SDIV(StatesElemNEW[nelem][nelem-ipISO][0].Pop)-1. >
						conv.EdenErrorAllowed)
					{
						conv.lgConvIoniz = false;
						sprintf( conv.chConvIoniz,"iso %3li%li",ipISO, nelem );
						conv.BadConvIoniz[0] = save_iso_grnd[ipISO][nelem]/SDIV(StatesElemNEW[nelem][nelem-ipISO][0].Pop);
						fixit();  // this looks like a bug, it Pop/Pop, unity for all Pop>SMALLFLOAT
						conv.BadConvIoniz[1] = StatesElemNEW[nelem][nelem-ipISO][0].Pop/SDIV(StatesElemNEW[nelem][nelem-ipISO][0].Pop);
					}
				}
				/*lint +e644 var possibly not init */
			}
		}
	}

	/* this counts how many times ConvBase has been called in this iteration,
	 * located here at bottom of routine so that number is false on first
	 * call, set to 0 in when iteration starts - used to create itr/zn 
	 * number in printout often used to tell whether this is the very 
	 * first attempt at solution in this iteration */
	++conv.nTotalIoniz;

	/* first sweep is over if this is not search phase */
	if( !conv.lgSearch )
		conv.lgFirstSweepThisZone = false;

	/* this was set with STOP NTOTALIONIZ command - only a debugging aid 
	 * by default is zero so false */
	if(StopCalc.nTotalIonizStop && conv.nTotalIoniz> StopCalc.nTotalIonizStop )
	{
		/* non-zero return indicates abort condition */
		fprintf(ioQQQ , "ABORT flag set since STOP nTotalIoniz was set and reached.\n");
		return 1;
	}

	if( save.lgTraceConvergeBase )
	{
		if( iteration > 1 && save.lgTraceConvergeBaseHash )
		{
			static int iter_punch=-1;
			if( iteration !=iter_punch )
				fprintf( save.ipTraceConvergeBase, "%s\n",save.chHashString );
			iter_punch = iteration;
		}

		fprintf( save.ipTraceConvergeBase, 
		"%li\t%.4e\t%.4e\t%.4e\n",
		nzone , thermal.htot , thermal.ctot , dense.eden  );
	}

	return 0;
}
