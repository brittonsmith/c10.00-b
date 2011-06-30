/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtComment analyze model, generating comments on its features */
/*chkCaHeps check whether CaII K and H epsilon overlap */
/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
#include "cddefines.h"
#include "physconst.h"
#include "cddrive.h"
#include "lines_service.h"
#include "iso.h"
#include "continuum.h"
#include "stopcalc.h"
#include "hyperfine.h"
#include "dense.h"
#include "grainvar.h"
#include "version.h"
#include "rt.h"
#include "he.h"
#include "ionbal.h"
#include "taulines.h"
#include "hydrogenic.h"
#include "lines.h"
#include "trace.h"
#include "hcmap.h"
#include "hmi.h"
#include "save.h"
#include "h2.h"
#include "conv.h"
#include "dynamics.h"
#include "opacity.h"
#include "geometry.h"
#include "elementnames.h"
#include "ca.h"
#include "broke.h"
#include "pressure.h"
#include "mole.h"
#include "atoms.h"
#include "abund.h"
#include "rfield.h"
#include "colden.h"
#include "phycon.h"
#include "timesc.h"
#include "hextra.h"
#include "radius.h"
#include "iterations.h"
#include "fudgec.h"
#include "called.h"
#include "magnetic.h"
#include "wind.h"
#include "secondaries.h"
#include "struc.h"
#include "oxy.h"
#include "input.h"
#include "thermal.h"
#include "atmdat.h"
#include "warnings.h"

/*chkCaHeps check whether CaII K and H epsilon overlap */
STATIC void chkCaHeps(double *totwid);

/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
STATIC void prt_smooth_predictions(void);

void PrtComment(void)
{
	char chLbl[11], 
	  chLine[INPUT_LINE_LENGTH];

	bool lgAbort_flag,
	  lgThick,
	  lgLots_of_moles;

	long int dum1,
	  dum2,
	  i, 
	  imas, 
	  ipLo ,
	  ipHi ,
	  ipISO,
	  nelem, 
	  isav, 
	  j, 
	  nc, 
	  nline, 
	  nn, 
	  nneg, 
	  ns, 
	  nw;

	double big_ion_jump, 
	  absint ,
	  aj, 
	  alpha, 
	  big, 
	  bigm, 
	  comfrc, 
	  differ, 
	  error, 
	  flur, 
	  freqn, 
	  rate, 
	  ratio, 
	  rel, 
	  small, 
	  tauneg, 
	  ts ,
	  HBeta, 
	  relfl ,
	  relhm,  
	  fedest,
	  GetHeat, 
	  SumNeg, 
	  c4363, 
	  t4363, 
	  totwid;

	double VolComputed , VolExpected , ConComputed , ConExpected;

	bool lgLotsSolids;

	DEBUG_ENTRY( "PrtComment()" );

	if( 0 && lgAbort )
	{
		return;
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PrtComment called.\n" );
	}

	/* 
	 * enter all comments cautions warnings and surprises into a large
	 * stack of statements
	 * at end of this routine is master call to printing routines
	 */
	iterations.lgIterAgain = false;

	/* initialize the line saver */
	wcnint();

	if( t_version::Inst().nBetaVer > 0 )
	{
		sprintf( chLine, 
			"  !This is beta test version %ld and is intended for testing only.", 
		  t_version::Inst().nBetaVer );
		bangin(chLine);
	}

	/* this flag set by call to fixit routine,
	 * to show that parts of the code need repair. 
	 * lgRelease is true if this is release version */
	if( broke.lgFixit && !t_version::Inst().lgRelease )
	{
		sprintf( chLine, "  !The code needs to be fixed - search for fixit()." );
		bangin(chLine);
	}

	/* this flag set by call to CodeReview routine,
	* to show that parts of the code need to be reviewed. 
	* lgRelease is true if this is release version */
	if( broke.lgCheckit  && !t_version::Inst().lgRelease )
	{
		sprintf( chLine, "  !New code needs to be reviewed - search for CodeReview()." );
		bangin(chLine);
	}

	/* say why calculation stopped */
	if( conv.lgBadStop )
	{
		/* this case stop probably was not intended */
		sprintf( warnings.chRgcln[0], " W-Calculation stopped because %s Iteration%3ld of %ld", 
		  StopCalc.chReasonStop, iteration, iterations.itermx + 1 );
		sprintf( chLine, " W-Calculation stopped because %s", 
		  StopCalc.chReasonStop );
		warnin(chLine);
		sprintf( chLine, " W-This was not intended." );
		warnin(chLine);
	}
	else
	{
		/* for iterate to convergence, print reason why it was not converged on 3rd and higher iterations */
		if( (conv.lgAutoIt && iteration != iterations.itermx + 1) && 
		  iteration > 2 )
		{
			sprintf( warnings.chRgcln[0], 
				"   Calculation stopped because %s Iteration %ld of %ld, not converged due to %s", 
			  StopCalc.chReasonStop, 
			  iteration, 
			  iterations.itermx + 1, 
			  conv.chNotConverged  );
		}
		else
		{
			sprintf( warnings.chRgcln[0], 
				"   Calculation stopped because %s Iteration %ld of %ld", 
			  StopCalc.chReasonStop, iteration, iterations.itermx + 1 );
		}
	}

	/* check whether stopped because default number of zones hit,
	 * not intended?? */
	if( (!geometry.lgZoneSet) && geometry.lgZoneTrp )
	{
		conv.lgBadStop = true;
		sprintf( chLine, 
			" W-Calculation stopped because default number of zones reached.  Was this intended???" );
		warnin(chLine);
		sprintf( chLine, 
			" W-Default limit can be increased while retaining this check with the SET NEND command." );
		warnin(chLine);
	}

	/* check whether stopped because zones too thin - only happens when glob set
	 * and depth + dr = depth
	 * not intended */
	if( radius.lgDrMinUsed || radius.lgdR2Small )
	{
		conv.lgBadStop = true;
		sprintf( chLine, 
			" W-Calculation stopped zone thickness became too thin.   This was not intended." );
		warnin(chLine);
		sprintf( chLine, 
			" W-The most likely reason was an uncontrolled oscillation." );
		warnin(chLine);
		ShowMe();
	}

	if( radius.lgdR2Small )
	{
		sprintf( chLine, 
			" W-This happened because the globule scale became very small relative to the depth." );
		warnin(chLine);
		sprintf( chLine, 
			" W-This problem is described in Hazy." );
		warnin(chLine);
	}

	/* possible that last zone does not have stored temp - if so
	 * then make it up - this happens for some bad stops */
	ASSERT( nzone < struc.nzlim );

	if( struc.testr[nzone-1] == 0. )
		struc.testr[nzone-1] = struc.testr[nzone-2];

	if( struc.ednstr[nzone-1] == 0. )
		struc.ednstr[nzone-1] = struc.ednstr[nzone-2];

	/* give indication of geometry */
	rel = radius.depth/radius.rinner;
	if( rel < 0.1 )
	{
		sprintf( warnings.chRgcln[1], "   The geometry is plane-parallel." );
	}
	else if( rel >= 0.1 && rel < 3. )
	{
		sprintf( warnings.chRgcln[1], "   The geometry is a thick shell." );
	}
	else
	{
		sprintf( warnings.chRgcln[1], "   The geometry is spherical." );
	}

	/* levels of warnings: Warning   (possibly major problems)
	 *                     Caution   (not likely to invalidate the results)
	 *                     [      !] surprise, but not a problem
	 *                     [nothing] interesting note
	 */
	/* this flag set by call to routine broken ( ); 
	 * and show that the code is broken. */
	if( broke.lgBroke )
	{
		sprintf( chLine, " W-The code is broken - search for broken()." );
		warnin(chLine);
	}

	/* incorrect electron density detected */
	if( dense.lgEdenBad )
	{
		if( dense.nzEdenBad == nzone )
		{
			sprintf( chLine, " C-The assumed electron density was incorrect for the last zone." );
			caunin(chLine);
			sprintf( chLine, " C-Did a temperature discontinuity occur??" );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, " W-The assumed electron density was incorrect during the calculation.  This is bad." );
			warnin(chLine);
			ShowMe();
		}
	}

	if( lgAbort )
	{
		sprintf( chLine, " W-The calculation aborted.  Something REALLY went wrong!" );
		warnin(chLine);
	}

	/* thermal map was done but results were not ok */
	if( hcmap.lgMapDone && !hcmap.lgMapOK )
	{
		sprintf( chLine, "  !The thermal map had changes in slope - check map output." );
		bangin(chLine);
	}

	/* first is greater than zero if fudge factors were entered, second is
	 * true if fudge ever evaluated, even to see if fudge factors are in place */
	if( fudgec.nfudge > 0 || fudgec.lgFudgeUsed )
	{
		sprintf( chLine, "  !Fudge factors were used or were checked.  Why?" );
		bangin(chLine);
	}

	if( dense.gas_phase[ipHYDROGEN] > 1.1e13 )
	{
		if( dense.gas_phase[ipHYDROGEN] > 1e15 )
		{
			sprintf( chLine, " C-Density greater than 10**15, heavy elements are very uncertain." );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, " C-Density greater than 10**13" );
			caunin(chLine);
		}
	}

	/* HBeta is used later in the code to check on line intensities */
	if( cdLine("Pump",4861.36f,&relfl,&absint)<=0 )
	{
		fprintf( ioQQQ, " PROBLEM Did not find Pump H-beta, set to unity\n" );
		relfl = 1.;
		absint = 1.;
	}

	/* now find total Hbeta */
	/* >>chng from "totl" Hbeta which was a special entry, to "H  1" Hbeta, which 
	 * is the general case */
	if( cdLine( "H  1",4861.36f,&HBeta,&absint)<=0 )
	{
		fprintf( ioQQQ, " NOTE Did not find H  1 H-beta - set intensity to unity, "
			"will not check on importance of H 1 pumping.\n" );
		HBeta = 1.;
		absint = 1.;
	}
	else 
	{
		/* check on continuum pumping of Balmer lines */
		if( HBeta>SMALLFLOAT )
		{
			flur = relfl/HBeta;
			if( flur > 0.1 )
			{
				sprintf( chLine, "  !Continuum fluorescent production of H-beta was very important." );
				bangin(chLine);
			}
			else if(flur > 0.01 )
			{
				sprintf( chLine, "   Continuum fluorescent production of H-beta was significant." );
				notein(chLine);
			}
		}
	}

	// iterate to convergence - status of this iteration
	if( iteration > 1 && conv.lgAutoIt && iteration<iterations.itermx)
	{
		sprintf( chLine , "   Iteration not converged because %s.",
			conv.chNotConverged );
		notein(chLine);
	}

	/* check if there were problems with initial wind velocity */
	if( wind.lgBallistic() && ((!wind.lgWindOK) || wind.windv < 1e6) )
	{
		sprintf( chLine, " C-Wind velocity below sonic point; solution is not valid." );
		caunin(chLine);
	}

	/* now confirm that mass flux here is correct */
	if( !wind.lgStatic() )
	{
		rel = wind.emdot/(wind.windv*dense.gas_phase[ipHYDROGEN])/radius.r1r0sq;
		if( fabs(1.-rel)> 0.02 )
		{
			sprintf( chLine, " C-Wind mass flux error is %g%%",fabs(1.-rel)*100. );
			caunin(chLine);
			fprintf(ioQQQ,"DEBUG emdot\t%.3e\t%.3e\t%.3e\t%.3e\n",
				wind.emdot , wind.windv*dense.gas_phase[ipHYDROGEN],wind.windv,dense.gas_phase[ipHYDROGEN]);
		}
	}

	/* check that we didn't overrun zone scale */
	if( nzone >= struc.nzlim )
	{
		TotalInsanity();
	}

	// check on energy conservation
	if( !lgConserveEnergy() )
	{
		ShowMe();
		lgAbort = true;
		fprintf(ioQQQ,"\n\n Enable per zone energy conservation check by setting "
				"CHECK_ENERGY_EVERY_ZONE=true in cloudy.cpp, recompile, then rerun.\n");
	}

	/* comments having to do with cosmic rays */
	/* comment if cosmic rays and magnetic field both present */
	if( hextra.cryden*magnetic.lgB > 0. )
	{
		sprintf( chLine, 
			"  !Magnetic field & cosmic rays both present.  Their interactions are not treated." );
		bangin(chLine);
	}

	/* comment if cosmic rays are not included and stop temp has been lowered to go into neutral gas */
	if( hextra.cryden== 0. && StopCalc.TempLoStopZone < phycon.TEMP_STOP_DEFAULT)
	{
		sprintf( chLine, 
			"  !Background cosmic rays are not included - is this physical?  It affects the chemistry." );
		bangin(chLine);
	}

	/* check whether cosmic rays on, but model thick to them */
	if( hextra.cryden > 0. && (colden.colden[ipCOL_H0]/10. + colden.colden[ipCOL_Hp]) > 1e23 )
	{
		sprintf( chLine, 
			" C-Model is thick to cosmic rays, which are on." );
		caunin(chLine);
	}

	/* was ionization rate less than cosmic ray ionization rate in ISM? */
	if( hextra.cryden == 0. && iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s] < 1e-17 )
	{
		sprintf( chLine, 
			"  !Ionization rate fell below background cosmic ray ionization rate.  Should this be added too?" );
		bangin(chLine);
		sprintf( chLine, 
			"  !   Use the COSMIC RAY BACKGROUND command." );
		bangin(chLine);
	}

	/* PrtComment if test code is in place */
	if( lgTestCodeCalled )
	{
		sprintf( chLine, "  !Test code is in place." );
		bangin(chLine);
	}

	/* lgComUndr set to .true. if Compton cooling rate underflows to 0 */
	if( rfield.lgComUndr )
	{
		sprintf( chLine, 
			"  !Compton cooling rate underflows to zero.  Is this important?" );
		bangin(chLine);
	}

	/* make note if input stream contained an underscore, which was converted into a space */
	if( input.lgUnderscoreFound )
	{
		sprintf( chLine, 
			"  !Some input lines contained underscores, these were changed to spaces." );
		bangin(chLine);
	}

	/* make note if input stream contained a left or right bracket, which was converted into a space */
	if( input.lgBracketFound )
	{
		sprintf( chLine, 
			"  !Some input lines contained [ or ], these were changed to spaces." );
		bangin(chLine);
	}

	/* lgHionRad set to .true. if no hydrogen ionizing radiation */
	if( rfield.lgHionRad )
	{
		sprintf( chLine, 
			"  !There is no hydrogen-ionizing radiation.  Was this intended?" );
		bangin(chLine);
	}

	/* check whether certain zones were thermally unstable */
	if( thermal.nUnstable > 0 )
	{
		sprintf( chLine, 
			"   Derivative of net cooling negative and so possibly thermally unstable in%4ld zones.", 
		  thermal.nUnstable );
		notein(chLine);
	}

	/* generate a bang if a large fraction of the zones were unstable */
	if( nzone > 1 && 
		(realnum)(thermal.nUnstable)/(realnum)(nzone) > 0.25 )
	{
		sprintf( chLine, 
			"  !A large fraction of the zones were possibly thermally unstable,%4ld out of%4ld", 
		  thermal.nUnstable, nzone );
		bangin(chLine);
	}

	/* comment if negative coolants were ever significant */
	if( thermal.CoolHeatMax > 0.2 )
	{
		sprintf( chLine, 
			"  !Negative cooling reached %6.1f%% of the local heating, due to %4.4s %.1f", 
		  thermal.CoolHeatMax*100., thermal.chCoolHeatMax, thermal.wlCoolHeatMax );
		bangin(chLine);
	}
	else if( thermal.CoolHeatMax > 0.05 )
	{
		sprintf( chLine, 
			"   Negative cooling reached %6.1f%% of the local heating, due to %4.4s %.2f", 
		  thermal.CoolHeatMax*100., thermal.chCoolHeatMax, thermal.wlCoolHeatMax );
		notein(chLine);
	}

	/* check if advection heating was important */
	if( dynamics.HeatMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Advection heating reached %.2f%% of the local heating.", 
		  dynamics.HeatMax*100. );
		bangin(chLine);
	}
	else if( dynamics.HeatMax > 0.005 )
	{
		sprintf( chLine, 
			"   Advection heating reached %.2f%% of the local heating.", 
		  dynamics.HeatMax*100. );
		notein(chLine);
	}

	/* check if advection cooling was important */
	if( dynamics.CoolMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Advection cooling reached %.2f%% of the local cooling.", 
		  dynamics.CoolMax*100. );
		bangin(chLine);
	}
	else if( dynamics.CoolMax > 0.005 )
	{
		sprintf( chLine, 
			"   Advection cooling reached %.2f%% of the local heating.", 
		  dynamics.CoolMax*100. );
		notein(chLine);
	}

	/* >>chng 06 mar 22, add this comment
	 * check if time dependent ionization front being done with too large a U */
	if( dynamics.lgTimeDependentStatic && dynamics.lgRecom )
	{
		if( rfield.uh > 1. )
		{
			sprintf( chLine, 
				" W-Time dependent ionization front cannot now handle strong-R cases - the ionization parameter is too large." );
			warnin(chLine);
		}
		else if( rfield.uh > 0.1 )
		{
			sprintf( chLine, 
				" C-Time dependent ionization front cannot now handle strong-R cases - the ionization parameter is too large." );
			caunin(chLine);
		}
	}

	/* check if thermal ionization of ground state of hydrogen was important */
	if( hydro.HCollIonMax > 0.10 )
	{
		sprintf( chLine, 
			"  !Thermal collisional ionization of H reached %.2f%% of the local ionization rate.", 
		  hydro.HCollIonMax*100. );
		bangin(chLine);
	}
	else if( hydro.HCollIonMax > 0.005 )
	{
		sprintf( chLine, 
			"   Thermal collisional ionization of H reached %.2f%% of the local ionization rate.", 
		  hydro.HCollIonMax*100. );
		notein(chLine);
	}

	/* check if lookup table for Hummer & Storey case B was exceeded */
	if( !atmdat.lgHCaseBOK[1][ipHYDROGEN]  )
	{
		sprintf( chLine, 
			"   Te-ne bounds of Case B lookup table exceeded, H I Case B line intensities set to zero." );
		notein(chLine);
	}
	if( !atmdat.lgHCaseBOK[1][ipHELIUM]  )
	{
		sprintf( chLine, 
			"   Te-ne bounds of Case B lookup table exceeded, He II Case B line intensities set to zero." );
		notein(chLine);
	}

	if( dense.EdenMax>1e8  )
	{
		sprintf( chLine, 
			"  !The high electron density makes the Nussbaumer/Storey CNO recombination predictions unreliable." );
		bangin(chLine);
	}

	/* check if secondary ionization of hydrogen was important */
	if( secondaries.SecHIonMax > 0.10 )
	{
		sprintf( chLine, 
			"  !Suprathermal collisional ionization of H reached %.2f%% of the local H ionization rate.", 
		  secondaries.SecHIonMax*100. );
		bangin(chLine);
	}
	else if( secondaries.SecHIonMax > 0.005 )
	{
		sprintf( chLine, 
			"   Suprathermal collisional ionization of H reached %.2f%% of the local H ionization rate.", 
		  secondaries.SecHIonMax*100. );
		notein(chLine);
	}

	/* check if H2 vib-deexcitation heating was important */
	if( hmi.HeatH2DexcMax > 0.05 )
	{
		sprintf( chLine, 
			"  !H2 vib deexec heating reached %.2f%% of the local heating.", 
		  hmi.HeatH2DexcMax*100. );
		bangin(chLine);
	}
	else if( hmi.HeatH2DexcMax > 0.005 )
	{
		sprintf( chLine, 
			"   H2 vib deexec heating reached %.2f%% of the local heating.", 
		  hmi.HeatH2DexcMax*100. );
		notein(chLine);
	}

	/* check if H2 vib-deexcitation heating was important */
	if( hmi.CoolH2DexcMax > 0.05 )
	{
		sprintf( chLine, 
			"  !H2 deexec cooling reached %.2f%% of the local heating.", 
		  hmi.CoolH2DexcMax*100. );
		bangin(chLine);
	}
	else if( hmi.CoolH2DexcMax > 0.005 )
	{
		sprintf( chLine, 
			"   H2 deexec cooling reached %.2f%% of the local heating.", 
		  hmi.CoolH2DexcMax*100. );
		notein(chLine);
	}

	/* check if charge transfer ionization of hydrogen was important */
	if( atmdat.HIonFracMax > 0.10 )
	{
		sprintf( chLine, 
			"  !Charge transfer ionization of H reached %.1f%% of the local H ionization rate.", 
		  atmdat.HIonFracMax*100. );
		bangin(chLine);
	}
	else if( atmdat.HIonFracMax > 0.005 )
	{
		sprintf( chLine, 
			"   Charge transfer ionization of H reached %.2f%% of the local H ionization rate.", 
		  atmdat.HIonFracMax*100. );
		notein(chLine);
	}

	/* check if charge transfer heating cooling was important */
	if( atmdat.HCharHeatMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Charge transfer heating reached %.2f%% of the local heating.", 
		  atmdat.HCharHeatMax*100. );
		bangin(chLine);
	}
	else if( atmdat.HCharHeatMax > 0.005 )
	{
		sprintf( chLine, 
			"   Charge transfer heating reached %.2f%% of the local heating.", 
		  atmdat.HCharHeatMax*100. );
		notein(chLine);
	}

	if( atmdat.HCharCoolMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Charge transfer cooling reached %.2f%% of the local heating.", 
		  atmdat.HCharCoolMax*100. );
		bangin(chLine);
	}
	else if( atmdat.HCharCoolMax > 0.005 )
	{
		sprintf( chLine, 
			"   Charge transfer cooling reached %.2f%% of the local heating.", 
		  atmdat.HCharCoolMax*100. );
		notein(chLine);
	}

	/* check whether photo from up level of Mg2 2798 ever important */
	if( atoms.xMg2Max > 0.1 )
	{
		sprintf( chLine, 
			"  !Photoionization of upper level of Mg II 2798 reached %.1f%% of the total Mg+ photo rate.", 
		  atoms.xMg2Max*100. );
		bangin(chLine);
	}
	else if( atoms.xMg2Max > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of upper level of Mg II 2798 reached %.1f%% of the total Mg+ photo rate.", 
		  atoms.xMg2Max*100. );
		notein(chLine);
	}

	/* check whether photo from up level of [O I] 6300 ever important */
	if( oxy.poimax > 0.1 )
	{
		sprintf( chLine, 
			"  !Photoionization of upper levels of [O I] reached %.1f%% of the total O destruction rate.", 
		  oxy.poimax*100. );
		bangin(chLine);
	}
	else if( oxy.poimax > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of upper levels of [O I] reached %.1f%% of the total O destruction rate.", 
		  oxy.poimax*100. );
		notein(chLine);
	}

	/* check whether photo from up level of [O III] 5007 ever important */
	if( (oxy.poiii2Max + oxy.poiii3Max) > 0.1 )
	{
		sprintf( chLine, 
			"  !Photoionization of upper levels of [O III] reached %.1f%% of the total O++ photo rate.", 
		  (oxy.poiii2Max + oxy.poiii3Max)*100. );
		bangin(chLine);
	}
	else if( (oxy.poiii2Max + oxy.poiii3Max) > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of upper levels of [O III] reached %.1f%% of the total O++ photo rate.", 
		  (oxy.poiii2Max + oxy.poiii3Max)*100. );
		notein(chLine);
	}

	/* check whether photoionization of He 2trip S was important */
	if( he.frac_he0dest_23S > 0.1 )
	{
		sprintf( chLine, 
			"  !Destruction of He 2TriS reached %.1f%% of the total He0 dest rate"
			" at zone %li, %.1f%% of that was photoionization.", 
		  he.frac_he0dest_23S*100, 
		  he.nzone, 
		  he.frac_he0dest_23S_photo*100.  );
		bangin(chLine);
	}
	else if( he.frac_he0dest_23S > 0.01 )
	{
		sprintf( chLine, 
			"   Destruction of He 2TriS reached %.1f%% of the total He0 dest rate"
			" at zone %li, %.1f%% of that was photoionization.", 
		  he.frac_he0dest_23S*100, 
		  he.nzone, 
		  he.frac_he0dest_23S_photo*100.  );
		notein(chLine);
	}

	/* check for critical density for l-mixing */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		if( !iso.lgCritDensLMix[ipISO] && dense.lgElmtOn[ipISO] )
		{
			sprintf( chLine,
				"   The density is too low to l-mix the lowest %s I collapsed level. "
				" More resolved levels are needed for accurate line ratios.",
				elementnames.chElementSym[ipISO]);
			notein(chLine);
		}
	}

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* report continuum lowering for xx-like iso-sequence. */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( iso.lgLevelsLowered[ipISO][nelem] && dense.lgElmtOn[nelem] )
			{
				sprintf( chLine, "  !Continuum was lowered into model %s-like %s due to high density.  Highest n is %li",
					elementnames.chElementSym[ipISO],
					elementnames.chElementSym[nelem],
					iso.n_HighestResolved_local[ipISO][nelem]+iso.nCollapsed_local[ipISO][nelem]);
				bangin(chLine);
			}
			else if( iso.lgLevelsEverLowered[ipISO][nelem] && dense.lgElmtOn[nelem] )
			{
				sprintf( chLine, "  !Continuum was lowered into model %s-like %s due to high density at SOME point but NOT at the last zone.",
					elementnames.chElementSym[ipISO],
					elementnames.chElementNameShort[nelem]);
				bangin(chLine);
			}
		}
	}

	/* frequency array may not have been defined for all energies */
	if( !rfield.lgMMok )
	{
		sprintf( chLine, 
			" C-Continuum not defined in extreme infrared - Compton scat, grain heating, not treated properly?" );
		caunin(chLine);
	}

	if( !rfield.lgHPhtOK )
	{
		sprintf( chLine, 
			" C-Continuum not defined at photon energies which ionize excited states of H, important for H- and ff heating." );
		caunin(chLine);
	}

	if( !rfield.lgXRayOK )
	{
		sprintf( chLine, 
			" C-Continuum not defined at X-Ray energies - Compton scattering and Auger ionization wrong?" );
		caunin(chLine);
	}

	if( !rfield.lgGamrOK )
	{
		sprintf( chLine, 
			" C-Continuum not defined at gamma-ray energies - pair production and Compton scattering OK?" );
		caunin(chLine);
	}

	if( continuum.lgCon0 )
	{
		sprintf( chLine, " C-Continuum zero at some energies." );
		caunin(chLine);
	}

	if( continuum.lgCoStarInterpolationCaution )
	{
		sprintf( chLine , " C-CoStarInterpolate interpolated between non-adjoining tracks, this may not be valid." );
		caunin(chLine);
	}

	if( rfield.lgOcc1Hi )
	{
		sprintf( chLine, 
			"  !The continuum occupation number at 1 Ryd is greater than unity." );
		bangin(chLine);
	}

	/* this flag set true it set dr forced first zone to be too big */
	if( radius.lgDR2Big )
	{
		sprintf( chLine, 
			" C-The thickness of the first zone was set larger than optimal by a SET DR command." );
		caunin(chLine);
		/* this is case where did one zone of specified thickness - but it 
		 * was too large */
		if( nzone<=1 )
			sprintf( chLine, 
			" C-Consider using the STOP THICKNESS command instead." );
		caunin(chLine);
	}

	/* check whether non-col excitation of 4363 was important */
	if( cdLine("TOTL",4363,&t4363,&absint)<=0 )
	{
		fprintf( ioQQQ, " PrtComment could not find total O III 4363 with cdLine.\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( cdLine("Coll",4363,&c4363,&absint)<=0 )
	{
		fprintf( ioQQQ, " PrtComment could not find collisional O III 4363 with cdLine.\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* only print this comment if 4363 is significant and density low */
	if( HBeta > 1e-35 )
	{
		/* >>chng 02 feb 27, lower ratio from -4 to -5, and raise density from 7 to 8 */
		if( t4363/HBeta > 1e-5 && dense.gas_phase[ipHYDROGEN] < 1e8 )
		{
			ratio = (t4363 - c4363)/t4363;
			if( ratio > 0.01 )
			{
				sprintf( chLine, 
					"  !Non-collisional excitation of [O III] 4363 reached %.2f%% of the total.", 
				  ratio*100. );
				bangin(chLine);
			}
			else if( ratio > 0.001 )
			{
				sprintf( chLine, 
					"   Non-collisional excitation of [O III] 4363 reached %.2f%% of the total.", 
				  ratio*100. );
				notein(chLine);
			}
		}
	}

	/* check for plasma shielding */
	if( rfield.lgPlasNu )
	{
		sprintf( chLine, 
			"  !The largest plasma frequency was %.2e Ryd = %.2e micron  The continuum is set to 0 below this.", 
		  rfield.plsfrqmax,
		  /* wavelength in microns */
		  RYDLAM/rfield.plsfrqmax/1e4);
		bangin(chLine);
	}

	if( rfield.occmax > 0.1 )
	{
		if( rfield.occmnu > 1e-4 )
		{
			sprintf( chLine, 
				"  !The largest continuum occupation number was %.3e at %.3e Ryd.", 
			  rfield.occmax, rfield.occmnu );
			bangin(chLine);
		}
		else
		{
			/* not surprising if occupation number bigger than 1 around 1e-5 Ryd,
			 * since this is the case for 3K background */
			sprintf( chLine, 
				"   The largest continuum occupation number was %.3e at %.3e Ryd.", 
			  rfield.occmax, rfield.occmnu );
			notein(chLine);
		}
	}

	if( rfield.occmax > 1e4 && rfield.occ1nu > 0. )
	{
		/* occ1nu is energy (ryd) where continuum occupation number falls below 1 */
		if( rfield.occ1nu < 0.0912 )
		{
			sprintf( chLine, 
				"   The continuum occupation number fell below 1 at %.3e microns.", 
			  0.0912/rfield.occ1nu );
			notein(chLine);
		}
		else if( rfield.occ1nu < 1. )
		{
			sprintf( chLine, 
				"   The continuum occupation number fell  below 1 at %.3e Angstroms.", 
			  912./rfield.occ1nu );
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The continuum occupation number fell  below 1 at %.3e Ryd.", 
			  rfield.occ1nu );
			notein(chLine);
		}
	}

	if( rfield.tbrmax > 1e3 )
	{
		sprintf( chLine, 
			"  !The largest continuum brightness temperature was %.3eK at %.3e Ryd.", 
		  rfield.tbrmax, rfield.tbrmnu );
		bangin(chLine);
	}

	if( rfield.tbrmax > 1e4 )
	{
		/* tbr4nu is energy (ryd) where continuum bright temp falls < 1e4 */
		if( rfield.tbr4nu < 0.0912 )
		{
			sprintf( chLine, 
				"   The continuum brightness temperature fell below 10,000K at %.3e microns.", 
			  0.0912/rfield.tbr4nu );
			notein(chLine);
		}
		else if( rfield.tbr4nu < 1. )
		{
			sprintf( chLine, 
				"   The continuum brightness temperature fell below 10,000K at %.3e Angstroms.", 
			  912./rfield.tbr4nu );
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The continuum brightness temperature fell below 10,000K at %.3e Ryd.", 
			  rfield.tbr4nu );
			notein(chLine);
		}
	}

	/* turbulence AND constant pressure do not make sense */
	if( DoppVel.TurbVel > 0. && strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		sprintf( chLine, 
			"  !Both constant pressure and turbulence makes no physical sense?" );
		bangin(chLine);
	}

	/* filling factor AND constant pressure do not make sense */
	if( geometry.FillFac < 1. && strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		sprintf( chLine, 
			"  !Both constant pressure and a filling factor makes no physical sense?" );
		bangin(chLine);
	}

	/* grains and solar abundances do not make sense */
	if( gv.lgDustOn() && abund.lgAbnSolar )
	{
		sprintf( chLine, 
			"  !Grains are present, but the gas phase abundances were left at the solar default.  This is not physical." );
		bangin(chLine);
	}

	/* check if depletion command set but no grains, another silly thing to do */
	if( abund.lgDepln && !gv.lgDustOn() )
	{
		sprintf( chLine, 
			"  !Grains are not present, but the gas phase abundances were depleted.  This is not physical." );
		bangin(chLine);
	}

	if( gv.lgDustOn() )
	{
		long nBin=0L, nFail=0L;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			if( gv.bin[nd]->QHeatFailures > 0L )
			{
				++nBin;
				nFail += gv.bin[nd]->QHeatFailures;
			}
		}
		if( nFail > 0 )
		{
			sprintf( chLine,
				 "  !The grain quantum heating treatment failed to converge %ld time(s) in %ld bin(s).", nFail, nBin );
			bangin(chLine);
		}
	}

#if 0
	/* check if PAHs were present in the ionized region */
	/* >>chng 05 jan 01, disabled this code now that PAH's have varying abundances by default, PvH */
	/** \todo	2	this statement needs to be reinstated with a better test for presence in the H II region */
	if( gv.lgDustOn() )
	{
		bool lgPAHsPresent_and_constant = false;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			lgPAHsPresent_and_constant = lgPAHsPresent_and_constant || 
				/* it is ok to have PAHs in the ionized region if the abundances vary */
				(gv.bin[nd]->lgPAHsInIonizedRegion /* && !gv.bin[nd]-> lgDustVary */);
		}
		if( lgPAHsPresent_and_constant )
		{
			sprintf( chLine,
				 " C-PAH's were present in the ionized region, this has never been observed in H II regions." );
			caunin(chLine);
		}
	}
#endif

	/* constant temperature greater than continuum energy density temperature */
	if( thermal.lgTemperatureConstant && thermal.ConstTemp*1.0001 < phycon.TEnerDen )
	{
		sprintf( chLine, 
			" C-The continuum energy density temperature (%g K)"
			" is greater than the gas kinetic temperature (%g K).",
			phycon.TEnerDen , thermal.ConstTemp);
		caunin(chLine);
		sprintf( chLine, " C-This is unphysical." );
		caunin(chLine);
	}

	/* remark that grains not present but energy density was low */
	if( !gv.lgDustOn() && phycon.TEnerDen < 800. )
	{
		sprintf( chLine, 
			"   Grains were not present but might survive in this environment (energy density temperature was %.2eK)", 
		  phycon.TEnerDen );
		notein(chLine);
	}

	/* call routine that will check age of cloud */
	AgeCheck();

	/* check on Ca H and H-epsilon overlapping
	 * need to do this since need to refer to lines arrays */
	chkCaHeps(&totwid);
	if( totwid > 121. )
	{
		sprintf( chLine, "   H-eps and Ca H overlap." );
		notein(chLine);
	}

	/* warning that something was turned off */
	if( !phycon.lgPhysOK )
	{
		sprintf( chLine, "  !A physical process has been disabled." );
		bangin(chLine);
	}

	/* check on lifetimes of [O III] against photoionization, only for low den */
	if( dense.gas_phase[ipHYDROGEN] < 1e8 )
	{
		if( oxy.r5007Max > 0.0263f )
		{
			sprintf( chLine, 
				"  !Photoionization of upper level of [O III] 5007 reached %.2e%% of the radiative lifetime.", 
			  oxy.r5007Max*100. );
			bangin(chLine);
		}
		else if( oxy.r5007Max > 0.0263f/10.f )
		{
			sprintf( chLine, 
				"   Photoionization of upper level of [O III] 5007 reached %.2e%% of the radiative lifetime.", 
			  oxy.r5007Max*100. );
			notein(chLine);
		}
		if( oxy.r4363Max > 1.78f )
		{
			sprintf( chLine, 
				"  !Photoionization of upper level of [O III] 4363 reached %.2e%% of the radiative lifetime.", 
			  oxy.r4363Max*100. );
			bangin(chLine);
		}
		else if( oxy.r4363Max > 1.78f/10.f )
		{
			sprintf( chLine, 
				"   Photoionization of upper level of [O III] 4363 reached %.2e%% of the radiative lifetime.", 
			  oxy.r4363Max*100. );
			notein(chLine);
		}
	}

	/* check whether total heating and cooling matched
	 * >>chng 97 mar 28, added GrossHeat, heat in terms normally heat-cool */
	error = fabs(thermal.power-thermal.totcol)/SDIV((thermal.power + thermal.totcol)/2.);
	if( thermal.lgTemperatureConstant )
	{
		if( error > 0.05 )
		{
			sprintf( chLine, 
				"  !Heating - cooling mismatch =%5.1f%%. Caused by constant temperature assumption. ", 
			  error*100. );
			bangin(chLine);
		}
	}

	else
	{
		if( error > 0.05 && error < 0.2 )
		{
			sprintf( chLine, " C-Heating - cooling mismatch =%.1f%%. What\'s wrong?", 
			  error*100. );
			caunin(chLine);
		}
		else if( error >= 0.2 )
		{
			sprintf( chLine, " W-Heating - cooling mismatch =%.2e%%. What\'s wrong????", 
			  error*100. );
			warnin(chLine);
		}
	}

	/* say if Ly-alpha photo of Ca+ excited levels was important */
	if( ca.Ca2RmLya > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of Ca+ 2D level by Ly-alpha reached %6.1f%% of the total rate out.", 
		  ca.Ca2RmLya*100. );
		notein(chLine);
	}

	/* check if Lya alpha ever hotter than gas */
	if( hydro.nLyaHot > 0 )
	{
		if( hydro.TLyaMax/hydro.TeLyaMax > 1.05 )
		{
			sprintf( chLine, 
				"  !The excitation temp of Lya exceeded the electron temp, largest value was %.2eK (gas temp there was %.2eK, zone%4ld)", 
			  hydro.TLyaMax, hydro.TeLyaMax, hydro.nZTLaMax );
			bangin(chLine);
		}
	}

	/* check if line absorption heating was important */

	/* get all negative lines, check if line absorption significant heat source
	 * this is used in "final" for energy budget print out */
	if( cdLine("Line",0,&SumNeg,&absint)<=0 )
	{
		fprintf( ioQQQ, " did not get sumneg cdLine\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* this is total heating */
	if( cdLine("TotH",0,&GetHeat,&absint)<=0 )
	{
		fprintf( ioQQQ, " did not get GetHeat cdLine\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( GetHeat > 0. )
	{
		SumNeg /= GetHeat;
		if( SumNeg > 0.1 )
		{
			sprintf( chLine, 
				"  !Line absorption heating reached %.2f%% of the global heating.", 
			  SumNeg*100. );
			bangin(chLine);
		}
		else if( SumNeg > 0.01 )
		{
			sprintf( chLine, 
				"   Line absorption heating reached %.2f%% of the global heating.", 
			  SumNeg*100. );
			notein(chLine);
		}
	}

	/* this is check of extra lines added with g-bar */
	if( input.lgSetNoBuffering )
	{
		sprintf( chLine, 
			"  !NO BUFFERING command was entered - this increases exec time by LARGE amounts.");
		bangin(chLine);
	}

	/* this is check of extra lines added with g-bar */
	if( thermal.GBarMax > 0.1 )
	{
		ASSERT( thermal.ipMaxExtra > 0 );
		strcpy( chLbl, chLineLbl(&TauLine2[thermal.ipMaxExtra-1]) );

		sprintf( chLine, 
			"  !G-bar cooling lines reached %.2f%% of the local cooling.  Line=%.10s", 
		  thermal.GBarMax*100., chLbl );
		bangin(chLine);
	}

	else if( thermal.GBarMax > 0.01 )
	{
		strcpy( chLbl, chLineLbl(&TauLine2[thermal.ipMaxExtra-1]) );

		sprintf( chLine, 
			"   G-bar cooling lines reached %.2f%% of the local cooling.  Line=%.10s", 
		  thermal.GBarMax*100., chLbl );
		notein(chLine);
	}

	/* this is check of hyperfine structure lines*/
	if( hyperfine.cooling_max > 0.1 )
	{
		sprintf( chLine, 
			"  !Hyperfine structure line cooling reached %.2f%% of the local cooling.", 
		  hyperfine.cooling_max*100.);
		bangin(chLine);
	}

	else if( hyperfine.cooling_max > 0.01 )
	{
		sprintf( chLine, 
			"   Hyperfine structure line cooling reached %.2f%% of the local cooling.", 
		  hyperfine.cooling_max*100. );
		notein(chLine);
	}

	/* line absorption heating reached more than 10% of local heating?
	 * HeatLineMax is largest heating(1,23)/htot */
	if( thermal.HeatLineMax > 0.1 )
	{
		if( thermal.levlmax == 1 )
		{
			/* main block of lines */
			/* >>chng 01 may 05, removed chGetLbl routine, which was here,
			 * replaced with chLineLbl routine and address of TauLines 
			 * should be no change in functionality */
			strcpy( chLbl, chLineLbl(&TauLines[thermal.ipHeatlmax] ) );
		}
		else if( thermal.levlmax == 2 )
		{
			/* level 2 lines */
			strcpy( chLbl, chLineLbl(&TauLine2[thermal.ipHeatlmax]) );
		}
		else if( thermal.levlmax == 3 )
		{
			/* hyperfine lines */
			strcpy( chLbl, chLineLbl(&HFLines[thermal.ipHeatlmax]) );
		}
		else if( thermal.levlmax == 4 )
		{
			/* 3rd-party database lines */
			strcpy( chLbl, chLineLbl(dBaseLines[thermal.ipHeatlmax].tran) );
		}
		else
		{
			fprintf( ioQQQ, " PROBLEM DISASTER PrtComment has insane levlmax,=%5ld\n", 
			  thermal.levlmax );
		}
		sprintf( chLine, 
			"  !Line absorption heating reached %.2f%% of the local heating - largest by level%2ld line %.10s", 
		  thermal.HeatLineMax*100., thermal.levlmax, chLbl );
		bangin(chLine);
	}

	else if( thermal.HeatLineMax > 0.01 )
	{
		sprintf( chLine, 
			"   Line absorption heating reached %.2f%% of the local heating.", 
		  thermal.HeatLineMax*100. );
		notein(chLine);
	}

	if( ionbal.CompHeating_Max > 0.05 )
	{
		sprintf( chLine, 
			"  !Bound Compton heating reached %.2f%% of the local heating.", 
		  ionbal.CompHeating_Max*100. );
		bangin(chLine);
	}
	else if( ionbal.CompHeating_Max > 0.01 )
	{
		sprintf( chLine, 
			"   Bound Compton heating reached %.2f%% of the local heating.", 
		  ionbal.CompHeating_Max*100. );
		notein(chLine);
	}

	/* check whether any lines in the iso sequences mased */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				long int nmax = iso.numLevels_local[ipISO][nelem];

				/* minus one here is to exclude highest level */
				for( ipHi=1; ipHi < nmax - 1; ++ipHi )
				{
					for( ipLo=0; ipLo < ipHi; ++ipLo )
					{
						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
							continue;

						/* did the line mase */
						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn < -0.1 )
						{
							sprintf( chLine, 
								"  !Some iso-structure lines mased: %s-like %s, line %li-%li had optical depth %.2e", 
								elementnames.chElementSym[ipISO],
								elementnames.chElementNameShort[nelem],
								ipHi , ipLo ,
								Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn );
							bangin(chLine);
						}
					}
				}
			}
		}
	}

	if( dense.gas_phase[ipHYDROGEN] < 1e7 )
	{
		/* check on IR fine structure lines - not necessary if dense since will be in LTE */
		lgThick = false;
		tauneg = 0.;
		alpha = 0.;
		/* loop from 3, since 0 is dummy, 1 and 2 are spin-flip transitions of H and He */
		for( i=3; i <= nLevel1; i++ )
		{
			/* define IR as anything longward of 1 micron */
			if( TauLines[i].EnergyWN < 10000. )
			{
				if( TauLines[i].Emis->TauIn > 1. )
				{
					lgThick = true;
					alpha = MAX2(alpha,(double)TauLines[i].Emis->TauIn);
				}
				else if( TauLines[i].Emis->TauIn < (realnum)tauneg )
				{
					tauneg = TauLines[i].Emis->TauIn;
					strcpy( chLbl, chLineLbl(&TauLines[i]) );
				}
			}
		}
		/* now print results, were any fine structure lines optically thick? */
		if( lgThick )
		{
			sprintf( chLine, 
				"  !Some infrared fine structure lines are optically thick:  largest tau was %.2e", 
			  alpha );
			bangin(chLine);
		}
		/* did any fine structure lines mase? */
		if( tauneg < -0.01 )
		{
			sprintf( chLine, 
				"  !Some fine structure lines mased: line %s had optical depth %.2e", 
			  chLbl, tauneg );
			bangin(chLine);
		}
	}

	/* were any other lines masing? */
	tauneg = 0.;
	alpha = 0.;
	for( i=1; i <= nLevel1; i++ )
	{
		/* define UV as anything shortward of 1 micron */
		if( TauLines[i].EnergyWN >= 10000. )
		{
			if( TauLines[i].Emis->TauIn < (realnum)tauneg )
			{
				tauneg = TauLines[i].Emis->TauIn;
				strcpy( chLbl, chLineLbl(&TauLines[i]) );
			}
		}
	}

	/* did any level1 lines mase? */
	if( tauneg < -0.01 )
	{
		sprintf( chLine, 
			"  !Some level1 lines mased: most negative ion and tau were: %s %.2e", 
		  chLbl, tauneg );
		bangin(chLine);
	}

	/* this is check that at least a second iteration was done with sphere static,
	 * the test is overridden with the (OK) option on the sphere static command,
	 * which sets geometry.lgStaticNoIt true */
	if( geometry.lgStatic && iterations.lgLastIt && (iteration == 1) && 
		!geometry.lgStaticNoIt)
	{
		sprintf( chLine, " C-I must iterate when SPHERE STATIC is set." );
		caunin(chLine);
		iterations.lgIterAgain = true;
	}

	/* caution if continuum is punched but only one iteration performed */
	if( save.lgPunContinuum && iteration == 1 && iterations.lgLastIt)
	{
		sprintf( chLine, " C-I must iterate when save continuum output is done." );
		caunin(chLine);
		iterations.lgIterAgain = true;
	}

	/** \todo	2	extend to all iso and elem */
	/* how important was induced two photon?? */
	if( iso.TwoNu_induc_dn_max[ipH_LIKE][ipHYDROGEN] > 1. )
	{
		sprintf( chLine, "  !Rate of induced H 2-photon emission reached %.2e s^-1", 
		  iso.TwoNu_induc_dn_max[ipH_LIKE][ipHYDROGEN] );
		bangin(chLine);
	}

	else if( iso.TwoNu_induc_dn_max[ipH_LIKE][ipHYDROGEN] > 0.01 )
	{
		sprintf( chLine, "   Rate of induced H 2-photon emission reached %.2e s^-1", 
		  iso.TwoNu_induc_dn_max[ipH_LIKE][ipHYDROGEN] );
		notein(chLine);
	}

	/* how important was induced recombination? */
	if( hydro.FracInd > 0.01 )
	{
		sprintf( chLine, 
			"   Induced recombination was %5.1f%% of the total for H level%3ld", 
		  hydro.FracInd*100., hydro.ndclev );
		notein(chLine);
	}

	if( hydro.fbul > 0.01 )
	{
		sprintf( chLine, 
			"   Stimulated emission was%6.1f%% of the total for H transition%3ld -%3ld", 
		  hydro.fbul*100., hydro.nbul + 1, hydro.nbul );
		notein(chLine);
	}

	/* check whether Fe II destruction of La was important - entry into lines stack 
	 * is in prt_lines_hydro.c */
	if( cdLine("Fe 2",1216,&fedest,&absint)<=0 )
	{
		fprintf( ioQQQ, " Did not find Fe II Lya\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* find total Lya for comparison */
	if( cdLine("TOTL",1215.68f,&relhm,&absint)<=0 )
	{
		fprintf( ioQQQ, " Did not find Lya\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( relhm > 0. )
	{
		ratio = fedest/(fedest + relhm);
		if( ratio > 0.1 )
		{
			sprintf( chLine, "  !Fe II destruction of Ly-a removed %.1f%% of the line.", 
			  ratio *100.);
			bangin(chLine);
		}
		else if( ratio > 0.01 )
		{
			sprintf( chLine, "   Fe II destruction of Ly-a removed %.1f%% of the line.", 
			  ratio );
			notein(chLine);
		}
	}

	if( cdLine("H-CT",6563,&relhm,&absint)<=0 )
	{
		fprintf( ioQQQ, " Comment did not find H-CT H-alpha\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( HBeta > 0. )
	{
		if( relhm/HBeta > 0.01 )
		{
			sprintf( chLine, 
				"  !Mutual neutralization production of H-alpha was significant." );
			bangin(chLine);
		}
	}

	/* note about very high population in H n=2 rel to ground, set in hydrogenic */
	if( hydro.lgHiPop2 )
	{
		sprintf( chLine, 
			"   The population of H n=2 reached %.2e relative to the ground state.", 
		  hydro.pop2mx );
		notein(chLine);
	}

	/* check where diffuse emission error */
	for( ipISO=ipH_LIKE; ipISO<=ipHE_LIKE; ++ipISO )
	{
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			if( iso.CaseBCheck[ipISO][nelem] > 1.5 )
			{
				sprintf( chLine, 
					"   Ratio of computed diffuse emission to case B reached %g for iso %li element %li",
					iso.CaseBCheck[ipISO][nelem] , ipISO , nelem+1 );
				notein(chLine);
			}
		}
	}

	/* check whether electrons were relativistic */
	if( thermal.thist > 1e9 )
	{
		/* >>chng 06 feb 19, from 5e9 K for warning to 1e10K.  add test case at 1e10K
		 * and don't want warning in test suite.  nothing is wrong at this temp - eeff
		 * is in correctly for relativistic temps and will eventually dominate cooling */
		if( thermal.thist > 1.0001e10 )
		{
			sprintf( chLine, " W-Electrons were relativistic; High TE=%.2e", 
			  thermal.thist );
			warnin(chLine);
		}
		else
		{
			sprintf( chLine, " C-Electrons were mildly relativistic; High TE=%.2e", 
			  thermal.thist );
			caunin(chLine);
		}
	}

	/* check on timescale for photoerosion of elements */
	rate = timesc.TimeErode*2e-26;
	if( rate > 1e-35 )
	{
		/*  2E-26 is roughly cross section for photoerosion
		 *  see 
		 * >>refer	all	photoerode	Boyd, R., & Ferland, G.J. ApJ, 318, L21. */
		ts = (1./rate)/3e7;
		if( ts < 1e3 )
		{
			sprintf( chLine, "  !Timescale-photoerosion of Fe=%.2e yr", 
			  ts );
			bangin(chLine);
		}
		else if( ts < 1e9 )
		{
			sprintf( chLine, "   Timescale-photoerosion of Fe=%.2e yr", 
			  ts );
			notein(chLine);
		}
	}

	/* check whether Compton heating was significant */
	comfrc = rfield.comtot/SDIV(thermal.power);
	if( comfrc > 0.01 )
	{
		sprintf( chLine, "   Compton heating was %5.1f%% of the total.", 
		  comfrc*100. );
		notein(chLine);
	}

	/* check on relative importance of induced Compton heating */
	if( comfrc > 0.01 && rfield.cinrat > 0.05 )
	{
		sprintf( chLine, 
			"  !Induced Compton heating was %.2e of the total Compton heating.", 
		  rfield.cinrat );
		bangin(chLine);
	}

	/* check whether equilibrium timescales are short rel to Hubble time */
	if( timesc.tcmptn > 5e17 )
	{
		if( comfrc > 0.05 )
		{
			sprintf( chLine, 
				" C-Compton cooling is significant and the equilibrium timescale (%.2e s) is longer than the Hubble time.", 
			  timesc.tcmptn );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   Compton cooling equilibrium timescale (%.2e s) is longer than Hubble time.", 
			  timesc.tcmptn );
			notein(chLine);
		}
	}

	if( timesc.time_therm_long > 5e17 )
	{
		sprintf( chLine, 
			" C-Thermal equilibrium timescale, %.2e s, longer than Hubble time; this cloud is not time-steady.", 
		  timesc.time_therm_long );
		caunin(chLine);
	}

	/* check whether model large relative to Jeans length
	 * DMEAN is mean density (gm per cc)
	 * mean temp is weighted by mass density */
	if( log10(radius.depth) > colden.rjnmin )
	{
		/* AJMIN is minimum Jeans mass, log in grams */
		aj = pow(10.,colden.ajmmin - log10(SOLAR_MASS));
		if( strcmp(dense.chDenseLaw,"CPRE") == 0 )
		{
			sprintf( chLine, 
				" C-Cloud thicker than smallest Jeans length=%8.2ecm; stability problems? (smallest Jeans mass=%8.2eMo)", 
			  pow((realnum)10.f,colden.rjnmin), aj );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   Cloud thicker than smallest Jeans length=%8.2ecm; stability problems? (smallest Jeans mass=%8.2eMo)", 
			  pow((realnum)10.f,colden.rjnmin), aj );
			notein(chLine);
		}
	}

	/* check whether grains too hot to survive */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( gv.bin[nd]->TeGrainMax > gv.bin[nd]->Tsublimat )
		{
			sprintf( chLine, 
				" W-Maximum temperature of grain%-12.12s was %.2eK, above its sublimation temperature, %.2eK.", 
			  gv.bin[nd]->chDstLab, gv.bin[nd]->TeGrainMax, 
			  gv.bin[nd]->Tsublimat );
			warnin(chLine);
		}
		else if( gv.bin[nd]->TeGrainMax > gv.bin[nd]->Tsublimat* 0.9 )
		{
			sprintf( chLine, 
				" C-Maximum temperature of grain%-12.12s was %.2eK, near its sublimation temperature, %.2eK.", 
			  gv.bin[nd]->chDstLab, gv.bin[nd]->TeGrainMax, 
			  gv.bin[nd]->Tsublimat );
			caunin(chLine);
		}
	}

	if( gv.lgNegGrnDrg )
	{
		sprintf( chLine, "  !Grain drag force <0." );
		bangin(chLine);
	}

	/* largest relative number of electrons donated by grains */
	if( gv.GrnElecDonateMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Grains donated %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecDonateMax*100. );
		bangin(chLine);
	}
	else if( gv.GrnElecDonateMax > 0.005 )
	{
		sprintf( chLine, 
			"   Grains donated %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecDonateMax*100. );
		notein(chLine);
	}

	/* largest relative number of electrons on grain surface */
	if( gv.GrnElecHoldMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Grains contained %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecHoldMax*100. );
		bangin(chLine);
	}
	else if( gv.GrnElecHoldMax > 0.005 )
	{
		sprintf( chLine, 
			"   Grains contained %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecHoldMax*100. );
		notein(chLine);
	}

	/* is photoelectric heating of gas by photoionization of grains important */
	if( gv.dphmax > 0.5 )
	{
		sprintf( chLine, 
			"  !Local grain-gas photoelectric heating rate reached %5.1f%% of the total.", 
		  gv.dphmax*100. );
		bangin(chLine);
	}
	else if( gv.dphmax > 0.05 )
	{
		sprintf( chLine, 
			"   Local grain-gas photoelectric heating rate reached %5.1f%% of the total.", 
		  gv.dphmax*100. );
		notein(chLine);
	}

	if( gv.TotalDustHeat/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, 
			"   Global grain photoelectric heating of gas was%5.1f%% of the total.", 
		  gv.TotalDustHeat/thermal.power*100. );
		notein(chLine);
		if( gv.TotalDustHeat/thermal.power > 0.25 )
		{
			sprintf( chLine, 
				"  !Grain photoelectric heating is VERY important." );
			bangin(chLine);
		}
	}

	/* grain-gas collisional cooling of gas */
	if( gv.dclmax > 0.05 )
	{
		sprintf( chLine, 
			"   Local grain-gas cooling of gas rate reached %5.1f%% of the total.", 
		  gv.dclmax*100. );
		notein(chLine);
	}

	/* check how H2 chemistry network performed */
	if( h2.renorm_max > 1.05 )
	{
		if( h2.renorm_max > 1.2 )
		{
			sprintf( chLine, 
				"  !The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_max);
			bangin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_max);
			notein(chLine);
		}
	}
	if( h2.renorm_min < 0.95 )
	{
		if( h2.renorm_min < 0.8 )
		{
			sprintf( chLine, 
				"  !The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_min);
			bangin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_min);
			notein(chLine);
		}
	}

	/* check whether photodissociation of H_2^+ molecular ion was important */
	if( hmi.h2pmax > 0.10 )
	{
		sprintf( chLine, 
			"  !The local H2+ photodissociation heating rate reached %5.1f%% of the total heating.", 
		  hmi.h2pmax*100. );
		bangin(chLine);
	}

	else if( hmi.h2pmax > 0.01 )
	{
		sprintf( chLine, 
			"   The local H2+ photodissociation heating rate reached %.1f%% of the total heating.", 
		  hmi.h2pmax*100. );
		notein(chLine);
	}

	/* check whether photodissociation of molecular hydrogen (H2)was important */
	if( hmi.h2dfrc > 0.1 )
	{
		sprintf( chLine, 
			"  !The local H2 photodissociation heating rate reached %.1f%% of the total heating.", 
		  hmi.h2dfrc*100. );
		bangin(chLine);
	}
	else if( hmi.h2dfrc > 0.01 )
	{
		sprintf( chLine, 
			"   The local H2 photodissociation heating rate reached %.1f%% of the total heating.", 
		  hmi.h2dfrc*100. );
		notein(chLine);
	}

	/* check whether cooling by molecular hydrogen (H2) was important */
	if( hmi.h2line_cool_frac > 0.1 )
	{
		sprintf( chLine, 
			"  !The local H2 cooling rate reached %.1f%% of the local cooling.", 
		  hmi.h2line_cool_frac*100. );
		bangin(chLine);
	}
	else if( hmi.h2line_cool_frac > 0.01 )
	{
		sprintf( chLine, 
			"   The local H2 cooling rate reached %.1f%% of the local cooling.", 
		  hmi.h2line_cool_frac*100. );
		notein(chLine);
	}

	if( hmi.h2dtot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, 
			"   Global H2 photodissociation heating of gas was %.1f%% of the total heating.", 
		  hmi.h2dtot/thermal.power*100. );
		notein(chLine);
		if( hmi.h2dtot/thermal.power > 0.25 )
		{
			sprintf( chLine, "   H2 photodissociation heating is VERY important." );
			notein(chLine);
		}
	}

	/* check whether photodissociation of carbon monoxide (co) was important */
	if( co.codfrc > 0.25 )
	{
		sprintf( chLine, 
			"  !Local CO photodissociation heating rate reached %.1f%% of the total.", 
		  co.codfrc*100. );
		bangin(chLine);
	}
	else if( co.codfrc > 0.05 )
	{
		sprintf( chLine, 
			"   Local CO photodissociation heating rate reached %.1f%% of the total.", 
		  co.codfrc*100. );
		notein(chLine);
	}

	if( co.codtot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, 
			"   Global CO photodissociation heating of gas was %.1f%% of the total.", 
		  co.codtot/thermal.power*100. );
		notein(chLine);
		if( co.codtot/thermal.power > 0.25 )
		{
			sprintf( chLine, "   CO photodissociation heating is VERY important." );
			notein(chLine);
		}
	}

	if( thermal.lgEdnGTcm )
	{
		sprintf( chLine, 
			"   Energy density of radiation field was greater than the Compton temperature. Is this physical?" );
		notein(chLine);
	}

	/* was cooling due to induced recombination important? */
	if( hydro.cintot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, "   Induced recombination cooling was %.1f%% of the total.", 
		  hydro.cintot/thermal.power*100. );
		notein(chLine);
	}

	/* check whether free-free heating was significant */
	if( thermal.FreeFreeTotHeat/SDIV(thermal.power) > 0.1 )
	{
		sprintf( chLine, "  !Free-free heating was %.1f%% of the total.", 
		  thermal.FreeFreeTotHeat/thermal.power*100. );
		bangin(chLine);
	}
	else if( thermal.FreeFreeTotHeat/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, "   Free-free heating was %.1f%% of the total.", 
		  thermal.FreeFreeTotHeat/thermal.power*100. );
		notein(chLine);
	}

	/* was heating due to H- absorption important? */
	if( hmi.hmitot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, "   H- absorption heating was %.1f%% of the total.", 
		  hmi.hmitot/SDIV(thermal.power)*100. );
		notein(chLine);
	}

	/* water destruction rate was zero */
	if( co.lgH2Ozer )
	{
		sprintf( chLine, "   Water destruction rate zero." );
		notein(chLine);
	}

	/* numerical instability in matrix inversion routine */
	if( atoms.nNegOI > 0 )
	{
		sprintf( chLine, " C-O I negative level populations %ld times.", 
		  atoms.nNegOI );
		caunin(chLine);
	}

	/* check for negative optical depths,
	 * optical depth in excited state helium lines */
	small = 0.;
	imas = 0;
	isav = 0;
	j = 0;
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* >>chng 06 aug 28, from numLevels_max to _local. */
			for( ipLo=ipH2p; ipLo < (iso.numLevels_local[ipH_LIKE][nelem] - 1); ipLo++ )
			{
				/* >>chng 06 aug 28, from numLevels_max to _local. */
				for( ipHi=ipLo + 1; ipHi < iso.numLevels_local[ipH_LIKE][nelem]; ipHi++ )
				{
					if( Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
						continue;

					if( Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->TauIn < (realnum)small )
					{
						small = Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->TauIn;
						imas = ipHi;
						j = ipLo;
						isav = nelem;
					}
				}
			}
		}
	}

	if( small < -0.05 )
	{
		sprintf( chLine, 
			"  !Some hydrogenic lines mased, species was %2s%2ld, smallest tau was %.2e, transition %li-%li", 
			elementnames.chElementSym[isav], 
			isav+1,small, imas , j );
		bangin(chLine);
	}

	/* check for negative opacities */
	if( opac.lgOpacNeg )
	{
		sprintf( chLine, "  !Some opacities were negative - the SET NEGOPC command will save which ones." );
		bangin(chLine);
	}

	/* now check continua */
	small = 0.;
	imas = 0;
	isav = 0;
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* >>chng 06 aug 28, from numLevels_max to _local. */
			for( i=0; i < iso.numLevels_local[ipH_LIKE][nelem]; i++ )
			{
				if( opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][nelem][i]-1] < -0.001 )
				{
					small = MIN2(small,(double)opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][nelem][i]-1]);
					imas = i;
					isav = nelem;
				}
			}
		}
	}

	if( small < -0.05 )
	{
		sprintf( chLine, "  !Some hydrogenic (%2s%2ld) continua optical depths were negative; smallest=%.2e level=%3ld", 
			elementnames.chElementSym[isav], 
			isav+1,
		  small, imas );
		bangin(chLine);
	}

	/* check whether any continuum optical depths are negative */
	nneg = 0;
	tauneg = 0.;
	freqn = 0.;
	for( i=0; i < rfield.nflux; i++ )
	{
		if( opac.TauAbsGeo[0][i] < -0.001 )
		{
			nneg += 1;
			/* only remember the smallest freq, and most neg optical depth */
			if( nneg == 1 )
				freqn = rfield.anu[i];
			tauneg = MIN2(tauneg,(double)opac.TauAbsGeo[0][i]);
		}
	}

	if( nneg > 0 )
	{
		sprintf( chLine, "  !Some continuous optical depths <0.  The lowest freq was %.3e Ryd, and a total of%4ld", 
		  freqn, nneg );
		bangin(chLine);
		sprintf( chLine, "  !The smallest optical depth was %.2e", 
		  tauneg );
		bangin(chLine);
	}

	/* say if Balmer continuum optically thick */
	if( opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]-1] > 0.05 )
	{
		sprintf( chLine, "   The Balmer continuum optical depth was %.2e.", 
		  opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]-1] );
		notein(chLine);
	}

	/* was correction for stimulated emission significant? */
	if( opac.stimax[0] > 0.02 && opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1] > 0.2 )
	{
		sprintf( chLine, "   The Lyman continuum stimulated emission correction to optical depths reached %.2e.", 
		  opac.stimax[0] );
		notein(chLine);
	}
	else if( opac.stimax[1] > 0.02 && opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2]-1] > 0.1 )
	{
		sprintf( chLine, "   The Balmer continuum stimulated emission correction to optical depths reached %.2e.", 
		  opac.stimax[1] );
		notein(chLine);
	}

	/* say if Paschen continuum optically thick */
	if( opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][3]-1] > 0.2 )
	{
		sprintf( chLine, 
			"   The Paschen continuum optical depth was %.2e.", 
		  opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][3]-1] );
		notein(chLine);
	}

	/* some comments about near IR total optical depth */
	if( opac.TauAbsGeo[0][0] > 1. )
	{
		sprintf( chLine, 
			"   The continuum optical depth at the lowest energy considered (%.3e Ryd) was %.3e.", 
		  rfield.anu[0], opac.TauAbsGeo[0][0] );
		notein(chLine);
	}

	/* comment if optical depth to Rayleigh scattering is big
	 * cs from VAL 76 */
	if( colden.colden[ipCOL_H0]*7e-24 > 0.01 )
	{
		sprintf( chLine, 
			"   The optical depth to Rayleigh scattering at 1300A is %.2e", 
		  colden.colden[ipCOL_H0]*6.71e-24 );
		notein(chLine);
	}

	if( colden.colden[ipCOL_H2p]*7e-18 > 0.1 )
	{
		sprintf( chLine, 
			"  !The optical depth to the H2+ molecular ion is %.2e", 
		  colden.colden[ipCOL_H2p]*7e-18 );
		bangin(chLine);
	}
	else if( colden.colden[ipCOL_H2p]*7e-18 > 0.01 )
	{
		sprintf( chLine, 
			"   The optical depth to the H2+ molecular ion is %.2e", 
		  colden.colden[ipCOL_H2p]*7e-18 );
		notein(chLine);
	}

	/* warn if optically thick to H- absorption */
	if( opac.thmin > 0.1 )
	{
		sprintf( chLine, 
			"  !Optical depth to negative hydrogen ion is %.2e", 
		  opac.thmin );
		bangin(chLine);
	}
	else if( opac.thmin > 0.01 )
	{
		sprintf( chLine, 
			"   Optical depth to negative hydrogen ion is %.2e", 
		  opac.thmin );
		notein(chLine);
	}

	/* say if 3-body recombination coefficient function outside range of validity
	 * tripped if te/z**2 < 100 or approx 10**13: => effect >50% of radiative
	 * other integers defined in source for da */
	if( ionbal.ifail > 0 && ionbal.ifail <= 10 )
	{
		sprintf( chLine, 
			"   3 body recombination coefficient outside range %ld", ionbal.ifail );
		notein(chLine);
	}
	else if( ionbal.ifail > 10 )
	{
		sprintf( chLine, 
			" C-3 body recombination coefficient outside range %ld", ionbal.ifail );
		caunin(chLine);
	}

	/* check whether energy density less than background */
	if( phycon.TEnerDen < 2.6 )
	{
		sprintf( chLine, 
			"  !Incident radiation field energy density is less than 2.7K.  Add background with CMB command." );
		bangin(chLine);
	}

	/* check whether CMB set at all */
	if( !rfield.lgCMB_set )
	{
		sprintf( chLine, 
			"  !The CMB was not included.  This is added with the CMB command." );
		bangin(chLine);
	}

	/* incident radiation field is less than background Habing ISM field */
	if( rfield.lgHabing )
	{
		sprintf( chLine, 
			"  !The intensity of the incident radiation field is less than 10 times the Habing diffuse ISM field.  Is this OK?" );
		bangin(chLine);
		sprintf( chLine, 
			"  !   Consider adding diffuse ISM emission with TABLE ISM command." );
		bangin(chLine);
	}

	/* some things dealing with molecules, or molecule formation */

	/* if C/O > 1 then chemistry will be carbon dominated rather than oxygen dominated */
	if( dense.lgElmtOn[ipOXYGEN] && dense.lgElmtOn[ipCARBON] )
	{
		if( dense.gas_phase[ipCARBON]/dense.gas_phase[ipOXYGEN] > 1. )
		{
			sprintf( chLine, "  !The C/O abundance ratio, %.1f, is greater than unity.  The chemistry will be carbon dominated.", 
				dense.gas_phase[ipCARBON]/dense.gas_phase[ipOXYGEN] );
			bangin(chLine);
		}
	}

       // Remove this bit as redundant with more generalized below?
#if    1
	/* more than 10% of H is in the H2 molecule */
	if( hmi.BiggestH2 > 0.1 )
	{
		sprintf( chLine, "  !The fraction of %s in %s reached %.1f%% at some point in the cloud.", 
			"H ",
			"H2   ",
			hmi.BiggestH2*100. );
		bangin(chLine);
	}
	else if( hmi.BiggestH2>0.01 )
	{
		sprintf( chLine, "   The fraction of %s in %s reached %.2f%% at some point in the cloud.", 
			"H ",
			"H2   ",
			hmi.BiggestH2*100. );
		notein(chLine);
	}
	else if( hmi.BiggestH2 > 1e-3 )
	{
		sprintf( chLine, "   The fraction of %s in %s reached %.3f%% at some point in the cloud.", 
			"H ",
			"H2   ",
			hmi.BiggestH2*100. );
		notein(chLine);
	}
#endif

	lgLots_of_moles = false;
	lgLotsSolids = false;
	/* largest fraction in any heavy element molecule */
	for( i=0; i<mole.num_comole_calc; ++i )
	{
		if( COmole[i]->n_nuclei <= 1)
			continue;

		if( COmole[i]->xMoleFracMax > 0.1 )
		{
			sprintf( chLine, "  !The fraction of %s in %s reached %.1f%% at some point in the cloud.", 
				elementnames.chElementSym[COmole[i]->nelem_hevmol],
				COmole[i]->label,
				COmole[i]->xMoleFracMax*100. );
			bangin(chLine);
			lgLots_of_moles = true;
			/* check whether molecules are on grains */
			if( !COmole[i]->lgGas_Phase )
				lgLotsSolids = true;
		}
		else if( COmole[i]->xMoleFracMax>0.01 )
		{
			sprintf( chLine, "   The fraction of %s in %s reached %.2f%% at some point in the cloud.", 
				elementnames.chElementSym[COmole[i]->nelem_hevmol],
				COmole[i]->label,
				COmole[i]->xMoleFracMax*100. );
			notein(chLine);
			lgLots_of_moles = true;
			/* check whether molecules are on grains */
			if( !COmole[i]->lgGas_Phase )
				lgLotsSolids = true;
		}
		else if( COmole[i]->xMoleFracMax > 1e-3 )
		{
			sprintf( chLine, "   The fraction of %s in %s reached %.3f%% at some point in the cloud.", 
				elementnames.chElementSym[COmole[i]->nelem_hevmol],
				COmole[i]->label,
				COmole[i]->xMoleFracMax*100. );
			notein(chLine);
			/* check whether molecules are on grains */
			if( !COmole[i]->lgGas_Phase )
				lgLotsSolids = true;
		}
	}

	/* generate comment if molecular fraction was significant but some heavy elements are turned off */
	if( lgLots_of_moles )
	{
		/* find all elements that are turned off */
		for(nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			/* >>chng 05 dec 23, add mole.lgElem_in_chemistry */
			if( !dense.lgElmtOn[nelem]  && mole.lgElem_in_chemistry[nelem] )
			{
				/* this triggers if element turned off but it is part of co chem net */
				sprintf( chLine, 
					" C-Molecules are important, but %s, part of the chemistry network, is turned off.", 
					elementnames.chElementName[nelem] );
				caunin(chLine);
			}
#			if 0
				/* this element has been turned off - now check if part of chemistry */
				for( i=NUM_HEAVY_MOLEC+NUM_ELEMENTS; i<NUM_COMOLE_CALC; ++i )
				{
					if( nelem==COmole[i].nelem_hevmol )
					{
						/* this triggers if element turned off but it is part of co chem net */
						sprintf( chLine, 
							" C-Molecules are important, but %s, part of the chemistry network, is turned off.", 
							elementnames.chElementName[nelem] );
						caunin(chLine);
					}
				}
			}
#			endif
		}
	}

	/* say if lots of molecules on grains,
	 * molecules with labels that *GR */
	if( lgLotsSolids ) 
	{
		sprintf( chLine, "  !A significant amount of molecules condensed onto grain surfaces." );
		bangin(chLine);
		sprintf( chLine, "  !These are the molecular species with \"grn\" above." );
		bangin(chLine);
	}

	/* bremsstrahlung optical depth */
	if( rfield.EnergyBremsThin > 0.09 )
	{
		sprintf( chLine, "  !The cloud is optically thick at optical wavelengths, extending to %.3e Ryd =%.3eA", 
		  rfield.EnergyBremsThin, RYDLAM/rfield.EnergyBremsThin );
		bangin(chLine);
	}
	else if( rfield.EnergyBremsThin > 0.009 )
	{
		sprintf( chLine, "   The continuum of the computed structure may be optically thick in the near infrared." );
		notein(chLine);
	}

	/* did model run away to very large radius? */
	if( radius.Radius > 1e23 && radius.Radius/radius.rinner > 10. )
	{
		sprintf( chLine, "   Is an outer radius of %.2e reasonable?", 
		  radius.Radius );
		notein(chLine);
	}

	/* following set true in RT_line_one_tauinc if maser capped at tau = -1 */
	if( rt.lgMaserCapHit )
	{
		sprintf( chLine, "   Laser maser optical depths capped in RT_line_one_tauinc." );
		notein(chLine);
	}

	/* following set true in adius_next if maser cap set dr */
	if( rt.lgMaserSetDR )
	{
		sprintf( chLine, "  !Line maser set zone thickness in some zones." );
		bangin(chLine);
	}

	/* lgPradCap is true if radiation pressure was capped on first iteration
	 * also check that this is a constant total pressure model */
	if( (pressure.lgPradCap && (strcmp(dense.chDenseLaw,"CPRE") == 0)) && 
	  pressure.lgPres_radiation_ON )
	{
		sprintf( chLine, "   Radiation pressure kept below gas pressure on this iteration." );
		notein(chLine);
	}

	if( pressure.RadBetaMax > 0.25 )
	{
		if( pressure.ipPradMax_line == 0 )
		{
			sprintf( chLine, 
				"  !The ratio of radiation to gas pressure reached %.2e at zone %li.  Caused by Lyman alpha.", 
			  pressure.RadBetaMax,
			  pressure.ipPradMax_nzone);
			bangin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"  !The ratio of radiation to gas pressure reached %.2e at zone %li.  "
				"Caused by line number %ld, label %s", 
			  pressure.RadBetaMax, 
			  pressure.ipPradMax_nzone,
			  pressure.ipPradMax_line,
			  pressure.chLineRadPres );
			bangin(chLine);
		}
	}

	else if( pressure.RadBetaMax > 0.025 )
	{
		if( pressure.ipPradMax_line == 0 )
		{
			sprintf( chLine, 
				"   The ratio of radiation to gas pressure reached %.2e at zone %li.  Caused by Lyman alpha.", 
			  pressure.RadBetaMax,
			  pressure.ipPradMax_nzone);
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The ratio of radiation to gas pressure reached %.2e at zone %li.  "
				"Caused by line number %ld, label %s", 
			  pressure.RadBetaMax, 
			  pressure.ipPradMax_nzone,
			  pressure.ipPradMax_line,
			  pressure.chLineRadPres );
			notein(chLine);
		}
	}

	if( opac.telec >= 5. )
	{
		sprintf( chLine, " W-The model is optically thick to electron "
			"scattering; tau=%.2e  Cloudy is NOT intended for this regime.", 
		  opac.telec );
		warnin(chLine);
	}
	else if( opac.telec > 2.0  )
	{
		sprintf( chLine, " C-The model is moderately optically thick to electron scattering; tau=%.1f", 
		  opac.telec );
		caunin(chLine);
	}
	else if( opac.telec > 0.1  )
	{
		sprintf( chLine, "  !The model has modest optical depth to electron scattering; tau=%.2f", 
		  opac.telec );
		bangin(chLine);
	}
	else if( opac.telec > 0.01 )
	{
		sprintf( chLine, "   The optical depth to electron scattering is %.3f", 
		  opac.telec );
		notein(chLine);
	}

	/* optical depth to 21 cm */
	if( HFLines[0].Emis->TauIn > 0.5 )
	{
		sprintf( chLine, "  !The optical depth in the H I 21 cm line is %.2e",HFLines[0].Emis->TauIn );
		bangin(chLine);
	}

	/* comment if level2 lines are off - they are used to pump excited states
	 * of ground term by UV light */
	if( nWindLine==0 )
	{
		/* generate comment */
		sprintf( chLine, "  !The level2 lines are disabled.  UV pumping of excited levels within ground terms is not treated." );
		bangin(chLine);
	}

	/* check on optical depth convergence of all hydrogenic lines */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			if( Transitions[ipH_LIKE][nelem][ipH3p][ipH2s].Emis->TauIn > 0.2 )
			{
				differ = fabs(1.-Transitions[ipH_LIKE][nelem][ipH3p][ipH2s].Emis->TauIn*
				  rt.DoubleTau/Transitions[ipH_LIKE][nelem][ipH3p][ipH2s].Emis->TauTot)*100.;

				/* check whether H-alpha optical depth changed by much on last iteration
				 * no tolerance can be finer than autocv, the tolerance on the
				 * iterate to convergence command.  It is 15% */
				if( ((iterations.lgLastIt && Transitions[ipH_LIKE][nelem][ipH3p][ipH2s].Emis->TauIn > 0.8) && 
					differ > 20.) && wind.lgStatic() )
				{
					sprintf( chLine, 
						" C-This is the last iteration and %2s%2ld Bal(a) optical depth"
						" changed by%6.1f%% (was %.2e). Try another iteration.", 
					  elementnames.chElementSym[nelem], 
					  nelem+1, differ, 
					  Transitions[ipH_LIKE][nelem][ipH3p][ipH2s].Emis->TauTot );
					caunin(chLine);
					iterations.lgIterAgain = true;
				}

				/* only check on Lya convergence if Balmer lines are thick */
				if( Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn > 0. )
				{
					differ = fabs(1.-Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn*
					  rt.DoubleTau/Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauTot)*100.;

					/* check whether Lya optical depth changed on last iteration
					 * no tolerance can be finer than autocv, the tolerance on the
					 * iterate to convergence command.  It is 15% */
					if( ((iterations.lgLastIt && Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn > 0.8) && 
						differ > 25.) && wind.lgStatic() )
					{
						sprintf( chLine, 
							" C-This is the last iteration and %2s%2ld Ly(a) optical depth"
							" changed by%7.0f%% (was %.2e). Try another iteration.", 
						elementnames.chElementSym[nelem], 
						  nelem+1,differ, Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauTot );
						caunin(chLine);
						iterations.lgIterAgain = true;
					}
				}
			}
		}
	}

	/* check whether sphere was set if dr/r large */
	if( radius.Radius/radius.rinner > 2. && !geometry.lgSphere )
	{
		sprintf( chLine, " C-R(out)/R(in)=%.2e and SPHERE was not set.", 
		  radius.Radius/radius.rinner );
		caunin(chLine);
	}

	/* check if thin in hydrogen or helium continua, but assumed to be thick */
	if( iterations.lgLastIt && !opac.lgCaseB )
	{

		/* check if thin in Lyman continuum, and assumed thick */
		if( rfield.nflux > iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s] )
		{
			if( opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1] < 2. && 
				opac.TauAbsGeo[1][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1] > 2. )
			{
				sprintf( chLine, " C-The H Lyman continuum is thin, and I assumed"
					" that it was thick.  Try another iteration." );
				caunin(chLine);
				iterations.lgIterAgain = true;
			}
		}

		/* check on the He+ ionizing continuum */
		if( rfield.nflux > iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s] && dense.lgElmtOn[ipHELIUM] )
		{
			if( (opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s]-1] < 2. && 
				 opac.TauAbsGeo[1][iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s]-1] > 2.) )
			{
				sprintf( chLine, 
					" C-The He II continuum is thin and I assumed that it was thick."
					"  Try another iteration." );
				caunin(chLine);
				iterations.lgIterAgain = true;
			}
		}

		if( rfield.nflux > iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0] && dense.lgElmtOn[ipHELIUM] )
		{ 
			if( (opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0]-1] < 2. && 
				 opac.TauAbsGeo[1][iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0]-1] > 2.) )
			{
				sprintf( chLine, 
					" C-The He I continuum is thin and I assumed that it was thick."
					"  Try another iteration." );
				caunin(chLine);
				iterations.lgIterAgain = true;
			}
		}
	}

	/* check whether column density changed by much on this iteration */
	if( iteration > 1 )
	{
		if( colden.colden_old[ipCOL_HTOT] <= 0. )
		{
			fprintf( ioQQQ, " colden_old is insane in PrtComment.\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		differ = fabs(1.-colden.colden[ipCOL_HTOT]/
			colden.colden_old[ipCOL_HTOT]);

		if( differ > 0.1 && differ <= 0.3 )
		{
			sprintf( chLine, 
				"   The H column density changed by %.2e%% between this and previous iteration.", 
			  differ*100. );
			notein(chLine);
		}

		else if( differ > 0.3 )
		{
			if( iterations.lgLastIt )
			{
				sprintf( chLine, 
					" C-The H column density changed by %.2e%% and this is the last iteration.  What happened?", 
				  differ*100. );
				caunin(chLine);
			}
			else
			{
				sprintf( chLine, 
					"  !The H column density changed by %.2e%%  What happened?", 
				  differ*100. );
				bangin(chLine);
			}
		}

		/* check on H2 column density, but only if significant fraction of H is molecular */
		if( (colden.colden[ipCOL_H2g]+colden.colden[ipCOL_H2s])/SDIV(colden.colden[ipCOL_HTOT]) > 1e-5 )
		{
			differ = fabs(1.-colden.colden[ipCOL_H2g]/
				SDIV(colden.colden_old[ipCOL_H2g]));

			if( differ > 0.1 && differ <= 0.3 )
			{
				sprintf( chLine, 
					"   The H2 column density changed by %.2e%% between this and previous iteration.", 
				differ*100. );
				notein(chLine);
			}

			else if( differ > 0.3 )
			{
				if( iterations.lgLastIt )
				{
					sprintf( chLine, 
						" C-The H2 column density changed by %.2e%% and this is the last iteration.  What happened?", 
					differ*100. );
					caunin(chLine);
				}
				else
				{
					sprintf( chLine, 
						"  !The H2 column density changed by %.2e%%  What happened?", 
					differ*100. );
					bangin(chLine);
				}
			}
		}
	}

	/* say if rad pressure caused by la and la optical depth changed too much */
	differ = fabs(1.-Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauIn/
	  SDIV(Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauTot))*100.;

	if( iterations.lgLastIt && (pressure.RadBetaMax > 0.1) && 
		(differ > 50.) && (pressure.ipPradMax_line == 1) && (pressure.lgPres_radiation_ON) && 
		wind.lgStatic() )
	{
		sprintf( chLine, " C-This is the last iteration, radiation pressure was significant, and the L-a optical depth changed by %7.2f%% (was %.2e)", 
		  differ, Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauTot );
		caunin(chLine);
	}

	/* caution that 21 cm spin temperature is incorrect when Lya optical depth
	 * scale is overrun */
	if( iterations.lgLastIt &&
		( Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauTot * 1.02 -
		Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauIn ) < 0. )
	{
		sprintf( chLine, " C-The Lya optical depth scale was overrun and this is the last iteration - Tspin(21 cm) is not valid." );
		caunin(chLine);
		sprintf( chLine, " C-Another iteration is needed for Tspin(21 cm) to be valid." );
		caunin(chLine);
	}

	/* say if la rad pressure capped by thermalization length */
	if( pressure.lgPradDen )
	{
		sprintf( chLine, "   Line radiation pressure capped by thermalization length." );
		notein(chLine);
	}

	/* print te failures */
	nline = MIN2(conv.nTeFail,10);
	if( conv.nTeFail != 0 )
	{
		long int _o;
		if( conv.failmx < 0.1 )
		{
			_o = sprintf( chLine, "   There were %ld minor temperature failures.  zones:", 
				conv.nTeFail );
			/* don't know how many zones we will save, there are nline,
			 * hence this use of pointer arith */
			for( i=0; i < nline; i++ )
			{
				_o += sprintf( chLine+_o, " %ld", conv.ifailz[i] );
			}
			notein(chLine);
		}
		else
		{
			_o = sprintf( chLine, 
				"  !There were %ld temperature failures, and some were large. The largest was %.1f%%.  What happened?", 
			  conv.nTeFail, conv.failmx*100. );
			bangin(chLine);

			/* don't know how many zones we will save, there are nline,
			 * hence this use of pointer arith */
			_o = sprintf( chLine , "  !The zones were" );
			for( i=0; i < nline; i++ )
			{
				_o += sprintf( chLine+_o, " %ld", conv.ifailz[i] );
			}
			bangin(chLine);

			if( struc.testr[0] > 8e4 && phycon.te < 6e5 )
			{
				sprintf( chLine, "  !I think they may have been caused by the change from hot to nebular gas phase.  The physics of this is unclear." );
				bangin(chLine);
			}
		}
	}

	/* check for temperature jumps */
	big_ion_jump = 0.;
	j = 0;
	for( i=1; i < nzone; i++ )
	{
		big = fabs(1.-struc.testr[i-1]/struc.testr[i]);
		if( big > big_ion_jump )
		{
			j = i;
			big_ion_jump = big;
		}
	}

	if( big_ion_jump > 0.2 )
	{
		/* this is a sanity check, but only do it if jump detected */
		if( j < 1 )
		{
			fprintf( ioQQQ, " j too small big jump check\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		if( big_ion_jump > 0.4 )
		{
			sprintf( chLine, " C-A temperature discontinuity occurred at zone %ld from %.2eK to %.2eK.", 
			  j, struc.testr[j-1], struc.testr[j] );
			caunin(chLine);
			/* check if the second temperature is between 100 and 1000K */
			/* >>chng 05 nov 07, test second not first temperature since second
			 * will be lower of the two */
			/*if( struc.testr[j-1] < 1000. && struc.testr[j-1]>100. )*/
			if( struc.testr[j]>100. && struc.testr[j] < 1000. )
			{
				sprintf( chLine, " C-This was probably due to a thermal front." );
				caunin(chLine);
			}
		}
		else if( big_ion_jump > 0.2 )
		{
			sprintf( chLine, "  !A temperature discontinuity occurred at zone %ld from %.2eK to %.2eK.", 
			  j, struc.testr[j-1], struc.testr[j] );
			bangin(chLine);
			/* check if the second temperature is between 100 and 1000K */
			/* >>chng 05 nov 07, test second not first temperature since second
			 * will be lower of the two */
			/*if( struc.testr[j-1] < 1000. && struc.testr[j-1]>100. )*/
			if( struc.testr[j]>100. && struc.testr[j] < 1000. )
			{
				sprintf( chLine, "  !This was probably due to a thermal front." );
				bangin(chLine);
			}
		}
	}

	/* check for largest error in local electron density */
	if( fabs(conv.BigEdenError) > conv.EdenErrorAllowed )
	{
		/* this only produces a warning if not the very last zone */
		if( fabs(conv.BigEdenError) > conv.EdenErrorAllowed*20. && dense.nzEdenBad != 
		  nzone )
		{
			sprintf( chLine, " W-The local error in the electron density reached %.1f%% at zone %ld", 
			  conv.BigEdenError*100, dense.nzEdenBad );
			warnin(chLine);
		}
		else if( fabs(conv.BigEdenError) > conv.EdenErrorAllowed*5. )
		{
			sprintf( chLine, " C-The local error in the electron density reached %.1f%% at zone %ld", 
			  conv.BigEdenError*100, dense.nzEdenBad );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, "   The local error in the electron density reached %.1f%% at zone %ld", 
			  conv.BigEdenError*100, dense.nzEdenBad );
			notein(chLine);
		}
	}

	/* check for temperature oscillations or fluctuations*/
	big_ion_jump = 0.;
	j = 0;
	for( i=1; i < (nzone - 1); i++ )
	{
		big = fabs( (struc.testr[i-1] - struc.testr[i])/struc.testr[i] );
		bigm = fabs( (struc.testr[i] - struc.testr[i+1])/struc.testr[i] );

		/* this is sign of change in temperature, we are looking for change in sign */
		rel = ( (struc.testr[i-1] - struc.testr[i])/struc.testr[i])*
			( (struc.testr[i] - struc.testr[i+1])/struc.testr[i] );

		if( rel < 0. && MIN2( bigm , big ) > big_ion_jump )
		{
			j = i;
			big_ion_jump = MIN2( bigm , big );
		}
	}

	if( big_ion_jump > 0.1 )
	{
		/* only do sanity check if jump detected */
		if( j < 1 )
		{
			fprintf( ioQQQ, " j too small bigjump2 check\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		if( big_ion_jump > 0.3 )
		{
			sprintf( chLine, 
				" C-A temperature oscillation occurred at zone%4ld by%5.0f%% from %.2e to %.2e to %.2e", 
			  j, big_ion_jump*100., struc.testr[j-1], struc.testr[j], struc.testr[j+1] );
			caunin(chLine);
		}
		else if( big_ion_jump > 0.1 )
		{
			sprintf( chLine, 
				"  !A temperature oscillation occurred at zone%4ld by%5.0f%% from %.2e to %.2e to %.2e", 
			  j, big_ion_jump*100., struc.testr[j-1], struc.testr[j], struc.testr[j+1] );
			bangin(chLine);
		}
	}

	/* check for eden oscillations */
	if( strcmp(dense.chDenseLaw,"CDEN") == 0 )
	{
		j = 0;
		big_ion_jump = 0.;
		for( i=1; i < (nzone - 1); i++ )
		{
			big = (struc.ednstr[i-1] - struc.ednstr[i])/struc.ednstr[i];
			if( fabs(big) < conv.EdenErrorAllowed )
				big = 0.;
			bigm = (struc.ednstr[i] - struc.ednstr[i+1])/struc.ednstr[i];
			if( fabs(bigm) < conv.EdenErrorAllowed )
				bigm = 0.;
			if( big*bigm < 0. && 
				fabs(struc.ednstr[i-1]-struc.ednstr[i])/struc.ednstr[i] > big_ion_jump )
			{
				j = i;
				big_ion_jump = fabs(struc.ednstr[i-1]-struc.ednstr[i])/
				  struc.ednstr[i];
			}
		}

		/* only check on j if there was a big jump detected, number must be
		 * smallest jump */
		if( big_ion_jump > conv.EdenErrorAllowed*3. )
		{
			if( j < 1 )
			{
				fprintf( ioQQQ, " j too small bigjump3 check\n" );
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}

			if( big_ion_jump > conv.EdenErrorAllowed*10. )
			{
				sprintf( chLine, " C-An electron density oscillation occurred at zone%4ld by%5.0f%% from %.2e to %.2e to %.2e", 
				  j, big_ion_jump*100., struc.ednstr[j-1], struc.ednstr[j], 
				  struc.ednstr[j+1] );
				caunin(chLine);
			}
			else if( big_ion_jump > conv.EdenErrorAllowed*3. )
			{
				sprintf( chLine, "  !An electron density oscillation occurred at zone%4ld by%5.0f%% from %.2e to %.2e to %.2e", 
				  j, big_ion_jump*100., struc.ednstr[j-1], struc.ednstr[j], 
				  struc.ednstr[j+1] );
				bangin(chLine);
			}
		}
	}

	/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
	/* >>chng 03 dec 05, add this test */
	prt_smooth_predictions();

	/**********************************************************
	 * check that the volume integrates out ok                *
	 **********************************************************/

	/* this was the number 1 fed through the line integrators,
	 * the number 1e-10 is sent to linadd in lineset1 as follows:*/
	/*linadd( 1.e-10 , 1 , "Unit" , 'i' );*/
	i = cdLine( "Unit" , 1 , &rate , &absint );
	ASSERT( i> 0 );

	/* this is now the linear vol, rel to inner radius */
	VolComputed = LineSv[i].SumLine[0] /  1e-10;

	/* spherical or plane parallel case? */
	if( radius.Radius/radius.rinner > 1.0001 )
	{
		/* spherical case, 
		 * geometry.iEmissPower is usually 2,
		 * and can be reset to 1 (long slit) or 0 (beam) with 
		 * slit and beam options on aperture */
		VolExpected = geometry.covaper*geometry.FillFac*radius.rinner/(geometry.iEmissPower+1)*
			( powi( radius.Radius/radius.rinner,geometry.iEmissPower+1 ) - 1. );
	}
	else
	{
		/* plane parallel case */
		/* next can be zero for very thin model, depth is always positive */
		VolExpected = geometry.covaper*geometry.FillFac*(radius.depth-DEPTH_OFFSET);
	}

	/* now get the relative difference between computed and expected volumes */
	error = fabs(VolComputed - VolExpected)/SDIV(VolExpected);

	/* we need to ignore this test if filling factor changes with radius, or
	 * cylinder geometry in place */
	if( radius.lgCylnOn || geometry.filpow!=0. )
	{
		error = 0.;
	}

	/* how large is relative error? */
	if( error > 0.001 && !lgAbort )
	{
		sprintf( chLine, 
			" W-PrtComment insanity - Line unit integration did not verify \n");
		warnin(chLine);
		fprintf( ioQQQ,
			" PROBLEM PrtComment insanity - Line unit integration did not verify \n");
		fprintf( ioQQQ,
			" expected, derived vols were %g %g \n",
			VolExpected , VolComputed );
		fprintf( ioQQQ,
			" relative difference is %g, ratio is %g.\n",error,VolComputed/VolExpected);
		TotalInsanity();
	}

	/* next do same thing for fake continuum point propagated in highest energy cell, plus 1 
	 *  = 
	 * the variable rfield.ConEmitLocal[rfield.nflux]
	 * are set to 
	 * the number 1.e-10f * Dilution in RT_diffuse.  this is the outward
	 * local emissivity, per unit vol.  It is then added to the outward beams
	 * by the rest of the code, and then checked here.
	 *
	 * insanity will be detected if diffuse emission is thrown into the outward beam
	 * in MadeDiffuse.  this happens if the range of ionization encompasses the full
	 * continuum array, up to nflux.  */
	ConComputed = rfield.ConInterOut[rfield.nflux]/ 1e-10;
	/* correct for fraction that went out, as set in ZoneStart,
	 * this should now be the volume of the emitting region */
	ConComputed /= ( (1. + geometry.covrt)/2. );

	/* we expect this to add up to the integral of unity over r^-2 */
	if( radius.Radius/radius.rinner < 1.001 )
	{
		/* plane parallel case, no dilution, use thickness */
		ConExpected = (radius.depth-DEPTH_OFFSET)*geometry.FillFac;
	}
	else
	{
		/* spherical case */
		ConExpected = radius.rinner*geometry.FillFac * (1. - radius.rinner/radius.Radius );
	}
	/* this is impossible */
	ASSERT( ConExpected > 0. );

	/* now get the relative difference between computed and expected volumes */
	error = fabs(ConComputed - ConExpected)/ConExpected;

	/* we need to ignore this test if filling factor changes with radius, or
	 * cylinder geometry in place */
	if( radius.lgCylnOn || geometry.filpow!=0. )
	{
		error = 0.;
	}

	/* \todo 2 - These "volumes" seem to be too small by a factor of two.
	 * rfield.ConInterOut[rfield.nflux] (hence ConComputed) and ConExpected  
	 * should be greater by a factor of 2 if comparison is really of "volume"
	 * of 1/cc pencil. */

	/* how large is relative error? */
	if( error > 0.001 && !lgAbort)
	{
		sprintf( chLine, 
			" W-PrtComment insanity - Continuum unit integration did not verify \n");
		warnin(chLine);
		fprintf( ioQQQ," PROBLEM PrtComment insanity - Continuum unit integration did not verify \n");
		fprintf( ioQQQ," exact vol= %g, derived vol= %g relative difference is %g \n",
			ConExpected , ConComputed ,error);
		fprintf( ioQQQ," ConInterOut= %g,  \n",
			rfield.ConInterOut[rfield.nflux]);
		TotalInsanity();
	}

	/* final printout of warnings, cautions, notes, */
	cdNwcns(&lgAbort_flag,&nw,&nc,&nn,&ns,&i,&j,&dum1,&dum2);

	warnings.lgWarngs = ( nw > 0 );
	warnings.lgCautns = ( nc > 0 );

	if( called.lgTalk )
	{
		/* print the title of the calculation */
		fprintf( ioQQQ, "   %s\n", input.chTitle  );
		/* say why the calculation stopped, and indicate the geometry*/
		cdReasonGeo(ioQQQ);
		/* print all warnings */
		cdWarnings(ioQQQ);
		/* all cautions */
		cdCautions(ioQQQ);
		/* surprises, beginning with a ! */
		cdSurprises(ioQQQ);
		/* notes about the calculations */
		cdNotes(ioQQQ);
	}

	/* option to print warnings on special io */
	if( lgPrnErr )
	{
		cdWarnings(ioPrnErr);
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PrtComment returns.\n" );
	}
	return;
}


/*chkCaHeps check whether CaII K and H epsilon overlap */
STATIC void chkCaHeps(double *totwid)
{
	double conca, 
		conalog ,
	  conhe;

	DEBUG_ENTRY( "chkCaHeps()" );

	*totwid = 0.;
	/* pumping of CaH overlapping with Hy epsilon, 6-2 of H */
	if( iso.n_HighestResolved_local[ipH_LIKE][ipHYDROGEN] + 
		iso.nCollapsed_local[ipH_LIKE][ipHYDROGEN] >= 6 )
	{
		/* this is 6P */
		long ip6p = iso.QuantumNumbers2Index[ipH_LIKE][ipHYDROGEN][6][1][2];

		if( TauLines[ipT3969].Emis->TauIn > 0. && 
			Transitions[ipH_LIKE][ipHYDROGEN][ip6p][ipH2s].Emis->TauIn >   0. )
		{
			/* casts to double here are to prevent FPE */
			conca = pow(6.1e-5* TauLines[ipT3969].Emis->TauIn,0.5);
			conalog = log((double)TauLines[ipT3969].Emis->TauIn);
			conalog = sqrt(MAX2(1., conalog));
			conca = MAX2(conalog,conca);

			conalog = log((double)Transitions[ipH_LIKE][ipHYDROGEN][ip6p][ipH2s].Emis->TauIn);
			conalog = sqrt(MAX2(1.,conalog));
			conhe = pow(1.7e-6*Transitions[ipH_LIKE][ipHYDROGEN][ip6p][ipH2s].Emis->TauIn,0.5);
			conhe = MAX2(conalog, conhe);

			*totwid = 10.*conhe + 1.6*conca;
		}
	}
	return;
}

/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
STATIC void prt_smooth_predictions(void)
{
	long int i,
		nzone_oscillation,
		nzone_ion_jump,
		nzone_den_jump,
		nelem,
		ion;
	double BigOscillation ,
		big_ion_jump,
		big_jump,
		rel,
		big,
		bigm;

	char chLine[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "prt_smooth_predictions()" );

	/* check for ionization oscillations or fluctuations and or jumps */
	nzone_oscillation = 0;
	nzone_ion_jump = 0;

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( ion=0; ion<=nelem+1; ++ion) 
			{
				BigOscillation = 0.;
				big_ion_jump = -15.;
				/* >>chng 05 mar 12, add -lgAbort, since in some bad aborts current zone never evaluated */
				for( i=1; i < (nzone - 1-(int)lgAbort); i++ )
				{

					/* only do check if all ions are positive */
					if( struc.xIonDense[nelem][ion][i-1]/struc.gas_phase[nelem][i-1]>struc.dr_ionfrac_limit &&
						struc.xIonDense[nelem][ion][i  ]/struc.gas_phase[nelem][i  ]>struc.dr_ionfrac_limit &&
						struc.xIonDense[nelem][ion][i+1]/struc.gas_phase[nelem][i+1]>struc.dr_ionfrac_limit )
					{

						/* this is check for oscillations */
						big = fabs( (struc.xIonDense[nelem][ion][i-1] - struc.xIonDense[nelem][ion][i])/struc.xIonDense[nelem][ion][i] );
						bigm = fabs( (struc.xIonDense[nelem][ion][i]  - struc.xIonDense[nelem][ion][i+1])/struc.xIonDense[nelem][ion][i] );

						/* this is sign of change in ionization, we are looking for change in sign */
						rel = ( (struc.xIonDense[nelem][ion][i-1] - struc.xIonDense[nelem][ion][i]  )/struc.xIonDense[nelem][ion][i])*
							  ( (struc.xIonDense[nelem][ion][i]   - struc.xIonDense[nelem][ion][i+1])/struc.xIonDense[nelem][ion][i] );

						if( rel < 0. && MIN2( bigm , big ) > BigOscillation )
						{
							nzone_oscillation = i;
							BigOscillation = MIN2( bigm , big );
						}

						/* check whether we tripped over an ionization front - a major source
						 * of instability in a complete linearization code like this one */
						/* neg sign picks up only increases in ionization */
						rel = -log10( (struc.xIonDense[nelem][ion][i]/struc.gas_phase[nelem][i]) / 
							(struc.xIonDense[nelem][ion][i+1]/struc.gas_phase[nelem][i+1] ) );
						/* only do significant stages of ionization */
						if( rel > big_ion_jump )
						{
							big_ion_jump = rel;
							nzone_ion_jump = i;
						}
					}
				}
				/* end loop over zones, 
				 * check whether this ion and element underwent fluctuations or jump */

				if( BigOscillation > 0.2 )
				{
					/* only do sanity check if jump detected */
					if( nzone_oscillation < 1 )
					{
						fprintf( ioQQQ, " nzone_oscillation too small bigjump2 check\n" );
						ShowMe();
						cdEXIT(EXIT_FAILURE);
					}
					if( BigOscillation > 3. )
					{
						sprintf( chLine, 
							" W-An ionization oscillation occurred at zone %ld, elem %.2s%2li, by %.0f%% from %.2e to %.2e to %.2e", 
							nzone_oscillation, 
							elementnames.chElementSym[nelem], ion+1,
							BigOscillation*100., 
							struc.xIonDense[nelem][ion][nzone_oscillation-1]/struc.gas_phase[nelem][nzone_oscillation-1], 
							struc.xIonDense[nelem][ion][nzone_oscillation]/struc.gas_phase[nelem][nzone_oscillation], 
							struc.xIonDense[nelem][ion][nzone_oscillation+1]/struc.gas_phase[nelem][nzone_oscillation+1] );
							warnin(chLine);
					}

					else if( BigOscillation > 0.7 )
					{
						sprintf( chLine, 
							" C-An ionization oscillation occurred at zone %ld, elem %.2s%2li, by %.0f%% from %.2e to %.2e to %.2e", 
							nzone_oscillation, 
							elementnames.chElementSym[nelem], ion+1,
							BigOscillation*100., 
							struc.xIonDense[nelem][ion][nzone_oscillation-1]/struc.gas_phase[nelem][nzone_oscillation-1], 
							struc.xIonDense[nelem][ion][nzone_oscillation]/struc.gas_phase[nelem][nzone_oscillation], 
							struc.xIonDense[nelem][ion][nzone_oscillation+1]/struc.gas_phase[nelem][nzone_oscillation+1] );
							caunin(chLine);
					}
					else if( BigOscillation > 0.2 )
					{
						sprintf( chLine, 
							"  !An ionization oscillation occurred at zone %ld, elem %.2s%2li, by %.0f%% from %.2e to %.2e to %.2e", 
							nzone_oscillation, 
							elementnames.chElementSym[nelem], ion+1,
							BigOscillation*100., 
							struc.xIonDense[nelem][ion][nzone_oscillation-1]/struc.gas_phase[nelem][nzone_oscillation-1], 
							struc.xIonDense[nelem][ion][nzone_oscillation]/struc.gas_phase[nelem][nzone_oscillation], 
							struc.xIonDense[nelem][ion][nzone_oscillation+1]/struc.gas_phase[nelem][nzone_oscillation+1] );
							bangin(chLine);
					}
				}

				/* big_ion_jump was a log above, convert to linear quantity */
				/* if no jump occurred then big_ion_jump is small and nzone_ion_jump is 0 */
				big_ion_jump = pow(10., big_ion_jump );
				if( big_ion_jump > 1.5 && nzone_ion_jump > 0 )
				{
					if( big_ion_jump > 10. )
					{
						sprintf( chLine, 
							" C-An ionization jump occurred at zone %ld, elem %.2s%2li, by %.0f%% from %.2e to %.2e to %.2e", 
							nzone_ion_jump, 
							elementnames.chElementSym[nelem], ion+1,
							big_ion_jump*100., 
							struc.xIonDense[nelem][ion][nzone_ion_jump-1]/struc.gas_phase[nelem][nzone_ion_jump-1], 
							struc.xIonDense[nelem][ion][nzone_ion_jump]/struc.gas_phase[nelem][nzone_ion_jump], 
							struc.xIonDense[nelem][ion][nzone_ion_jump+1]/struc.gas_phase[nelem][nzone_ion_jump+1] );
							caunin(chLine);
					}
					else
					{
						sprintf( chLine, 
							"  !An ionization jump occurred at zone %ld, elem %.2s%2li, by %.0f%% from %.2e to %.2e to %.2e", 
							nzone_ion_jump, 
							elementnames.chElementSym[nelem], ion+1,
							big_ion_jump*100., 
							struc.xIonDense[nelem][ion][nzone_ion_jump-1]/struc.gas_phase[nelem][nzone_ion_jump-1], 
							struc.xIonDense[nelem][ion][nzone_ion_jump]/struc.gas_phase[nelem][nzone_ion_jump], 
							struc.xIonDense[nelem][ion][nzone_ion_jump+1]/struc.gas_phase[nelem][nzone_ion_jump+1] );
							bangin(chLine);
					}
				}
			}
		}
	}

	big_jump = -15;
	nzone_den_jump = 0;

	/* >>chng 05 mar 12, add -lgAbort, since in some bad aborts current zone never evaluated */
	for( i=1; i < (nzone - 1 - (int)lgAbort); i++ )
	{
		/* this first check is on how the total hydrogen density has changed */
		rel = fabs(log10( struc.gas_phase[ipHYDROGEN][i] / 
			struc.gas_phase[ipHYDROGEN][i+1] ) );
		/* only do significant stages of ionization */
		if( rel > big_jump )
		{
			big_jump = rel;
			nzone_den_jump = i;
		}
	}

	/* check how stable density was */
	big_jump = pow( 10., big_jump );
	if( big_jump > 1.2 )
	{
		if( big_jump > 3. )
		{
			sprintf( chLine, 
				" C-The H density jumped at by %.0f%% at zone %ld, from %.2e to %.2e to %.2e", 
				big_jump*100., 
				nzone_den_jump, 
				struc.gas_phase[ipHYDROGEN][nzone_den_jump-1], 
				struc.gas_phase[ipHYDROGEN][nzone_den_jump], 
				struc.gas_phase[ipHYDROGEN][nzone_den_jump+1] );
				caunin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"  !An H density jump occurred at zone %ld, by %.0f%% from %.2e to %.2e to %.2e", 
				nzone_den_jump, 
				big_jump*100., 
				struc.gas_phase[ipHYDROGEN][nzone_den_jump-1], 
				struc.gas_phase[ipHYDROGEN][nzone_den_jump], 
				struc.gas_phase[ipHYDROGEN][nzone_den_jump+1] );
				bangin(chLine);
		}
	}

	/* now do check on smoothness of radiation pressure */
	big_jump = -15;
	nzone_den_jump = 0;

	/* loop starts on zone 3 since dramatic fall in radiation pressure across first
	 * few zones is normal behavior */
	/* >>chng 05 mar 12, add -lgAbort, since in some bad aborts current zone never evaluated */
	for( i=3; i < (nzone - 2 - (int)lgAbort); i++ )
	{
		/* this first check is on how the total hydrogen density has changed */
		rel = fabs(log10( SDIV(struc.pres_radiation_lines_curr[i]) / 
			SDIV(0.5*(struc.pres_radiation_lines_curr[i-1]+struc.pres_radiation_lines_curr[i+1])) ) );
		/* only do significant stages of ionization */
		if( rel > big_jump )
		{
			big_jump = rel;
			nzone_den_jump = i;
		}
	}
	/* note that changing log big_jump to linear takes place in next branch */

	/* check how stable radiation pressure was, but only if significant */
	if( pressure.RadBetaMax > 0.01 )
	{
		big_jump = pow( 10., big_jump );
		if( big_jump > 1.2 )
		{
			/* only make it a caution is pressure jumped, and we were trying
			* to do a constant pressure model */
			if( big_jump > 3. && strcmp(dense.chDenseLaw,"CPRE") == 0)
			{
				sprintf( chLine, 
					" C-The radiation pressure jumped by %.0f%% at zone %ld, from %.2e to %.2e to %.2e", 
					big_jump*100., 
					nzone_den_jump, 
					struc.pres_radiation_lines_curr[nzone_den_jump-1], 
					struc.pres_radiation_lines_curr[nzone_den_jump], 
					struc.pres_radiation_lines_curr[nzone_den_jump+1] );
					caunin(chLine);
			}
			else
			{
				sprintf( chLine, 
					"  !The radiation pressure jumped by %.0f%% at zone %ld, from %.2e to %.2e to %.2e", 
					big_jump*100., 
					nzone_den_jump, 
					struc.pres_radiation_lines_curr[nzone_den_jump-1], 
					struc.pres_radiation_lines_curr[nzone_den_jump], 
					struc.pres_radiation_lines_curr[nzone_den_jump+1] );
					bangin(chLine);
			}
		}
	}

	/* these will be used to check on continuity */
	phycon.BigJumpTe = 0.;
	phycon.BigJumpne = 0.;
	phycon.BigJumpH2 = 0.;
	phycon.BigJumpCO = 0.;

	for( i=1; i < (nzone - 1 - (int)lgAbort); i++ )
	{
		/* check on how much temperature has changed */
		rel = fabs(log10( struc.testr[i] / struc.testr[i+1] ) );
		if( rel > phycon.BigJumpTe )
		{
			phycon.BigJumpTe = (realnum)rel;
		}

		/* check on how much electron density has changed */
		rel = fabs(log10( struc.ednstr[i] / struc.ednstr[i+1] ) );
		if( rel > phycon.BigJumpne )
		{
			phycon.BigJumpne = (realnum)rel;
		}

		/* check on how much H2 density has changed */
		if( (struc.H2_molec[ipMH2g][i]+struc.H2_molec[ipMH2s][i])>SMALLFLOAT &&
			(struc.H2_molec[ipMH2g][i+1]+struc.H2_molec[ipMH2s][i+1]) > SMALLFLOAT 
			/* only do this test if H2 abund is significant */
			&& (struc.H2_molec[ipMH2g][i]+struc.H2_molec[ipMH2s][i])/struc.gas_phase[ipHYDROGEN][i]>1e-3)
		{
			rel = fabs(log10( (struc.H2_molec[ipMH2g][i]+struc.H2_molec[ipMH2s][i]) / 
				SDIV(struc.H2_molec[ipMH2g][i+1]+struc.H2_molec[ipMH2s][i+1]) ) );
			if( rel > phycon.BigJumpH2 )
			{
				phycon.BigJumpH2 = (realnum)rel;
			}
		}

		{ 
			int ipCO = findspecies("CO")->index;
			/* check on how much CO density has changed */
			if( struc.CO_molec[ipCO][i]>SMALLFLOAT &&
					struc.CO_molec[ipCO][i+1]>SMALLFLOAT &&
					struc.CO_molec[ipCO][i]/SDIV(struc.gas_phase[ipCARBON][i])>1e-3 )
			{
				rel = fabs(log10( struc.CO_molec[ipCO][i] / struc.CO_molec[ipCO][i+1] ) );
				if( rel > phycon.BigJumpCO )
				{
					phycon.BigJumpCO = (realnum)rel;
				}
			}
		}
	}

	/* convert to linear change - subtract 1 to make it the residual difference */
	if(phycon.BigJumpTe>0. )
		phycon.BigJumpTe = (realnum)pow( 10., (double)phycon.BigJumpTe ) -1.f;

	if(phycon.BigJumpne>0. )
		phycon.BigJumpne = (realnum)pow( 10., (double)phycon.BigJumpne )-1.f;

	if(phycon.BigJumpH2>0. )
		phycon.BigJumpH2 = (realnum)pow( 10., (double)phycon.BigJumpH2 )-1.f;

	if(phycon.BigJumpCO>0. )
		phycon.BigJumpCO = (realnum)pow( 10., (double)phycon.BigJumpCO )-1.f;
	/*fprintf(ioQQQ,"DEBUG continuity large change %.2e %.2e %.2e %.2e \n",
		phycon.BigJumpTe , phycon.BigJumpne , phycon.BigJumpH2 , phycon.BigJumpCO );*/

	if( phycon.BigJumpTe > 0.3 )
	{
		sprintf( chLine, 
			" C-The temperature varied by %.1f%% between two zones", 
			phycon.BigJumpTe*100.);
			caunin(chLine);
	}
	else if( phycon.BigJumpTe > 0.1 )
	{
		sprintf( chLine, 
			"  !The temperature varied by %.1f%% between two zones", 
			phycon.BigJumpTe*100.);
			bangin(chLine);
	}

	if( phycon.BigJumpne > 0.3 )
	{
		sprintf( chLine, 
			" C-The electron density varied by %.1f%% between two zones", 
			phycon.BigJumpne*100.);
			caunin(chLine);
	}
	else if( phycon.BigJumpne > 0.1 )
	{
		sprintf( chLine, 
			"  !The electron density varied by %.1f%% between two zones", 
			phycon.BigJumpne*100.);
			bangin(chLine);
	}

	if( phycon.BigJumpH2 > 0.8 )
	{
		sprintf( chLine, 
			" C-The H2 density varied by %.1f%% between two zones", 
			phycon.BigJumpH2*100.);
			caunin(chLine);
	}
	else if( phycon.BigJumpH2 > 0.1 )
	{
		sprintf( chLine, 
			"  !The H2 density varied by %.1f%% between two zones", 
			phycon.BigJumpH2*100.);
			bangin(chLine);
	}

	if( phycon.BigJumpCO > 0.8 )
	{
		sprintf( chLine, 
			" C-The CO density varied by %.1f%% between two zones", 
			phycon.BigJumpCO*100.);
			caunin(chLine);
	}
	else if( phycon.BigJumpCO > 0.2 )
	{
		sprintf( chLine, 
			"  !The CO density varied by %.1f%% between two zones", 
			phycon.BigJumpCO*100.);
			bangin(chLine);
	}
	return;
}
