/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* iso_solve main routine to call iso_level and determine iso level balances */
/* HydroRenorm - renormalize H so that it agrees with the chemistry */
/* AGN_He1_CS routine to save table needed for AGN3 - collision strengths of HeI */
#include "cddefines.h"
#include "atmdat.h"
#include "conv.h"
#include "dense.h"
#include "opacity.h"
#include "elementnames.h"
#include "h2.h"
#include "helike.h"
#include "helike_cs.h"
#include "hmi.h"
#include "hydrogenic.h"
#include "ionbal.h"
#include "iso.h"
#include "phycon.h"
#include "secondaries.h"
#include "taulines.h"
#include "thermal.h"
#include "trace.h"

valarray<double> PumpSave;

void iso_drive( void )
{
	DEBUG_ENTRY( "iso_drive()" );

	iso_charge_transfer_update();

	/* update and solve iso levels */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
		for( long nelem=ipISO; nelem<LIMELM; nelem++ )
		{
			iso_update_rates( ipISO, nelem );
			iso_solve( ipISO, nelem );
		}

	return;
}

/* iso_update_rates routine to set up iso rates, level balance is done elsewhere */
void iso_update_rates(long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_update_rates()" );

	/* 
	 * This is an option to account for intrinsic absorption or emission line the lyman 
	 * lines.  This is done without changing the incident coarse continuum.  
	 * Check if continuum pumping of H Lyman lines to be multiplied by scale factor
	 * hydro.xLymanPumpingScaleFactor is set with atom h-like Lyman pumping scale command 
	 * Lya pump rate is always updated - this is simplest way to finese intrinsic absorption
	 * or emission 
	 */
	if( ipISO==ipH_LIKE && nelem==ipHYDROGEN && hydro.xLymanPumpingScaleFactor!=1.f )
	{
		if( PumpSave.size() != (unsigned int)iso.numLevels_local[ipISO][nelem] )
			PumpSave.resize( iso.numLevels_local[ipISO][nelem] );
		for( long ipHi=1; ipHi < iso.numLevels_local[ipH_LIKE][ipHYDROGEN]; ++ipHi )
		{
			PumpSave[ipHi] = Transitions[ipISO][ipHYDROGEN][ipHi][0].Emis->pump;
			Transitions[ipISO][ipHYDROGEN][ipHi][0].Emis->pump *= hydro.xLymanPumpingScaleFactor;
		}
	}

	/* do not consider elements that have been turned off */
	if( dense.lgElmtOn[nelem] )
	{
		/* note that nelem scale is totally on c not physical scale, so 0 is h */
		/* evaluate balance if ionization reaches this high */
		if( ((dense.IonHigh[nelem] >= nelem - ipISO) &&
			(dense.IonLow[nelem] <= nelem - ipISO))  || !conv.nTotalIoniz )
		{
			/* truncate atom if physical conditions limit the maximum principal quantum number of a
			 * bound electron to a number less than the malloc'd size */
			if( iso.lgContinuumLoweringEnabled[ipISO] )
				iso_continuum_lower( ipISO, nelem );

			/* evaluate collisional rates */
			iso_collide( ipISO,  nelem );

			/* evaluate photoionization rates */
			iso_photo( ipISO , nelem );

			/* evaluate recombination rates */
			iso_radiative_recomb( ipISO , nelem );

			/* Generate Gaussian errors if turned on. */
			if( iso.lgRandErrGen[ipISO] && nzone==0 && !iso.lgErrGenDone[ipISO][nelem] )
			{
				iso_error_generation(ipISO, nelem );
			}

			if( opac.lgRedoStatic )
			{
				if( nelem<=ipHELIUM )
				{
					iso_collapsed_bnl_set( ipISO, nelem );

					//iso_collapsed_bnl_print( ipISO, nelem );

					iso_collapsed_Aul_update( ipISO, nelem );

					iso_collapsed_lifetimes_update( ipISO, nelem );
				}

				iso_cascade( ipISO, nelem );

				iso_radiative_recomb_effective( ipISO, nelem );
			}

			/* evaluate state specific creation and destruction processes,
			 * also define iso.xIonSimple */
			iso_ionize_recombine( ipISO , nelem );
		}
	}

	return;
}

void iso_solve(long ipISO, long nelem)
{
	DEBUG_ENTRY( "iso_solve()" );

	/* do not consider elements that have been turned off */
	if( dense.lgElmtOn[nelem] )
	{
		/* option to force ionization */
		if( dense.lgSetIoniz[nelem] )
		{
			dense.xIonDense[nelem][nelem+1-ipISO] = dense.SetIoniz[nelem][nelem+1-ipISO]*dense.gas_phase[nelem];
			dense.xIonDense[nelem][nelem-ipISO] = dense.SetIoniz[nelem][nelem-ipISO]*dense.gas_phase[nelem];
		}

		/* note that nelem scale is totally on c not physical scale, so 0 is h */
		/* evaluate balance if ionization reaches this high */
		if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
			(dense.IonLow[nelem] <= nelem - ipISO) )
		{
			/* solve for the level populations */
			iso_level( ipISO , nelem );

			/* this contains a bunch of trace statements and HydroRenorm.  
			 * Will eventually come over to iso treatment. */
			if( ipISO == ipH_LIKE )
				HydroLevel(nelem);
		}
		else
		{
			iso.xIonSimple[ipISO][nelem] = 0.;

			/* zero it out since no population*/
			StatesElemNEW[nelem][nelem-ipISO][0].Pop = 0.;
			for( long ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop = 0.;
				for( long ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
						continue;

					/* population of lower level rel to ion, corrected for stim em */
					Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc =  0.;
				}
			}
		}

		/* now evaluate departure coefficients */
		iso_departure_coefficients( ipISO, nelem );

		ASSERT( Transitions[ipISO][nelem][iso.nLyaLevel[ipISO]][0].Lo->Pop == StatesElemNEW[nelem][nelem-ipISO][0].Pop );
	}

	return;
}

void IonHydro( void )
{
	double coltot,
		gamtot,
		sum, 
		error;
	bool lgH_chem_conv;
	int loop_H_chem;
	double solomon_assumed;

	DEBUG_ENTRY( "IonHydro()" );

	/* ============================================================================== */
	/* rest is for hydrogen only */

	/* this block appears redundant and could be removed? */
	/* >> 02 nov 21 rjrw -- xIonDense is used in hmole but only set in ion_solver, 
	 * so we need this at least for the first iteration. */

	/* do molecular balance 
	 * hmovh1 will be ratio of molecular to atomic hydrogen
	 * HIonFrac is fraction of H that is ionized, ratio of ion to atom */

	lgH_chem_conv = false;
	loop_H_chem = 0;
	while( loop_H_chem < 5 && !lgH_chem_conv )
	{
		/* save Solomon rate that was assumed in the above - want to make sure that assumed
		 * Solomon rate is stable, and does not change after call to large H2 molecule */
		if( hmi.H2_rate_destroy < SMALLFLOAT )
			solomon_assumed = 1.;
		else
			solomon_assumed = hmi.H2_Solomon_dissoc_rate_used_H2g/hmi.H2_rate_destroy;

		/* now do chem, this will reset hmi.H2_Solomon_dissoc_rate_used */
		hmole();

		/* now do H2 line RT and level populations */
		/* the large H2 molecule */
		H2_RTMake();
		H2_LevelPops();

		/* check that Solomon rate from big molecule is equal to
		 * rate used in H chem */
		lgH_chem_conv = true;
		/* check whether Solomon rate has changed */
		if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			/* use H2g rather than total for solomon rate */
			if( fabs( solomon_assumed - hmi.H2_Solomon_dissoc_rate_BigH2_H2g/SDIV(hmi.H2_rate_destroy) ) > 
				conv.EdenErrorAllowed/5.)
			{
				lgH_chem_conv = false;
			}
		}
		++loop_H_chem;
	}

	{
		/*@-redef@*/
		/* often the H- route is the most efficient formation mechanism for H2,
		 * will be through rate called ratach
		 * this debug print statement is to trace h2 oscillations */
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if(DEBUG_LOC )
		{
			fprintf(ioQQQ,"DEBUG  \t%.2f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
				fnzone,
				hmi.H2_total ,
				hmi.Hmolec[ipMH2g],
				dense.xIonDense[ipHYDROGEN][0],
				dense.xIonDense[ipHYDROGEN][1],
				hmi.H2_Solomon_dissoc_rate_used_H2g,
				hmi.H2_Solomon_dissoc_rate_BD96_H2g,
				hmi.H2_Solomon_dissoc_rate_TH85_H2g);
		}
	}
	/* >>chng 01 may 09, add option to force abundance, with element name ioniz cmmnd */
	if( dense.lgSetIoniz[ipHYDROGEN] )
	{
		dense.xIonDense[ipHYDROGEN][1] = dense.SetIoniz[ipHYDROGEN][1]*dense.gas_phase[ipHYDROGEN];
		dense.xIonDense[ipHYDROGEN][0] = dense.SetIoniz[ipHYDROGEN][0]*dense.gas_phase[ipHYDROGEN];

		/* initialize ground state pop too */
		StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop = dense.xIonDense[ipHYDROGEN][0];
	}
	else
	{
		/* 
		 * >> chng 03 jan 15 rjrw:- terms are now in ion_solver, to allow for
		 * molecular sources and sinks of H and H+.  ion_solver renormalizes
		 * to keep the total H abundance correct -- only the molecular
		 * network is allowed to change this. 
		 */
		ion_solver( ipHYDROGEN , false );
	}

	/* confirm that species still add up correctly */
	/* this exposes a "leak" that occurs somewhere, almost certainly in hmole
	 * if DEBUG_LOC is set true in the following there will be comments printing
	 * due to a loss of some protons */

	sum = dense.xIonDense[ipHYDROGEN][0] + dense.xIonDense[ipHYDROGEN][1];
	for(long mol=0;mol<N_H_MOLEC;mol++) 
	{
		sum += hmi.Hmolec[mol]*hmi.nProton[mol];
	}
	/* do not double count H0 and H+ */
	sum -=  hmi.Hmolec[ipMH]+hmi.Hmolec[ipMHp];

	/* >>chng 06 jun 30, remove H in CO network, as much as a few percent can be
	 * in the form of OH, H2O, etc 
	 * >>chng 06 jul 03, undo this change so function is as before - need to
	 * think about co protons */
	error = ( dense.gas_phase[ipHYDROGEN] - sum )/dense.gas_phase[ipHYDROGEN];
	/*error = ( dense.gas_phase[ipHYDROGEN] - sum - dense.H_sum_in_CO)/dense.gas_phase[ipHYDROGEN];*/
	{
		/*@-redef@*/
		/* often the H- route is the most efficient formation mechanism for H2,
		 * will be through rate called ratach
		 * this debug print statement is to trace h2 oscillations */
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if(DEBUG_LOC && (fabs(error) > 1e-4) )
			fprintf(ioQQQ,"PROBLEM hydrogenic zone %li hden %.4e, sum %.4e (h-s)/h %.3e \n", nzone, dense.gas_phase[ipHYDROGEN] , sum , 
				error );
	}

	fixit(); /* this is called in HydroLevel above, is it needed in both places? */
	/* >>hcng 05 mar 24,
	 * renormalize the populations and emission of H atom to agree with chemistry */
	HydroRenorm();

	/* remember the ratio of pops of 2p to 1s for possible printout in prtComment
	 * and to obtain Lya excitation temps.  the pop of ground is not defined if
	 * NO ionization at all since these pops are relative to ion */
	/* >>chng 99 jun 03, added MAX2 to protect against totally neutral gas */
	if( StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop/MAX2(SMALLDOUBLE,StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop) > 0.1 &&
		StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop > SMALLDOUBLE )
	{
		hydro.lgHiPop2 = true;
		hydro.pop2mx = (realnum)MAX2(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop/StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop,
		  hydro.pop2mx);
	}

	gamtot = iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s] + secondaries.csupra[ipHYDROGEN][0];

	coltot = iso.ColIoniz[ipH_LIKE][ipHYDROGEN][ipH1s] + 
	  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Coll.ColUL*4.*
	  iso.Boltzmann[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s];

	/* if ground state destruction rate is significant, recall different dest procceses */
	if( iso.RateLevel2Cont[ipH_LIKE][ipHYDROGEN][ipH1s] > SMALLFLOAT )
	{
		hydro.H_ion_frac_photo = 
			(realnum)(iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]/iso.RateLevel2Cont[ipH_LIKE][ipHYDROGEN][ipH1s] );

		/* fraction of ionizations of H from ground, due to thermal collisions */
		hydro.H_ion_frac_collis = 
			(realnum)(iso.ColIoniz[ipH_LIKE][ipHYDROGEN][ipH1s]*dense.eden/iso.RateLevel2Cont[ipH_LIKE][ipHYDROGEN][ipH1s]);

		/* this flag is used in ConvBase to decide whether we
		 * really need to converge the secondary ionization rates */
		secondaries.sec2total = 
			(realnum)(secondaries.csupra[ipHYDROGEN][0] / iso.RateLevel2Cont[ipH_LIKE][ipHYDROGEN][ipH1s]);

		/* frac of ionizations due to ct */
		atmdat.HIonFrac = atmdat.HCharExcIonTotal / iso.RateLevel2Cont[ipH_LIKE][ipHYDROGEN][ipH1s];
	}
	else
	{
		hydro.H_ion_frac_collis = 0.;
		hydro.H_ion_frac_photo = 0.;
		secondaries.sec2total = 0.;
		atmdat.HIonFrac = 0.;
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "       Hydrogenic return %.2f ",fnzone);
		fprintf(ioQQQ,"H0:%.3e ", dense.xIonDense[ipHYDROGEN][0]);
		fprintf(ioQQQ,"H+:%.3e ", dense.xIonDense[ipHYDROGEN][1]);
		fprintf(ioQQQ,"H2:%.3e ", hmi.H2_total);
		fprintf(ioQQQ,"H-:%.3e ", hmi.Hmolec[ipMHm]);
		fprintf(ioQQQ,"ne:%.3e ", dense.eden);
		fprintf( ioQQQ, " REC, COL, GAMT= ");
		/* recomb rate coef, cm^3 s-1 */
		fprintf(ioQQQ,"%.2e ", iso.RadRec_effec[ipH_LIKE][ipHYDROGEN] );
		fprintf(ioQQQ,"%.2e ", coltot);
		fprintf(ioQQQ,"%.2e ", gamtot);
		fprintf( ioQQQ, " CSUP=");
		PrintE82( ioQQQ, secondaries.csupra[ipHYDROGEN][0]);
		fprintf( ioQQQ, "\n");
	}

	// reset the pump rates
	if( hydro.xLymanPumpingScaleFactor!=1.f )
	{
		for( long ipHi=1; ipHi < iso.numLevels_local[ipH_LIKE][ipHYDROGEN]; ++ipHi )
		{
			Transitions[ipH_LIKE][ipHYDROGEN][ipHi][0].Emis->pump = PumpSave[ipHi];
		}
	}

	return;
}

/* HydroRenorm - renormalize H so that it agrees with the chemistry */
void HydroRenorm( void )
{
	double sum_atom_iso , renorm;

	DEBUG_ENTRY( "HydroRenorm()" );

	/*>>chng 04 mar 23, add this renorm */
	/* renormalize the state specific populations, so that they
	 * add up to the results that came from ion_solver 
	 * units at first is sum div by H+ density, since Pop2Ion defn this way */
	sum_atom_iso = 0.;
	/* >> chng 06 aug 31, from numLevels_max to _local */
	for( long level=ipH1s; level < iso.numLevels_local[ipH_LIKE][ipHYDROGEN]; level++ )
	{
		/* cm-3 - total population in iso solved model */
		sum_atom_iso += StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][level].Pop;
	}

	/* >>chng 04 may 25, humunculus sum_atom_iso is zero */
	if( sum_atom_iso > SMALLFLOAT )
	{
		renorm = dense.xIonDense[ipHYDROGEN][0] / sum_atom_iso;
	}
	else
	{
		renorm = 0.;
	}
	ASSERT( renorm < BIGFLOAT );
	/*fprintf(ioQQQ,"DEBUG renorm\t%.2f\t%.3e\n",fnzone, renorm);*/
	/* renormalize populations from iso-model atom so that they agree with ion solver & chemistry */
	StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop *= renorm;
	/*fprintf(ioQQQ,"DEBUG h \t%.3e hydrogenic renorm %.3e\n", 
		StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop ,
		StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop/renorm );*/

	/* >> chng 06 aug 31, from numLevels_max to _local */
	for( long ipHi=ipH2s; ipHi < iso.numLevels_local[ipH_LIKE][ipHYDROGEN]; ipHi++ )
	{
		StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipHi].Pop *= renorm;

		for( long ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			if( Transitions[ipH_LIKE][ipHYDROGEN][ipHi][ipLo].Emis->Aul <= iso.SmallA )
				continue;

			/* population of lower level rel to ion, corrected for stim em */
			Transitions[ipH_LIKE][ipHYDROGEN][ipHi][ipLo].Emis->PopOpc *= renorm;
		}
	}

	return;
}

void iso_departure_coefficients( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_departure_coefficients()" );
		
	for( long level=0; level < iso.numLevels_local[ipISO][nelem]; level++ )
	{
		double denom = dense.xIonDense[nelem][nelem+1-ipISO]*
			iso.PopLTE[ipISO][nelem][level]*dense.eden;

		if( iso.PopLTE[ipISO][nelem][level] > 0. && denom > SMALLFLOAT )
			iso.DepartCoef[ipISO][nelem][level] = safe_div( 
			  StatesElemNEW[nelem][nelem-ipISO][level].Pop, denom );
		else
			iso.DepartCoef[ipISO][nelem][level] = 0.;
	}

	for( long level=iso.numLevels_local[ipISO][nelem]; level < iso.numLevels_max[ipISO][nelem]; level++ )
		iso.DepartCoef[ipISO][nelem][level] = 0.;

	return;
}

/*iso_prt_pops print out iso sequence populations or departure coefficients */
void iso_prt_pops( long ipISO, long nelem, bool lgPrtDeparCoef )
{
	long int in, il, is, i, ipLo, nResolved, ipFirstCollapsed=LONG_MIN;
	char chPrtType[2][12]={"populations","departure"};
	/* first dimension is multiplicity */
	char chSpin[3][9]= {"singlets", "doublets", "triplets"};

#define ITEM_TO_PRINT(A_)	( lgPrtDeparCoef ? iso.DepartCoef[ipISO][nelem][A_] : StatesElemNEW[nelem][nelem-ipISO][A_].Pop )

	DEBUG_ENTRY( "iso_prt_pops()" );

	ASSERT( ipISO < NISO );

	for( is = 1; is<=3; ++is)
	{
		if( ipISO == ipH_LIKE && is != 2 )
			continue;
		else if( ipISO == ipHE_LIKE && is != 1 && is != 3 )
			continue;

		ipFirstCollapsed= iso.numLevels_local[ipISO][nelem]-iso.nCollapsed_local[ipISO][nelem];
		nResolved = StatesElemNEW[nelem][nelem-ipISO][ipFirstCollapsed-1].n;
		ASSERT( nResolved == iso.n_HighestResolved_local[ipISO][nelem] );
		ASSERT(nResolved > 0 );

		/* give element number and spin */
		fprintf(ioQQQ," %s %s  %s %s\n",
			iso.chISO[ipISO],
			elementnames.chElementSym[nelem],
			chSpin[is-1],
			chPrtType[lgPrtDeparCoef]);

		/* header with the l states */
		fprintf(ioQQQ," n\\l=>    ");
		for( i =0; i < nResolved; ++i)
		{
			fprintf(ioQQQ,"%2ld         ",i);
		}
		fprintf(ioQQQ,"\n");

		/* loop over prin quant numbers, one per line, with l across */
		for( in = 1; in <= nResolved; ++in)
		{
			if( is==3 && in==1 )
				continue;

			fprintf(ioQQQ," %2ld      ",in);

			for( il = 0; il < in; ++il)
			{
				if( ipISO==ipHE_LIKE && (in==2) && (il==1) && (is==3) )
				{
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipHe2p3P0) );
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipHe2p3P1) );
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipHe2p3P2) );
				}
				else
				{
					ipLo = iso.QuantumNumbers2Index[ipISO][nelem][in][il][is];
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipLo) );
				}
			}
			fprintf(ioQQQ,"\n");
		}
	}
	/* above loop was over spin, now do collapsed levels, no spin or ang momen */
	for( il = ipFirstCollapsed; il < iso.numLevels_local[ipISO][nelem]; ++il)
	{
		in = StatesElemNEW[nelem][nelem-ipISO][il].n;
		/* prin quan number of collapsed levels */
		fprintf(ioQQQ," %2ld      ",in);
		fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(il) );
		fprintf(ioQQQ,"\n");
	}
	return;
}

/* routine to save table needed for AGN3 - collision strengths of HeI */
void AGN_He1_CS( FILE *ioPun )
{

	long int i;

	/* list of temperatures where cs will be printed */
	const int NTE = 5;
	double TeList[NTE] = {6000.,10000.,15000.,20000.,25000.};
	double TempSave;

	DEBUG_ENTRY( "AGN_He1_CS()" );

	/* put on a header */
	fprintf(ioPun, "Te\t2 3s 33s\n");

	/* Restore the original temp when this routine done.	*/
	TempSave = phycon.te;

	for( i=0; i<NTE; ++i )
	{
		TempChange(TeList[i] , false);

		fprintf(ioPun , "%.0f\t", 
			TeList[i] );
		fprintf(ioPun , "%.2f\t", 
			HeCSInterp( 1 , ipHe3s3S , ipHe2s3S, ipELECTRON ) );
		fprintf(ioPun , "%.2f\t", 
			HeCSInterp( 1 , ipHe3p3P , ipHe2s3S, ipELECTRON ) );
		fprintf(ioPun , "%.2f\t", 
			HeCSInterp( 1 , ipHe3d3D , ipHe2s3S, ipELECTRON ) );
		fprintf(ioPun , "%.3f\t", 
			HeCSInterp( 1 , ipHe3d1D , ipHe2s3S, ipELECTRON ) );
		/*fprintf(ioPun , "%.1f\t%.1f\t%.1f\n", */
		fprintf(ioPun , "%.1f\n", 
			HeCSInterp( 1 , ipHe2p3P0 , ipHe2s3S, ipELECTRON ) +
			HeCSInterp( 1 , ipHe2p3P1 , ipHe2s3S, ipELECTRON ) +
			HeCSInterp( 1 , ipHe2p3P2 , ipHe2s3S, ipELECTRON ));
	}

	/* no need to force update since didn't do above	*/
	TempChange(TempSave , false);
	return;
}
