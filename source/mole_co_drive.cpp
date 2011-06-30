/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lgMolecAver average old and new molecular equilibrium balance from CO_solve */
/*CO_drive - public routine, calls CO_solve to converge molecules */
#include "cddefines.h"
#include "taulines.h"
#include "dense.h"
#include "hmi.h"
#include "conv.h"
#include "phycon.h"
#include "trace.h"
#include "thermal.h"
#include "called.h"
#include "hcmap.h"
#include "coolheavy.h"
#include "mole.h"
#include "mole_co_priv.h"

/* says whether the co network is currently set to zero */
static bool lgMoleZeroed=true;

/* the limit to H2/H where we will solve for CO */
static double h2lim;

/* this is the relative change in OH, which is used as a convergence
 * criteria in lgMolecAver */
static const double COTOLER_MOLAV = 0.10;

/* flag to control debug statement, if 0 then none, just loops when 1, more when 2 */
static const int CODEBUG = 0;
/* this is the maximum number of loops CO_drive will go over while
 * trying to converge the molecular network */
static const int LUPMAX_CODRIV = 100;

/*lgMolecAver average old and new molecular equilibrium balance from CO_solve,
 * returns true if converged */
STATIC bool lgMolecAver(
	/* job == "SETUP", set things up during first call,
	 * job == "AVER" take average of arrays */
	const char chJob[10] ,
	char *chReason );

/*CO_drive main driver for heavy molecular equilibrium routines      */
void CO_drive(void)
{
	bool lgNegPop, 
	  lgZerPop,
	  lgFirstTry;
	long int i;
	bool lgConverged;
	/* count how many times through loop */
	long int loop;
	long int nelem , ion;


	/* this will give reason CO not converged */
	char chReason[100];
	/*static bool lgMoleZeroed=true;*/

	DEBUG_ENTRY( "CO_drive()" );

	/* h2lim is smallest ratio of H2/HDEN for which we will
	 * even try to invert heavy element network */
	if( hmi.lgNoH2Mole )
	{
		/* H2 network is turned off, ignore this limit */
		h2lim = 0.;
	}
	else
	{
		/* this limit is important since the co molecular network is first turned
		 * on when this is hit.  It's conservation law will only ever include the
		 * initial O, O+, and C, C+, so must not start before nearly all
		 * C and O is in these forms */
		h2lim = 1e-15;
		/* >>chng 05 jul 15, from 1e-15 to 1e-12, see whether results are stable,
		 * this does include CO in the H+ zone in orion_hii_pdr
		 * a problem is that, during search phase, the first temp is 4000K and the
		 * H2 abundance is larger than it will be at the illuminated face.  try to
		 * avoid turning H2 on too soon */
		h2lim = 1e-12;
	}

	/* do we want to try to calculate the CO? 
	 * >>chng 03 sep 07, add first test, so that we do not turn off
	 * heavy element network if it has ever been turned on
	 * >>chng 04 apr 23, change logic so never done when temp > 20000K */
	/* make series of tests setting co.lgCODoCalc for simplicity */
	co.lgCODoCalc = true;
	/* too hot - don't bother Te > 2e4K */
	if( phycon.te > 20000. )
		co.lgCODoCalc = false;

	/* molecules have been turned off */
	if( co.lgNoCOMole )
		co.lgCODoCalc = false;

	/* a critical element has been turned off */
	if( dense.lgElmtOn[ipCARBON]*dense.lgElmtOn[ipOXYGEN] ==  0 )
		co.lgCODoCalc = false;

	/* >>chng 04 jul 20, do not attempt co if no O or C atoms present,
	 * since co net needs their atomic data */
	if( dense.IonLow[ipCARBON]>0 || dense.IonLow[ipOXYGEN]>0 )
		co.lgCODoCalc = false;

	/* do not do if H2/Htot < h2lim which is set to 1e-12 by default above */
	if( iteration!=co.iteration_co && 
		hmi.H2_total/dense.gas_phase[ipHYDROGEN] < h2lim )
		co.lgCODoCalc = false;

	/* >>chng 06 sep 10, do not call CO network if H+/Htot is greater than 0.5 -
	 * this change suggested by Nick Abel after cooling plasma developed CO
	 * convergence problems when H+/H = 0.97 */
	if( dense.xIonDense[ipHYDROGEN][1] / dense.gas_phase[ipHYDROGEN] > 0.5 )
		co.lgCODoCalc = false;

	if( dense.lgElmtOn[ipSILICON] )
	{
		/* >>chng 04 jun 28, add test on atomic Si, some H II region
		 * models crashed due to matrix failure when co network turned
		 * on with atomic Si at very small abundance */
		if( dense.xIonDense[ipSILICON][0]/dense.gas_phase[ipSILICON] < 1e-15 )
			co.lgCODoCalc = false;
	}

	/* this branch we will not do the CO equilibrium, set some variables that
	 * would be calculated by the chemistry, zero others, and return */
	if( !co.lgCODoCalc )
	{
		/* these are heavy - heavy charge transfer rate that will still be needed
		 * for recombination of Si+, S+, and C+ */
		struct COmole_rate_s *rate;

		mole.xMoleChTrRate[ipSILICON][1][0] = 0.;
		rate = CO_findrate_s("Si+,Fe=>Si,Fe+");
		mole.xMoleChTrRate[ipSILICON][1][0] += (realnum)(rate->rk*rate->rate_species[1]->hevmol);
		rate = CO_findrate_s("Si+,Mg=>Si,Mg+");
		mole.xMoleChTrRate[ipSILICON][1][0] += (realnum)(rate->rk*rate->rate_species[1]->hevmol);

		mole.xMoleChTrRate[ipSULPHUR][1][0] = 0.;
		rate = CO_findrate_s("S+,Fe=>S,Fe+");
		mole.xMoleChTrRate[ipSULPHUR][1][0] += (realnum)(rate->rk*rate->rate_species[1]->hevmol);
		rate = CO_findrate_s("S+,Mg=>S,Mg+");
		mole.xMoleChTrRate[ipSULPHUR][1][0] += (realnum)(rate->rk*rate->rate_species[1]->hevmol);

		mole.xMoleChTrRate[ipCARBON][1][0] = 0.;
		rate = CO_findrate_s("C+,Fe=>C,Fe+");
		mole.xMoleChTrRate[ipCARBON][1][0] += (realnum)(rate->rk*rate->rate_species[1]->hevmol);
		rate = CO_findrate_s("C+,Mg=>C,Mg+");
		mole.xMoleChTrRate[ipCARBON][1][0] += (realnum)(rate->rk*rate->rate_species[1]->hevmol);

#if 0
		mole.xMoleChTrRate[ipSILICON][1][0] = (realnum)(dense.xIonDense[ipMAGNESIUM][0]*
			COmole_rate[ipR_SiP_Mg_Si_MgP].rk + 
			dense.xIonDense[ipIRON][0]*COmole_rate[ipR_SiP_Fe_Si_FeP].rk);

		mole.xMoleChTrRate[ipSULPHUR][1][0] = (realnum)(dense.xIonDense[ipMAGNESIUM][0]*
			COmole_rate[ipR_SP_Mg_S_MgP].rk + 
			dense.xIonDense[ipIRON][0]*COmole_rate[ipR_SP_Fe_S_FeP].rk);

		mole.xMoleChTrRate[ipCARBON][1][0] = (realnum)(dense.xIonDense[ipMAGNESIUM][0]*
			COmole_rate[ipR_CP_Mg_C_MgP].rk + 
			dense.xIonDense[ipIRON][0]*COmole_rate[ipR_CP_Fe_C_FeP].rk);
#endif

		/* do we need to zero out the arrays and vars? */
		if( !lgMoleZeroed )
		{
			lgMoleZeroed = true;
			for( i=0; i<mole.num_comole_calc; ++i )
			{
				COmole[i]->hevmol = 0.;
			}
			/* heating due to CO photodissociation */
			thermal.heating[0][9] = 0.;
			co.CODissHeat = 0.;

			for( nelem=ipLITHIUM; nelem<LIMELM; ++nelem )
			{
				for( ion=0; ion<nelem+2; ++ion )
				{
					mole.source[nelem][ion] = 0.;
					mole.sink[nelem][ion] = 0.;
				}
			}
		}

		if( trace.nTrConvg>=4 )
		{
			/* trace ionization  */
			fprintf( ioQQQ, 
				"    CO_drive4 do not evaluate CO chemistry.\n");
		}
		return;
	}

	CO_update_species_cache();  /* Update densities of species controlled outside the network */

	/* Update the rate constants -- final generic version */
	CO_update_rks();

	/* this flag is used to stop calculation due to negative abundances */
	loop = 0;
	lgNegPop = false;
	lgConverged = lgMolecAver("SETUP",chReason);
	do
	{

		if( trace.nTrConvg>=5 )
		{
			/* trace ionization  */
			fprintf( ioQQQ, 
				"      CO_drive5 CO chemistry loop %li chReason %s.\n" , loop, chReason );
		}

		/*oh_h_o_h2ld = findspecies("OH")->hevmol;*/
		CO_solve(&lgNegPop,&lgZerPop );
		/* this one takes average first, then updates molecular densities in co.hevmol,
		 * returns true if change in OH density is less than macro COTOLER_MOLAV set above*/
		lgConverged = lgMolecAver("AVER",chReason);
		if( CODEBUG )
			fprintf(ioQQQ,"codrivbug %.2f %li Neg?:%c\tZero?:%c\tOH new\t%.3e\tCO\t%.3e\tTe\t%.3e\tH2/H\t%.2e\n",
				fnzone ,
				loop ,
				TorF(lgNegPop) , 
				TorF(lgZerPop) , 
				findspecies("OH")->hevmol ,
				findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]),
				phycon.te,
				hmi.H2_total/dense.gas_phase[ipHYDROGEN]);

		++loop;
	}   while (
			  /* these are a series of conditions to stop this loop - 
			   * has the limit to the number of loops been reached? */
			  (loop < LUPMAX_CODRIV) && !lgConverged  );

	/* say that we have found a solution before */
	lgMoleZeroed = false;

	/* check whether we have converged */
	if( loop == LUPMAX_CODRIV)
	{
		conv.lgConvIoniz = false;
		strcpy( conv.chConvIoniz, "CO con1" );
		if( CODEBUG )
			fprintf(ioQQQ,"CO_drive not converged\n");
	}

	/* this flag says whether this is the first zone we are trying
	 * to converge the molecules - there will be problems during initial
	 * search so do not print anything in this case */
	lgFirstTry = (nzone==co.co_nzone && iteration==co.iteration_co);

	/* did the molecule network have negative pops? */
	if( lgNegPop )
	{
		if( conv.lgSearch && (hmi.H2_total/dense.gas_phase[ipHYDROGEN] <h2lim) &&
			(iteration==1) )
		{
			/* we are in search phase during the first iteration, 
			 * and the H2/H ratio has fallen
			 * substantially below the threshold for solving the CO network.
			 * Turn it off */
			/* >> chng 07 jan 10 rjrw: this was CO_Init(), but the comment suggests
			 * it should really be CO_zero */
			CO_zero();
		}
		else
		{
			if( called.lgTalk && !lgFirstTry )
			{
				conv.lgConvPops = false;
				fprintf( ioQQQ, 
						" PROBLEM CO_drive failed to converge1 due to negative pops, zone=%.2f,  CO/C=%.2e  OH/C=%.2e H2/H=%.2e\n", 
						fnzone, 
						findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]),
						findspecies("OH")->hevmol/SDIV(dense.gas_phase[ipCARBON]),
						hmi.H2_total/dense.gas_phase[ipHYDROGEN]);
				ConvFail( "pops" , "CO" );
			}
		}
	}

	/* this test, hit upper limit to number of iterations */
	else if( loop == LUPMAX_CODRIV )
	{
		/* do not print this if we are in the first zone where molecules are turned on */
		if( called.lgTalk && !lgFirstTry )
		{
			conv.lgConvPops = false;
			fprintf( ioQQQ, 
				" PROBLEM CO_drive failed to converge2 in %li iter, zone=%.2f, CO/C=%.2e negpop?%1c reason:%s\n", 
					 loop,
					 fnzone, 
					 findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]), 
					 TorF(lgNegPop) ,
					 chReason);
			ConvFail( "pops" , "CO" );
		}
	}

	/* make sure we have not used more than all the CO */
	ASSERT(conv.lgSearch || findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]) <= 1.01 ||
		/* this is true when loop did not converge co */
		loop == LUPMAX_CODRIV );
	/*fprintf(ioQQQ,"ratioo o\t%c\t%.2f\t%f\n", TorF(conv.lgSearch),
		fnzone , co.hevmol[ipCO]/dense.gas_phase[ipCARBON] );*/

	/* these are the elements whose converge we check */
	/* this is a self consistency check, for converged solution */
	/* >>chng 04 dec 02, this test is turn  back on - don't know why it was turned off */
	if(0) {
		double sum_ion , sum_mol;
		sum_ion = dense.xIonDense[ipCARBON][0] + dense.xIonDense[ipCARBON][1];
		sum_mol = findspecies("C")->hevmol + findspecies("C+")->hevmol;
		if( fabs(sum_ion-sum_mol)/SDIV(sum_mol) > 1e-2 )
		{
			/*fprintf(ioQQQ,
				"non conservation element %li zone %.2f sum_ion %.3e sum_mol %.3e\n",
				5, fnzone, sum_ion , sum_mol);*/
			conv.lgConvIoniz = false;
			strcpy( conv.chConvIoniz, "CO con2" );
			conv.BadConvIoniz[0] = sum_ion;
			conv.BadConvIoniz[1] = sum_mol;
			if( CODEBUG )
				fprintf(ioQQQ,"CO_drive not converged\n");
		}
		sum_ion = dense.xIonDense[ipOXYGEN][0] + dense.xIonDense[ipOXYGEN][1];
		sum_mol = findspecies("O")->hevmol + findspecies("O+")->hevmol;
		if( fabs(sum_ion-sum_mol)/SDIV(sum_mol) > 1e-2 )
		{
			/*fprintf(ioQQQ,
				"non conservation element %li zone %.2f sum_ion %.3e sum_mol %.3e\n",
				7, fnzone, sum_ion , sum_mol);*/
			conv.lgConvIoniz = false;
			strcpy( conv.chConvIoniz, "CO con3" );
			conv.BadConvIoniz[0] = sum_ion;
			conv.BadConvIoniz[1] = sum_mol;
			if( CODEBUG )
				fprintf(ioQQQ,"CO_drive not converged\n");
		}
	}
	return;
}

STATIC bool lgMolecAver(
	/* job == "SETUP", set things up during first call,
	 * job == "AVER" take average of arrays */
	const char chJob[10] ,
	char *chReason )
{
	long int i;
	bool lgReturn;
	/* this will be used to rescale old saved abundances in constant pressure model */
	static realnum hden_save_prev_call;
	double conv_fac;
	struct chem_element_s *element;

	DEBUG_ENTRY( "lgMolecAver()" );
	const bool DEBUG_MOLECAVER = false;

	/* this will become reason not converged */
	strcpy( chReason , "none given" );

	/* when called with SETUP, just about to enter solver loop */
	if( strcmp( "SETUP" , chJob ) == 0 )
	{
		static realnum hden_save_prev_iter;

		/* >>chng 04 apr 29, use PDR solution, scaled to density of neutral carbon, as first guess*/
		/* CO_Init set this to -2 when code initialized, so negative
		 * number shows very first call in this model */
		/* >>chng 05 jul 15, do not init if we have never evaluated CO in this iteration
		 * and we are below limit where it should be evaluated */
		if( iteration!=co.iteration_co && 
			hmi.H2_total/dense.gas_phase[ipHYDROGEN] < h2lim )
		{

			lgReturn = true;
			return lgReturn;
		}

		else if( co.iteration_co < 0 || lgMoleZeroed )
		{

			/* very first attempt at ever obtaining a solution -
			 * called one time per model since co.iteration_co set negative
			 * when cdInit called */

			/* >>chng 05 jun 24, during map phase do not reset molecules to zero
			 * since we probably have a better estimate right now */
			if( !hcmap.lgMapBeingDone || lgMoleZeroed )
			{
				for( i=0; i < mole.num_comole_calc; i++ )
				{
					COmole[i]->hevmol = dense.xIonDense[ipCARBON][0]*COmole[i]->pdr_mole_co;
					COmole[i]->co_save = COmole[i]->hevmol;
				}
			}

			/* we should have a neutral carbon solution at this point */
			ASSERT( dense.xIonDense[ipCARBON][0]>SMALLFLOAT );

			/* first guess of the elements and charges come from ionization solvers */
			for(i=0;i<mole.num_elements;i++) {
				element = chem_element[i];
				COmole[element->ipMl]->hevmol = dense.xIonDense[element->ipCl][0];
				COmole[element->ipMlP]->hevmol = dense.xIonDense[element->ipCl][1];
			}

			/* set iteration flag */
			co.iteration_co = iteration;
			co.co_nzone = nzone;

			/* for constant pressure models when molecules are reset on second
			 * and higher iterations, total density will be different, so
			 * must take this into account when abundances are reset */
			hden_save_prev_iter = dense.gas_phase[ipHYDROGEN];
			hden_save_prev_call = dense.gas_phase[ipHYDROGEN];

			if( DEBUG_MOLECAVER )
				fprintf(ioQQQ," DEBUG lgMolecAver iteration %li zone %li zeroing since very first call. H2/H=%.2e\n", 
				iteration,nzone,hmi.H2_total/dense.gas_phase[ipHYDROGEN]);
		}
		else if( iteration > co.iteration_co )
		{
			realnum ratio = dense.gas_phase[ipHYDROGEN] / hden_save_prev_iter;

			/* this is first call on new iteration, reset
			 * to first initial abundances from previous iteration */
			for( i=0; i < mole.num_comole_calc; i++ )
			{
				COmole[i]->hevmol = COmole[i]->hev_reinit*ratio;
				COmole[i]->co_save = COmole[i]->hev_reinit*ratio;
			}

			co.iteration_co = iteration;
			co.co_nzone = nzone;

			if( DEBUG_MOLECAVER )
				fprintf(ioQQQ," DEBUG lgMolecAver iteration %li zone %li resetting since first call on new iteration. H2/H=%.2e\n", 
				iteration,
				nzone,
				hmi.H2_total/dense.gas_phase[ipHYDROGEN]);
		}
		else if( iteration == co.iteration_co && nzone==co.co_nzone+1 )
		{
			/* this branch, second zone with solution, so save previous
			 * zone's solution to reset things in next iteration */
			for( i=0; i < mole.num_comole_calc; i++ )
			{
				COmole[i]->hev_reinit = COmole[i]->hevmol;
			}

			co.co_nzone = -2;
			hden_save_prev_iter = dense.gas_phase[ipHYDROGEN];
			hden_save_prev_call = dense.gas_phase[ipHYDROGEN];

			if( DEBUG_MOLECAVER )
				fprintf(ioQQQ,"DEBUG lgMolecAver iteration %li zone %li second zone on new iteration, saving reset.\n", iteration,nzone);
		}

		/* didn't do anything, but say converged */
		lgReturn = true;
	}

	else if( strcmp( "AVER" , chJob ) == 0 )
	{
		/* >>chng 04 jan 22, co.hevmol is current value, we want to save old
		 * value, which is in co.co_save */
		/*realnum oldoh = findspecies("OH")->hevmol;*/
		realnum oldoh = findspecies("OH")->co_save;
		lgReturn = true;

		/* get new numbers - take average of old and new */
		/* >>chng 03 sep 16, only use damper for molecular species,
		 * not ion/atoms */
		for( i=0; i < mole.num_comole_calc; i++ )
		{
			/* parameter to mix old and new,
			 * damper is fraction of old solution to take */
			realnum damper = 0.2f;

			/* >>chng 04 sep 11, correct for den chng in var den models */
			/* the ratio of densities is unity for constant density model, but in const pres or other
			 * models, account for change in entire density scale */
			COmole[i]->co_save *= dense.gas_phase[ipHYDROGEN] / hden_save_prev_call;

			/* if co.co_save is zero then don't take average of it and 
			 * new solution */
			if( COmole[i]->co_save< SMALLFLOAT )
			{
				COmole[i]->co_save = COmole[i]->hevmol;
			}
			else
			{
				/* >>chng 04 jan 24, add logic to check how different old and
				 * new solutions are.  if quite different then take average
				 * in the log */
				double ratio = (double)COmole[i]->hevmol/(double)COmole[i]->co_save;
				ASSERT( COmole[i]->co_save>0. && COmole[i]->hevmol>=0. );
				if( fabs(ratio-1.) < 1.5 ||
					fabs(ratio-1.) > 0.66 )
				{
					/* this branch moderate changes */
					COmole[i]->hevmol = COmole[i]->hevmol*(1.f - damper) + 
						COmole[i]->co_save*damper;
				}
				else
				{
					/* this branch - large changes take average in the log */
					COmole[i]->hevmol = (realnum)pow(10. , 
						(log10(COmole[i]->hevmol) + log10(COmole[i]->co_save))/2. );
				}
				COmole[i]->co_save = COmole[i]->hevmol;
			}
		}
		/* remember the current H density */
		hden_save_prev_call = dense.gas_phase[ipHYDROGEN];

		/* >>chng 05 jul 14, add logic to not be as fine a tolerance when
		 * only trace amounts of molecules are present */
		if( oldoh/SDIV(dense.gas_phase[ipOXYGEN]) < 1e-10  )
			conv_fac = 3.;
		else
			conv_fac = 1.;

		if(fabs(oldoh/SDIV(findspecies("OH")->hevmol)-1.) < COTOLER_MOLAV*conv_fac )
		{
			lgReturn = true;
		}
		else
		{
			/* say we are not converged */
			lgReturn = false;
			sprintf( chReason , "OH change from %.3e to %.3e",
				oldoh , 
				findspecies("OH")->hevmol);
		}
	}
	else
		TotalInsanity();

	/* also not converged if co exceeds c */
	if( findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]) > 1.+FLT_EPSILON )
	{
		lgReturn = false;
		sprintf( chReason , "CO/C >1, value is %e",
			findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]) );
	}

	if( (trace.lgTrace && trace.lgTr_CO_Mole) || DEBUG_MOLECAVER )
	{
		fprintf( ioQQQ, " DEBUG lgMolecAver converged? %c" ,TorF(lgReturn) );
		fprintf( ioQQQ, " CO/C=%9.1e", findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]) );
		fprintf( ioQQQ, " OH/O=%9.1e", findspecies("OH")->hevmol/SDIV(dense.gas_phase[ipOXYGEN]) );
		fprintf( ioQQQ, "\n" );
	}

	return lgReturn;
}
