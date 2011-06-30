/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_solver solve the bi-diagonal matrix for ionization balance */
#include "cddefines.h"
#include "yield.h"
#include "prt.h"
#include "continuum.h"
#include "iso.h"
#include "dynamics.h"
#include "grainvar.h"
#include "hmi.h"
#include "mole.h"
#include "thermal.h"
#include "thirdparty.h"
#include "conv.h"
#include "secondaries.h"
#include "phycon.h"
#include "atmdat.h"
#include "heavy.h"
#include "elementnames.h"
#include "dense.h"
#include "radius.h"
#include "ionbal.h"
#include "taulines.h"
#include "trace.h"

bool lgOH_ChargeTransferDominant(  void );

STATIC bool lgTrivialSolution( long nelem, double abund_total );

STATIC void find_solution( long nelem, long ion_range, valarray<double> &xmat, valarray<double> &source);

STATIC void fill_array( long int nelem, long ion_range, valarray<double> &xmat, valarray<double> &source, valarray<double> &auger, double *abund_total );

STATIC void combine_arrays( valarray<double> &xmat, const valarray<double> &xmat1, const valarray<double> &xmat2, long ion_range1, long ion_range2 );

STATIC double get_total_abundance_ions( long int nelem );

STATIC bool lgHomogeneousSource( long nelem, long ion_low, long ion_range, valarray<double> &xmat, valarray<double> &source, double abund_total );

STATIC void store_new_densities( long nelem, long ion_range, long ion_low, double *source, long lgHomogeneous, double abund_total, bool *lgNegPop ) ;

STATIC void PrintRates( long nelem, bool lgNegPop, double abund_total, valarray<double> &auger, bool lgPrintIt );

void solveions(double *ion, double *rec, double *snk, double *src,
			   long int nlev, long int nmax);

	/* this will be used to address the 2d arrays */
#	undef MAT
#	define MAT(M_,I_,J_)	((M_)[(I_)*(ion_range)+(J_)])

#	undef MAT1
#	define MAT1(M_,I_,J_)	((M_)[(I_)*(ion_range1)+(J_)])

#	undef MAT2
#	define MAT2(M_,I_,J_)	((M_)[(I_)*(ion_range2)+(J_)])

void ion_solver( long int nelem, bool lgPrintIt )
{
	double abund_total;
	long ion_low = dense.IonLow[nelem];
	bool lgNegPop = false;

	DEBUG_ENTRY( "ion_solver()" );

	abund_total = get_total_abundance_ions( nelem );

	if( lgTrivialSolution( nelem, abund_total ) )
	   return;

	long ion_range = dense.IonHigh[nelem]-dense.IonLow[nelem]+1;

	valarray<double> xmat(ion_range*ion_range);
	valarray<double> source(ion_range);
	valarray<double> auger(LIMELM+1);

	// fill xmat and source with appropriate terms
	fill_array( nelem, ion_range, xmat, source, auger, &abund_total );

	// decide if matrix is homogeneous
	bool lgHomogeneous = lgHomogeneousSource( nelem, ion_low, ion_range, xmat, source, abund_total );

	// Now find the solution
	find_solution( nelem, ion_range, xmat, source);

	// save the results in the global density variables
	store_new_densities( nelem, ion_range, ion_low, &source[0], lgHomogeneous, abund_total, &lgNegPop );

	if( prt.lgPrtArry[nelem] || lgPrintIt )
		PrintRates( nelem, lgNegPop, abund_total, auger, lgPrintIt );
	
	return;
}

// ion_solver is overloaded - this version looks for a simultaneous solution to the ionization balance of two elements 
void ion_solver( long int nelem1, long int nelem2, bool lgPrintIt)
{
	bool lgHomogeneous = true;
	bool lgNegPop1 = false;
	bool lgNegPop2 = false;

	DEBUG_ENTRY( "ion_solver()" );
	
	long ion_low1 = dense.IonLow[nelem1];
	long ion_low2 = dense.IonLow[nelem2];
	long ion_range1 = dense.IonHigh[nelem1]-dense.IonLow[nelem1]+1;
	long ion_range2 = dense.IonHigh[nelem2]-dense.IonLow[nelem2]+1;
	long ion_range = ion_range1 + ion_range2;
	double abund_total1 = get_total_abundance_ions( nelem1 );
	double abund_total2 = get_total_abundance_ions( nelem2 );
	valarray<double> xmat(ion_range*ion_range);
	valarray<double> xmat1(ion_range1*ion_range1);
	valarray<double> xmat2(ion_range2*ion_range2);
	valarray<double> source(ion_range);
	valarray<double> auger1(LIMELM+1);
	valarray<double> auger2(LIMELM+1);	

#define ENABLE_SIMULTANEOUS_SOLUTION 0
	
	if( lgTrivialSolution( nelem1, abund_total1 ) )
	{
		// nelem2 has only an ancillary role in solving nelem1.  We do not need to solve nelem2 balance if nelem1 is trivial
		return;
	}
	else if( lgTrivialSolution( nelem2, abund_total2 ) || !ENABLE_SIMULTANEOUS_SOLUTION )
	{
		// if nelem2 has a trivial solution, we still need to solve nelem1, just do it by itself
		ion_solver( nelem1, lgPrintIt );
		return;
	}
					
	/* don't do simultaneous if one or both does not include the neutral stage. */
	if( dense.IonLow[nelem1]>0 || dense.IonLow[nelem2]>0 || nelem1!=ipOXYGEN || nelem2!=ipHYDROGEN )
	{
		fprintf( ioQQQ, "This routine is currently intended to do only O and H when both have significant neutral fractions.\n" );
		fprintf( ioQQQ, "It should be generalized further for other cases.  Exiting.  Sorry.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	ASSERT( nelem1 != nelem2 );
	
	// fill active element's array 
	fill_array( nelem1, ion_range1, xmat1, source, auger1, &abund_total1 );
	// fill passive element's array 
	fill_array( nelem2, ion_range2, xmat2, source, auger2, &abund_total2 );
		
	combine_arrays( xmat, xmat1, xmat2, ion_range1, ion_range2 );
	
	// force homogeneity in simultaneous solutions
	for( long i=0; i< ion_range; ++i )
		source[i] = 0.;
		
	// replace neutral stages with abundance total
	for( long i=0; i< ion_range1; ++i )
		MAT( xmat, i, 0 ) = 1.;
		
	for( long i=ion_range1; i< ion_range; ++i )
		MAT( xmat, i, ion_range1 ) = 1.;
		
	source[0] = dense.xIonDense[nelem1][0] + dense.xIonDense[nelem1][1];
	source[ion_range1] = dense.xIonDense[nelem2][0] + dense.xIonDense[nelem2][1];
		
	// declare homogeneous.
	lgHomogeneous = true;
	
	// Now find the solution
	find_solution( nelem1, ion_range, xmat, source);

	// save the results in the global density variables
	store_new_densities( nelem1, ion_range1, ion_low1, &source[0], lgHomogeneous, abund_total1, &lgNegPop1 );
	store_new_densities( nelem2, ion_range2, ion_low2, &source[ion_range1], lgHomogeneous, abund_total2, &lgNegPop2 );
	ASSERT( lgNegPop2 == false );
	
	if( prt.lgPrtArry[nelem1] || lgPrintIt )
		PrintRates( nelem1, lgNegPop1, abund_total1, auger1, lgPrintIt );
	
	return;
}

STATIC bool lgTrivialSolution( long nelem, double abund_total )
{
	/* return if IonHigh is zero, since no ionization at all */
	if( dense.IonHigh[nelem] == 0 )
	{
		/* set the atom to the total gas phase abundance */
		dense.xIonDense[nelem][0] = (realnum)abund_total;
		return true;
	}
	else if( dense.lgSetIoniz[nelem] )
	{
		/* option to force ionization distribution with element name ioniz */
		for( long ion=0; ion<nelem+2; ++ion )
			dense.xIonDense[nelem][ion] = dense.SetIoniz[nelem][ion]*(realnum)abund_total;
		return true;
	}
	else
		return false;
}

STATIC void find_solution( long nelem, long ion_range, valarray<double> &xmat, valarray<double> &source)
{
	int32 nerror;
	valarray<double> xmatsave(ion_range*ion_range);
	valarray<double> sourcesave(ion_range);
	valarray<int32> ipiv(ion_range);

	DEBUG_ENTRY("find_solution()");
	
	// save copy of xmat before continuing.
	for( long i=0; i< ion_range; ++i )
	{
		sourcesave[i] = source[i];
		for( long j=0; j< ion_range; ++j )
		{
			MAT( xmatsave, i, j ) = MAT( xmat, i, j );
		}
	}
	
	
	if (0)
	{
		// this is just a placeholder for the solver on newmole
		
	}
	else
	{
		/* this is the default solver - now get new solution */
		nerror = 0;
		/* Use general matrix solver */
		getrf_wrapper(ion_range, ion_range, &xmat[0], ion_range, &ipiv[0], &nerror);
		if( nerror != 0 )
		{
			fprintf( ioQQQ, 
					" DISASTER ion_solver: dgetrf finds singular or ill-conditioned matrix nelem=%li %s ion_range=%li, nConv %li IonLow %li IonHi %li\n",
					nelem , 
					elementnames.chElementSym[nelem],
					ion_range,
					conv.nTotalIoniz ,
                                dense.IonLow[nelem], dense.IonHigh[nelem]);
			fprintf( ioQQQ, " xmat follows\n");
			for( long i=0; i<ion_range; ++i )
			{
				for( long j=0;j<ion_range;j++ )
				{
					fprintf(ioQQQ,"%e\t",MAT(xmatsave,j,i));
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ,"source follows\n");
			for( long i=0; i<ion_range;i++ )
			{
				fprintf(ioQQQ,"%e\t",sourcesave[i]);
			}
			fprintf(ioQQQ,"\n");
			cdEXIT(EXIT_FAILURE);
		}
		getrs_wrapper('N', ion_range, 1, &xmat[0], ion_range, &ipiv[0], &source[0], ion_range, &nerror);
		if( nerror != 0 )
		{
			fprintf( ioQQQ, " DISASTER ion_solver: dgetrs finds singular or ill-conditioned matrix nelem=%li ionrange=%li\n",
					nelem , ion_range );
			cdEXIT(EXIT_FAILURE);
		}
	}
	
	{
		/* this is to debug following failed assert */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nelem == ipHYDROGEN )
		{
			fprintf(ioQQQ,"debuggg ion_solver1 %ld\t%.2f\t%.4e\t%.4e\tIon\t%.3e\tRec\t%.3e\n", 
					nelem,
					fnzone,
					phycon.te,
					dense.eden,
					ionbal.RateIonizTot(nelem,0) , 
					ionbal.RateRecomTot[nelem][0]);
			fprintf(ioQQQ," Msrc %.3e %.3e\n", mole.source[nelem][0], mole.source[nelem][1]);
			fprintf(ioQQQ," Msnk %.3e %.3e\n", mole.sink[nelem][0], mole.sink[nelem][1]);
			fprintf(ioQQQ," Poprat %.3e nomol %.3e\n",source[1]/source[0],
					ionbal.RateIonizTot(nelem,0)/ionbal.RateRecomTot[nelem][0]);
		}
	}
	
	for( long i=0; i<ion_range; i++ )
	{
		ASSERT( !isnan( source[i] ) );
		ASSERT( source[i] < MAX_DENSITY );
	}
		
	return;
}

bool lgOH_ChargeTransferDominant(  void )
{
#define THRESHOLD 0.75

	bool lgDominant = (
		(atmdat.HCharExcRecTo[ipOXYGEN][0]*dense.xIonDense[ipOXYGEN][1]/(ionbal.RateIonizTot(ipHYDROGEN,0) + atmdat.HCharExcIonTotal + mole.sink[ipHYDROGEN][0]) > THRESHOLD ) ||
		(atmdat.HCharExcIonOf[ipOXYGEN][0]*dense.xIonDense[ipOXYGEN][0]/(ionbal.RateRecomTot[ipHYDROGEN][0] + atmdat.HCharExcRecTotal ) > THRESHOLD ) );

	return lgDominant;
}

STATIC void combine_arrays( valarray<double> &xmat, const valarray<double> &xmat1, const valarray<double> &xmat2, long ion_range1, long ion_range2 )
{
	DEBUG_ENTRY( "combine_arrays()" );

	long ion_range = ion_range1 + ion_range2;

	for( long i=0; i<ion_range1; i++ )
		for( long j=0; j<ion_range1; j++ )
			MAT( xmat, i, j ) = MAT1( xmat1, i, j );

	for( long i=0; i<ion_range2; i++ )
		for( long j=0; j<ion_range2; j++ )
			MAT( xmat, i+ion_range1, j+ion_range1 ) = MAT2( xmat2, i, j );
	
#if 0
	bool lgNoDice = false;
	for( long i=0; i<ion_range1; i++ )
	{
		if( !fp_equal( MAT( xmat, i, 0), 1.0, 1 ) )
		{
			lgNoDice = true;
			break;
		}
	}
	
	if( !lgNoDice )
	{
		for( long i=ion_range1; i<ion_range; i++ )
			MAT( xmat, i, 0 ) = 1.0;
	}
#endif

	return;
}

STATIC void store_new_densities( long nelem, long ion_range, long ion_low, double *source, long lgHomogeneous, double abund_total, bool *lgNegPop ) 
{
	double renormnew;
	double dennew;

	DEBUG_ENTRY( "store_new_densities()" );

	ASSERT( nelem >= 0 );
	ASSERT( nelem < LIMELM );
	ASSERT( ion_range <= nelem + 2 );
	ASSERT( ion_low >= 0 );
	ASSERT( ion_low <= nelem + 1 );

#if 0
#	define RJRW 0
	if( RJRW && 0 )
	{
		/* verify that the rates are sensible */
		double test;
		for(long i=0; i<ion_range; i++) {
			test = 0.;
			for(long j=0; j<ion_range; j++) {
				test = test+source[j]*MAT(xmatsave,j,i);
			}
			fprintf(ioQQQ,"%e\t",test);
		}
		fprintf(ioQQQ,"\n");

		test = 0.;
		fprintf(ioQQQ," ion %li abundance %.3e\n",nelem,abund_total);
		for( long ion=dense.IonLow[nelem]; ion < dense.IonHigh[nelem]; ion++ )
		{
			if( ionbal.RateRecomTot[nelem][ion] != 0 && source[ion-ion_low] != 0 )
				fprintf(ioQQQ," %li %.3e %.3e : %.3e\n",ion,source[ion-ion_low],
								source[ion-ion_low+1]/source[ion-ion_low],
								ionbal.RateIonizTot(nelem,ion)/ionbal.RateRecomTot[nelem][ion]);
			else
				fprintf(ioQQQ," %li %.3e [One ratio infinity]\n",ion,source[ion-ion_low]);
			test += source[ion-ion_low];
		}
		test += source[dense.IonHigh[nelem]-ion_low];
	}
#endif

	/* 
	 * >> chng 03 jan 15 rjrw:- terms are now included for
	 * molecular sources and sinks of H and H+.
	 *
	 * When the network is not in equilibrium, this will lead to a
	 * change in the derived abundance of H and H+ when after the
	 * matrix solution -- the difference between `renorm' and 1. is a
	 * measure of the quality of the solution (it will be 1. if the
	 * rate of transfer into Ho/H+ balances the rate of transfer
	 * out, for the consistent relative abundances).
	 *
	 * We therefore renormalize to keep the total H abundance
	 * correct -- only the molecular network is allowed to change
	 * this.
	 *
	 * To do this, only the ion abundances are corrected, as the
	 * molecular abundances may depend on several different
	 * conserved species.
	 *
	 */
	if( lgHomogeneous )
	{
		dense.xIonDense[nelem][ion_low] = (realnum)abund_total;
		for( long i=1;i < ion_range; i++ )
		{
			dense.xIonDense[nelem][i+ion_low] = 0.;
		}
	}

	renormnew = 1.;
	/* >>chng 06 mar 17, comment out test on old full depth - keep old solution if overrun scale */
	if(iteration > dynamics.n_initial_relax+1 && dynamics.lgAdvection &&
		 dynamics.Rate != 0.0 && 
			nelem == ipHYDROGEN && hmi.lgNoH2Mole)
	{
		/* The normalization out of the matrix solution is correct and
		 * should be retained if: dynamics is on and the total
		 * abundance of HI & H+ isn't being controlled by the
		 * molecular network */
		renormnew = 1.;
	}
	else
	{
		dennew = 0.;
		double sum_dense = 0.;

		/* find total population to renorm - also here check that negative pops did not occur */
		for( long i=0;i < ion_range; i++ )
		{
			long ion = i+ion_low;
			sum_dense += dense.xIonDense[nelem][ion];
			dennew += source[i];
		} 

		if( dennew > 0.)
		{
			renormnew = sum_dense / dennew;
			/** \todo	2	renorm should == 1 when the molecules and
			 * ionization are in equilibrium.  Should monitor
			 * this figure of merit in calling routine.
			 * */
		}
		else
		{
			renormnew = 1.;
		}
	}
	/* check not negative, should be +ve, can be zero if species has become totally molecular.
	 * this happens for hydrogen if no cosmic rays, or cr ion set very low */
	if( renormnew < 0)
	{
		fprintf(ioQQQ,"PROBLEM impossible value of renorm \n");
	}
	ASSERT( renormnew>=0 );

	/* source[i] contains new solution for ionization populations
	 * save resulting abundances into main ionization density array, 
	 * while checking whether any negative level populations occurred */
	*lgNegPop = false;
	for( long i=0; i < ion_range; i++ )
	{
		long ion = i+ion_low;

		if( source[i] < 0. )
		{
			/* >>chng 04 dec 04, put test on neg abund here, don't print uncles value is very -ve */
			/* >>chng 06 feb 28, from -1e-10 to -1e-9, sim func_t10 had several negative
			 * during initial search, due to extremely high ionization */
			/* >>chng 06 mar 11, from 1e-9 to 2e-9 make many struc elements floats from double */
			if( source[i]<-2e-9 )
				fprintf(ioQQQ,
				" PROBLEM negative ion population in ion_solver, nelem=%li, %s ion=%li val=%.3e Search?%c zone=%li iteration=%li\n",
				nelem , 
				elementnames.chElementSym[nelem],
				ion , 
				source[i] , 
				TorF(conv.lgSearch) ,
				nzone ,
				iteration );
			source[i] = 0.;
			/* if this is one of the iso seq model atoms then must also zero out pops */
			if( ion > nelem-NISO && ion < nelem + 1 )
			{
				long int ipISO = nelem - ion;
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso.numLevels_max[ipISO][nelem]; level++ )
					StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
			}
		}

		/* use new solution */
		dense.xIonDense[nelem][ion] = (realnum)(source[i]*renormnew);
		if( dense.xIonDense[nelem][ion] >= MAX_DENSITY )
		{
			fprintf( ioQQQ, "PROBLEM DISASTER Huge density in ion_solver, nelem %ld ion %ld source %e renormnew %e\n",
				nelem, ion, source[i], renormnew );
		}
		ASSERT( dense.xIonDense[nelem][ion] < MAX_DENSITY );
	}
	
	fixit(); // this should only be done if trimming is not disabled?

	/* Zero levels with abundances < 1e-25 which which will suffer numerical noise */
	while( dense.IonHigh[nelem] > dense.IonLow[nelem] && 
		dense.xIonDense[nelem][dense.IonHigh[nelem]] < 1e-25*abund_total )
	{
		ASSERT( dense.xIonDense[nelem][dense.IonHigh[nelem]] >= 0. );
		/* zero out abundance and heating due to stage of ionization we are about to zero out */
		dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
		thermal.heating[nelem][dense.IonHigh[nelem]-1] = 0.;
		/* decrement counter */
		--dense.IonHigh[nelem];
	}

	// sanity check, either offset stages of low and high ionization
	ASSERT( (dense.IonLow[nelem] <= dense.IonHigh[nelem]) ||
		// both totally neutral
		(dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 ) ||
		// both fully stripped
		( dense.IonLow[nelem]==nelem+1 && dense.IonHigh[nelem]==nelem+1 ) );

	fixit();  //this routine does not ever set *lgNegPop true.  Think it does on newmole though.
	//return *lgNegPop;
	return;
}

STATIC double get_total_abundance_ions( long int nelem )
{
	double abund_total = 0.;

	DEBUG_ENTRY( "get_total_abundance_ions()" );

	ASSERT( nelem >= 0 );
	ASSERT( nelem < LIMELM );

	/* H is special because its abundance spills into three routines -
	 * the ion/atom solver (this routine), the H-mole solvers (hmole), and
	 * the heavy molecule solver.  xmolecules only includes the heavy mole
	 * part for H only.  So the difference between gas_phase and xmolecules
	 * includes the H2 part of the molecular network.  This branch does
	 * this special H case, then the general case (heavy elements are
	 * not part of the H2 molecular network) */

	/* >>chng 01 dec 07, define abund_total, total atom and ion abundance 
	 * here, removing molecules */
	if( nelem == ipHYDROGEN )
	{
		/* Hydrogen is a special case since hmole does not include the
		 * atom/molecules - hmole collapses its network into H = H0 + H+
		 * and then forces the solution determined here for partitioning
		 * between these two */
		abund_total = dense.xIonDense[nelem][0] + dense.xIonDense[nelem][1];
	}
	else
	{
		abund_total = SDIV( dense.gas_phase[nelem] -  dense.xMolecules[nelem] );
	}


	/* protect against case where all gas phase abundances are in molecules, use the
	 * atomic and first ion density from the molecule solver 
	 * >>chng 04 aug 15, NA change from 10 0000 to 10 pre-coef on
	 * FLT_EPSILON for stability in PDR 
	 * the factor 10.*FLT_EPSILON also appears in mole_h_step in total H 
	 * conservation */
	if( fabs( dense.gas_phase[nelem] -  dense.xMolecules[nelem])/SDIV(dense.gas_phase[nelem] ) <
		10.*FLT_EPSILON )
	{
		double renorm;
		/* >>chng 04 jul 31, add logic to conserve nuclei in fully molecular limit;
		 * in first calls, when searching for solution, we may be VERY far 
		 * off, and sum of first ion and atom density may be far larger 
		 * than difference between total gas and molecular densities,
		 * since they reflect the previous evaluation of the solution.  Do 
		 * renorm to cover this case */
		/* first form sum of all atoms and ions */
		realnum sum = 0.;
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
			sum += dense.xIonDense[nelem][ion];
		/* now renorm to this sum - this should be unity, and is not if we have
		 * now conserved particles, due to molecular fraction changing */
		renorm = dense.gas_phase[nelem] / SDIV(sum + dense.xMolecules[nelem] );

		abund_total = renorm * sum;
	}

	/* negative residual density */
	if( abund_total < 0. )
	{
		/* >>chng 05 dec 16, do not abort when net atomic/ ionic abundance 
		 * is negative due to chem net having too much of a species - this 
		 * happens naturally when ices become well formed, but the code will 
		 * converge away from it.  Rather, we set the ionization
		 * convergence flag to "not converge" and soldier on 
		 * if negative populations do not go away, we will eventually 
		 * terminate due to convergence failures */
		/* print error if not search */
		if(!conv.lgSearch )
			fprintf(ioQQQ,
				"PROBLEM ion_solver - neg net atomic abundance zero for nelem= %li, rel val= %.2e conv.nTotalIoniz=%li, fixed\n",
				nelem,
				fabs(abund_total) / SDIV(dense.xMolecules[nelem]),
				conv.nTotalIoniz );
		/* fix up is to use half the positive abundance, assuming chem is 
		 * trying to get too much of this species */
		abund_total = -abund_total/2.;

		/* say that ionization is not converged - do not abort - but if 
		 * cannot converge away from negative solution, this will become a 
		 * convergence failure abort */
		conv.lgConvIoniz = false;
		strcpy( conv.chConvIoniz, "neg ion" );
	}

	ASSERT( abund_total < MAX_DENSITY );

	return abund_total;
}

STATIC void fill_array( long int nelem, long ion_range, valarray<double> &xmat, valarray<double> &source, valarray<double> &auger, double *abund_total )
{
	long int limit, 
	  IonProduced;
	double rateone;
	long ion_low;

	valarray<double> sink(ion_range);
	valarray<int32> ipiv(ion_range);

	DEBUG_ENTRY( "fill_array()" );

	/* this is on the c scale, so H is 0 */
	ASSERT( nelem >= 0);
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonHigh[nelem] >= 0 );

	/* >>chng 04 nov 22,
	 * if gas nearly all atomic/ionic do not let source - sink terms from 
	 * molecular network force system of balance equations to become 
	 * inhomogeneous
	 * what constitutes a source or sink IS DIFFERENT for hydrogen and the rest 
	 * the H solution must couple with hmole - and its defn of source and sink.  For instance, oxygen charge
	 * transfer goes into source and sink terms for hydrogen.  So we never hose source and sink for H.
	 * for the heavy elements, which couple onto comole, mole.source and sink represent terms that
	 * remove atoms and ions from the ionization ladder.  their presence makes the system of
	 * equations inhomogeneous.  we don't want to do this when comole has a trivial effect on
	 * the ionization balance since the matrix becomes unstable */
	/* >>chng 04 dec 06, limit from 0.01 to 1e-10 as per NA suggestion */
	if( nelem>ipHYDROGEN && dense.xMolecules[nelem]/SDIV(dense.gas_phase[nelem]) < 1e-10 )
	{
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
		{
			mole.source[nelem][ion] = 0.;
			mole.sink[nelem][ion] = 0.;
		}
	}

	/* impossible for HIonFrac[nelem] to be zero if IonHigh(nelem)=nelem+1
	 * HIonFrac(nelem) is stripped to hydrogen */
	/* >>chng 01 oct 30, to assert */
	//ASSERT( (dense.IonHigh[nelem] < nelem + 1) || dense.xIonDense[nelem][nelem+1-ipH_LIKE] > 0. );

	/* zero out the ionization and recombination rates that we will modify here,
	 * but not the iso-electronic stages which are done elsewhere,
	 * the nelem stage of ionization is he-like,
	 * the nelem+1 stage of ionization is h-like */

	/* loop over stages of ionization that we solve for here, 
	 * up through and including one less than nelem-NISO,
	 * never actually do highest NISO stages of ionization since they
	 * come from the ionization ratio from the next lower stage */
	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);

	/* the full range of ionization - this is number of ionization stages */
	ASSERT( ion_range <= nelem+2 );

	ion_low = dense.IonLow[nelem];

	for( long i=0; i<ion_range;i++ )
	{
		source[i] = 0.;
	}

	/* zero-out loop comes before main loop since there are off-diagonal
	 * elements in the main ionization loop, due to multi-electron processes,
	 * TotIonizRate and TotRecom were already set in h-like and he-like solvers 
	 * other recombination rates were already set by routines responsible for them */
	for( long ion_from=0; ion_from <= limit; ion_from++ )
	{
		for( long ion_to=0; ion_to < nelem+2; ion_to++ )
		{
			ionbal.RateIoniz[nelem][ion_from][ion_to] = 0.;
		}
	}

	/* auger is used only for debug printout - it is special because with multi-electron
	 * Auger ejection, very high stages of ionization can be produced, well beyond
	 * the highest stage considered here.  so we malloc to the limit */
	for( long i=0; i< LIMELM+1; ++i )
	{
		auger[i] = 0.;
	}

	/* zero out xmat */
	for( long i=0; i< ion_range; ++i )
	{
		for( long j=0; j< ion_range; ++j )
		{
			MAT( xmat, i, j ) = 0.;
		}
	}

	{
		/* this sets up a fake ionization balance problem, with a trivial solution,
		 * for debugging the ionization solver */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nelem==ipCARBON )
		{
			broken();/* within optional debug print statement */
			dense.IonLow[nelem] = 0;
			dense.IonHigh[nelem] = 3;
			*abund_total = 1.;
			/* make up ionization and recombination rates */
			for( long ion=dense.IonLow[nelem]; ion <= limit; ion++ )
			{
				double fac=1;
				if(ion)
					fac = 1e-10;
				ionbal.RateRecomTot[nelem][ion] = 100.;
				for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
				{
					/* direct photoionization of this shell */
					ionbal.PhotoRate_Shell[nelem][ion][ns][0] = fac;
				}
			}
		}
	}

	/* Now put in all recombination and ionization terms from CO_mole() that 
	 * come from molecular reactions. this traces molecular process that 
	 * change ionization stages with this ladder - but do not remove from 
	 * the ladder */
	for( long ion_to=dense.IonLow[nelem]; ion_to <= dense.IonHigh[nelem]; ion_to++ )
	{
		for( long ion_from=dense.IonLow[nelem]; ion_from <= dense.IonHigh[nelem]; ++ion_from )
		{
			/* do not do ion onto itself */
			if( ion_to != ion_from )
			{
				/* this is the rate coefficient for charge transfer from ion to ion_to */
				/*rateone = gv.GrainChTrRate[nelem][ion_from][ion_to];*/
				rateone = mole.xMoleChTrRate[nelem][ion_from][ion_to] * atmdat.lgCTOn;
				// ASSERT( rateone >= 0. );
				MAT( xmat, ion_from-ion_low, ion_from-ion_low ) -= rateone;
				MAT( xmat, ion_from-ion_low, ion_to-ion_low ) += rateone;
			}
		}
	}

	/* now get actual arrays of ionization and recombination processes,
	 * but only for the ions that are done as two-level systems */
	/* in two-stage system, atom + first ion, limit is zero but must
	 * include gv.GrainChTrRate[nelem][1][0] */
	/* grain charge transfer */
	if( gv.lgDustOn() && ionbal.lgGrainIonRecom && gv.lgGrainPhysicsOn )
	{
		long int low;
		/* do not double count this process for atoms that are in the co network - we use
		 * a net recombination coefficient derived from the co solution, this includes grain ct */
		/* >>chng 05 dec 23, add mole.lgElem_in_chemistry */
		if( mole.lgElem_in_chemistry[nelem] )
		/*if( nelem==ipHYDROGEN ||nelem==ipCARBON ||nelem== ipOXYGEN ||nelem==ipSILICON ||nelem==ipSULPHUR 
			 ||nelem==ipNITROGEN ||nelem==ipCHLORINE )*/
		{
			 low = MAX2(1, dense.IonLow[nelem] );
		}
		else
			low = dense.IonLow[nelem];

		for( long ion_to=low; ion_to <= dense.IonHigh[nelem]; ion_to++ )
		{
			for( long ion_from=dense.IonLow[nelem]; ion_from <= dense.IonHigh[nelem]; ++ion_from )
			{
				/* do not do ion onto itself */
				if( ion_to != ion_from )
				{
					/* this is the rate coefficient for charge transfer from ion to ion_to
					 * both grain charge transfer ionization and recombination */
					rateone = gv.GrainChTrRate[nelem][ion_from][ion_to];
					MAT( xmat, ion_from-ion_low, ion_from-ion_low ) -= rateone;
					MAT( xmat, ion_from-ion_low, ion_to-ion_low ) += rateone;
				}
			}
		}
	}

	for( long ion=dense.IonLow[nelem]; ion <= limit; ion++ )
	{
		/* thermal & secondary collisional ionization */
		rateone = ionbal.CollIonRate_Ground[nelem][ion][0] +
			secondaries.csupra[nelem][ion] +
			/* inner shell ionization by UTA lines */
			ionbal.UTA_ionize_rate[nelem][ion];
		ionbal.RateIoniz[nelem][ion][ion+1] += rateone;

		if( ion+1-ion_low < ion_range )
		{
			/* UTA ionization */
			/* depopulation processes enter with negative sign */
			MAT( xmat, ion-ion_low, ion-ion_low ) -= rateone;
			MAT( xmat, ion-ion_low, ion+1-ion_low ) += rateone;

			/* total recombination rate */
			/* loss of next higher stage due to recombination to this ion stage */
			MAT( xmat, ion+1-ion_low, ion+1-ion_low ) -= ionbal.RateRecomTot[nelem][ion];
			MAT( xmat, ion+1-ion_low, ion-ion_low ) += ionbal.RateRecomTot[nelem][ion];
		}

		/* loop over all atomic sub-shells to include photoionization */
		for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
		{
			/* this is the primary ionization rate - add to diagonal element,
			 * test on ion stage is so that we don't include ionization from the very highest
			 * ionization stage to even higher - since those even higher stages are not considered
			 * this would appear as a sink - but populations of this highest level is ensured to
			 * be nearly trivial and neglecting it production of even higher ionization OK */
			/* >>chng 04 nov 29 RJRW, include following in this branch so only
			 * evaluated when below ions done with iso-sequence */
			if( ion+1-ion_low < ion_range )
			{
				/* this will be redistributed into charge states in following loop */ 

				/* t_yield::Inst().nelec_eject(nelem,ion,ns) is total number of electrons that can
				 * possibly be freed 
				 * loop over nej, the number of electrons ejected including the primary,
				 * nej = 1 is primary, nej > 1 includes primary plus Auger 
				 * t_yield::Inst().elec_eject_frac is prob of nej electrons */
				for( long nej=1; nej <= t_yield::Inst().nelec_eject(nelem,ion,ns); nej++ )
				{
					/* this is the ion that is produced by this ejection,
					 * limited by highest possible stage of ionization -
					 * do not want to ignore ionization that go beyond this */
					IonProduced = MIN2(ion+nej,dense.IonHigh[nelem]);
					rateone = ionbal.PhotoRate_Shell[nelem][ion][ns][0]*
						t_yield::Inst().elec_eject_frac(nelem,ion,ns,nej-1);

					/* direct photoionization of this shell */
					ionbal.RateIoniz[nelem][ion][IonProduced] += rateone;

					/* >>chng 04 sep 06, above had included factor of nej to get rate, but
					 * actually want events into particular ion */
					/* number of electrons ejected
					 *(double)nej; */
					/* it goes into this charge state - recall upper cap due to ion stage trimming 
					 * note that compensating loss term on diagonal was done before this
					 * loop, since frac_elec_eject adds to unity */
					MAT( xmat, ion-ion_low, ion-ion_low ) -= rateone;
					MAT( xmat, ion-ion_low, IonProduced-ion_low ) += rateone;

					/* only used for possible printout - multiple electron Auger rate  -
					 * do not count one-electron as Auger */
					if( nej>1 )
						auger[IonProduced-1] += rateone;
				}
			}
		}

		/* this is charge transfer ionization of this species by hydrogen and helium */
		rateone = 
			atmdat.HeCharExcIonOf[nelem][ion]*dense.xIonDense[ipHELIUM][1]+ 
			atmdat.HCharExcIonOf[nelem][ion]*dense.xIonDense[ipHYDROGEN][1];
		ionbal.RateIoniz[nelem][ion][ion+1] += rateone;
		if( ion+1-ion_low < ion_range )
		{
			MAT( xmat, ion-ion_low, ion-ion_low ) -= rateone;
			MAT( xmat, ion-ion_low, ion+1-ion_low ) += rateone;
		}
	}

	/* after this loop, ionbal.RateIonizTot and ionbal.RateRecomTot have been defined for the
	 * stages of ionization that are done with simple */
	/* begin loop at first species that is treated with full model atom
	 * but possible that lowest stage of ionization is higher than this, hence MAX3 */
	for( long ion= MAX3(0,limit+1,ion_low); ion<=dense.IonHigh[nelem]; ion++ )
	{
		ASSERT( ion>=0 && ion<nelem+2 );
		/* use total ionization/recombination rates for species done with ISO solver */
		if( ion+1-ion_low < ion_range )
		{
			fixit();  //This should also include multiple-electron ionizations

			/* depopulation processes enter with negative sign */
			MAT( xmat, ion-ion_low, ion-ion_low ) -= ionbal.RateIoniz[nelem][ion][ion+1];
			MAT( xmat, ion-ion_low, ion+1-ion_low ) += ionbal.RateIoniz[nelem][ion][ion+1];

			/* loss of next higher ion due to recombination to this ion stage */
			MAT( xmat, ion+1-ion_low, ion+1-ion_low ) -= ionbal.RateRecomTot[nelem][ion];
			MAT( xmat, ion+1-ion_low, ion-ion_low ) += ionbal.RateRecomTot[nelem][ion];
		}
	}

	/*>>chng 06 nov 28 only include source from molecules if we have an estimated first 
	 * solution - first test is that we have called mole at least twice,
	 * second test is that we are on a later iteration */
	if( conv.nTotalIoniz > 1 || iteration > 1 )
	{

		for( long i=0; i<ion_range;i++ )
		{
			long ion = i+ion_low;

			/* these are the external source and sink terms */
			/* source first */
			/* need negative sign to get positive pops */
			source[i] -= mole.source[nelem][ion];
			MAT( xmat, i, i ) -= mole.sink[nelem][ion];
		}
	}

	return;
}

STATIC bool lgHomogeneousSource( long nelem, long ion_low, long ion_range, valarray<double> &xmat, valarray<double> &source, double abund_total )
{
	bool lgHomogeneous = true;
	double totsrc;

	DEBUG_ENTRY( "lgHomogeneousSource()" );

	/* this will be sum of source and sink terms, will be used to decide if
	 * matrix is singular */
	totsrc = 0.;
	for( long i=0; i<ion_range;i++ )
		totsrc -= source[i];

	/* matrix is not homogeneous if source is non-zero */
	if( totsrc != 0. )
		lgHomogeneous = false;

	fixit();  // dynamics rates need to be moved into fill_array?
	/* chng 03 jan 13 rjrw, add in dynamics if required here,
	 * - only do advection if we have not overrun the radius scale */
	if( iteration > dynamics.n_initial_relax+1 && dynamics.lgAdvection && !dynamics.lgEquilibrium
		&& dynamics.Rate != 0. )
	{
		for( long i=0; i<ion_range;i++ )
		{			 
			long ion = i+ion_low;
			MAT( xmat, i, i ) -= dynamics.Rate;
			source[i] -= dynamics.Source[nelem][ion];
			/* fprintf(ioQQQ," %li %li %.3e (%.3e %.3e)\n",i,i,MAT(*xmat,i,i),
			  dynamics.Rate, dynamics.Source[nelem][ion]);*/
		}
		lgHomogeneous = false;
	}

	/* >>chng 06 nov 21, for very high ionization parameter sims the H molecular
	* fraction can become so small that atom = atom + molecule.  In this case we
	* will not count system as an inhomogeneous case since linear algebra package
	* will fail.  If molecules are very small, count as homogeneous.
	* Note that we have already added sink terms to the main matrix and they
	* will not be removed, a possible source of error, but they must not
	* have been significant, given that the molecular fraction is so small */
	if( !lgHomogeneous && ion_range==2 )
	{
		/* solve algebraically */
		double a = MAT( xmat, 0, 0 ), 
			b = MAT( xmat, 1, 0 ) ,
			c = MAT( xmat, 0, 1 ) ,
			d = MAT( xmat, 1, 1 );
		b = SDIV(b);
		d = SDIV(d);
		double ratio1 = a/b , ratio2 = c/d , fratio1=fabs(a/b),fratio2=fabs(c/d);
		if( fabs(ratio1-ratio2)/MAX2(fratio1,fratio2) <DBL_EPSILON*1e4 )
		{
			//abund_total = SDIV( dense.gas_phase[nelem] -  dense.xMolecules[nelem] );
			lgHomogeneous = true;
		}
	}
#	if 0
	if( nelem==ipHYDROGEN && 
		fabs(dense.xMolecules[nelem]) / SDIV(dense.gas_phase[ipHYDROGEN]) <DBL_EPSILON*100. )
	{
		abund_total = SDIV( dense.gas_phase[nelem] -  dense.xMolecules[nelem] );
		lgHomogeneous = true;
	}
#	endif

	/* this is true if no source terms 
	 * we will use total population and species conservation to replace one
	 * set of balance equations since overdetermined */
	if( lgHomogeneous  )
	{
		double rate_ioniz=1., rate_recomb /*, scale = 0.*/;
		/* Simple estimate of most abundant ion */
		long jmax = 0;
		for( long i=0; i<ion_range-1;i++)
		{ 
			long ion = i+ion_low;
			rate_ioniz *= ionbal.RateIonizTot(nelem,ion);
			rate_recomb = ionbal.RateRecomTot[nelem][ion];
			/* trips if ion rate zero, so ll the gas will be more neutral than this */
			if( rate_ioniz == 0)
				break;
			/* rec rate is zero */
			if( rate_recomb <= 0.) 
				break;

			rate_ioniz /= rate_recomb;
			if( rate_ioniz > 1.) 
			{
				/* this is peak ionization stage */
				jmax = i;
				rate_ioniz = 1.;
			}
		}

		/* replace its matrix elements with population sum */
		for( long i=0; i<ion_range;i++ )
		{
			MAT(xmat,i,jmax) = 1.;
		}
		source[jmax] = abund_total;
	}

#if 0
	if( false && nelem == ipHYDROGEN && dynamics.lgAdvection&& iteration>1 ) 
	{
		fprintf(ioQQQ,"DEBUGG Rate %.2f %.3e \n",fnzone,dynamics.Rate);
		fprintf(ioQQQ," %.3e %.3e\n", ionbal.RateIonizTot(nelem,0), ionbal.RateIonizTot(nelem,1) );
		fprintf(ioQQQ," %.3e %.3e\n", ionbal.RateRecomTot[nelem][0], ionbal.RateRecomTot[nelem][1]);
		fprintf(ioQQQ," %.3e %.3e %.3e\n\n", dynamics.Source[nelem][0], dynamics.Source[nelem][1], dynamics.Source[nelem][2]);
	}

	{
		/* option to print matrix */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone==1 && lgPrintIt )
		{
			fprintf( ioQQQ, 
				" DEBUG ion_solver: nelem=%li ion_range=%li, limit=%li, nConv %li xmat follows\n",
				nelem , ion_range,limit  , conv.nTotalIoniz );
			if( lgHomogeneous )
				fprintf(ioQQQ , "Homogeneous \n");
			for( long i=0; i<ion_range; ++i )
			{
				for( long j=0;j<ion_range;j++ )
				{
					fprintf(ioQQQ,"%e\t",MAT(xmat,i,j));
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ,"source follows\n");
			for( long i=0; i<ion_range;i++ )
			{
				fprintf(ioQQQ,"%e\t",source[i]);
			}
			fprintf(ioQQQ,"\n");
		}
	}
#endif

	return lgHomogeneous;
}

STATIC void PrintRates( long nelem, bool lgNegPop, double abund_total, valarray<double> &auger, bool lgPrintIt )
{
	DEBUG_ENTRY( "PrintRates()" );

	long ion;

	/* this should not happen */
	if( lgNegPop )
	{
		fprintf( ioQQQ, " PROBLEM Negative population found for abundance of ionization stage of element %4.4s, ZONE=%4ld\n", 
		  elementnames.chElementNameShort[nelem], nzone );

		fprintf( ioQQQ, " Populations were" );
		for( ion=1; ion <= dense.IonHigh[nelem]+1; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", dense.xIonDense[nelem][ion-1] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " destroy vector =" );
		for( ion=1; ion <= dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", ionbal.RateIonizTot(nelem,ion-1) );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " CTHeavy  vector =" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", atmdat.HeCharExcIonOf[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " HCharExcIonOf vtr=" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", atmdat.HCharExcIonOf[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " CollidRate  vtr=" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", ionbal.CollIonRate_Ground[nelem][ion][0] );
		}
		fprintf( ioQQQ, "\n" );

		/* photo rates per subshell */
		fprintf( ioQQQ, " photo rates per subshell, ion\n" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%3ld", ion );
			for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
			{
				fprintf( ioQQQ, "%9.1e", ionbal.PhotoRate_Shell[nelem][ion][ns][0] );
			}
			fprintf( ioQQQ, "\n" );
		}

		/* now check out creation vector */
		fprintf( ioQQQ, " create  vector =" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", ionbal.RateRecomTot[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );
	}

	/* option to print ionization and recombination arrays
	 * prt flag set with print array print arrays command */
	if( lgPrintIt || prt.lgPrtArry[nelem] || lgNegPop )
	{
		/* say who we are, what we are doing .... */
		fprintf( ioQQQ, 
			"\n %s ion_solver ion/rec rt [s-1] %s nz%.2f Te%.4e ne%.4e Tot abun:%.3e ion abun%.2e mole%.2e\n", 
			elementnames.chElementSym[nelem],
			elementnames.chElementName[nelem],
			fnzone,
			phycon.te , 
			dense.eden,
			dense.gas_phase[nelem],
			abund_total ,
			dense.xMolecules[nelem] );
		/* total ionization rate, all processes */
		fprintf( ioQQQ, " %s Ioniz total " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", ionbal.RateIonizTot(nelem,ion) );
		}
		fprintf( ioQQQ, "\n" );

		/* sinks from the chemistry network */
		fprintf(ioQQQ," %s molecule snk",
			elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf(ioQQQ," %9.2e", mole.sink[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		if( dynamics.lgAdvection )
		{
			/* sinks from the dynamics */
			fprintf(ioQQQ," %s dynamics snk",
				elementnames.chElementSym[nelem]);
			for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
			{
				fprintf(ioQQQ," %9.2e", dynamics.Rate );
			}
			fprintf( ioQQQ, "\n" );
		}

		/* sum of all creation processes */
		fprintf( ioQQQ, " %s Recom total " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", ionbal.RateRecomTot[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* collisional ionization */
		fprintf( ioQQQ, " %s Coll ioniz  " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double ColIoniz = ionbal.CollIonRate_Ground[nelem][ion][0];
			if( ion > nelem - NISO )
			{
				long ipISO = nelem-ion;
				ASSERT( ipISO >=0 && ipISO < NISO );
				ColIoniz *= StatesElemNEW[nelem][nelem-ipISO][0].Pop;
				if( dense.xIonDense[nelem][nelem-ipISO] > 0. )
				{
					for( long ipLevel=1; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
					{
						ColIoniz += iso.ColIoniz[ipISO][nelem][ipLevel] * dense.EdenHCorr * StatesElemNEW[nelem][nelem-ipISO][ipLevel].Pop;
					}
					ColIoniz /= dense.xIonDense[nelem][nelem-ipISO];
				}
			}
			fprintf( ioQQQ, " %9.2e", ColIoniz );
		}
		fprintf( ioQQQ, "\n" );		

		/* UTA ionization */
		fprintf( ioQQQ, " %s UTA ioniz   " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", ionbal.UTA_ionize_rate[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* photo ionization */
		fprintf( ioQQQ, " %s Photoion snk" ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double PhotIoniz = 0.;
			for( long ipShell = 0; ipShell < Heavy.nsShells[nelem][ion]; ipShell++ )
				PhotIoniz += ionbal.PhotoRate_Shell[nelem][ion][ipShell][0];
			
			// still don't have the total if stage is in one of the iso-sequences
			if( ion > nelem - NISO )
			{
				long ipISO = nelem-ion;
				ASSERT( ipISO >=0 && ipISO < NISO );
				if( dense.xIonDense[nelem][nelem-ipISO]>0  )
				{
					PhotIoniz *= StatesElemNEW[nelem][nelem-ipISO][0].Pop;
					for( long ipLevel=1; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
					{
						PhotIoniz += iso.gamnc[ipISO][nelem][ipLevel] * StatesElemNEW[nelem][nelem-ipISO][ipLevel].Pop;
					}
					PhotIoniz /= dense.xIonDense[nelem][nelem-ipISO];
				}
			}
			fprintf( ioQQQ, " %9.2e", PhotIoniz );
		}
		fprintf( ioQQQ, "\n" );		

		/* photoionization source (source of this stage due to ionization of next stage) */
		fprintf( ioQQQ, " %s Photoion src" ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double source = 0.;
			if( ion>0 )
				source = ionbal.RateIoniz[nelem][ion-1][ion];
				
			fprintf( ioQQQ, " %9.2e", source );
		}
		fprintf( ioQQQ, "\n" );
				
		/* auger production (of this stage due to auger ionization of another) */
		fprintf( ioQQQ, " %s Auger src   " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double source = 0.;
			if( ion>0 )
				source = auger[ion-1];
			
			fprintf( ioQQQ, " %9.2e", source );
		}
		fprintf( ioQQQ, "\n" );

		/* secondary ionization */
		fprintf( ioQQQ, " %s Secon ioniz " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", 
				secondaries.csupra[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* grain ionization - not total rate but should be dominant process */
		fprintf( ioQQQ, " %s ion on grn  "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", gv.GrainChTrRate[nelem][ion][ion+1] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge transfer ionization from chemistry */
		fprintf( ioQQQ, " %s ion xfr mol "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", mole.xMoleChTrRate[nelem][ion][ion+1] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge exchange ionization */
		fprintf( ioQQQ, " %s chr trn ion " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			/* sum has units s-1 */
			double sum;

			if( nelem==ipHELIUM && ion==0 )
			{
				sum = atmdat.HeCharExcIonTotal;
			}
			else if( nelem==ipHYDROGEN && ion==0 )
			{
				sum = atmdat.HCharExcIonTotal;
			}
			else
				sum = atmdat.HeCharExcIonOf[nelem][ion] * dense.xIonDense[ipHELIUM][1]+ 
				atmdat.HCharExcIonOf[nelem][ion] * dense.xIonDense[ipHYDROGEN][1];
			
			fprintf( ioQQQ, " %9.2e", sum );
		}
		fprintf( ioQQQ, "\n" );

		/* radiative recombination */
		fprintf( ioQQQ, " %s radiati rec "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", dense.eden*ionbal.RR_rate_coef_used[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* DR recombination */
		/* DR rates are only defined for Li-like and lower, so never access higher stages */
		fprintf( ioQQQ, " %s dielect rec "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < min(nelem-1,dense.IonHigh[nelem]); ion++ )
		{
			fprintf( ioQQQ, " %9.2e", dense.eden*ionbal.DR_rate_coef_used[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* old DR recombination */
		fprintf( ioQQQ, " %s dr old  rec "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < min(nelem-1,dense.IonHigh[nelem]); ion++ )
		{
			fprintf( ioQQQ, " %9.2e", dense.eden*ionbal.DR_old_rate_coef[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* Badnell DR recombination */
		fprintf( ioQQQ, " %s drBadnel rec"  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < min(nelem-1,dense.IonHigh[nelem]); ion++ )
		{
			fprintf( ioQQQ, " %9.2e", dense.eden*ionbal.DR_Badnell_rate_coef[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* Cota rate */
		fprintf( ioQQQ, " %s CotaRate rec"  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", dense.eden*ionbal.CotaRate[ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* grain recombination - not total but from next higher ion, should
		 * be dominant */
		fprintf( ioQQQ, " %s rec on grn  "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", gv.GrainChTrRate[nelem][ion+1][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge transfer recombination from chemistry */
		fprintf( ioQQQ, " %s rec xfr mol "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", mole.xMoleChTrRate[nelem][ion+1][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge exchange recombination */
		fprintf( ioQQQ, " %s chr trn rec "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double sum;
			
			if( nelem==ipHELIUM && ion==0 )
			{
				sum = atmdat.HeCharExcRecTotal;
			}
			else if( nelem==ipHYDROGEN && ion==0 )
			{
				sum = atmdat.HCharExcRecTotal;
			}
			else
				sum = atmdat.HCharExcRecTo[nelem][ion] * dense.xIonDense[ipHYDROGEN][0] +
				atmdat.HeCharExcRecTo[nelem][ion] * dense.xIonDense[ipHELIUM][0];
				
			fprintf( ioQQQ, " %9.2e", sum );
		}
		fprintf( ioQQQ, "\n" );

		// NB NB NB -- these are all in units cm-3 s-1.  They are NOT comparable to the above rates in s-1.
		{
			fprintf(ioQQQ," %s src cm-3 s-1\n",
					elementnames.chElementSym[nelem]);
			
			/* sources from the chemistry network */
			fprintf(ioQQQ," %s molecule src",
					elementnames.chElementSym[nelem]);
			for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
			{
				fprintf(ioQQQ," %9.2e", mole.source[nelem][ion] );
			}
			fprintf( ioQQQ, "\n" );
			
			if( dynamics.lgAdvection )
			{
				/* source from the dynamcs */
				fprintf(ioQQQ," %s dynamics src",
						elementnames.chElementSym[nelem]);
				for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
				{
					fprintf(ioQQQ," %9.2e", dynamics.Source[nelem][ion] );
				}
				fprintf( ioQQQ, "\n" );
			}
		
		}
		
		/* the "new" abundances the resulted from the previous ratio */
		fprintf( ioQQQ, " %s Abun [cm-3] " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion <= dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, " %9.2e", dense.xIonDense[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );
	}

	if( lgNegPop )
	{
		ContNegative();
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	return;
}

/* 

	 Solve an ionization level system with specified ionization and
	 recombination rates between neighboring ions, and additional sink
	 and source terms.  The sink array is overwritten, and the results
	 appear in the source array.

	 Written in matrix form, the algorithm is equivalent to the
	 tridiagonal algorithm in Numerical Recipes applied to:

	 / i_0+a_0     -r_0          .           .    .  \ / x_0 \   / s_0 \
	 |  -i_0    i_1+a_1+r_0    -r_1          .    .  | | x_1 |   | s_1 |
	 |    .        -i_1      i_2+a_2+r_1   -r_2   .  | | x_2 |   | s_2 |
	 |    .          .       (etc....)               | | ... | = | ... |
	 \    .          .          .                    / \     /   \     /

	 where i, r are the ionization and recombination rates, s is the
	 source rate and a is the sink rate.

	 This matrix is diagonally dominant only when the sink terms are
	 large -- the alternative method coded here prevents rounding error
	 in the diagonal terms disturbing the solution when this is not the
	 case.

*/

/* solveions tridiagonal solver but optimized for structure of balance matrix */
void solveions(double *ion, double *rec, double *snk, double *src,
	       long int nlev, long int nmax)
{
	double kap, bet;
	long int i;

	DEBUG_ENTRY( "solveions()" );

	if(nmax != -1) 
	{
		/* Singular case */
		src[nmax] = 1.;
		for(i=nmax;i<nlev-1;i++)
			src[i+1] = src[i]*ion[i]/rec[i];
		for(i=nmax-1;i>=0;i--)
			src[i] = src[i+1]*rec[i]/ion[i];
	} 
	else 
	{
		i = 0;
		kap = snk[0];    
		for(i=0;i<nlev-1;i++) 
		{
			bet = ion[i]+kap;
			if(bet == 0.)
			{
				fprintf(ioQQQ,"Ionization solver error\n");
				cdEXIT(EXIT_FAILURE);
			}
			bet = 1./bet;
			src[i] *= bet;
			src[i+1] += ion[i]*src[i];
			snk[i] = bet*rec[i];
			kap = kap*snk[i]+snk[i+1];
		}
		bet = kap;
		if(bet == 0.)
		{
			fprintf(ioQQQ,"Ionization solver error\n");
			cdEXIT(EXIT_FAILURE);
		}
		src[i] /= bet;

		for(i=nlev-2;i>=0;i--)
		{
			src[i] += snk[i]*src[i+1];
		}
	}
}

void ion_wrapper( long nelem )
{

	DEBUG_ENTRY( "ion_wrapper()" );

	ASSERT( nelem >= ipHYDROGEN );
	ASSERT( nelem < LIMELM );

	switch( nelem )
	{
	case ipHYDROGEN:
		IonHydro();
		break;
	case ipHELIUM:
		IonHelium();
		break;
	case ipCARBON:
		IonCarbo();
		break;
	case ipOXYGEN:
		IonOxyge();
		break;
	case ipNITROGEN:
		IonNitro();
		break;
	case ipSILICON:
		IonSilic();
		break;
	case ipSULPHUR:
		IonSulph();
		break;
	case ipCHLORINE:
		IonChlor();
		break;
	case ipLITHIUM:
		IonLithi();
		break;
	case ipBERYLLIUM:
		IonBeryl();
		break;
	case ipBORON:
		IonBoron();
		break;
	case ipFLUORINE:
		IonFluor();
		break;
	case ipNEON:
		IonNeon();
		break;
	case ipSODIUM:
		IonSodiu();
		break;
	case ipMAGNESIUM:
		IonMagne();
		break;
	case ipALUMINIUM:
		IonAlumi();
		break;
	case ipPHOSPHORUS:
		IonPhosi();
		break;
	case ipARGON:
		IonArgon();
		break;
	case ipPOTASSIUM:
		IonPotas();
		break;
	case ipCALCIUM:
		IonCalci();
		break;
	case ipSCANDIUM:
		IonScand();
		break;
	case ipTITANIUM:
		IonTitan();
		break;
	case ipVANADIUM:
		IonVanad();
		break;
	case ipCHROMIUM:
		IonChrom();
		break;
	case ipMANGANESE:
		IonManga();
		break;
	case ipIRON:
		IonIron();
		break;
	case ipCOBALT:
		IonCobal();
		break;
	case ipNICKEL:
		IonNicke();
		break;
	case ipCOPPER:
		IonCoppe();
		break;
	case ipZINC:
		IonZinc();
		break;
	default:
		TotalInsanity();
	}
	
	if( trace.lgTrace && dense.lgElmtOn[nelem] && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     ion_wrapper returns; %s fracs = ", elementnames.chElementSym[nelem] );
		for( long ion = 0; ion<=nelem+1; ion++ )
			fprintf( ioQQQ,"%10.3e ", dense.xIonDense[nelem][ion]/dense.gas_phase[nelem] );
		fprintf( ioQQQ, "\n" );
	}	

	return;
}
