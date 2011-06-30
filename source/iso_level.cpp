/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_level solve for iso-sequence level populations */
#include "cddefines.h"
#include "atmdat.h"
#include "continuum.h"
#include "conv.h"
#include "dense.h"
#include "dynamics.h"
#include "elementnames.h"
#include "grainvar.h"
#include "he.h"
#include "helike.h"
#include "hmi.h"
#include "hydrogenic.h"
#include "ionbal.h"
#include "iso.h"
#include "mole.h"
#include "phycon.h"
#include "physconst.h"
#include "rfield.h"
#include "secondaries.h"
#include "taulines.h"
#include "thirdparty.h"
#include "trace.h"

/* solve for level populations  */
void iso_level( const long int ipISO, const long int nelem)
{
	long int ipHi,
		ipLo,
		i,
		level,
		level_error;

  const long int numlevels_local = iso.numLevels_local[ipISO][nelem];

	double BigError;

	int32 nerror;
	double HighestPColOut = 0.,
		TotalPop;
	bool lgNegPop=false;
	valarray<int32> ipiv(numlevels_local); 
	/* this block of variables will be obtained and freed here */
	valarray<double>
		creation(numlevels_local),
		error(numlevels_local),
		work(numlevels_local),
		Save_creation(numlevels_local);
	double source=0.,
		sink=0.;
	valarray<double> PopPerN(iso.n_HighestResolved_local[ipISO][nelem]+1);

	multi_arr<double,2,C_TYPE> z, SaveZ;

	DEBUG_ENTRY( "iso_level()" );

	/* check that we were called with valid charge */
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	/* now do the 2D array */
	z.alloc(numlevels_local,numlevels_local);

	/* fill in recombination vector - values were set in iso_ionize_recombine.cpp */
	for( level=0; level < numlevels_local; level++ )
	{
		ASSERT( dense.xIonDense[nelem][nelem+1-ipISO] >= 0.f );
		/* total recombination from once more ionized [cm-3 s-1] */
		creation[level] = iso.RateCont2Level[ipISO][nelem][level] * dense.xIonDense[nelem][nelem+1-ipISO];
	}

	/* these two collision rates must be the same or we are in big trouble,
	 * since used interchangeably */
	ASSERT( ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0]< SMALLFLOAT ||
		fabs( (iso.ColIoniz[ipISO][nelem][0]* dense.EdenHCorr) /
		SDIV(ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0] ) - 1.) < 0.001 );

	/* which case atom to solve??? */
	if( dense.xIonDense[nelem][nelem+1-ipISO] < SMALLFLOAT || iso.xIonSimple[ipISO][nelem] < 1e-35 )
	{
		/* don't bother if no ionizing radiation */
		strcpy( iso.chTypeAtomUsed[ipISO][nelem], "zero " );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_level iso=%2ld nelem=%2ld simple II/I=%10.2e so not doing equilibrium, doing %s.\n", 
				ipISO, nelem, iso.xIonSimple[ipISO][nelem], iso.chTypeAtomUsed[ipISO][nelem] );
		}

		/* total ionization will just be the ground state */
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = iso.RateLevel2Cont[ipISO][nelem][0];
		lgNegPop = false;
		StatesElemNEW[nelem][nelem-ipISO][0].Pop = dense.xIonDense[nelem][nelem-ipISO];

		for( long n=1; n < numlevels_local; n++ )
		{
			StatesElemNEW[nelem][nelem-ipISO][n].Pop =  0.;
		}
		iso.qTot2S[ipISO][nelem] = 0.;
	}
	else
	{
		/* this branch is main solver, full level populations 
		 * assert since this code must change if NISO ever increased */
		ASSERT( NISO == 2 );
		long SpinUsed[NISO] = {2,1};
		long indexNmaxP =
			iso.QuantumNumbers2Index[ipISO][nelem][ iso.n_HighestResolved_local[ipISO][nelem] ][1][SpinUsed[ipISO]];

		strcpy( iso.chTypeAtomUsed[ipISO][nelem], "Popul" );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_level iso=%2ld nelem=%2ld doing regular matrix inversion, %s\n", 
				ipISO, nelem, iso.chTypeAtomUsed[ipISO][nelem] );
		}

		multi_arr<quantumState,3>::const_iterator StElm = StatesElemNEW.begin(nelem,nelem-ipISO);

		/* master balance equation, use when significant population */
		for( level=0; level < numlevels_local; level++ )
		{
			/* all process depopulating level and placing into the continuum
			 * this does NOT include grain charge transfer ionization, added below */
			z[level][level] = iso.RateLevel2Cont[ipISO][nelem][level];
			
			if( level == 1 )
				/* >>chng 05 dec 21, rm eden to make into rate coefficient */
				iso.qTot2S[ipISO][nelem] = iso.ColIoniz[ipISO][nelem][level];

			multi_arr<transition,4>::const_iterator Trans = Transitions.begin(ipISO,nelem,level);
			md4ci Boltz = iso.Boltzmann.begin(ipISO,nelem,level);

			/* all processes populating level from below */
			for( ipLo=0; ipLo < level; ipLo++ )
			{
				double coll_down = Trans[ipLo].Coll.ColUL * dense.EdenHCorr;

				double RadDecay = MAX2( iso.SmallA, Trans[ipLo].Emis->Aul*
							(Trans[ipLo].Emis->Pesc + 
							 Trans[ipLo].Emis->Pelec_esc + 
							 Trans[ipLo].Emis->Pdest)*
							KILL_BELOW_PLASMA(Trans[ipLo].EnergyWN*WAVNRYD) );

				double pump = MAX2( iso.SmallA, Trans[ipLo].Emis->pump *
						    KILL_BELOW_PLASMA(Trans[ipLo].EnergyWN*WAVNRYD) );

				if( iso.lgRandErrGen[ipISO] )
				{
					coll_down *= iso.ErrorFactor[ipISO][nelem][level][ipLo][IPCOLLIS];
					RadDecay *= iso.ErrorFactor[ipISO][nelem][level][ipLo][IPRAD];
					pump *= iso.ErrorFactor[ipISO][nelem][level][ipLo][IPRAD];
				}

				double coll_up = coll_down * 
					(double)StElm[level].g/
					(double)StElm[ipLo].g*
					Boltz[ipLo];

				z[ipLo][ipLo] += coll_up + pump ;
				z[ipLo][level] = - ( coll_up + pump );

				double pump_down = pump *
					(double)StElm[ipLo].g/
					(double)StElm[level].g;

				z[level][level] += RadDecay + pump_down + coll_down;
				z[level][ipLo] = - (RadDecay + pump_down + coll_down);

				if( level == indexNmaxP )
				{
					HighestPColOut += coll_down;
				}				 
				if( ipLo == indexNmaxP )
				{
					HighestPColOut += coll_up;
				}

				/* find total collisions out of 2^3S to singlets. */
				if( (level == 1) && (ipLo==0) )
				{
					iso.qTot2S[ipISO][nelem] += coll_down/dense.EdenHCorr;
				}
				if( (ipLo == 1) && (ipISO==ipH_LIKE || StElm[level].S==1) )
				{
					iso.qTot2S[ipISO][nelem] += coll_up/dense.EdenHCorr;
				}
			}
		}

		if( ipISO == nelem )
		{
			/* iso.lgCritDensLMix[ipISO] is a flag used to print warning if density is
			* too low for first collapsed level to be l-mixed.  Check is if l-mixing collisions
			* out of highest resolved singlet P are greater than sum of transition probs out.	*/
			if( HighestPColOut < 1./StatesElemNEW[nelem][nelem-ipISO][indexNmaxP].lifetime )
			{
				iso.lgCritDensLMix[ipISO] = false;
			}
		}

		/** \todo 2 the indices for the two-photon rates must be changed for further iso sequences. */  
		ASSERT( ipISO <= ipHE_LIKE );

		/* induced two photon emission - special because upward and downward are
		 * not related by ratio of statistical weights */
		/* iso.lgInd2nu_On is controlled with SET IND2 ON/OFF command */
		z[1+ipISO][0] -= iso.TwoNu_induc_dn[ipISO][nelem]*iso.lgInd2nu_On;
		z[0][1+ipISO] -= iso.TwoNu_induc_up[ipISO][nelem]*iso.lgInd2nu_On;

		/* rates out of 1s, and out of 2s */
		z[0][0] += iso.TwoNu_induc_up[ipISO][nelem]*iso.lgInd2nu_On;
		z[1+ipISO][1+ipISO] += iso.TwoNu_induc_dn[ipISO][nelem]*iso.lgInd2nu_On;

		/* grain charge transfer recombination and ionization to ALL other stages */
		for( long ion=0; ion<=nelem+1; ++ion )
		{
			if( ion!=nelem-ipISO )
			{
				source += gv.GrainChTrRate[nelem][ion][nelem-ipISO] *
					dense.xIonDense[nelem][ion];
				sink += gv.GrainChTrRate[nelem][nelem-ipISO][ion];
#if 0
				source += mole.xMoleChTrRate[nelem][ion][nelem-ipISO] * 
					dense.xIonDense[nelem][ion] * atmdat.lgCTOn;
				sink += mole.xMoleChTrRate[nelem][nelem-ipISO][ion] * atmdat.lgCTOn;
#endif
			}
		}
		
#if	0
		/* >>chng 02 Sep 06 rjrw -- all elements have these terms */
		/*>>>chng 02 oct 01, only include if lgAdvection is set */
		if( iteration > dynamics.n_initial_relax+1 && dynamics.lgAdvection && 
				dynamics.Rate != 0.0 &&
			!dynamics.lgEquilibrium && dynamics.lgISO[ipISO])
		{
			/* add in advection - these terms normally zero */
			source += dynamics.Source[nelem][nelem-ipISO];
			/* >>chng 02 Sep 06 rjrw -- advective term not recombination */
			sink += dynamics.Rate;
		}
#else
		/*>>>chng 02 oct 01, only include if lgAdvection is set */
		if( iteration > dynamics.n_initial_relax+1 && dynamics.lgAdvection &&
				dynamics.Rate != 0.0 &&
			!dynamics.lgEquilibrium && dynamics.lgISO[ipISO])
		{
			for( level=0; level < numlevels_local; level++ )
			{
				creation[level] += dynamics.StatesElemNEW[nelem][nelem-ipISO][level];
				z[level][level] += dynamics.Rate;
			}
		}
#endif

		/* ionization from/recombination from lower ionization stages */
		for(long ion_from=dense.IonLow[nelem]; ion_from < MIN2( dense.IonHigh[nelem], nelem-ipISO ) ; ion_from++ )
		{
			if( ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] >= 0. )
			{
				double damper = sexp( dense.gas_phase[nelem] / 1e10 );
				
				/* ionization from lower ionization stages, cm-3 s-1 */
				source += ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] * dense.xIonDense[nelem][ion_from] * damper;
				/* recombination to next lower ionization stage, s-1 */
				if( ion_from == nelem-1-ipISO )
					sink += ionbal.RateRecomTot[nelem][ion_from] * damper;
			}
		}

		ASSERT( source >= 0.f );
		creation[0] += source;
		for( level=0; level < numlevels_local; level++ )
		{
			z[level][level] += sink;
		}

		/* >>chng 04 nov 30, atom XX-like collisions off turns this off too */
		if( secondaries.Hx12[ipISO][nelem][iso.nLyaLevel[ipISO]] * iso.lgColl_excite[ipISO] > 0. )
		{
			/* now add on supra thermal excitation */
			for( level=1; level < numlevels_local; level++ )
			{
				double RateUp , RateDown;

				RateUp = secondaries.Hx12[ipISO][nelem][level];
				RateDown = RateUp * (double)StatesElemNEW[nelem][nelem-ipISO][ipH1s].g /
					(double)StatesElemNEW[nelem][nelem-ipISO][level].g;

				/* total rate out of lower level */
				z[ipH1s][ipH1s] += RateUp;

				/* rate from the upper level to ground */
				z[level][ipH1s] -= RateDown;

				/* rate from ground to upper level */
				z[ipH1s][level] -= RateUp;

				z[level][level] += RateDown;  
			}
		}

		/* =================================================================== 
		 *
		 * at this point all matrix elements have been established 
		 *
		 * ==================================================================== */
		/* save matrix, this allocates SaveZ */
		SaveZ = z;

		for( ipLo=0; ipLo < numlevels_local; ipLo++ )
			Save_creation[ipLo] = creation[ipLo];

		if( trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "  pop level     others => (iso_level)\n" );
			for( long n=0; n < numlevels_local; n++ )
			{
				fprintf( ioQQQ, "  %s %s %2ld", iso.chISO[ipISO], elementnames.chElementNameShort[nelem], n );
				for( long j=0; j < numlevels_local; j++ )
				{
					fprintf( ioQQQ,"\t%.9e", z[j][n] );
				}
				fprintf( ioQQQ, "\n" );
			}
			fprintf(ioQQQ," recomb          ");
			for( long n=0; n < numlevels_local; n++ )
			{
				fprintf( ioQQQ,"\t%.9e", creation[n] );
			}
			fprintf( ioQQQ, "\n" );
			fprintf(ioQQQ," recomb ct %.2e co %.2e hectr %.2e hctr %.2e\n",
				atmdat.HeCharExcRecTotal,
				findspecies("CO")->hevmol ,
				atmdat.HeCharExcRecTo[nelem][nelem-ipISO]*dense.xIonDense[ipHELIUM][0] ,
				atmdat.HCharExcRecTo[nelem][nelem-ipISO]*dense.xIonDense[ipHYDROGEN][0] );
		}

		nerror = 0;

		getrf_wrapper(numlevels_local,numlevels_local,
			      z.data(),numlevels_local,&ipiv[0],&nerror);

		getrs_wrapper('N',numlevels_local,1,z.data(),numlevels_local,&ipiv[0],
			      &creation[0],numlevels_local,&nerror);

		if( nerror != 0 )
		{
			fprintf( ioQQQ, " iso_level: dgetrs finds singular or ill-conditioned matrix\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* check whether solution is valid */
		/* >>chng 06 aug 28, both of these from numLevels_max to _local. */
		for( level=ipH1s; level < numlevels_local; level++ )
		{
			double qn = 0., qx = 0.;
			error[level] = 0.;
			for( long n=ipH1s; n < numlevels_local; n++ )
			{
				double q = SaveZ[n][level]*creation[n];
				
				/* remember the largest size of element in sum to div by below */
				if ( q > qx )
					qx = q;
				else if (q < qn)
					qn = q;

				error[level] += q;
			}
			
			if (-qn > qx)
				qx = -qn;

			if( qx > 0. )
			{
				error[level] = (error[level] - Save_creation[level])/qx;
			}
			else
			{
				error[level] = 0.;
			}
		}

		/* remember largest residual in matrix inversion */
		BigError = -1.;
		level_error = -1;
		/* >>chng 06 aug 28, from numLevels_max to _local. */
		for( level=ipH1s; level < numlevels_local; level++ )
		{
			double abserror;
			abserror = fabs( error[level]);
			/* this will be the largest residual in the matrix inversion */
			if( abserror > BigError )
			{
				BigError = abserror;
				level_error = level;
			}
		}

		/* matrix inversion should be nearly as good as the accuracy of a double,
		 * but demand that it is better than epsilon for a float */
		if( BigError > FLT_EPSILON ) 
		{
			if( !conv.lgSearch )
				fprintf(ioQQQ," PROBLEM" );

			fprintf(ioQQQ,
				" iso_level zone %.2f - largest residual in iso=%li %s nelem=%li matrix inversion is %g "
				"level was %li Search?%c \n", 
				fnzone,
				ipISO,
				elementnames.chElementName[nelem],
				nelem , 
				BigError , 
				level_error,
				TorF(conv.lgSearch) );
		}

		/* put level populations into master array */
		for( level=0; level < numlevels_local; level++ )
		{
			StatesElemNEW[nelem][nelem-ipISO][level].Pop = creation[level];

			//ASSERT( StatesElemNEW[nelem][nelem-ipISO][level].Pop < dense.gas_phase[nelem] );

			/* check for negative level populations */
			if( StatesElemNEW[nelem][nelem-ipISO][level].Pop < 0. )
				lgNegPop = true;

			if( StatesElemNEW[nelem][nelem-ipISO][level].Pop <= 0 && !conv.lgSearch )
			{
				fprintf(ioQQQ,
					"PROBLEM non-positive level pop for iso = %li, nelem = "
					"%li = %s, level=%li val=%.3e nTotalIoniz %li\n", 
					ipISO,
					nelem , 
					elementnames.chElementSym[nelem],
					level,
					StatesElemNEW[nelem][nelem-ipISO][level].Pop ,
					conv.nTotalIoniz);
			}

		}

		/* zero populations of unused levels. */
		for( level=numlevels_local; level < iso.numLevels_max[ipISO][nelem]; level++ )
		{
			StatesElemNEW[nelem][nelem-ipISO][level].Pop = 0.;
			/* >>chng 06 jul 25, no need to zero this out, fix limit to 3-body heating elsewhere. */
			/* iso.PopLTE[ipISO][nelem][level] = 0.; */
		}

		/* TotalPop is sum of level pops */
		TotalPop = 0.;
		/* create sum of populations */
		for( level=0; level < numlevels_local; level++ )
			TotalPop += StatesElemNEW[nelem][nelem-ipISO][level].Pop;

		/* option to force ionization */
		if( dense.lgSetIoniz[nelem] )
		{
			lgNegPop = false;
			for( level=0; level < numlevels_local; level++ )
			{
				StatesElemNEW[nelem][nelem-ipISO][level].Pop *=
					dense.xIonDense[nelem][nelem-ipISO]/TotalPop;
				lgNegPop = lgNegPop || (StatesElemNEW[nelem][nelem-ipISO][level].Pop < 0.);
			}
			TotalPop = dense.xIonDense[nelem][nelem-ipISO];
		}

		ASSERT( TotalPop >= 0. );

		/* this is total ionization rate, s-1, of this species referenced to
		 * the total abundance */
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = 0.;
		for( level=0; level < numlevels_local; level++ )
		{
			/* sum of all ionization processes from this atom to ion, cm-3 s-1 now,
			 * but is divided by TotalPop below to become s-1 */
			ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] += 
				StatesElemNEW[nelem][nelem-ipISO][level].Pop * iso.RateLevel2Cont[ipISO][nelem][level];
		}

		if( ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] > BIGDOUBLE )
		{
			fprintf( ioQQQ, "DISASTER RateIonizTot for Z=%li, ion %li is larger than BIGDOUBLE.  This is a big problem.",
				nelem+1, nelem-ipISO);
			cdEXIT(EXIT_FAILURE);
		}
		else
			ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] /= MAX2(SMALLDOUBLE , TotalPop);
	}

	/* all solvers end up here */

	if( ionbal.RateRecomTot[nelem][nelem-ipISO] > 0. )
		iso.xIonSimple[ipISO][nelem] = ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1]/ionbal.RateRecomTot[nelem][nelem-ipISO];
	else
		iso.xIonSimple[ipISO][nelem] = 0.;

	ASSERT( ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] >= 0. );

	/* check on the sum of the populations */
	if( lgNegPop || dense.xIonDense[nelem][nelem-ipISO] < 0. )
	{
		fprintf( ioQQQ, 
			" DISASTER iso_level finds negative ion fraction for iso=%2ld nelem=%2ld "\
			"%s using solver %s, IonFrac=%.3e simple=%.3e TE=%.3e ZONE=%4ld\n", 
			ipISO,
			nelem,
			elementnames.chElementSym[nelem],
			iso.chTypeAtomUsed[ipISO][nelem],
			dense.xIonDense[nelem][nelem+1-ipISO]/SDIV(dense.xIonDense[nelem][nelem-ipISO]),
			iso.xIonSimple[ipISO][nelem],
			phycon.te,
			nzone );

		fprintf( ioQQQ, " level pop are:\n" );
		for( i=0; i < numlevels_local; i++ )
		{
			fprintf( ioQQQ,PrintEfmt("%8.1e", StatesElemNEW[nelem][nelem-ipISO][i].Pop ));
			if( (i!=0) && !(i%10) ) fprintf( ioQQQ,"\n" );
		}
		fprintf( ioQQQ, "\n" );
		ContNegative();
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( ipISO == ipHE_LIKE && nelem==ipHELIUM && nzone>0 )
	{
		/* find fraction of He0 destructions due to photoionization of 2^3S */
		double ratio;
		double rateOutOf2TripS = StatesElemNEW[nelem][nelem-ipISO][ipHe2s3S].Pop * iso.RateLevel2Cont[ipISO][nelem][ipHe2s3S];
		if( rateOutOf2TripS > SMALLFLOAT )
		{
			ratio = rateOutOf2TripS /
				( rateOutOf2TripS + StatesElemNEW[nelem][nelem-ipISO][ipHe1s1S].Pop*ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] );
		}
		else
		{
			ratio = 0.;
		}
		if( ratio > he.frac_he0dest_23S )
		{
			/* remember zone where this happended and fraction, and frac due to photoionization */
			he.nzone = nzone;
			he.frac_he0dest_23S = ratio;
			rateOutOf2TripS = StatesElemNEW[nelem][nelem-ipISO][ipHe2s3S].Pop *iso.gamnc[ipISO][nelem][ipHe2s3S];
			if( rateOutOf2TripS > SMALLFLOAT )
			{
				he.frac_he0dest_23S_photo = rateOutOf2TripS /
					( rateOutOf2TripS + StatesElemNEW[nelem][nelem-ipISO][ipHe1s1S].Pop*ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] );
			}
			else
			{
				he.frac_he0dest_23S_photo = 0.;
			}
		}
	}

	for( ipHi=1; ipHi<numlevels_local; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
				continue;

			/* population of lower level, corrected for stimulated emission */
			Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc = 
				StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop - StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop*
				StatesElemNEW[nelem][nelem-ipISO][ipLo].g/StatesElemNEW[nelem][nelem-ipISO][ipHi].g;

			/* don't allow masers from collapsed levels */
			if( N_(ipHi) > iso.n_HighestResolved_local[ipISO][nelem] )
				Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc = StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop;
		}
	}

	return;
}
