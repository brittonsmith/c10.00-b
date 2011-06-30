/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_levelN compute an arbitrary N level atom */
#include "cddefines.h"
#include "physconst.h"
#include "thermal.h"
#include "conv.h"
#include "phycon.h"
#include "dense.h"
#include "trace.h"
#include "thirdparty.h"
#include "atoms.h"

void atom_levelN(
	/* nlev is the number of levels to compute*/ 
	long int nLevelCalled,
	/* ABUND is total abundance of species, used for nth equation
	 * if balance equations are homogeneous */
	realnum abund, 
	/* G(nlev) is stat weight of levels */
	const double g[], 
	/* EX(nlev) is excitation potential of levels, deg K or wavenumbers
	 * 0 for lowest level, all are energy rel to ground NOT d(ENER) */
	const double ex[], 
	/* this is 'K' for ex[] as Kelvin deg, is 'w' for wavenumbers */
	char chExUnits,
	/* populations [cm-3] of each level as deduced here */
	double pops[], 
	/* departure coefficient, derived below */
	double depart[],
	/* net transition rate, A * esc prob, s-1 */
	double ***AulEscp, 
	/* col str from up to low */
	double ***col_str, 
	/* AulDest[ihi][ilo] is destruction rate, trans from ihi to ilo, A * dest prob,
	 * asserts confirm that [ihi][ilo] is zero */
	double ***AulDest, 
	/* AulPump[ilo][ihi] is pumping rate from lower to upper level (s^-1), (hi,lo) must be zero  */
	double ***AulPump, 
	/* collision rates (s^-1), evaluated here and returned for cooling by calling function,
	 * unless following flag is true.  If true then calling function has already filled
	 * in these rates.  CollRate[i][j] is rate from i to j */
	double ***CollRate,
	/* this is an additional creation rate from continuum, normally zero, units cm-3 s-1 */
	const double source[] ,
	/* this is an additional destruction rate to continuum, normally zero, units s-1 */
	const double sink[] ,
	/* flag saying whether CollRate already done, or we need to do it here,
	 * this is stored in data)[ihi][ilo] as either downward rate or collis strength*/
	bool lgCollRateDone,
	/* total cooling and its derivative, set here but nothing done with it*/
	double *cooltl, 
	double *coolder, 
	/* string used to identify calling program in case of error */
	const char *chLabel, 
	/* nNegPop flag indicating what we have done
	 * positive if negative populations occurred
	 * zero if normal calculation done
	 * negative if too cold, matrix not solved, since highest level had little excitation
	 * (for some atoms other routine will be called in this case) */
	int *nNegPop,
	/* true if populations are zero, either due to zero abundance of very low temperature */
	bool *lgZeroPop ,
	/* option to print debug information */
	bool lgDeBug )
{
	bool lgHomogeneous;

	long int level, 
	  ihi, 
	  ilo, 
	  j; 

	int32 ner;

	double cool1,
	  TeInverse,
	  TeConvFac,
	  sum;

	DEBUG_ENTRY( "atom_levelN()" );

	*nNegPop = -1;

	/* >>chng 05 dec 14, units of ex[] can be Kelvin (old default) or wavenumbers */
	if( chExUnits=='K' )
	{
		/* ex[] is in temperature units - this will multiply ex[] to
		 * obtain Boltzmann factor */
		TeInverse = 1./phycon.te;
		/* this multiplies ex[] to obtain energy difference between levels */
		TeConvFac = 1.;
	}
	else if( chExUnits=='w' )
	{
		/* ex[] is in wavenumber units */
		TeInverse = 1./phycon.te_wn;
		TeConvFac = T1CM;
	}
	else
		TotalInsanity();

	long int nlev = nLevelCalled;
	// decrement number of levels until we have positive excitation rate,
	while( ( (dsexp((ex[nlev-1]-ex[0])*TeInverse) + (*AulPump)[0][nlev-1]) < SMALLFLOAT ) &&
			source[nlev-1]==0. && nlev > 1)
	{
		pops[nlev-1] = 0.;
		depart[nlev-1] = 0.;
		--nlev;
	}

	/* exit if zero abundance or all population in ground */
	ASSERT( abund>= 0. );
	if( abund == 0. || nlev==1 )
	{
		*cooltl = 0.;
		*coolder = 0.;
		/* says calc was ok */
		*nNegPop = 0;
		*lgZeroPop = true;

		pops[0] = abund;
		depart[0] = 1.;
		for( level=1; level < nlev; level++ )
		{
			pops[level] = 0.;
			depart[level] = 0.;
		}

		/* there are TWO abort returns in this sub,
		 * this one is for zero abundance */
		return;
	}

	// these are all automatically deallocated when they go out of scope
	auto_vec<int32> ipiv( new int32[nlev] );
	auto_vec<double> bvec( new double[nlev] );
	multi_arr<double,2,C_TYPE> amat(nlev,nlev); 
	multi_arr<double,2> excit(nLevelCalled,nLevelCalled);

#	ifndef NDEBUG
	/* excitation temperature of lowest level must be zero */
	ASSERT( ex[0] == 0. );

	for( ihi=1; ihi < nlev; ihi++ )
	{
		for( ilo=0; ilo < ihi; ilo++ )
		{
			/* following must be zero:
			 * AulDest[ilo][ihi] - so that spontaneous transitions only proceed from high energy to low
			 * AulEscp[ilo][ihi] - so that spontaneous transitions only proceed from high energy to low
			 * AulPump[ihi][ilo] - so that pumping only proceeds from low energy to high */
			ASSERT( (*AulDest)[ilo][ihi] == 0. );
			ASSERT( (*AulEscp)[ilo][ihi] == 0 );
			ASSERT( (*AulPump)[ihi][ilo] == 0. );

			ASSERT( (*AulEscp)[ihi][ilo] >= 0 );
			ASSERT( (*AulDest)[ihi][ilo] >= 0 );
			ASSERT( (*col_str)[ihi][ilo] >= 0 );
		}
	}
#	endif

	if( lgDeBug || (trace.lgTrace && trace.lgTrLevN) )
	{
		fprintf( ioQQQ, " atom_levelN trace printout for atom=%s with tot abund %e \n", chLabel, abund);
		fprintf( ioQQQ, " AulDest\n" );

		fprintf( ioQQQ, "  hi  lo" );

		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (*AulDest)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " A*esc\n" );
		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (*AulEscp)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " AulPump\n" );

		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (*AulPump)[ilo][ihi] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " coll str\n" );
		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (*col_str)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " coll rate\n" );
		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		if( lgCollRateDone )
		{
			for( ihi=1; ihi < nlev; ihi++ )
			{
				fprintf( ioQQQ, "%3ld", ihi );
				for( ilo=0; ilo < ihi; ilo++ )
				{
					fprintf( ioQQQ, "%10.2e", (*CollRate)[ihi][ilo] );
				}
				fprintf( ioQQQ, "\n" );
			}
		}
	}

	/* find set of Boltzmann factors
	 * evaluate full range of levels since calling routine will use these values
	 * */
	for( ilo=0; ilo < (nLevelCalled - 1); ilo++ )
	{
		for( ihi=ilo + 1; ihi < nLevelCalled; ihi++ )
		{
			/* >>chng 05 dec 14, option to have ex be either Kelvin or wavenumbers */
			excit[ilo][ihi] = dsexp((ex[ihi]-ex[ilo])*TeInverse);
		}
	}

	if( trace.lgTrace && trace.lgTrLevN )
	{
		fprintf( ioQQQ, " excit, te=%10.2e\n", phycon.te );
		fprintf( ioQQQ, "  hi  lo" );

		for( ilo=0; ilo < (nlev-1); ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", excit[ilo][ihi] );
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* we will predict populations */
	*lgZeroPop = false;

	/* already have excitation pumping, now get deexcitation */
	for( ilo=0; ilo < (nlev - 1); ilo++ )
	{
		for( ihi=ilo + 1; ihi < nlev; ihi++ )
		{
			/* (*AulPump)[low][ihi] is excitation rate, low -> ihi, due to external continuum,
			 * so derive rate from upper to lower */
			(*AulPump)[ihi][ilo] = (*AulPump)[ilo][ihi]*g[ilo]/g[ihi];
		}
	}

	/* evaluate collision rates from collision strengths, but only if calling
	 * routine has not already done this
	 * evaluate full range of levels since calling routine will use these rates
	 */
	if( !lgCollRateDone )
	{
		for( ilo=0; ilo < (nLevelCalled - 1); ilo++ )
		{
			for( ihi=ilo + 1; ihi < nLevelCalled; ihi++ )
			{
				/* this should be a collision strength */
				ASSERT( (*col_str)[ihi][ilo]>= 0. );
				/* this is deexcitation rate */
				(*CollRate)[ihi][ilo] = (*col_str)[ihi][ilo]/g[ihi]*dense.cdsqte;
				/* this is excitation rate */
				(*CollRate)[ilo][ihi] = (*CollRate)[ihi][ilo]*g[ihi]/g[ilo]*
					excit[ilo][ihi];
			}
		}
	}

	/* rate equations equal zero */
	amat.zero();

	/* following is column of vector - represents source terms from elsewhere,
	 * if this is zero then matrix is singular and must replace one row with
	 * population sum equation - if sum is non-zero then get total abundance
	 * from source and sink terms */
	sum = 0.;
	lgHomogeneous = false;
	for( level=0; level < nlev; level++ )
	{
		bvec[level] = source[level];
		sum += bvec[level];
	}
	if( sum==0. )
		lgHomogeneous = true;

	/* eqns for destruction of level
	 * AulEscp[iho][ilo] are A*esc p, CollRate is coll excit in direction
	 * AulPump[low][high] is excitation rate from low to high */
	for( level=0; level < nlev; level++ )
	{
		amat[level][level] = sink[level];

		/* leaving level to below */
		for( ilo=0; ilo < level; ilo++ )
		{
			amat[level][level] += (*CollRate)[level][ilo] + (*AulEscp)[level][ilo] +
				(*AulDest)[level][ilo] + (*AulPump)[level][ilo];
			/* >>chng 97 jul 31, added pumping down
			 * coming to level I from below */
			amat[ilo][level] = -(*CollRate)[ilo][level] - (*AulPump)[ilo][level];
		}

		/* leaving level to above */
		for( ihi=level + 1; ihi < nlev; ihi++ )
		{
			amat[level][level] += (*CollRate)[level][ihi] + (*AulPump)[level][ihi];
			/* coming to level from above */
			amat[ihi][level] = -(*CollRate)[ihi][level] - (*AulEscp)[ihi][level] -
				(*AulDest)[ihi][level] - (*AulPump)[ihi][level];
			/* >>chng 97 jul 31, added pumping down */
		}
	}

	/* homogeneous case, all source terms add up to zero, so use the population,
	 * system has no total population information,
	 * equation to replace redundant equation */
	if( lgHomogeneous )
	{
		for( level=0; level<nlev; ++level )
		{
			amat[level][0] = 1.0;
		}
		/* these add up to total abundance */
		bvec[0] = abund;
	}

	if( lgDeBug )
	{
		fprintf(ioQQQ," amat matrix follows:\n");
		for( level=0; level < nlev; level++ )
		{
			for( j=0; j < nlev; j++ )
			{
				fprintf(ioQQQ," %.4e" , amat[level][j]);
			}
			fprintf(ioQQQ,"\n");
		}
		if( sum==0. )
		{
			fprintf(ioQQQ," Sum creation zero so pop sum used\n");
		}
		else
		{
			fprintf(ioQQQ," Sum creation non-zero (%e), vector follows:\n",sum);
			for( j=0; j < nlev; j++ )
			{
				fprintf(ioQQQ," %.4e" , bvec[j] );
			}
			fprintf(ioQQQ,"\n");
		}
	}

	ner = 0;
	getrf_wrapper(nlev,nlev,amat.data(),nlev,ipiv.data(),&ner);
	/* usage DGETRS, 'N' = no transpose
		* order of matrix,
		* number of cols in bvec, =1
		* array
		* leading dim of array */
	getrs_wrapper('N',nlev,1,amat.data(),nlev,ipiv.data(),bvec.data(),nlev,&ner);

	if( ner != 0 )
	{
		fprintf( ioQQQ, " atom_levelN: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* set populations */
	for( level=0; level < nlev; level++ )
	{
		/* save bvec into populations */
		pops[level] = bvec[level];
	}

	/* now find total energy exchange rate, normally net cooling and its 
	 * temperature derivative */
	*cooltl = 0.;
	*coolder = 0.;
	for( ihi=1; ihi < nlev; ihi++ )
	{
		for( ilo=0; ilo < ihi; ilo++ )
		{
			/* this is now net cooling rate [K cm-3 s-1] */
			cool1 = (pops[ilo]*(*CollRate)[ilo][ihi] - pops[ihi]*(*CollRate)[ihi][ilo])*
			  (ex[ihi] - ex[ilo]);
			*cooltl += cool1;

			/* derivative wrt temperature - use Boltzmann factor relative to ground */
			/* >>chng 03 aug 28, use real cool1 */
			*coolder += cool1*( (ex[ihi] - ex[0])*thermal.tsq1 - thermal.halfte);
		}
	}
	/* convert from units of ex[] into ergs */
	/* >>chng 05 dec 14, ex[] may be K or wn, TeConvFac will take care of either case */
	*cooltl *= BOLTZMANN*TeConvFac;
	*coolder *= BOLTZMANN*TeConvFac;

	/* fill in departure coefficients */
	if( pops[0] > SMALLFLOAT && excit[0][nlev-1] > SMALLFLOAT )
	{
		/* >>chng 00 aug 10, loop had been from 1 and 0 was set to total abundance */
		depart[0] = 1.;
		for( ihi=1; ihi < nlev; ihi++ )
		{
			/* these are off by one - lowest index is zero */
			depart[ihi] = (pops[ihi]/pops[0])*(g[0]/g[ihi])/excit[0][ihi];
		}
	}

	else
	{
		/* >>chng 00 aug 10, loop had been from 1 and 0 was set to total abundance */
		for( ihi=0; ihi < nlev; ihi++ )
		{
			/* Boltzmann factor or abundance too small, set departure coefficient to zero */
			depart[ihi] = 0.;
		}
		depart[0] = 1.;
	}

	/* sanity check for valid solution - non negative populations */
	*nNegPop = 0;
	/* the limit we allow the fractional population to go below zero before announcing failure. */
	double poplimit = 1.0E-10;
	for( level=0; level < nlev; level++ )
	{
		if( pops[level] < 0. )
		{
			if( fabs(pops[level]/abund) > poplimit )
			{
				//nNegPop = 1 leads to a failure
				*nNegPop = 1;
			}
			else
			{
				pops[level] = SMALLFLOAT;
				//fprintf( ioQQQ, "\n PROBLEM Small negative populations were found in atom = %s . "
				//		"The problem was ignored and the negative populations were set to SMALLFLOAT",chLabel );
			}
		}
	}

	if( *nNegPop!=0 )
	{
		ASSERT( *nNegPop==1 );
		// *nNegPop 0 is valid solution, nNegPop> 1 negative populations found
		fprintf( ioQQQ, "\n PROBLEM atom_levelN found negative population, nNegPop=%i3, atom=%s lgSearch=%c\n",
				*nNegPop , chLabel , TorF( conv.lgSearch ) );

		for( level=0; level < nlev; level++ )
		{
			fprintf( ioQQQ, "%10.2e", pops[level] );
		}

		fprintf( ioQQQ, "\n" );
		for( level=0; level < nlev; level++ )
		{
			pops[level] = (double)MAX2(0.,pops[level]);
		}
	}

	if(  lgDeBug || (trace.lgTrace && trace.lgTrLevN) )
	{
		fprintf( ioQQQ, "\n atom_leveln absolute population   " );
		for( level=0; level < nlev; level++ )
		{
			fprintf( ioQQQ, " %10.2e", pops[level] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " departure coefficients" );
		for( level=0; level < nlev; level++ )
		{
			fprintf( ioQQQ, " %10.2e", depart[level] );
		}
		fprintf( ioQQQ, "\n\n" );
	}

#	ifndef NDEBUG
	/* these were reset to non zero values by the solver, but we will
	 * assert that they are zero (for safety) when routine reenters so must
	 * set to zero here, but only if asserts are in place */
	for( ihi=1; ihi < nlev; ihi++ )
	{
		for( ilo=0; ilo < ihi; ilo++ )
		{
			/* zero ots destruction rate */
			(*AulDest)[ilo][ihi] = 0.;
			/* both AulDest and AulPump (low, hi) are not used, should be zero */
			(*AulPump)[ihi][ilo] = 0.;
			(*AulEscp)[ilo][ihi] = 0;
		}
	}
#	endif

	// -1 return had meant too cool and no populations determined, no longer used
	// since we now decrement until populations can be determined
	ASSERT( *nNegPop>=0 );

	return;
}
