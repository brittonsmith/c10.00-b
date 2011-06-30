/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*gridXspec handles all grid calculations, called by griddo */
/*gridFunc */
/*GridGatherInCloudy - gathers spectra for each grid calculation to save for end */
#include "cddefines.h"
#include "save.h"
#include "warnings.h"
#include "optimize.h"
#include "cddrive.h"
#include "continuum.h"
#include "rfield.h"
#include "grid.h"
#include "ipoint.h"
#include "called.h"
#include "physconst.h"
#include "prt.h"
#include "mpi_utilities.h"

/*gridXspec handles all grid calculations, called by grid_do */
void gridXspec(realnum xc[], long int nInterpVars)
{
	long int i;

	DEBUG_ENTRY( "gridXspec()" );

	if( nInterpVars > LIMPAR )
	{
		fprintf( ioQQQ, "grid_do: too many parameters are varied, increase LIMPAR\n" );
		cdEXIT(EXIT_FAILURE);
	}

	optimize.nOptimiz = 0;
	grid.nintparm = nInterpVars;

	/* if this is changed there must be some change made to actually
	 * stuff the additional parameter information. */
	grid.naddparm = 0;

	ASSERT( grid.nintparm + grid.naddparm <= LIMPAR );

	grid.totNumModels = 1;
	/* >>chng 06 aug 21, allow the number of parameter values to be different for different parameters. */
	for( i=0; i<nInterpVars; i++ )
	{
		/* >>chng 06 sep 4, use grid variable instead of passing to routine. */
		grid.totNumModels *= grid.numParamValues[i];
	}
	/* grid.totNumModels = (long)pow((double)nParVals, (double)nInterpVars); */
	ASSERT( grid.totNumModels > 1 );

	grid.paramNames = (char**)MALLOC(sizeof(char*)*(unsigned)(nInterpVars+grid.naddparm) );
	grid.paramMethods = (long*)MALLOC(sizeof(long)*(unsigned)(nInterpVars+grid.naddparm) );
	grid.paramRange = (realnum**)MALLOC(sizeof(realnum*)*(unsigned)(nInterpVars+grid.naddparm) );
	grid.paramData = (realnum**)MALLOC(sizeof(realnum*)*(unsigned)(nInterpVars+grid.naddparm) );
	grid.interpParameters = (realnum**)MALLOC(sizeof(realnum*)*(unsigned)(grid.totNumModels) );

	for( i=0; i<nInterpVars+grid.naddparm; i++ )
	{
		grid.paramNames[i] = (char*)MALLOC(sizeof(char)*(unsigned)(12) );
		grid.paramRange[i] = (realnum*)MALLOC(sizeof(realnum)*(unsigned)(6) );
		grid.paramData[i] = (realnum*)MALLOC(sizeof(realnum)*(unsigned)(grid.numParamValues[i]) );

		sprintf( grid.paramNames[i],	"%s%ld", "PARAM", i+1 );
		/* Method is 0 for linear, 1 for logarithmic */
		grid.paramMethods[i] = grid.lgLinearIncrements[i] ? 0 : 1;
		/* Initial */
		grid.paramRange[i][0] = xc[i]+grid.paramIncrements[i]*(grid.numParamValues[i]-1.f)/2.f;
		/* Delta */
		grid.paramRange[i][1] = grid.paramIncrements[i]/10.f;
		/* Minimum */
		grid.paramRange[i][2] = xc[i];
		/* Bottom */
		grid.paramRange[i][3] = xc[i]+grid.paramIncrements[i]/10.f;
		/* Top */
		grid.paramRange[i][4] = xc[i]+grid.paramIncrements[i]*(grid.numParamValues[i]-1.f)-grid.paramIncrements[i]/10.f;
		/* Maximum */
		grid.paramRange[i][5] = xc[i]+grid.paramIncrements[i]*(grid.numParamValues[i]-1.f);

		for( long j=0; j<grid.numParamValues[i]; j++ )
		{
			grid.paramData[i][j] = xc[i]+grid.paramIncrements[i]*j;
		}
	}

	for( i=0; i<grid.totNumModels; i++ )
	{
		grid.interpParameters[i] = (realnum*)MALLOC(sizeof(realnum)*(unsigned)(nInterpVars) );
	}

	for( i=0; i< grid.totNumModels; i++ )
	{
		long j;
		realnum variableVector[LIMPAR];

		for( j=0; j<nInterpVars; j++ )
		{
			int index;
			long volumeOtherDimensions = 1;

			/* by "volume", we simply mean the product of the parameter dimensions 
			 * AFTER the present index, i.e.:
			 * first "volume" is product of grid.numParamValues[1]*grid.numParamValues[2]*....grid.numParamValues[n]
			 * second "volume" is product of grid.numParamValues[2]*grid.numParamValues[3]*....grid.numParamValues[n]
			 * last "volume" is unity.  */
			for( long k=j+1; k<nInterpVars; k++ )
			{
				volumeOtherDimensions *= grid.numParamValues[k];
			}

			/* For each successive parameter, the "volume" is less than the previous one.
			 * So the left-hand side of this modulus operation increases for each parameter,
			 * which causes the index of each parameter to be incremented more often than the
			 * index of the previous parameter.  Thus, the last dimension is incremented
			 * every time, then the second to last dimension is incremented with each repeat
			 * of the last dimension.  This repeats until, finally, the first dimension is
			 * incremented.  For example, the indices of the parameter vectors for a 2x2x3
			 * box would be ordered as such:
			 * [0][0][0]
			 * [0][0][1]
			 * [0][0][2]
			 * [0][1][0]
			 * [0][1][1]
			 * [0][1][2]
			 * [1][0][0]
			 * [1][0][1]
			 * [1][0][2]
			 * [1][1][0]
			 * [1][1][1]
			 * [1][1][2]
			 */
			index = (int)( (i/volumeOtherDimensions)%(grid.numParamValues[j]) );

			/* this prevents parameter incrementation for debugging purposes.  */
			if( grid.lgStrictRepeat )
				variableVector[j] = xc[j];
			else
				variableVector[j] = xc[j] + grid.paramIncrements[j]*index;

			grid.interpParameters[i][j] = variableVector[j];

			if( grid.lgLinearIncrements[j] && !optimize.lgOptimizeAsLinear[j] )
				variableVector[j] = log10(variableVector[j]);
		}

		for( j=nInterpVars; j<LIMPAR; j++ )
		{
			variableVector[j] = xc[j];
		}

		if( i == grid.totNumModels - 1 )
		{
			fixit(); // is this needed ??
			called.lgTalk = cpu.lgMPI_talk();
			called.lgTalkIsOK = cpu.lgMPI_talk();
			prt.lgFaintOn = true;
			grid.lgGridDone = true;
		}

		(void)optimize_func(variableVector,i);
	}
	return;
}

/*GridGatherInCloudy - gathers spectra for each grid calculation to save for end */
void GridGatherInCloudy( void )
{
	long i;

	DEBUG_ENTRY( "GridGatherInCloudy()" );

	ASSERT( grid.lgGrid );

	/* malloc some arrays if first call and save continuum energies. */
	if( grid.Energies.empty() )
	{
		long i1, i2;

		// this will happen if we have more MPI ranks than grid points
		// the highest ranks will not have executed any model
		if( rfield.nupper <= 0 )
			ContCreateMesh();

		if( grid.LoEnergy_keV == 0. )
			grid.ipLoEnergy = 0;
		else
			grid.ipLoEnergy = ipoint( grid.LoEnergy_keV * 1000. / EVRYD );

		if( grid.HiEnergy_keV == 0. || grid.HiEnergy_keV >= continuum.filbnd[continuum.nrange] )
			grid.ipHiEnergy = rfield.nupper - 1;
		else
			grid.ipHiEnergy = ipoint( grid.HiEnergy_keV * 1000. / EVRYD );

		grid.numEnergies = grid.ipHiEnergy - grid.ipLoEnergy + 1;
		ASSERT( grid.numEnergies > 0 );
		grid.Energies.resize( grid.numEnergies );
		grid.Spectra.reserve(NUM_OUTPUT_TYPES);
		for( i1=0; i1 < NUM_OUTPUT_TYPES; i1++ )
		{
			if( grid.lgOutputTypeOn[i1] )
			{
				grid.Spectra.reserve(i1,grid.totNumModels);
				for( i2=0; i2 < grid.totNumModels; i2++ )
				{
					grid.Spectra.reserve(i1,i2,grid.numEnergies);
				}
			}
		}
		grid.Spectra.alloc();
		// this needs to be zeroed out for MPI runs
		grid.Spectra.zero();

		for( i1=0; i1<grid.numEnergies; i1++ )
		{
			long j = grid.ipLoEnergy + i1;
			grid.Energies[i1] = rfield.AnuOrg[j];
		}
	}

	if( optimize.nOptimiz < grid.totNumModels )
	{
		ASSERT( optimize.nOptimiz >= 0 );

		for( i=0; i < NUM_OUTPUT_TYPES; i++ )
		{
			/* Grab spectrum for xspec printout 
			 * at this point nOptimiz has already been incremented for first model */
			if( grid.lgOutputTypeOn[i] )
				cdSPEC2( i, grid.numEnergies, grid.ipLoEnergy, grid.ipHiEnergy,
					 &grid.Spectra[i][optimize.nOptimiz][0]);
		}
	}
	else if( optimize.nOptimiz == grid.totNumModels )
	{
		ASSERT( cpu.lgMPI() );

		multi_arr<realnum,3> Spectra_Copy = grid.Spectra;

		// combine the grid.Spectra data from all ranks. This is done by adding up
		// the results from all ranks. All but one should be zero. This is needed
		// because we do not know which rank calculated which grid point...
		for( int i=0; i < NUM_OUTPUT_TYPES; ++i )
		{
			if( grid.lgOutputTypeOn[i] )
			{
				for( int j=0; j < grid.totNumModels; ++j )
				{
					MPI::COMM_WORLD.Reduce( &Spectra_Copy[i][j][0], &grid.Spectra[i][j][0], 
								grid.numEnergies, MPI::type(grid.Spectra[i][j][0]),
								MPI::SUM, 0 );
				}
			}
		}
		MPI::COMM_WORLD.Barrier();
	}
	else
	{
		TotalInsanity();
	}
	return;
}
