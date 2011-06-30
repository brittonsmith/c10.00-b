/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseGrid parse the grid command lines */
#include "cddefines.h"
#include "optimize.h"
#include "grid.h"
#include "parser.h"

/* ParseGrid - called from ParseCommands if GRID command found */
void ParseGrid(
	/* command line, which was changed to all caps in main parsing routine */
	Parser &p)
{
	DEBUG_ENTRY( "ParseGrid()" );

	/* RP fake optimizer to run a grid of calculations, also accepts
	 * keyword XSPEC */
	strcpy( optimize.chOptRtn, "XSPE" );
	grid.lgGrid = true;

	if( p.nMatch("REPE") )
	{
		/* just keep repeating, don't actually change the values in the grid.
		 * useful for debugging unintentional crosstalk */
		grid.lgStrictRepeat = true;
	}

	/* 06 aug 22, change to accept three parameters: lower and upper limit and number of points. */
	/* scan off range for the previously selected variable */
	if( optimize.nparm > 0 )
	{
		ASSERT( optimize.nparm <= LIMPAR );

		grid.paramLimits[optimize.nparm-1][0] = (realnum)p.FFmtRead();
		grid.paramLimits[optimize.nparm-1][1] = (realnum)p.FFmtRead();
		grid.paramIncrements[optimize.nparm-1] = (realnum)p.FFmtRead();
		grid.lgLinearIncrements[optimize.nparm-1] = p.nMatch("LINE") ? true : false ;

		if( grid.paramIncrements[optimize.nparm-1] <= 0. )
		{
			fprintf( ioQQQ," The increment (third parameter) must be a positive number.\n" );
			fprintf( ioQQQ," Sorry.\n" );
			cdEXIT( EXIT_FAILURE );
		}

		if( grid.paramIncrements[optimize.nparm-1] > fabs( grid.paramLimits[optimize.nparm-1][1] - grid.paramLimits[optimize.nparm-1][0] ) )
		{
			fprintf( ioQQQ," The increment (third parameter) must not be greater than the difference between the limits (first and second parameters).\n" );
			fprintf( ioQQQ," Sorry.\n" );
			cdEXIT( EXIT_FAILURE );
		}

		if( p.lgEOL() )
		{
			fprintf( ioQQQ," This command has changed since the definition given in Porter et al. 2006, PASP, 118, 920.\n" );
			fprintf( ioQQQ," The grid command now requires three parameters: lower limit, upper limit, and increment.\n" );
			fprintf( ioQQQ," The keywords RANGE and STEPS are no longer necessary.\n" );
			fprintf( ioQQQ," Sorry.\n" );
			cdEXIT( EXIT_FAILURE );
		}
		else
		{
			++optimize.nRangeSet;
		}

		/* Swap if second range parameter is less than the first. */
		if( grid.paramLimits[optimize.nparm-1][1] < grid.paramLimits[optimize.nparm-1][0] )
		{
			realnum temp = grid.paramLimits[optimize.nparm-1][0];
			grid.paramLimits[optimize.nparm-1][0] = grid.paramLimits[optimize.nparm-1][1];
			grid.paramLimits[optimize.nparm-1][1] = temp;
		}

		ASSERT( grid.paramLimits[optimize.nparm-1][1] - grid.paramLimits[optimize.nparm-1][0] > 0. );

		realnum ratio =
			(grid.paramLimits[optimize.nparm-1][1] - grid.paramLimits[optimize.nparm-1][0])/
			grid.paramIncrements[optimize.nparm-1];

		// this takes care of the blowup in the error due to cancellation in limits[1]-limits[0]
		// it assumes that limits[1]-limits[0] is accurate within 3*eps*(limits[0]+limits[1])/2
		// which should be a very conservative estimate...
		realnum feps = realnum(1.5)*
			(grid.paramLimits[optimize.nparm-1][1] + grid.paramLimits[optimize.nparm-1][0])/
			(grid.paramLimits[optimize.nparm-1][1] - grid.paramLimits[optimize.nparm-1][0]);
		long eps = max(nint(abs(feps)),3);

		// this will blow for pathologically narrow grid ranges
		ASSERT( eps <= INT_MAX );

		// take special care if step is integer fraction of max-min (which is nearly always the case)
		if( fp_equal( ratio, realnum(nint(ratio)), int(eps) ) )
			grid.numParamValues[optimize.nparm-1] = nint(ratio) + 1;
		else
			grid.numParamValues[optimize.nparm-1] = long(ratio) + 1;

		if( grid.numParamValues[optimize.nparm-1] < 2 )
			fprintf( ioQQQ, " NOTE must have at least two grid points\n" );

		grid.numParamValues[optimize.nparm-1] = MAX2( 2, grid.numParamValues[optimize.nparm-1] );

		// Create some buffer area in the allowed range of parameter values to prevent
		// accidentally going over the limit due to roundoff error. The buffer is 1/10th
		// of a step, so should still guard against doing too many steps due to bugs.
		realnum safety = 0.1f*grid.paramIncrements[optimize.nparm-1];

		if( grid.lgLinearIncrements[optimize.nparm-1] )
		{
			optimize.varang[optimize.nparm-1][0] = log10(grid.paramLimits[optimize.nparm-1][0]-safety);
			optimize.varang[optimize.nparm-1][1] = log10(grid.paramLimits[optimize.nparm-1][1]+safety);
		}
		else
		{
			optimize.varang[optimize.nparm-1][0] = grid.paramLimits[optimize.nparm-1][0]-safety;
			optimize.varang[optimize.nparm-1][1] = grid.paramLimits[optimize.nparm-1][1]+safety;
		}
	}

	return;
}
