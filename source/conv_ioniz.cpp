/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvIoniz called by ConvEdenIonz, it calls ConvBase until converged */
#include "cddefines.h"
#include "thermal.h"
#include "trace.h"
#include "conv.h"

/* this routine is called by ConvEdenIoniz, it calls ConvBase
 * until it converges or overruns the loop limit */
int ConvIoniz(void)
{
	DEBUG_ENTRY( "ConvIoniz()" );

	/* expand limit to number of calls to ConvBase during search phase */
	int LoopLimit = conv.lgSearch ? 30 : 20;

	/* do not go into the loop with first call to ionization,
	 * since results will be bogus - do it here */
	if( !conv.lgSearch && conv.nPres2Ioniz == 0 )
	{
		if( ConvBase(0) )
			return 1;
	}

	strcpy( conv.chConvIoniz, " NONE!!!!!" );
	/* this is ionization/electron density convergence loop
	 * keep calling ConvBase until lgIonDone is true */
	for( int i=0; i < LoopLimit; ++i )
	{
		/* compute the current ionization, ots rates, secondary ionization rates */
		if( ConvBase(i) )
			return 1;

		if( trace.nTrConvg >= 4 )
		{
			/* cooling has not been evaluated yet */
			fprintf( ioQQQ, "    ConvIoniz4 %d heat: %.2e cool: %.2e ", 
				 i, thermal.htot , thermal.ctot );

			/* this is flag saying whether or not ionization/eden has converged */
			if( conv.lgConvIoniz )
			{
				fprintf( ioQQQ, " ioniz converged\n" );
			}
			else
			{
				fprintf( ioQQQ, " ioniz no conv: %s old %.4e new %.4e OscilOTS %c\n", 
				  conv.chConvIoniz , 
				  conv.BadConvIoniz[0] ,
				  conv.BadConvIoniz[1] ,
				  TorF(conv.lgOscilOTS) );
			}
		}

		if( conv.lgConvIoniz || lgAbort )
			break;
	}

	/* if trace convergence is in operation and we did not converge, give warning */
	if( trace.nTrConvg >= 4 && !conv.lgConvIoniz )
		fprintf( ioQQQ, 
			 "    ConvIoniz4>>>>>>>>>>exit without converging after %i tries!!!!\n", LoopLimit );

	return 0;
}
