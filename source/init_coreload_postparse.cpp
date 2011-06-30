/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitCoreloadPostparse initialization at start, called from cloudy
* after parser one time per core load */
#include "cddefines.h" 
#include "monitor_results.h"
#include "dense.h"
#include "init.h" 
#include "iso.h"
#include "lines_service.h"
#include "taulines.h"

/*InitCoreloadPostparse initialization at start, called from cloudy
* after parser, one time per core load */
void InitCoreloadPostparse( void )
{

	static int nCalled = 0;

	DEBUG_ENTRY( "InitCoreloadPostparse()" );

	/* only do this once per coreload */
	if( nCalled > 0 )
	{
		return;
	}

	/* this is first call, increment the nCalled counter so we never do this again */
	++nCalled;

	MonitorResults.SumErrorCaseMonitor = 0.;
	MonitorResults.nSumErrorCaseMonitor = 0;

	StatesElemNEW.reserve( LIMELM );

	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		/* only grab core for elements that are turned on */
		if( nelem < 2 || dense.lgElmtOn[nelem] )
		{
			StatesElemNEW.reserve( nelem, nelem+1 );
			for( long ion=0; ion<nelem+1; ion++ )
			{
				long ipISO = nelem-ion;
				if( ipISO < NISO )
				{
					iso_update_num_levels( ipISO, nelem );
					ASSERT( iso.numLevels_max[ipISO][nelem] > 0 );
					StatesElemNEW.reserve( nelem, ion, iso.numLevels_max[ipISO][nelem] );
				}
				else
				{
					fixit();  // for now, point non-iso ions to NULL
					StatesElemNEW.reserve( nelem, ion, 1 );
				}
			}
		}
	}

	StatesElemNEW.alloc();

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( long ipLo=0; ipLo<iso.numLevels_max[ipISO][nelem]; ipLo++ )
				{
					StateJunk( &StatesElemNEW[nelem][nelem-ipISO][ipLo] );
				}
			}
		}
	}

	return;
}
