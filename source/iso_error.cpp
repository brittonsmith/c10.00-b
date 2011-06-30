/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeLikeError fills uncertainty arrays */
#include "cddefines.h" 
#include "iso.h"

/* This routine handles errors when that option is turned on (via the command
 * "atom he-like error generation" */
void iso_put_error(long int ipISO,
			  long int nelem,
			  long int ipHi,
			  long int ipLo,
			  long int whichData,
			  realnum errorOpt,
			  realnum errorPess)
{

	DEBUG_ENTRY( "iso_put_error()" );

	if( iso.lgRandErrGen[ipISO] )
	{
		/* whichData is either IPRAD, IPCOLLIS, or IPENERGY */
		ASSERT( whichData <= 2 );
		ASSERT( ipISO < NISO );
		ASSERT( nelem < LIMELM );
		ASSERT( ipHi <= iso.numLevels_max[ipISO][nelem] );
		ASSERT( ipLo <= iso.numLevels_max[ipISO][nelem] );
		ASSERT( errorOpt >= 0. );
		ASSERT( errorPess >= 0. );
		
		if( !iso.lgPessimisticErrors )
			iso.Error[ipISO][nelem][ipHi][ipLo][whichData] = errorOpt;
		else
			iso.Error[ipISO][nelem][ipHi][ipLo][whichData] = errorPess;
	}
	return;
}

void iso_error_generation( long ipISO, long nelem )
{
	long ipHi, ipLo, typeOfRate;

	DEBUG_ENTRY( "iso_error_generation()" );

	iso.ErrorFactor[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][iso.numLevels_max[ipISO][nelem]][IPRAD] =
		(realnum)MyGaussRand( iso.Error[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][iso.numLevels_max[ipISO][nelem]][IPRAD] );

	for( ipHi=1; ipHi<= iso.numLevels_max[ipISO][nelem]; ipHi++ )
	{
		/* >>chng 06 mar 15, the upper limit incorrectly went to numLevels_max */
		for( ipLo=0; ipLo< ipHi; ipLo++ )
		{
			for( typeOfRate=0; typeOfRate<=1; typeOfRate++ )
			{
				if( iso.Error[ipISO][nelem][ipHi][ipLo][typeOfRate] >= 0. )
				{
					iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][typeOfRate] =  
						(realnum)MyGaussRand( iso.Error[ipISO][nelem][ipHi][ipLo][typeOfRate] );
					ASSERT( iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][typeOfRate] > 0. );
				}
				else
				{
					iso.ErrorFactor[ipISO][nelem][ipHi][ipLo][typeOfRate] = 1.0f;
				}
			}
		}
	}

	/* set flag saying that error generation has been done.  */
	iso.lgErrGenDone[ipISO][nelem] = true;
	return;
}
