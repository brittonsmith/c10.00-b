/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*dense_tabden interpolate on table of points for density with dlaw table command, by K Volk */
#include "cddefines.h"
#include "dense.h"

double dense_tabden(double r0, 
  double depth)
{
	bool lgHit;
	long int j;
	double frac, 
	  tabden_v, 
	  x;

	DEBUG_ENTRY( "dense_tabden()" );
	/*interpolate on table of points for density with dlaw table command, by K Volk
	 *each line is log radius and H density per cc. */

	/*begin sanity check */
	if( r0 <= 0. || depth <= 0. )
	{
		fprintf( ioQQQ, " dense_tabden called with insane depth, radius, =%10.2e%10.2e\n", 
		  depth, r0 );
	}
	/*end sanity check */

	/* interpolate on radius or depth? */
	if( dense.lgDLWDepth )
	{
		/* depth key appeared = we want depth */
		x = log10(depth);
	}
	else
	{
		/* use radius */
		x = log10(r0);
	}

	/* set to impossible value, will crash if not reset */
	tabden_v = -DBL_MAX;

	if( x < dense.frad[0] || x >= dense.frad[dense.nvals-1] )
	{
		fprintf( ioQQQ, " requested radius outside range of dense_tabden\n" );
		fprintf( ioQQQ, " radius was%10.2e min, max=%10.2e%10.2e\n", 
		  x, dense.frad[0], dense.frad[dense.nvals-1] );
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		lgHit = false;
		j = 1;

		while( !lgHit && j <= dense.nvals - 1 )
		{
			if( dense.frad[j-1] <= (realnum)x && dense.frad[j] > (realnum)x )
			{
				frac = (x - dense.frad[j-1])/(dense.frad[j] - 
				  dense.frad[j-1]);
				tabden_v = dense.fhden[j-1] + frac*(dense.fhden[j] - 
				  dense.fhden[j-1]);
				lgHit = true;
			}
			j += 1;
		}

		if( !lgHit )
		{
			fprintf( ioQQQ, " radius outran dlaw table scale, requested=%6.2f largest=%6.2f\n", 
			  x, dense.frad[dense.nvals-1] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* got it, now return value, not log of density */
	tabden_v = pow(10.,tabden_v);
	return( tabden_v );
}
