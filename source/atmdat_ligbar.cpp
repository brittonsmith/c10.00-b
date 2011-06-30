/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ligbar obtain collision strength for any Li-sequence line */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "phycon.h"
#include "ligbar.h"
#include "transition.h"

void ligbar(long int ized, 
  transition * t2s2p, 
  transition * t2s3p, 
  double *cs2s2p, 
  double *cs2s3p)
{
	double a, 
	  b, 
	  c, 
	  excit, 
	  gbar;

	DEBUG_ENTRY( "ligbar()" );

	/* compute collision strength for Li-seq g-bar approx
	 * summarized in Cochrane and McWhirter Physica Scripta 28, 25.
	 * ized is nuclear charge, can be anything larger than 2
	 * Kirk Korista
	 *
	 * ized is nuclear charge, 6 for carbon, etc
	 * t2s2p is tau array for stronger memeber of 2s 2p multiplet,
	 * which is treated as two separate lines
	 * t2s3p is next transition up, treated as an averaged multiplet
	 *
	 * cs2s2p is the cs for the single 2s2p line that comes in
	 * cs2s3p is the multiplet cs for that transition, which is
	 * treated as  a multiplet average.  If t2s3p is ever separated
	 * (as t2s2p was) then cs2s3p will be the single line not the multiplet
	 *
	 * T2S2P, T2S3P are line information array, defined in block data
	 */

	/* no need to evaluate coll strength if population is zero */
	if(dense.xIonDense[ t2s2p->Hi->nelem -1 ][ t2s2p->Hi->IonStg-1 ] == 0)
	{
		*cs2s2p = 1.;
		*cs2s3p = 1.;
		return;
	}

	if( ized < 3 )
	{
		/* this is a sanity check */
		fprintf( ioQQQ, " LIGBAR called with insane charge, ized=%4ld\n", 
		  ized );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	else if( ized == 6 )
	{
		/* CIV 1549 */
		a = 0.292;
		b = 0.289;
		c = 2.67;
	}

	else if( ized == 7 )
	{
		/* NV 1240 */
		a = 0.387;
		b = 0.247;
		c = 3.93;
	}

	else if( ized == 8 )
	{
		/* OVI 1035 -- values interpolated */
		a = 0.40;
		b = 0.256;
		c = 4.12;
	}

	else if( ized == 10 )
	{
		/* NeVIII 774 */
		a = 0.426;
		b = 0.273;
		c = 4.50;
	}

	else if( ized == 12 )
	{
		/* Mg 10 615 -- these values are general */
		a = 0.45;
		b = 0.27;
		c = 5.0;
	}

	else if( ized == 18 )
	{
		/* Ar 16 365 */
		a = 0.311;
		b = 0.294;
		c = 6.65;
	}

	else if( ized == 26 )
	{
		/* Fe 24 213 */
		a = 0.435;
		b = 0.314;
		c = 6.92;
	}

	else
	{
		/* use general formula for all other cases */
		a = 0.6 - 1.5/((realnum)(ized) - 2.);
		b = 0.27;
		c = 5.;
	}

	/* evaluate expression in terms of coefficients
	 * tarray(ipLnBolt) = line energy in degrees kelvin */
	excit = t2s2p->EnergyK/phycon.te;

	/* excit = e1/(te * 1.380622e-16) */
	gbar = a + b*log(1./excit+c);

	/* tarray(ipLnGF) = gf; tarray(ipLnBolt) excit temp kelvin */
	/*
	 *cs2s2p = gbar*197.47*EVDEGK*t2s2p->Lo->gf/t2s2p->EnergyK;
	 */
	*cs2s2p = gbar*197.47*EVDEGK*t2s2p->Emis->gf/t2s2p->EnergyK;
	/* small correction factors to CMcW83 2s-2p fits:
	 * fits 3.57% too small compared to R-matrix calc. for Mg X.
	 * scaled all, initially, by this constant. Pradhan & Peng (1994)
	 * compilation cites a pc with Burgess, which further scales
	 * cs(C IV) by 1.0429 and cs(N V) by 0.9691, approximately. 
	 * The scaled cs(OVI) matched well with Burgess, so no further
	 * scaling was done for more highly ionized species. */

	if( ized == 6 )
	{
		*cs2s2p *= 1.08013;
	}

	else if( ized == 7 )
	{
		*cs2s2p *= 1.00370;
	}

	else
	{
		*cs2s2p *= 1.0357;
	}


	/* use general formula for 2s3p */
	a = -0.244;
	b = 0.25;
	c = 4.;

	/* excit = e2/(te * 1.380622e-16) */
	excit = t2s3p->EnergyK/phycon.te;
	gbar = a + b*log(1./excit+c);
	/* tarray(ipLnGF) = gf */
	/**cs2s3p = gbar*197.47*EVDEGK*t2s3p->Lo->gf/t2s3p->EnergyK;*/
	*cs2s3p = gbar*197.47*EVDEGK*t2s3p->Emis->gf/t2s3p->EnergyK;
	/* cs2s3p = gbar * 197.47*eVdegK *  GF2/(e2/1.60184e-12)
	 * */
	return;
}
