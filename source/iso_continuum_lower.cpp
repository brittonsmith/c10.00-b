/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_continuum_lower - limit max prin. quan. no. due to continuum lowering processes	*/
#include "cddefines.h"
#include "phycon.h"
#include "dense.h"
#include "conv.h"
#include "iso.h"
#include "hydrogenic.h"
#include "trace.h"

void iso_continuum_lower( long ipISO, long nelem )
{
	double a;
	long int np, nd, ns, nc;
	long eff_charge;

	/* size of rate matrices will be defined according to the n calculated here	*/

	ASSERT( dense.xNucleiTotal < MAX_DENSITY );
	ASSERT( nelem < LIMELM );
	/* this may change at a future date. */
	ASSERT( ipISO <= 1 );

	eff_charge = nelem + 1 - ipISO;

	/* Particle packing - the density here should be density of all nuclei in the plasma */
	/* This one is just nuclear charge, which is independent of iso, and always nelem+1. */
	a = sqrt( 1.8887E8 * (nelem+1.) / pow((double)dense.xNucleiTotal, 0.333) );
	ASSERT( a > 0. );
	if( a > (double)iso.n_HighestResolved_max[ipISO][nelem]+(double)iso.nCollapsed_max[ipISO][nelem] )
	{
		np = iso.n_HighestResolved_max[ipISO][nelem]+iso.nCollapsed_max[ipISO][nelem] + 1;
	}
	else
		np = (long)a;

	/* Debye shielding - the density here is electron density	*/
	/* This one depends on effective charge. */
	a = 2.6E7 * eff_charge * eff_charge * pow( phycon.te/dense.eden, 0.25);
	ASSERT( a > 0. );
	if( a > (double)iso.n_HighestResolved_max[ipISO][nelem]+(double)iso.nCollapsed_max[ipISO][nelem] )
	{
		nd = iso.n_HighestResolved_max[ipISO][nelem]+iso.nCollapsed_max[ipISO][nelem] + 1;
	}
	else
		nd = (long)a;

	/* Stark broadening - this should be the density of singly charged ions, 
	 * both positive and negative.  The sum of protons, electrons, and HeII should be
	 * good enough.	*/
	/* This one depends on effective charge. */
	a = 3171. * pow( (double)eff_charge, 0.8 ) * pow( dense.eden + (double)dense.xIonDense[ipHYDROGEN][1]
		+ (double)dense.xIonDense[ipHELIUM][1], -0.1333);
	ASSERT( a > 0. );
	if( a > (double)iso.n_HighestResolved_max[ipISO][nelem]+(double)iso.nCollapsed_max[ipISO][nelem] )
	{
		ns = iso.n_HighestResolved_max[ipISO][nelem]+iso.nCollapsed_max[ipISO][nelem] + 1;
	}
	else
		ns = (long)a;

	nc = MIN3(np, nd, ns);
	/* Don't allow continuum to be lowered below n=3. */
	nc = MAX2( nc, 3 );

	// on very first call we must define valid data up to top of all possible levels.
	// Rates are evaluated to local highest level when temperature changes.
	// Continuum can be raised during a constant temperature part of the solver calls,
	// as the ionization is determined.  If continuum is raised to include a level that has never
	// been defined, because electron density is lower, then the rate will be NaN
	// this insures that valid data of some sort will always exist.
	//
	//TODO - the number of levels should only be changed when we know that the
	// collisional data for all levels will be reevaluated.  The current setup
	// can cause grief since may get different results depending on how identical
	// conditions are approached.
	if( !conv.nTotalIoniz )
		nc = iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] + 1;

	if( nc < iso.n_HighestResolved_max[ipISO][nelem])
	{
		iso.lgLevelsLowered[ipISO][nelem] = true;
		iso.lgLevelsEverLowered[ipISO][nelem] = true;
		iso.lgMustReeval[ipISO][nelem] = true;
		iso.n_HighestResolved_local[ipISO][nelem] = nc;
		iso.nCollapsed_local[ipISO][nelem] = 0;
		iso.numLevels_local[ipISO][nelem] = iso_get_total_num_levels( ipISO, nc, 0 );
	}
	/* Here is the case where the critical n lies among the one or more collapsed levels */
	/* we just get rid of any that are too high. */
	else if( nc <= iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] )
	{
		iso.lgLevelsLowered[ipISO][nelem] = true;
		iso.lgLevelsEverLowered[ipISO][nelem] = true;
		iso.lgMustReeval[ipISO][nelem] = true;
		iso.n_HighestResolved_local[ipISO][nelem] = iso.n_HighestResolved_max[ipISO][nelem];
		iso.nCollapsed_local[ipISO][nelem] = nc - iso.n_HighestResolved_local[ipISO][nelem];
		iso.numLevels_local[ipISO][nelem] = 
			iso_get_total_num_levels( ipISO, iso.n_HighestResolved_max[ipISO][nelem], iso.nCollapsed_local[ipISO][nelem] );
	}
	/* This is usually where control will flow, because in most conditions the continuum will not be lowered.
	* Nothing changes in this case. */
	else
	{
		iso.numLevels_local[ipISO][nelem] = iso.numLevels_max[ipISO][nelem];
		iso.nCollapsed_local[ipISO][nelem] = iso.nCollapsed_max[ipISO][nelem];
		iso.n_HighestResolved_local[ipISO][nelem] = iso.n_HighestResolved_max[ipISO][nelem];
		
		/* if levels were lowered on last pass but are not now, must reeval */
		if( iso.lgLevelsLowered[ipISO][nelem] )
		{
			iso.lgMustReeval[ipISO][nelem] = true;
		}
		else
		{
			iso.lgMustReeval[ipISO][nelem] = false;
		}

		iso.lgLevelsLowered[ipISO][nelem] = false;
	}

	/* None of these can be greater than that which was originally malloc'd. */
	ASSERT( iso.numLevels_local[ipISO][nelem] <= iso.numLevels_max[ipISO][nelem] );
	ASSERT( iso.nCollapsed_local[ipISO][nelem] <= iso.nCollapsed_max[ipISO][nelem] );
	ASSERT( iso.n_HighestResolved_local[ipISO][nelem] <= iso.n_HighestResolved_max[ipISO][nelem] );

	/* Lyman lines can not be greater than original malloc or critical pqn. */
	iso.nLyman[ipISO] = MIN2( nc, iso.nLyman_malloc[ipISO]);

	if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug)  )
	{
		fprintf( ioQQQ,"     iso_continuum_lower: ipISO %li nelem %li nc %li numLevels %li nCollapsed %li n_HighestResolved %li \n",
			ipISO, 
			nelem,
			nc, 
			iso.numLevels_local[ipISO][nelem],
			iso.nCollapsed_local[ipISO][nelem],
			iso.n_HighestResolved_local[ipISO][nelem]
			);
	}

	return;
}
