/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_pop2 do level population for simple two level atom, no radiative transfer */
#include "cddefines.h"
#include "phycon.h"
#include "dense.h"
#include "atoms.h"

double atom_pop2(double omega, 
  double g1, 
  double g2, 
  double a21, 
  double bltz, 
  double abund)
{
	double boltz, 
	  popexc_v, 
	  q12, 
	  q21, 
	  r;

	DEBUG_ENTRY( "atom_pop2()" );

	/* result is density (cm-3) of excited state times a21
	 * result normalized to n1+n2=abund
	 * cdsqte is eden / sqrte * 8.629e-6
	 * */
	boltz = bltz*phycon.teinv;
	if( abund == 0. || boltz > 15. )
	{
		popexc_v = 0.;
		return( popexc_v );
	}

	/*begin sanity check */
	ASSERT( omega > 0. );

	q21 = dense.cdsqte*omega;
	q12 = q21/g1*exp(-boltz);
	q21 /= g2;
	r = (a21 + q21)/q12;
	popexc_v = abund*a21/(r + 1.);
	return( popexc_v );
}
