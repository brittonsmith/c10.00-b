/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_pop3 solve 3-level atom without radiative transfer, returns pops of level 2 and 3 */
#include "cddefines.h"
#include "phycon.h"
#include "dense.h"
#include "atoms.h"

/*atom_pop3 return value is population for 3-level atom, cm^-3 */
double atom_pop3(
	/* statictical weights of levels 1, 2, and 3 */
	double g1, double g2, double g3, 

	/* collision strengths between three levels */
	double o12, double o13, double o23, 

	/* transition probabilities between three levels */
	double a21, double a31, double a32, 

	/* excitation energy in Kelvin */
	double Tex12, double Tex23, 

	/* returned population of level 2, cm^-3 */
	realnum *pop2, 

	/* incoming total abundance of ion */
	double abund, 

	/* possible photodestruction of level 2, normally 0 */
	double gam2,

	/* excitation rates (s-1) due to "other" processes,
	 * these are not included in the energy exchange */
	double r12, 
	double r13 )
{
	double alf, 
	  b12, 
	  b13, 
	  b23, 
	  bet, 
	  c21, 
	  c23, 
	  c31, 
	  c32, 
	  ex, 
	  fac, 
	  pop3_v;

	DEBUG_ENTRY( "atom_pop3()" );

	/* computes level populations for 3 level atom, all col rad coupling
	 * results are populations of levels 2 and 3, (cm^-3) no A included */
	ex = Tex12/phycon.te;
	if( (abund <= 0.) || (ex > 20. && r12<SMALLFLOAT ) )
	{
		pop3_v = 0.;
		*pop2 = 0.;
		return( pop3_v );
	}

	/* confirm that these are sane values */
	ASSERT( g1>0. && g2>0. && g3>0. && o12>=0. && o13>=0. && o23>=0. && a21>=0. && a31>=0. && a32>=0. && 
		Tex12>=0. && Tex23>=0. );

	b12 = exp(-ex);
	b23 = exp(-Tex23/phycon.te);

	b13 = b12*b23;
	if( b13 == 0. && r12<SMALLFLOAT )
	{
		pop3_v = 0.;
		*pop2 = 0.;
		return( pop3_v );
	}

	/* these rates have units s-1 */
	atoms.c12 = dense.cdsqte*o12/g1*b12 + r12;
	atoms.c13 = dense.cdsqte*o13/g1*b13 + r13;
	c23 = dense.cdsqte*o23/g2*b23;
	c32 = dense.cdsqte*o23/g3;
	c31 = dense.cdsqte*o13/g3;
	c21 = dense.cdsqte*o12/g2;

	alf = a21 + c21 + c23 + gam2;
	bet = a31 + a32 + c31 + c32;
	*pop2 = (realnum)((atoms.c13/bet + atoms.c12/(c32 + a32))/(alf/(c32 + a32) - c23/bet));
	pop3_v = (atoms.c13 + *pop2*c23)/bet;

	/* renorm to 1+pop2+atom_pop3=1 */
	fac = abund/(1. + *pop2 + pop3_v);
	*pop2 *= (realnum)fac;
	pop3_v *= fac;

	return( pop3_v );
}
