/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atmdat_DielSupres derive scale factors for suppression of Burgess dielectronic recombination */
#include "cddefines.h"
#include "ionbal.h"
#include "dense.h"
#include "phycon.h"
#include "save.h"
#include "atmdat.h"

/* >>chng 07 feb 28 by Mitchell Martin, added new dr suppression routine */
/* This function computes the standard electron density-dependent
 * suppression factor of the DR rate coefficient of the C3+ ion 
 * at T = 10^5 K, based on Nigel Badnell's 1993 ApJ paper.
 * It is then scalable for other choices of ionic charge and temperature.
 */
//#define USENEW
#ifdef USENEW
STATIC double dr_suppress(
	/* This routine takes the following arguments:
	 * atomic_number = nuclear charge */
	long int atomic_number, 
	/*ionic_charge = ionic charge*/
	long int ionic_charge, 
	/*eden = electron density */
	double eden, 
	/*T = temperature (K)*/
	double T )
{


	/* fitting constants to compute nominal suppression factor function */
	const double A = 12.4479;	/* various fitting parameters follow */
	const double mu = 0.46665;
	const double w = 4.96916;
	const double y_c0 = 0.5498; 	/* log10 of the electron density fitting parameter for C3+ */

	/* Nigel's 1993 ApJ paper computed the DR rate as a function of log n_e
	 * at a temperature T = 100,000 K.
	 */
	const double T_0 = 1.e5;	/* the standard temperature in Nigel's original C3+ fit */
	const double q_0 = 3.;		/* the ionic charge of C3+, addressed in Nigel's paper */

	/* a fitting constant to compute the suppression factor corrected for an
	 * estimate of surviving DR based on the lowest dipole allowed core
	 * excitation energy
	 */
	const double c = 3.0;	/* smaller c means larger fraction will survive, and vice versa */

	double s, snew, y_c, E_c, E_c0, x, g_x;
	long int iso_sequence;

	eden = log10(eden);
	iso_sequence = atomic_number - ionic_charge;	/* the isoelectronic sequence number, iso_sequence=3 for Li-like, etc */

	y_c = y_c0 + log10( pow( (ionic_charge/q_0), 7. ) * sqrt( T/T_0 ) );

	/* here we compute the standard suppression factor function, s( n_e, T, ionic_charge ) */
	if( eden >= y_c )
	{
		s = (A/(PI*w)) * ( mu/( 1. + pow((eden-y_c)/w, 2.) ) + 
		    (1. - mu) * sqrt(PI*LN_TWO) * exp( -LN_TWO * 
		    pow((eden-y_c)/w, 2.) ) );
	}
	else
	{
		s = (A/(PI*w)) * ( mu + (1.- mu) * sqrt(PI*LN_TWO) );
	}

	/* Now we're going to modify this standard suppression factor curve to
	 * allow for the survival of some fraction of the total DR rate at
	 * generally lower temperatures T, when appropriate.
	 */

	/* Computational estimates of lowest dipole allowed core excitation
	 * energies for various iso-electronic sequences of recombining ion;
	 * these are fits to NIST statistical weighted energies
	 */
	if( iso_sequence == 3 )	/* Li-like ions */
	{
		E_c = 2.08338 + 19.1356*(ionic_charge/10.) + 0.974 *
		      pow( ionic_charge/10., 2. ) - 0.309032*pow( ionic_charge/10., 3. ) +
		      0.419951*pow( ionic_charge/10., 4. );
	}
	else if( iso_sequence == 4 ) /* Be-like ions */
	{
		E_c = 5.56828 + 34.6774*(ionic_charge/10.) + 1.005 *
		      pow( ionic_charge/10., 2. ) - 0.994177*pow( ionic_charge/10., 3. ) +
		      0.726053*pow( ionic_charge/10., 4. );
	}
	else if( iso_sequence == 7 ) /* N-like ions */
	{
		E_c = 10.88361 + 39.7851*(ionic_charge/10.) + 0.423 *
		      pow( ionic_charge/10., 2. ) - 0.310368*pow( ionic_charge/10., 3. ) +
		      0.937186*pow( ionic_charge/10., 4. );
	}
	else if( iso_sequence == 11 ) /* Na-like ions */
	{
		E_c = 2.17262 + 22.5038*(ionic_charge/10.) - 1.227*pow( ionic_charge/10., 2. ) +
		      0.801291*pow( ionic_charge/10., 3. ) +
		      0.0434168*pow( ionic_charge/10., 4. );
	}
	else if( iso_sequence == 1 || iso_sequence == 2 || iso_sequence == 10 )
	{
		/* set to a very large number to force suppression factor to 1.0
		 * for H, He, Ne-like ions */
		E_c = 1.e10;
	}
	else
	{
		/* specifically B, C, O, or F-like ions (iso_sequence = 5, 6, 8, 9) */
		E_c = 0.0;	/* forces suppression factor to s for all T */
		/* iso_sequence.B. ion sequences beyond Na-like, iso_sequence > 11, are currently not
 		 * treated */
	}

	/* the lowest dipole allowed energy for Li-like C3+, atomic_number = 6, iso_sequence = atomic_number-ionic_charge = 3 */
	E_c0 = 2.08338 + 19.1356*(q_0/10.) + 0.974 * 
	       pow( (q_0/10.), 2. ) - 0.309032 * 
	       pow( (q_0/10.), 3. ) + 0.419951 *
	       pow( (q_0/10.), 4. );

	/* and important factor that determines what survives */
	x = ( (E_c*EVDEGK)/(c*T) - (E_c0*EVDEGK)/(c*T_0) );

	if( x > 1 )
	{
		g_x = x;
	}
	else if( x >= 0 && x <= 1 )
	{
		g_x = x*x;
	}
	else
	{
		g_x = 0.0;
	}

	/* converting the standard curve to the revised one allowing for
	 * survival at lower energies
	 */
	snew = 1. + (s-1.)*exp(-g_x);

	return snew;
	ASSERT( snew >=0. && snew <= 1. );
}
#endif

void atmdat_DielSupres(void)
{
	long int i;

	DEBUG_ENTRY( "atmdat_DielSupres()" );

	/* dielectronic burgess recombination suppression, default is true */
	if( ionbal.lgSupDie[0] )
	{
		for( i=0; i < LIMELM; i++ )
		{
#			ifdef USENEW
			ionbal.DielSupprs[0][i] = (realnum)dr_suppress( i+1, 3, dense.eden, phycon.te );
#			else
			/* effective density for scaling from Davidson's plot
			 * first do temperature scaling, term in () is SQRT(te/15,000) */
			double effden = dense.eden/(phycon.sqrte/122.47);

			/* this is rough charge dependence, z^7 from Davidson */
			effden /= powi((realnum)(i+1)/3.,7);

			ionbal.DielSupprs[0][i] = (realnum)(1.-0.092*log10(effden));
			ionbal.DielSupprs[0][i] = (realnum)MIN2(1.,ionbal.DielSupprs[0][i]);
			ionbal.DielSupprs[0][i] = (realnum)MAX2(0.08,ionbal.DielSupprs[0][i]);
#			endif
		}
	}

	else
	{
		for( i=0; i < LIMELM; i++ )
		{
			ionbal.DielSupprs[0][i] = 1.;
		}
	}

	/* nussbaumer and storey recombination
	 * default is this to be false */
	if( ionbal.lgSupDie[1] )
	{
		for( i=0; i < LIMELM; i++ )
		{
			/* assume same factors as above */
			ionbal.DielSupprs[1][i] = ionbal.DielSupprs[0][i];
		}
	}
	else
	{
		for( i=0; i < LIMELM; i++ )
		{
			ionbal.DielSupprs[1][i] = 1.;
		}
	}

	/* option to save recombination coefficients, set with *save recombination 
	 * coefficients* command*/
	if( save.lgioRecom )
	{
		fprintf( save.ioRecom, " atmdat_DielSupres finds following dielectronic"
			" recom suppression factors.\n" );
		fprintf( save.ioRecom, "  Z    fac \n" );
		for( i=0; i < LIMELM; i++ )
		{
			fprintf( save.ioRecom, "%3ld %10.3e %10.3e\n", i+1, 
			  ionbal.DielSupprs[0][i], ionbal.DielSupprs[1][i] );
		}
		fprintf( save.ioRecom, "\n");
	}
	return;
}
