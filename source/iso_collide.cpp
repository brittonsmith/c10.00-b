/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atmdat.h"
#include "conv.h"
#include "dense.h"
#include "heavy.h"
#include "helike_cs.h"
#include "hydrogenic.h"
#include "hydro_vs_rates.h"
#include "ionbal.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "physconst.h"
#include "rfield.h"
#include "secondaries.h"
#include "trace.h"
#include "taulines.h"

/* These are masses relative to the proton mass of the electron, proton, he+, and alpha particle. */
static double ColliderMass[4] = {ELECTRON_MASS/PROTON_MASS, 1.0, 4.0, 4.0};

void iso_collisional_ionization( long ipISO, long nelem )
{
	ASSERT( ipISO < NISO );

	DEBUG_ENTRY( "iso_collisional_ionization()" );

	/* collisional ionization from ground */
	iso.ColIoniz[ipISO][nelem][0] = iso.lgColl_ionize[ipISO] *
		t_ADfA::Inst().coll_ion( nelem+1, 1+ipISO, phycon.te );

	iso_put_error(ipISO,nelem,iso.numLevels_max[ipISO][nelem],0,IPCOLLIS,0.20f,0.20f);
	
	for( long ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
	{
		if( nelem == ipISO )
		{
			/* use routine from Vriens and Smeets (1981). */
			/* >>refer	iso neutral	col.ion.	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940 */
			iso.ColIoniz[ipISO][nelem][ipHi] = hydro_vs_ioniz( iso.xIsoLevNIonRyd[ipISO][nelem][ipHi], phycon.te );
		}
		else
		{
			/* ions */
			/* use hydrogenic ionization rates for ions
			 * >>refer	iso ions	col.ion.	Allen 1973, Astro. Quan. for low Te.
			 * >>refer	iso ions	col.ion.	Sampson and Zhang 1988, ApJ, 335, 516 for High Te.
			 * */
			iso.ColIoniz[ipISO][nelem][ipHi] = 
				Hion_coll_ioniz_ratecoef( ipISO, nelem, N_(ipHi), iso.xIsoLevNIonRyd[ipISO][nelem][ipHi], phycon.te  );
		}
		
		/* iso.lgColl_ionize is option to turn off collisions, "atom XX-like collis off" comnd */
		/* always leave highest level coupled to continuum */
		if( ipHi < iso.numLevels_max[ipISO][nelem] - 1 )
			iso.ColIoniz[ipISO][nelem][ipHi] *= iso.lgColl_ionize[ipISO];
 	
		iso_put_error(ipISO,nelem,iso.numLevels_max[ipISO][nelem],ipHi,IPCOLLIS,0.20f,0.20f);
	}
	
	/* Here we arbitrarily scale the highest level ionization to account for the fact
	 * that, if the atom is not full size, this level should be interacting with higher
	 * levels and not just the continuum.  We did add on collisional excitation terms instead
	 * but this caused a problem at low temperatures because the collisional ionization was 
	 * a sum of terms with different Boltzmann factors, while PopLTE had just one Boltzmann
	 * factor.  The result was a collisional recombination that had residual exponentials of
	 * the form exp(x/kT), which blew up at small T.	*/
	if( !iso.lgLevelsLowered[ipISO][nelem] )
	{
		iso.ColIoniz[ipISO][nelem][iso.numLevels_max[ipISO][nelem]-1] *= 100.;
	}
	
	return;
}

void iso_suprathermal( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_suprathermal()" );

	/* check that we were called with valid parameters */
	ASSERT( ipISO < NISO );
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	/***********************************************************************
	 *                                                                     *
	 * get secondary excitation by suprathermal electrons                  *
	 *                                                                     *
	 ***********************************************************************/

	for( long i=1; i < iso.numLevels_max[ipISO][nelem]; i++ )
	{
		if( Transitions[ipISO][nelem][i][0].ipCont > 0 )
		{
			/* get secondaries for all iso lines by scaling LyA 
			 * excitation by ratio of cross section (oscillator strength/energy) 
			 * Born approximation or plane-wave approximation based on
			 *>>refer	HI	excitation	Shemansky, D.E., et al., 1985, ApJ, 296, 774 */
			secondaries.Hx12[ipISO][nelem][i] = secondaries.x12tot *
				(Transitions[ipISO][nelem][i][0].Emis->gf/
				Transitions[ipISO][nelem][i][0].EnergyWN) /
				(Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][0].Emis->gf/
				Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][0].EnergyWN);
		}
		else
			secondaries.Hx12[ipISO][nelem][i] = 0.;
	}

	return;
}

/*============================*/
/* evaluate collisional rates */
void iso_collide( long ipISO, long nelem )
{
	long ipHi, ipLo;
	double factor, ConvLTEPOP;

	/* this array stores the last temperature at which collision data were evaluated for
	 * each species of the isoelectronic sequence. */
	static double TeUsed[NISO][LIMELM]={
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0} };

	/* collision strengths are assumed to be roughly constant for small changes in temperature
	 * and are not recalculated as often as other data.  This array stores the last temperature
	 * at which collision strengths were evaluated for each species of the isoelectronic sequence. */
	static double TeUsedForCS[NISO][LIMELM]={
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0} };

	DEBUG_ENTRY( "iso_collide()" );

	/* check that we were called with valid parameters */
	ASSERT( ipISO < NISO );
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );
	
	/* skip most of this routine if temperature has not changed,
	 * the check on conv.nTotalIoniz is to make sure that we redo this
	 * on the very first call in a grid calc - it is 0 on the first call */
	if( fp_equal( TeUsed[ipISO][nelem], phycon.te ) && conv.nTotalIoniz && !iso.lgMustReeval[ipISO][nelem] )
	{
		ASSERT( Transitions[ipISO][nelem][ iso.nLyaLevel[ipISO] ][0].Coll.ColUL >= 0. );

		if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
		{
			fprintf( ioQQQ, 
				"     iso_collide called %s nelem %li - no reeval Boltz fac, LTE dens\n",
				iso.chISO[ipISO], nelem );
		}
	}
	else
	{
		TeUsed[ipISO][nelem] = phycon.te;

		if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
		{
			fprintf( ioQQQ, 
				"     iso_collide called %s nelem %li - will reeval Boltz fac, LTE dens\n",
				iso.chISO[ipISO], nelem );
		}
		
		/**********************************************************
		 *                                                        *
		 * Boltzmann factors for all levels,  and                 *
		 * collisional ionization and excitation                  *
		 *                                                        *
		 **********************************************************/

		for( ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
		{
			for( ipLo=0; ipLo<ipHi; ipLo++ )
			{
				iso.Boltzmann[ipISO][nelem][ipHi][ipLo] = 
					sexp( Transitions[ipISO][nelem][ipHi][ipLo].EnergyK / phycon.te );
			}
		}

		/* HION_LTE_POP	is planck^2 / (2 pi m_e k ), raised to 3/2 next */
		factor = HION_LTE_POP*dense.AtomicWeight[nelem]/
			(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);

		/* term in () is stat weight of electron * ion */
		ConvLTEPOP = pow(factor,1.5)/(2.*iso.stat_ion[ipISO])/phycon.te32;

		iso.lgPopLTE_OK[ipISO][nelem] = true;

		// this is the maximum value of iso.PopLTE (units cm^3) that corresponds
		// to the minimum positive density values.  A smaller density will be
		// regarded as zero, and the product PopLTE*n_e*n_Z+ will also be zero.
#define MAX_POP_LTE	(MAX_DENSITY/dense.density_low_limit/dense.density_low_limit)

		/* fully define Boltzmann factors to continuum for model levels */
		for( ipLo=0; ipLo<iso.numLevels_max[ipISO][nelem]; ipLo++ )
		{
			/* this Boltzmann factor is exp( +ioniz energy / Te ) */
			StatesElemNEW[nelem][nelem-ipISO][ipLo].ConBoltz =
				dsexp(iso.xIsoLevNIonRyd[ipISO][nelem][ipLo]/phycon.te_ryd);

			/***************************************
			 *                                     *
			 * LTE abundances for all levels       *
			 *                                     *
			 ***************************************/

			if( StatesElemNEW[nelem][nelem-ipISO][ipLo].ConBoltz > SMALLDOUBLE )
			{
				/* LTE population of given level. */
				iso.PopLTE[ipISO][nelem][ipLo] = 
					StatesElemNEW[nelem][nelem-ipISO][ipLo].g / StatesElemNEW[nelem][nelem-ipISO][ipLo].ConBoltz * ConvLTEPOP;
				ASSERT( iso.PopLTE[ipISO][nelem][ipLo] < BIGDOUBLE );
			}
			else
			{
				iso.PopLTE[ipISO][nelem][ipLo] = 0.;
			}

			iso.PopLTE[ipISO][nelem][ipLo] = MIN2( iso.PopLTE[ipISO][nelem][ipLo], MAX_POP_LTE );

			/* now check for any zeros - if present then matrix cannot be used */
			if( iso.PopLTE[ipISO][nelem][ipLo] <= 0. )
			{
				iso.lgPopLTE_OK[ipISO][nelem] = false;
			}
		}

		iso_collisional_ionization( ipISO, nelem );

		/***********************************************************
		 *                                                         *
		 * collision strengths for all lines in iso sequence       *
		 *                                                         *
		 ***********************************************************/

		if( iso.lgColl_excite[ipISO] &&
			( iso.lgRandErrGen[ipISO] ||
			! fp_equal( phycon.te, TeUsedForCS[ipISO][nelem] ) ) )
		{
			/* Update temperature at which collision strengths are evaluated. */
			TeUsedForCS[ipISO][nelem] = phycon.te;

			for( ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					for( long ipCollider = ipELECTRON; ipCollider <= ipALPHA; ipCollider++ )
					{
						double cs_temp = 0.;

						if( N_(ipHi) == N_(ipLo) && !iso.lgColl_l_mixing[ipISO] )
							cs_temp = 0.;
						else if( N_(ipHi)-N_(ipLo) > 2 && ipCollider > ipELECTRON )
							cs_temp = 0.;
						else if( ipISO == ipH_LIKE )
							cs_temp =  HydroCSInterp( nelem , ipHi , ipLo, ipCollider );
						else if( ipISO == ipHE_LIKE )
							cs_temp = HeCSInterp( nelem , ipHi , ipLo, ipCollider );
						else
							TotalInsanity();

						Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipCollider] = (realnum) cs_temp;

        					/* check for sanity */
                                                ASSERT( Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipELECTRON] >= 0. );
					}

					if( N_(ipHi) <= 5 && N_(ipLo) <= 2 )
						iso_put_error( ipISO, nelem, ipHi, ipLo, IPCOLLIS, 0.10f, 0.30f );
					else
						iso_put_error( ipISO, nelem, ipHi, ipLo, IPCOLLIS, 0.20f, 0.30f );

					// store electron collision strength in generic collision strength
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_str =
						Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipELECTRON];

				}
			}
		}

		/***********************************************************************
		 *                                                                     *
		 * collisional deexcitation for all lines                              *
		 *                                                                     *
		 ***********************************************************************/
		for( ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
		{
			for( ipLo=0; ipLo<ipHi; ipLo++ )
			{
				double reduced_mass_proton = dense.AtomicWeight[nelem]*ColliderMass[ipPROTON]/
					(dense.AtomicWeight[nelem]+ColliderMass[ipPROTON])*ATOMIC_MASS_UNIT;

				double reduced_mass_heplus = dense.AtomicWeight[nelem]*ColliderMass[ipHE_PLUS]/
					(dense.AtomicWeight[nelem]+ColliderMass[ipHE_PLUS])*ATOMIC_MASS_UNIT;

				/********************************************************
				 ********************************************************
				 * NB - the collision strengths for proton and helium  *
				 * ion impact are multiplied by the ratio of collider  *
				 * to electron densities to precorrect for getting     *
				 * multiplied by the electron density when being put   *
				 * into rate matrix later.  They are also multiplied   *
				 * by the pow(m_e/reduced mass)^1.5 to correct for the *
				 * fact that COLL_CONST contains the electron mass.    *
				 * Here, reduced mass means the reduced mass of the    *
				 * collider-target system.                             *
				 ********************************************************
				 ********************************************************/
				Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL = (realnum)(
					(
					/* due to electron impact */
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_str	+
					/* due to proton impact */
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipPROTON]*
					(realnum)(dense.xIonDense[ipHYDROGEN][1]/dense.EdenHCorr)*
					pow(ELECTRON_MASS/reduced_mass_proton, 1.5) +
					/* due to he+ impact */
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipHE_PLUS]*
					(realnum)(dense.xIonDense[ipHELIUM][1]/dense.EdenHCorr)*
					pow(ELECTRON_MASS/reduced_mass_heplus, 1.5) +
					/* due to he++ impact */
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipALPHA]*
					(realnum)(dense.xIonDense[ipHELIUM][2]/dense.EdenHCorr)*
					pow(ELECTRON_MASS/reduced_mass_heplus, 1.5) 
					) / phycon.sqrte*COLL_CONST/(double)StatesElemNEW[nelem][nelem-ipISO][ipHi].g );

				if( ipISO == ipH_LIKE )
				{
					if( N_(ipHi) > iso.n_HighestResolved_max[ipISO][nelem] &&
						N_(ipLo) <= iso.n_HighestResolved_max[ipISO][nelem] )
					{
						Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL *= 
							(8.f/3.f)*(log((realnum)N_(ipHi))+2.f);
					}
				}
				else if( ipISO == ipHE_LIKE )
				{
					fixit();
					/* This is intended to be a trick to get the correct collisional excitation from 
					 * collapsed levels to resolved levels.  Is it needed or does the stat weight used 
					 * above handle it automatically? If it is needed, is this correct? */
					if( N_(ipHi) > iso.n_HighestResolved_max[ipISO][nelem] &&
						N_(ipLo) <= iso.n_HighestResolved_max[ipISO][nelem] )
					{
						Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL *=
							(2.f/3.f)*(log((realnum)N_(ipHi))+2.f);
					}
				}
					
				Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL *= iso.lgColl_excite[ipISO];
			}
		}

		if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO]) && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_collide: %s Z=%li de-excitation rates coefficients\n", iso.chISO[ipISO], nelem + 1 );
			for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				fprintf( ioQQQ, " %li\t", ipHi );
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					fprintf( ioQQQ,PrintEfmt("%10.3e", Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL ));
				}
				fprintf( ioQQQ, "\n" );
			}

			fprintf( ioQQQ, "     iso_collide: %s Z=%li collisional ionization coefficients\n", iso.chISO[ipISO], nelem + 1 );
			for( ipHi=0; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				fprintf( ioQQQ,PrintEfmt("%10.3e",  iso.ColIoniz[ipISO][nelem][ipHi] ));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "     iso_collide: %s Z=%li continuum boltzmann factor\n", iso.chISO[ipISO], nelem + 1 );
			for( ipHi=0; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				fprintf( ioQQQ,PrintEfmt("%10.3e",  StatesElemNEW[nelem][nelem-ipISO][ipHi].ConBoltz ));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "     iso_collide: %s Z=%li continuum boltzmann factor\n", iso.chISO[ipISO], nelem + 1 );
			for( ipHi=0; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				fprintf( ioQQQ,PrintEfmt("%10.3e",  iso.PopLTE[ipISO][nelem][ipHi] ));
			}
			fprintf( ioQQQ, "\n" );
		}

		/* the case b hummer and storey option,
		 * this kills collisional excitation and ionization from n=1 and n=2 */
		if( opac.lgCaseB_HummerStorey )
		{
			for( ipLo=0; ipLo<iso.numLevels_max[ipISO][nelem]-1; ipLo++ )
			{
				if( N_(ipLo)>=3 )
					break;
				
				iso.ColIoniz[ipISO][nelem][ipLo] = 0.;

				for( ipHi=ipLo+1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					/* don't disable 2-2 collisions */
					if( N_(ipLo)==2 && N_(ipHi)==2 )
						continue;

					Transitions[ipISO][nelem][ipHi][ipLo].Coll.ColUL = 0.;
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_str = 0.;
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipPROTON] = 0.;
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipHE_PLUS] = 0.;
					Transitions[ipISO][nelem][ipHi][ipLo].Coll.col_stri[ipALPHA] = 0.;
				}
			}
		}
	}

	iso_suprathermal( ipISO, nelem );

	/* this must be reevaluated every time since eden can change when Te does not */
	/* save into main array - collisional ionization by thermal electrons */
	ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0] = 
		iso.ColIoniz[ipISO][nelem][0]*dense.EdenHCorr;

	/* cooling due to collisional ionization, which only includes thermal electrons */
	ionbal.CollIonRate_Ground[nelem][nelem-ipISO][1] = 
		ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0]*
		rfield.anu[Heavy.ipHeavy[nelem][nelem-ipISO]-1]*EN1RYD;

	return;
}
