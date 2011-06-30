/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_radiative_recomb find state specific creation and destruction rates for iso sequences */
/*iso_RRCoef_Te - evaluate radiative recombination coef at some temperature */
/*iso_recomb_check - called by SanityCheck to confirm that recombination coef are ok,
 * return value is relative error between new calculation of recom, and interp value */

#include "cddefines.h"
#include "dynamics.h"
#include "atmdat.h"
#include "conv.h"
#include "cosmology.h"
#include "dense.h"
#include "elementnames.h"
#include "helike.h"
#include "helike_recom.h"
#include "hydrogenic.h"
#include "ionbal.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "physconst.h" 
#include "prt.h"
#include "save.h"
#include "thermal.h"
#include "thirdparty.h"
#include "trace.h"
#include "rt.h"
#include "taulines.h"

/* this will save log of radiative recombination rate coefficients at N_ISO_TE_RECOMB temperatures.
 * there will be NumLevRecomb[ipISO][nelem] levels saved in RRCoef[nelem][level][temp] */
static double ****RRCoef/*[LIMELM][NumLevRecomb[ipISO][nelem]][N_ISO_TE_RECOMB]*/;
static long **NumLevRecomb;
static double ***TotalRecomb;	/*[ipISO][nelem][i]*/

/* the array of logs of temperatures at which RRCoef was defined */
static double TeRRCoef[N_ISO_TE_RECOMB];
static double kTRyd,global_EthRyd; 
static long int ipLev,globalZ,globalISO;
static long int globalN, globalL, globalS;

STATIC double TempInterp( double* TempArray , double* ValueArray, long NumElements );
STATIC double iso_recomb_integrand(double EE);
STATIC void iso_put_recomb_error( long ipISO, long nelem );

double iso_radrecomb_from_cross_section(long ipISO, double temp, long nelem, long ipLo)
{
	double alpha,RecomIntegral=0.,b,E1,E2,step,OldRecomIntegral,TotChangeLastFive;
	double change[5] = {0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "iso_radrecomb_from_cross_section()" );

	if( ipISO==ipH_LIKE && ipLo == 0 )
		return t_ADfA::Inst().H_rad_rec(nelem+1,ipLo,temp);

	global_EthRyd = iso.xIsoLevNIonRyd[ipISO][nelem][ipLo];

	/* Factors outside integral in Milne relation.	*/
	b = MILNE_CONST * StatesElemNEW[nelem][nelem-ipISO][ipLo].g * pow(temp,-1.5);

	if( ipISO==ipH_LIKE )
		b /= 2.;
	else if( ipISO==ipHE_LIKE ) 	
		b /= 4.;
		
	/* kT in Rydbergs.	*/
	kTRyd = temp / TE1RYD;
	globalISO = ipISO;
	globalZ = nelem;
	ipLev = ipLo;
	globalN = N_(ipLo);
	globalL = L_(ipLo);
	globalS = S_(ipLo);

	/* Begin integration.	*/
	/* First define characteristic step */
	E1 = global_EthRyd;

	if( ipISO==ipH_LIKE )
		step = MIN2( 0.125*kTRyd, 0.5*E1 );
	else if( ipISO==ipHE_LIKE ) 	
		step = MIN2( 0.25*kTRyd, 0.5*E1 );
	else
		TotalInsanity();
	
	E2 = E1 + step;
	/* Perform initial integration, from threshold to threshold + step.	*/
	RecomIntegral = qg32( E1, E2, iso_recomb_integrand);
	/* Repeat the integration, adding each new result to the total, 
	 * except that the step size is doubled every time, since values away from 
	 * threshold tend to fall off more slowly.	*/
	do
	{
		OldRecomIntegral = RecomIntegral;
		E1 = E2;
		step *= 1.25;
		E2 = E1 + step;
		RecomIntegral += qg32( E1, E2, iso_recomb_integrand);
		change[4] = change[3];
		change[3] = change[2];
		change[2] = change[1];
		change[1] = change[0];
		change[0] = (RecomIntegral - OldRecomIntegral)/RecomIntegral;
		TotChangeLastFive = change[0] + change[1] + change[2] + change[3] + change[4];
	/* Continue integration until the upper limit exceeds 100kTRyd, an arbitrary
	 * point at which the integrand component exp(electron energy/kT) is very small,
	 * making neglible cross-sections at photon energies beyond that point,
	 * OR when the last five steps resulted in less than a 1 percent change.	*/
	} while ( ((E2-global_EthRyd) < 100.*kTRyd) && ( TotChangeLastFive > 0.0001) );

	/* Calculate recombination coefficient.	*/
	alpha = b * RecomIntegral;

	alpha = MAX2( alpha, SMALLDOUBLE );

	return alpha;
}

/*iso_recomb_integrand, used in Milne relation for iso sequences - the energy is photon Rydbergs.	*/
STATIC double iso_recomb_integrand(double ERyd)
{
	double x1, temp;

	/* Milne relation integrand	*/
	x1 = ERyd * ERyd * exp(-1.0 * ( ERyd - global_EthRyd ) / kTRyd);
	temp = iso_cross_section( ERyd , global_EthRyd, globalN, globalL, globalS, globalZ, globalISO );
	x1 *= temp;

	return x1;
}

double iso_cross_section( double EgammaRyd , double EthRyd, long n, long l, long S, long globalZ, long globalISO )
{
	double cross_section;
	DEBUG_ENTRY( "iso_cross_section()" );

	if( globalISO == ipH_LIKE )
		cross_section = H_cross_section( EgammaRyd , EthRyd, n, l, globalZ );
	else if( globalISO == ipHE_LIKE )
		cross_section =  He_cross_section( EgammaRyd , EthRyd, n, l, S, globalZ );	
	else
		TotalInsanity();

	return cross_section;
}

/*=======================================================*/
/* iso_radiative_recomb get rad recomb rate coefficients for iso sequences */
void iso_radiative_recomb( 
						  long ipISO,
						  /* nelem on the c scale, He is 1 */
						  long nelem )
{
	long ipFirstCollapsed, LastN=0L, ThisN=1L, ipLevel; 
	double topoff, TotMinusExplicitResolved,
		TotRRThisN=0., TotRRLastN=0., Total_DR_Added=0.;
	double RecExplictLevels, TotalRadRecomb, RecCollapsed;
	static double TeUsed[NISO][LIMELM]={
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};

	DEBUG_ENTRY( "iso_radiative_recomb()" );

	iso_put_recomb_error( ipISO, nelem );

	/* evaluate recombination escape probability for all levels */

	/* define radiative recombination rates for all levels */ 
	/* this will be the sum of all levels explicitly included in the model atom */
	RecExplictLevels = 0.;

	/* number of resolved levels, so first collapsed level is [ipFirstCollapsed] */
	ipFirstCollapsed = iso.numLevels_local[ipISO][nelem]-iso.nCollapsed_local[ipISO][nelem]; 

	if( !fp_equal(phycon.te,TeUsed[ipISO][nelem]) || iso.lgMustReeval[ipISO][nelem] || !conv.nTotalIoniz || cosmology.lgDo )
	{
		TeUsed[ipISO][nelem] = phycon.te;

		for( ipLevel=0; ipLevel<ipFirstCollapsed; ++ipLevel )
		{
			/* this is radiative recombination rate coefficient */
			double RadRec;

			if( !iso.lgNoRecombInterp[ipISO] )
			{
				RadRec = iso_RRCoef_Te( ipISO, nelem , ipLevel );
			}
			else
			{
				RadRec = iso_radrecomb_from_cross_section( ipISO, phycon.te, nelem, ipLevel);
			}
			ASSERT( RadRec > 0. );

			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] = RadRec;

			ASSERT( iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] > 0. );

			RecExplictLevels += iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad];

			if( iso.lgDielRecom[ipISO] )
			{
				/* \todo 2 suppress these rates for continuum lowering using factors from Jordan (1969). */
				iso.DielecRecomb[ipISO][nelem][ipLevel] = iso_dielec_recomb_rate( ipISO, nelem, ipLevel );
				Total_DR_Added	+= iso.DielecRecomb[ipISO][nelem][ipLevel];
			}
		}

		/**************************************************/
		/***  Add on recombination to collapsed levels  ***/
		/**************************************************/
		RecCollapsed = 0.;
		for( ipLevel=ipFirstCollapsed; ipLevel<iso.numLevels_local[ipISO][nelem]; ++ipLevel )
		{
			/* use hydrogenic for collapsed levels */
			double RadRec = t_ADfA::Inst().H_rad_rec(nelem+1-ipISO, N_(ipLevel), phycon.te);

			/* this is radiative recombination rate coefficient */
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] = RadRec;

			/* This kills recombination into the collapsed level so that the forced
			 * statistical weighting can be bypassed */
			if( !iso.lgTopoff[ipISO] )
			{
				iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] *= 1E-10;
			}

			RecCollapsed += iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad];

			ASSERT( iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] > 0. );

			if( iso.lgDielRecom[ipISO] )
			{
				/* \todo 2 suppress these rates for continuum lowering using factors from Jordan (1969). */
				iso.DielecRecomb[ipISO][nelem][ipLevel] = iso_dielec_recomb_rate( ipISO, nelem, ipLevel );
				Total_DR_Added	+= iso.DielecRecomb[ipISO][nelem][ipLevel];
			}
		}

		/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
		for( ipLevel = 0; ipLevel<iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{
			if( N_(ipLevel) == ThisN )
			{
				TotRRThisN += iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad];
			}
			else
			{
				ASSERT( iso.lgRandErrGen[ipISO] || (TotRRThisN<TotRRLastN) || (ThisN<=2L) || (phycon.te>3E7) || (nelem!=ipHELIUM) || (ThisN==iso.n_HighestResolved_local[ipISO][nelem]+1) );
				LastN = ThisN;
				ThisN = N_(ipLevel);
				TotRRLastN = TotRRThisN;
				TotRRThisN = iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad];

				{
					/* Print the sum of recombination coefficients per n at current temp.	*/
					/*@-redef@*/
					enum {DEBUG_LOC=false};
					/*@+redef@*/
					static long RUNONCE = false;

					if( !RUNONCE && DEBUG_LOC )
					{
						static long FIRSTTIME = true;

						if( FIRSTTIME )
						{
							fprintf( ioQQQ,"Sum of Radiative Recombination at current iso, nelem, temp = %li %li %.2f\n", 
								ipISO, nelem, phycon.te);
							FIRSTTIME= false;
						}

						fprintf( ioQQQ,"%li\t%.2e\n",LastN,TotRRLastN );
					}
					RUNONCE = true;
				}
			}
		}

		/* Get total recombination into all levels, including those not explicitly considered. */
		if( iso.lgNoRecombInterp[ipISO] )
		{
			/* We are not interpolating, must calculate total now, as sum of resolved and collapsed levels... */
			TotalRadRecomb = RecCollapsed + RecExplictLevels;

			/* Plus additional levels out to a predefined limit... */
			for( long nLo = N_(ipFirstCollapsed) + iso.nCollapsed_max[ipISO][nelem]; nLo < NHYDRO_MAX_LEVEL; nLo++ )
			{
				TotalRadRecomb += t_ADfA::Inst().H_rad_rec(nelem+1-ipISO, nLo, phycon.te);
			}
			/* Plus a bunch more for good measure */
			for( long nLo = NHYDRO_MAX_LEVEL; nLo<=SumUpToThisN; nLo++ )
			{
				TotalRadRecomb += Recomb_Seaton59( nelem+1-ipISO, phycon.te, nLo );
			}
		}
		else
		{
			/* We are interpolating, and total was calculated here in iso_recomb_setup */
			TotalRadRecomb = iso_RRCoef_Te( ipISO, nelem, iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem] );
		}
		if(ipISO==0 && nelem==0 )
		{
			// insure rec coef will be always evaluated
			TeUsed[ipISO][nelem] = phycon.te*0.999;

			 if(  iteration > dynamics.n_initial_relax+2 && fudge(-1)>0 )
				TotalRadRecomb *= fudge(0);
		}

		/* If generating random error, apply one to total recombination */
		if( iso.lgRandErrGen[ipISO] )
		{
			/* Error for total recombination */
			iso_put_error(ipISO,nelem,iso.numLevels_max[ipISO][nelem],iso.numLevels_max[ipISO][nelem],IPRAD,0.0001f,0.0001f);
			iso.ErrorFactor[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][iso.numLevels_max[ipISO][nelem]][IPRAD] =
				(realnum)MyGaussRand( iso.Error[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][iso.numLevels_max[ipISO][nelem]][IPRAD] );

			/* this has to be from iso.numLevels_max instead of iso.numLevels_local because
			 * the error factors for rrc are always stored at iso.numLevels_max, regardless of
			 * continuum lowering effects. */
			TotalRadRecomb *= 
				iso.ErrorFactor[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][iso.numLevels_max[ipISO][nelem]][IPRAD];
		}

		/* this is case B recombination, sum without the ground included */
		iso.RadRec_caseB[ipISO][nelem] = TotalRadRecomb - iso.RadRecomb[ipISO][nelem][0][ipRecRad];
		ASSERT( iso.RadRec_caseB[ipISO][nelem] > 0.);

		/**********************************************************************/
		/***  Add topoff (excess) recombination	to top level.  This is only	***/
		/***  done if atom is not full size.                                ***/
		/**********************************************************************/
		if( !iso.lgLevelsLowered[ipISO][nelem] )
		{
			/* at this point we have RecExplictLevels, the sum of radiative recombination 
			 * to all levels explicitly included in the model atom the total 
			 * recombination rate.  The difference is the amount of "topoff" we will need to do */
			TotMinusExplicitResolved = TotalRadRecomb - RecExplictLevels;

			topoff = TotMinusExplicitResolved - RecCollapsed;

			/* the t_ADfA::Inst().rad_rec fits are too small at high temperatures, so this atom is
			 * better than the topoff.  Only a problem for helium itself, at high temperatures.
			 * complain if the negative topoff is not for this case */
			if( topoff < 0. && (nelem!=ipHELIUM || phycon.te < 1e5 ) &&
				fabs( topoff/TotalRadRecomb ) > 0.01 )
			{
				fprintf(ioQQQ," PROBLEM  negative RR topoff for %li, rel error was %.2e, temp was %.2f\n",  
					nelem, topoff/TotalRadRecomb, phycon.te );
			}

			if( !iso.lgTopoff[ipISO] )
				topoff *= 1E-20;

			/* We always have at least one collapsed level if continuum is not lowered.  Put topoff there.	*/
			iso.RadRecomb[ipISO][nelem][iso.numLevels_max[ipISO][nelem]-1][ipRecRad] += MAX2(0., topoff );

			/* check for negative DR topoff, but only if Total_DR_Added is not negligible compared with TotalRadRecomb */
			if( Total_DR_Added > TotalRadRecomb/100. )
			{
				if( Total_DR_Added / ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] > 1.02 )
				{
					fprintf(ioQQQ," PROBLEM  negative DR topoff for %li, rel error was %.2e, temp was %.2f\n",  
						nelem,
						Total_DR_Added / ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] - 1.0,
						phycon.te );
				}
			}

			if( iso.lgDielRecom[ipISO] )
			{
				/* \todo 2 suppress this total rate for continuum lowering using factors from Jordan (1969). */
				/* put extra DR in top level */
				iso.DielecRecomb[ipISO][nelem][iso.numLevels_max[ipISO][nelem]-1] = 
					ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] - Total_DR_Added;
			}
		}

	}

	/**************************************************************/
	/***  Stuff escape probabilities, and effective rad recomb  ***/
	/**************************************************************/

	/* total effective radiative recombination, initialize to zero */
	iso.RadRec_effec[ipISO][nelem] = 0.;

	for( ipLevel=0; ipLevel<iso.numLevels_local[ipISO][nelem]; ++ipLevel )
	{
		/* option for case b conditions, kill ground state recombination */
		if( opac.lgCaseB && ipLevel==0 )
		{
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] = 1e-10;
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc] = 1e-10;
		}
		else if( cosmology.lgDo && ipLevel==0 )
		{
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] = 1e-30;
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc] = 1e-30;
		}
		else
		{
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] = 
				RT_recom_effic(iso.ipIsoLevNIonCon[ipISO][nelem][ipLevel]);

			/* net escape prob includes dest by background opacity */
			iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc] = 
				MIN2( 1.,
				iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] +
				(1.-iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc])*
				iso.ConOpacRatio[ipISO][nelem][ipLevel] );
		}

		ASSERT( iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] >= 0. );
		ASSERT( iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc] >= 0. );
		ASSERT( iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc] <= 1. );

		/* sum of all effective rad rec */
		iso.RadRec_effec[ipISO][nelem] += iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad]*
		  iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc];
	}

	/* zero out escape probabilities of levels above numLevels_local */
	for( ipLevel=iso.numLevels_local[ipISO][nelem]; ipLevel<iso.numLevels_max[ipISO][nelem]; ++ipLevel )
	{
		/* this is escape probability */
		iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] = 0.; 
		/* net escape prob includes dest by background opacity */
		iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc] = 0.;
	}

	/* trace escape probabilities */
	if( trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO]) )
	{
		fprintf( ioQQQ, "       iso_radiative_recomb trace ipISO=%3ld Z=%3ld\n", ipISO, nelem );

		/* print continuum escape probability */
		fprintf( ioQQQ, "       iso_radiative_recomb recomb effic" );
		for( ipLevel=0; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso.RadRecomb[ipISO][nelem][ipLevel][ipRecEsc] ));
		}
		fprintf( ioQQQ, "\n" );

		/* net recombination efficiency factor, including background opacity*/
		fprintf( ioQQQ, "       iso_radiative_recomb recomb net effic" );
		for( ipLevel=0; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso.RadRecomb[ipISO][nelem][ipLevel][ipRecNetEsc]) );
		}
		fprintf( ioQQQ, "\n" );

		/* inward continuous optical depths */
		fprintf( ioQQQ, "       iso_radiative_recomb in optic dep" );
		for( ipLevel=0; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{	
			fprintf( ioQQQ,PrintEfmt("%10.3e",  opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipISO][nelem][ipLevel]-1] ));
		}
		fprintf( ioQQQ, "\n" );

		/* outward continuous optical depths*/
		fprintf( ioQQQ, "       iso_radiative_recomb out op depth" );
		for( ipLevel=0; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e",  opac.TauAbsGeo[1][iso.ipIsoLevNIonCon[ipISO][nelem][ipLevel]-1] ));
		}
		fprintf( ioQQQ, "\n" );

		/* print radiative recombination coefficients */
		fprintf( ioQQQ, "       iso_radiative_recomb rad rec coef " );
		for( ipLevel=0; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad]) );
		}
		fprintf( ioQQQ, "\n" );
	}

	if( trace.lgTrace && (trace.lgHBug||trace.lgHeBug) )
	{
		fprintf( ioQQQ, "     iso_radiative_recomb total rec coef" );
		fprintf( ioQQQ,PrintEfmt("%10.3e", iso.RadRec_effec[ipISO][nelem] ));
		fprintf( ioQQQ, " case A=" );
		fprintf( ioQQQ,PrintEfmt("%10.3e", 
			iso.RadRec_caseB[ipISO][nelem] + iso.RadRecomb[ipISO][nelem][ipH1s][ipRecRad] ) );
		fprintf( ioQQQ, " case B=");
		fprintf( ioQQQ,PrintEfmt("%10.3e", iso.RadRec_caseB[ipISO][nelem] ));
		fprintf( ioQQQ, "\n" );
	}

	/****************************/
	/***  begin sanity check  ***/
	/****************************/
	{
		bool lgOK = true;
		for( ipLevel=0; ipLevel < iso.numLevels_local[ipISO][nelem]; ipLevel++ )
		{
			if( iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] <= 0. )
			{
				fprintf( ioQQQ, 
					" PROBLEM iso_radiative_recomb non-positive recombination coefficient for ipISO=%3ld Z=%3ld lev n=%3ld rec=%11.2e te=%11.2e\n", 
				  ipISO, nelem, ipLevel, iso.RadRecomb[ipISO][nelem][ipLevel][ipRecRad] , phycon.te);
					lgOK = false;
			}
		}
		/* bail if we found problems */
		if( !lgOK )
		{
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
		/*end sanity check */
	}

	/* confirm that we have good rec coef at bottom and top of atom/ion */
	ASSERT( iso.RadRecomb[ipISO][nelem][0][ipRecRad] > 0. );
	ASSERT( iso.RadRecomb[ipISO][nelem][iso.numLevels_local[ipISO][nelem]-1][ipRecRad] > 0. );

	/* set true when save recombination coefficients command entered */
	if( save.lgioRecom )
	{
		/* this prints Z on physical, not C, scale */
		fprintf( save.ioRecom, "%s %s %2li ", 
			iso.chISO[ipISO], elementnames.chElementSym[nelem], nelem+1 );
		fprintf( save.ioRecom,PrintEfmt("%9.2e", iso.RadRec_effec[ipISO][nelem] ));
		fprintf( save.ioRecom, "\n" );
	}

	return;
}

STATIC void iso_put_recomb_error( long ipISO, long nelem )
{
	long level;
	
	/* optimistic and pessimistic errors for HeI recombination */
	realnum He_errorOpti[31] = {
		0.0000f,
		0.0009f,0.0037f,0.0003f,0.0003f,0.0003f,0.0018f,
		0.0009f,0.0050f,0.0007f,0.0003f,0.0001f,0.0007f,
		0.0045f,0.0071f,0.0005f,0.0005f,0.0004f,0.0005f,0.0004f,0.0009f,
		0.0045f,0.0071f,0.0005f,0.0005f,0.0004f,0.0005f,0.0004f,0.0005f,0.0004f,0.0009f };

	realnum He_errorPess[31] = {
		0.0100f,
		0.0100f,0.0060f,0.0080f,0.0080f,0.0080f,0.0200f,
		0.0200f,0.0200f,0.0200f,0.0600f,0.0600f,0.0080f,
		0.0200f,0.0200f,0.0070f,0.0100f,0.0100f,0.0020f,0.0030f,0.0070f,
		0.0300f,0.0300f,0.0100f,0.0200f,0.0200f,0.0200f,0.0200f,0.0010f,0.0004f,0.0090f };

	/* now put recombination errors into iso.Error[ipISO] array */
	for( long ipHi=0; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
	{
		if( ipISO==ipHE_LIKE && nelem==ipHELIUM )
		{
			if( ipHi >= iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem] )
				level = 21;
			else if( ipHi<=30 )
				level = ipHi;
			else
				level = iso.QuantumNumbers2Index[ipISO][nelem][5][ MIN2(L_(ipHi),4) ][S_(ipHi)];

			ASSERT( level <=30 );

			iso_put_error(ipISO,nelem,iso.numLevels_max[ipISO][nelem],ipHi,IPRAD, He_errorOpti[level], He_errorPess[level]);
		}
		else
			iso_put_error(ipISO,nelem,iso.numLevels_max[ipISO][nelem],ipHi,IPRAD,0.1f,0.1f);
	}

	return;
}

void iso_radiative_recomb_effective( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_radiative_recomb_effective()" );

	/* Find effective recombination coefficients */
	for( long ipHi=0; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
	{
		iso.RadEffec[ipISO][nelem][ipHi] = 0.;

		/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
		for( long ipHigher=ipHi; ipHigher < iso.numLevels_local[ipISO][nelem]; ipHigher++ )
		{
			ASSERT( iso.CascadeProb[ipISO][nelem][ipHigher][ipHi] >= 0. );
			ASSERT( iso.RadRecomb[ipISO][nelem][ipHigher][ipRecRad] >= 0. );

			iso.RadEffec[ipISO][nelem][ipHi] += iso.CascadeProb[ipISO][nelem][ipHigher][ipHi] *
				iso.RadRecomb[ipISO][nelem][ipHigher][ipRecRad];
		}
	}

	/**************************************************************/
	/***  option to print effective recombination coefficients  ***/
	/**************************************************************/
	{
		enum {DEBUG_LOC=false};

		if( DEBUG_LOC )
		{
			const int maxPrt=10;

			fprintf( ioQQQ,"Effective recombination, ipISO=%li, nelem=%li, Te = %e\n", ipISO, nelem, phycon.te );
			fprintf( ioQQQ, "N\tL\tS\tRadEffec\tLifetime\n" );

			for( long i=0; i<maxPrt; i++ )
			{
				fprintf( ioQQQ, "%li\t%li\t%li\t%e\t%e\n", N_(i), L_(i), S_(i),
					iso.RadEffec[ipISO][nelem][i],
					MAX2( 0., StatesElemNEW[nelem][nelem-ipISO][i].lifetime ) );
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* If we have the variable set, find errors in rad rates */
	if( iso.lgRandErrGen[ipISO] )
	{
		dprintf( ioQQQ, "ipHi\tipLo\tWL\tEmiss\tSigmaEmiss\tRadEffec\tSigRadEff\tBrRat\tSigBrRat\n" );

		/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
		for( long ipHi=0; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
		{
			iso.SigmaRadEffec[ipISO][nelem][ipHi] = 0.;

			/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
			for( long ipHigher=ipHi; ipHigher < iso.numLevels_local[ipISO][nelem]; ipHigher++ )
			{
				ASSERT( iso.Error[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][ipHigher][IPRAD] >= 0. );
				ASSERT( iso.SigmaCascadeProb[ipISO][nelem][ipHigher][ipHi] >= 0. );

				/* iso.RadRecomb has to appear here because iso.Error is only relative error */ 
				iso.SigmaRadEffec[ipISO][nelem][ipHi] += pow( iso.Error[ipISO][nelem][iso.numLevels_max[ipISO][nelem]][ipHigher][IPRAD] *
					iso.CascadeProb[ipISO][nelem][ipHigher][ipHi] * iso.RadRecomb[ipISO][nelem][ipHigher][ipRecRad], 2.) + 
					pow( iso.SigmaCascadeProb[ipISO][nelem][ipHigher][ipHi] * iso.RadRecomb[ipISO][nelem][ipHigher][ipRecRad], 2.);
			}

			ASSERT( iso.SigmaRadEffec[ipISO][nelem][ipHi] >= 0. );
			iso.SigmaRadEffec[ipISO][nelem][ipHi] = sqrt( iso.SigmaRadEffec[ipISO][nelem][ipHi] );

			for( long ipLo = 0; ipLo < ipHi; ipLo++ )
			{
				if( (( L_(ipLo) == L_(ipHi) + 1 ) || ( L_(ipHi) == L_(ipLo) + 1 )) )
				{	
					double EnergyInRydbergs = iso.xIsoLevNIonRyd[ipISO][nelem][ipLo] - iso.xIsoLevNIonRyd[ipISO][nelem][ipHi];
					realnum wavelength = (realnum)(RYDLAM/MAX2(1E-8,EnergyInRydbergs));
					double emissivity = iso.RadEffec[ipISO][nelem][ipHi] * iso.BranchRatio[ipISO][nelem][ipHi][ipLo] * EN1RYD * EnergyInRydbergs;
					double sigma_emiss = 0., SigmaBranchRatio = 0.;

					if( ( emissivity > 2.E-29 ) && ( wavelength < 1.E6 ) && (N_(ipHi)<=5) )
					{
						SigmaBranchRatio = iso.BranchRatio[ipISO][nelem][ipHi][ipLo] * sqrt( 
							pow( (double)iso.Error[ipISO][nelem][ipHi][ipLo][IPRAD], 2. ) +
							pow( iso.SigmaAtot[ipISO][nelem][ipHi]*StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime, 2.) );

						sigma_emiss =  EN1RYD * EnergyInRydbergs * sqrt( 
							pow( (double)iso.SigmaRadEffec[ipISO][nelem][ipHi] * iso.BranchRatio[ipISO][nelem][ipHi][ipLo], 2.) +
							pow( SigmaBranchRatio * iso.RadEffec[ipISO][nelem][ipHi], 2.) );

						/* \todo 2 make this a trace option */
						dprintf( ioQQQ, "%li\t%li\t", ipHi, ipLo );
						prt_wl( ioQQQ, wavelength );
						fprintf( ioQQQ, "\t%e\t%e\t%e\t%e\t%e\t%e\n", 
							emissivity,
							sigma_emiss,
							iso.RadEffec[ipISO][nelem][ipHi],
							iso.SigmaRadEffec[ipISO][nelem][ipHi],
							iso.BranchRatio[ipISO][nelem][ipHi][ipLo],
							SigmaBranchRatio);
					}
				}
			}
		}
	}

	return;
}
/*iso_RRCoef_Te evaluated radiative recombination coef at some temperature */
double iso_RRCoef_Te( long ipISO, long nelem , long n )
{
	double rate = 0.;

	DEBUG_ENTRY( "iso_RRCoef_Te()" );

	ASSERT( !iso.lgNoRecombInterp[ipISO] );

	/* if n is equal to the number of levels, return the total recomb, else recomb for given level.	*/
	if( n == iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem] )
	{
		rate = TempInterp( TeRRCoef, TotalRecomb[ipISO][nelem], N_ISO_TE_RECOMB );
	}
	else
	{
		rate = TempInterp( TeRRCoef, RRCoef[ipISO][nelem][n], N_ISO_TE_RECOMB );
	}

	/* that was the log, now make linear */
	rate = pow( 10. , rate );

	return rate;
}

/*iso_recomb_check called by SanityCheck to confirm that recombination coef are ok
 * return value is relative error between new calculation of recom, and interp value */
double iso_recomb_check( long ipISO, long nelem, long level, double temperature )
{
	double RecombRelError ,
		RecombInterp,
		RecombCalc,
		SaveTemp;

	DEBUG_ENTRY( "iso_recomb_check()" );

	/* first set temp to desired value */
	SaveTemp = phycon.te;
	// uses overloaded version of TempChange that does not check
	// on floor since this is just a sanity check
	TempChange(temperature);

	/* actually calculate the recombination coefficient from the Milne relation,
	 * normally only done due to compile he-like command */
	RecombCalc = iso_radrecomb_from_cross_section( ipISO, temperature , nelem , level );

	/* interpolate the recombination coefficient, this is the usual method */
	RecombInterp = iso_RRCoef_Te( ipISO, nelem , level );

	/* reset temp */
	TempChange(SaveTemp);

	RecombRelError = ( RecombInterp - RecombCalc ) / MAX2( RecombInterp , RecombCalc );

	return RecombRelError;
}

/* malloc space needed for iso recombination tables */
void iso_recomb_malloc( void )
{
	DEBUG_ENTRY( "iso_recomb_malloc()" );

	NumLevRecomb = (long **)MALLOC(sizeof(long *)*(unsigned)NISO );
	TotalRecomb = (double ***)MALLOC(sizeof(double **)*(unsigned)NISO );
	RRCoef = (double ****)MALLOC(sizeof(double ***)*(unsigned)NISO );

	for( long ipISO=0; ipISO<NISO; ipISO++ )
	{
		TotalRecomb[ipISO] = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
		RRCoef[ipISO] = (double ***)MALLOC(sizeof(double **)*(unsigned)LIMELM );
		/* The number of recombination coefficients to be read from file for each element.	*/
		NumLevRecomb[ipISO] = (long*)MALLOC(sizeof(long)*(unsigned)LIMELM );

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			long int MaxLevels, maxN;

			TotalRecomb[ipISO][nelem] = (double*)MALLOC(sizeof(double)*(unsigned)N_ISO_TE_RECOMB );

			if( nelem == ipISO )
				maxN = RREC_MAXN;
			else
				maxN = LIKE_RREC_MAXN( nelem );

			NumLevRecomb[ipISO][nelem] = iso_get_total_num_levels( ipISO, maxN, 0 );

			if( nelem == ipISO || dense.lgElmtOn[nelem] )
			{
				/* must always have at least NumLevRecomb[ipISO][nelem] levels since that is number 
				* that will be read in from he rec data file, but possible to require more */
				MaxLevels = MAX2( NumLevRecomb[ipISO][nelem] , iso.numLevels_max[ipISO][nelem] );

				/* always define this */
				/* >>chng 02 jan 24, RRCoef will be iso.numLevels_max[ipISO][nelem] levels, not iso.numLevels_max,
				* code will stop if more than this is requested */
				RRCoef[ipISO][nelem] = (double**)MALLOC(sizeof(double*)*(unsigned)(MaxLevels) );

				for( long ipLo=0; ipLo < MaxLevels;++ipLo )
				{
					RRCoef[ipISO][nelem][ipLo] = (double*)MALLOC(sizeof(double)*(unsigned)N_ISO_TE_RECOMB );
				}
			}
		}
	}

	for(long i = 0; i < N_ISO_TE_RECOMB; i++)
	{
		/* this is the vector of temperatures */
		TeRRCoef[i] = 0.25*(i);
	}

	/* >>chng 06 jun 06, NP, assert thrown at T == 1e10 K, just bump the 
	 * high temperature end slightly.  */
	TeRRCoef[N_ISO_TE_RECOMB-1] += 0.01f;

	return;
}

void iso_recomb_auxiliary_free( void )
{
	DEBUG_ENTRY( "iso_recomb_auxiliary_free()" );

	for( long i = 0; i< NISO; i++ )
	{
		free( NumLevRecomb[i] );
	}
	free( NumLevRecomb );

	return;
}

void iso_recomb_setup( long ipISO )
{
	double RadRecombReturn;
	long int i, i1, i2, i3, i4, i5;
	long int ipLo, nelem;

	const int chLine_LENGTH = 1000;
	char chLine[chLine_LENGTH];
	/* this must be longer than data path string, set in path.h/cpu.cpp */
	const char* chFilename[NISO] = 
		{ "h_iso_recomb.dat", "he_iso_recomb.dat" };

	FILE *ioDATA;
	bool lgEOL;

	DEBUG_ENTRY( "iso_recomb_setup()" );

	/* if we are compiling the recombination data file, we must interpolate in temperature */
	if( iso.lgCompileRecomb[ipISO] )
	{
		iso.lgNoRecombInterp[ipISO] = false;
	}

	if( !iso.lgNoRecombInterp[ipISO] )
	{
		/******************************************************************/
		/**  Establish radiative recombination rate coefficients - RRC	***/
		/******************************************************************/
		/* This flag says we are not compiling the data file	*/
		if( !iso.lgCompileRecomb[ipISO] )
		{
			if( trace.lgTrace )
				fprintf( ioQQQ," iso_recomb_setup opening %s:", chFilename[ipISO] );

			/* Now try to read from file...*/
			try
			{
				ioDATA = open_data( chFilename[ipISO], "r" );
			}
			catch( cloudy_exit )
			{
				fprintf( ioQQQ, " Defaulting to on-the-fly computation, ");
				fprintf( ioQQQ, " but the code runs much faster if you compile %s!\n", chFilename[ipISO]);
				ioDATA = NULL;
			}
			if( ioDATA == NULL )
			{
				/* Do on the fly computation of R.R. Coef's instead.	*/
				for( nelem = ipISO; nelem < LIMELM; nelem++ )
				{
					if( dense.lgElmtOn[nelem] )
					{
						/* Zero out the recombination sum array.	*/
						for(i = 0; i < N_ISO_TE_RECOMB; i++)
						{
							TotalRecomb[ipISO][nelem][i] = 0.;
						}

						/* NumLevRecomb[ipISO][nelem] corresponds to n = 40 for H and He and 20 for ions, at present	*/
						/* There is no need to fill in values for collapsed levels, because we do not need to
						* interpolate for a given temperature, just calculate it directly with a hydrogenic routine.	*/
						for( ipLo=0; ipLo < iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem]; ipLo++ )
						{
							/* loop over temperatures to produce array of recombination coefficients	*/
							for(i = 0; i < N_ISO_TE_RECOMB; i++)
							{
								/* Store log of recombination coefficients, in N_ISO_TE_RECOMB half dec steps */
								RadRecombReturn = iso_radrecomb_from_cross_section( ipISO, pow( 10.,TeRRCoef[i] ) ,nelem,ipLo);
								TotalRecomb[ipISO][nelem][i] += RadRecombReturn;
								RRCoef[ipISO][nelem][ipLo][i] = log10(RadRecombReturn);
							}
						}
						for(i = 0; i < N_ISO_TE_RECOMB; i++)
						{
							for( i1 = iso.n_HighestResolved_max[ipISO][nelem]+1; i1< NHYDRO_MAX_LEVEL; i1++ )
							{
								TotalRecomb[ipISO][nelem][i] += t_ADfA::Inst().H_rad_rec(nelem+1-ipISO,i1, pow(10.,TeRRCoef[i]));
							}
							for( i1 = NHYDRO_MAX_LEVEL; i1<=SumUpToThisN; i1++ )
							{
								TotalRecomb[ipISO][nelem][i] += Recomb_Seaton59( nelem+1-ipISO, pow(10.,TeRRCoef[i]), i1 );
							}
							TotalRecomb[ipISO][nelem][i] = log10( TotalRecomb[ipISO][nelem][i] );
						}
					}
				}
			}
			/* Data file is present and readable...begin read.	*/
			else 
			{
				/* check that magic number is ok */
				if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				{
					fprintf( ioQQQ, " iso_recomb_setup could not read first line of %s.\n", chFilename[ipISO]);
					cdEXIT(EXIT_FAILURE);
				}
				i = 1;
				i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				i3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				i4 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				if( i1 !=RECOMBMAGIC || i2 !=NumLevRecomb[ipISO][ipISO] || i3 !=NumLevRecomb[ipISO][ipISO+1] || i4 !=N_ISO_TE_RECOMB )
				{
					fprintf( ioQQQ, 
						" iso_recomb_setup: the version of %s is not the current version.\n", chFilename[ipISO] );
					fprintf( ioQQQ, 
						" iso_recomb_setup: I expected to find the numbers  %i %li %li %i and got %li %li %li %li instead.\n" ,
						RECOMBMAGIC ,
						NumLevRecomb[ipISO][ipISO],
						NumLevRecomb[ipISO][ipISO+1],
						N_ISO_TE_RECOMB,
						i1 , i2 , i3, i4 );
					fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
					fprintf( ioQQQ, 
						" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
					cdEXIT(EXIT_FAILURE);
				}

				i5 = 1;
				/* now read in the data */
				for( nelem = ipISO; nelem < LIMELM; nelem++ )
				{
					for( ipLo=0; ipLo <= NumLevRecomb[ipISO][nelem]; ipLo++ )
					{
						i5++;
						/* get next line image */
						if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
						{
							fprintf( ioQQQ, " iso_recomb_setup could not read line %li of %s.\n", i5,
								chFilename[ipISO] );
							cdEXIT(EXIT_FAILURE);
						}
						/* each line starts with element and level number */
						i3 = 1;
						i1 = (long)FFmtRead(chLine,&i3,INPUT_LINE_LENGTH,&lgEOL);
						i2 = (long)FFmtRead(chLine,&i3,INPUT_LINE_LENGTH,&lgEOL);
						/* check that these numbers are correct */
						if( i1!=nelem || i2!=ipLo )
						{
							fprintf( ioQQQ, " iso_recomb_setup detected insanity in %s.\n", chFilename[ipISO]);
							fprintf( ioQQQ, 
								" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
							cdEXIT(EXIT_FAILURE);
						}

						/* loop over temperatures to produce array of recombination coefficients	*/
						for(i = 0; i < N_ISO_TE_RECOMB; i++)
						{
							double ThisCoef = FFmtRead(chLine,&i3,chLine_LENGTH,&lgEOL);

							if( nelem == ipISO || dense.lgElmtOn[nelem] )
							{
								/* The last line for each element is the total recombination for each temp.	*/
								if( ipLo == NumLevRecomb[ipISO][nelem] )
								{
									TotalRecomb[ipISO][nelem][i] = ThisCoef;
								}
								else
									RRCoef[ipISO][nelem][ipLo][i] = ThisCoef;
							}

							if( lgEOL )
							{
								fprintf( ioQQQ, " iso_recomb_setup detected insanity in %s.\n", chFilename[ipISO]);
								fprintf( ioQQQ, 
									" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
								cdEXIT(EXIT_FAILURE);
							}
						}
					}

					/* following loop only executed if we need more levels than are
						* stored in the recom coef data set
						* do not do collapsed levels since will use H-like recom there */
					if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
					{
						for( ipLo=NumLevRecomb[ipISO][nelem]; ipLo < iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem]; ipLo++ )
						{
								for(i = 0; i < N_ISO_TE_RECOMB; i++)
								{
									/* Store log of recombination coefficients, in N_ISO_TE_RECOMB half dec steps */
									RRCoef[ipISO][nelem][ipLo][i] = log10(iso_radrecomb_from_cross_section( ipISO, pow( 10.,TeRRCoef[i] ) ,nelem,ipLo));
								}
							}
						}
				}

				/* check that ending magic number is ok */
				if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				{
					fprintf( ioQQQ, " iso_recomb_setup could not read last line of %s.\n", chFilename[ipISO]);
					cdEXIT(EXIT_FAILURE);
				}
				i = 1;
				i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				i3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				i4 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

				if( i1 !=RECOMBMAGIC || i2 !=NumLevRecomb[ipISO][ipISO] || i3 !=NumLevRecomb[ipISO][ipISO+1] || i4 !=N_ISO_TE_RECOMB )
				{
					fprintf( ioQQQ, 
						" iso_recomb_setup: the version of %s is not the current version.\n", chFilename[ipISO] );
					fprintf( ioQQQ, 
						" iso_recomb_setup: I expected to find the numbers  %i %li %li %i and got %li %li %li %li instead.\n" ,
						RECOMBMAGIC ,
						NumLevRecomb[ipISO][ipISO],
						NumLevRecomb[ipISO][ipISO+1],
						N_ISO_TE_RECOMB,
						i1 , i2 , i3, i4 );
					fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
					fprintf( ioQQQ, 
						" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
					cdEXIT(EXIT_FAILURE);
				}

				/* close the data file */
				fclose( ioDATA );
			}
		}
		/* We are compiling the he_iso_recomb.dat data file.	*/
		else if( iso.lgCompileRecomb[ipISO] )
		{
			/* option to create table of recombination coefficients,
				* executed with the compile he-like command */
			FILE *ioRECOMB;

			ASSERT( !iso.lgNoRecombInterp[ipISO] );

			/*RECOMBMAGIC the magic number for the table of recombination coefficients */
			/*NumLevRecomb[ipISO][nelem] the number of levels that will be done */
			/* create recombination coefficients  */
			ioRECOMB = open_data( chFilename[ipISO], "w", AS_LOCAL_ONLY );
			fprintf(ioRECOMB,"%i\t%li\t%li\t%i\t%s isoelectronic sequence recomb data, created by COMPile RECOmb COEFficient H-LIke [or HE-Like] command, with %li %s levels, %li ion levels, and %i temperatures.\n",
				RECOMBMAGIC ,
				NumLevRecomb[ipISO][ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB,
				iso.chISO[ipISO],
				NumLevRecomb[ipISO][ipISO],
				elementnames.chElementSym[ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB );

			for( nelem = ipISO; nelem < LIMELM; nelem++ )
			{
				/* this must pass since compile xx-like command reset numLevels to the macro */
				ASSERT( NumLevRecomb[ipISO][nelem] <= iso.numLevels_max[ipISO][nelem] );

				/* Zero out the recombination sum array.	*/
				for(i = 0; i < N_ISO_TE_RECOMB; i++)
				{
					TotalRecomb[ipISO][nelem][i] = 0.;
				}

				for( ipLo=ipHe1s1S; ipLo < NumLevRecomb[ipISO][nelem]; ipLo++ )
				{
					fprintf(ioRECOMB, "%li\t%li", nelem, ipLo );
					/* loop over temperatures to produce array of recombination coefficients	*/
					for(i = 0; i < N_ISO_TE_RECOMB; i++)
					{
						/* Store log of recombination coefficients, in N_ISO_TE_RECOMB half dec steps */
						RadRecombReturn = iso_radrecomb_from_cross_section( ipISO, pow( 10.,TeRRCoef[i] ) ,nelem,ipLo);
						TotalRecomb[ipISO][nelem][i] += RadRecombReturn;
						RRCoef[ipISO][nelem][ipLo][i] = log10(RadRecombReturn);
						fprintf(ioRECOMB, "\t%f", RRCoef[ipISO][nelem][ipLo][i] );
					}
					fprintf(ioRECOMB, "\n" );
				}

				/* Store one additional line in XX_iso_recomb.dat that gives the total recombination,
				 * as computed by the sum so far, plus levels up to NHYDRO_MAX_LEVEL using Verner's fits,
				 * plus levels up to SumUpToThisN using Seaton 59, for each element and each temperature.	*/
				fprintf(ioRECOMB, "%li\t%li", nelem, NumLevRecomb[ipISO][nelem] );
				for(i = 0; i < N_ISO_TE_RECOMB; i++)
				{
					for( i1 = ( (nelem == ipISO) ? (RREC_MAXN + 1) : (LIKE_RREC_MAXN( nelem ) + 1) ); i1< NHYDRO_MAX_LEVEL; i1++ )
					{
						TotalRecomb[ipISO][nelem][i] += t_ADfA::Inst().H_rad_rec(nelem+1-ipISO,i1, pow(10.,TeRRCoef[i]));
					}
					for( i1 = NHYDRO_MAX_LEVEL; i1<=SumUpToThisN; i1++ )
					{
						TotalRecomb[ipISO][nelem][i] += Recomb_Seaton59( nelem+1-ipISO, pow(10.,TeRRCoef[i]), i1 );
					}
					fprintf(ioRECOMB, "\t%f", log10( TotalRecomb[ipISO][nelem][i] ) );
				}
				fprintf(ioRECOMB, "\n" );
			}
			/* end the file with the same information */
			fprintf(ioRECOMB,"%i\t%li\t%li\t%i\t%s isoelectronic sequence recomb data, created by COMPile RECOmb COEFficient [H-LIke/HE-Like] command, with %li %s levels, %li ion levels, and %i temperatures.\n",
				RECOMBMAGIC ,
				NumLevRecomb[ipISO][ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB,
				iso.chISO[ipISO],
				NumLevRecomb[ipISO][ipISO],
				elementnames.chElementSym[ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB );

			fclose( ioRECOMB );

			fprintf( ioQQQ, "iso_recomb_setup: compilation complete, %s created.\n", chFilename[ipISO] );
			fprintf( ioQQQ, "The compilation is completed successfully.\n");
			cdEXIT(EXIT_SUCCESS);
		}
	}

	return;
}

double iso_dielec_recomb_rate( long ipISO, long nelem, long ipLo )
{
	double rate;
	long ipTe, i;
	double TeDRCoef[NUM_DR_TEMPS];
	const double Te_over_Z1_Squared[NUM_DR_TEMPS] = {
		1.00000,	1.30103,	1.69897,	2.00000,	2.30103,	2.69897,	3.00000,
		3.30103,	3.69897,	4.00000,	4.30103,	4.69897,	5.00000,	5.30103,
		5.69897,	6.00000,	6.30103,	6.69897,	7.00000 };

	DEBUG_ENTRY( "iso_dielec_recomb_rate()" );

	/* currently only two iso sequences and only he-like is applicable. */
	ASSERT( ipISO == ipHE_LIKE );
	ASSERT( ipLo >= 0 );

	/* temperature grid is nelem^2 * constant temperature grid above. */
	for( i=0; i<NUM_DR_TEMPS; i++ )
	{
		TeDRCoef[i] = Te_over_Z1_Squared[i] + 2. * log10( (double) nelem );
	}

	if( ipLo == ipHe1s1S )
	{
		rate = 0.;
	}
	else if( ipLo<iso.numLevels_max[ipISO][nelem] )
	{
		if( phycon.alogte <= TeDRCoef[0] )
		{
			/* Take lowest tabulated value for low temperature end. */
			rate = iso.DielecRecombVsTemp[ipISO][nelem][ipLo][0];
		}
		else if( phycon.alogte >= TeDRCoef[NUM_DR_TEMPS-1] )
		{
			/* use T^-1.5 extrapolation at high temperatures. */
			rate = iso.DielecRecombVsTemp[ipISO][nelem][ipLo][NUM_DR_TEMPS-1] * 
				pow( 10., 1.5* (TeDRCoef[NUM_DR_TEMPS-1] - phycon.alogte ) ) ;
		}
		else
		{
			/* find temperature in tabulated values.  */
			ipTe = hunt_bisect( TeDRCoef, NUM_DR_TEMPS, phycon.alogte );			

			ASSERT( (ipTe >=0) && (ipTe < NUM_DR_TEMPS-1)  );

			if( iso.DielecRecombVsTemp[ipISO][nelem][ipLo][ipTe+1] == 0. )
				rate = 0.;
			else if( iso.DielecRecombVsTemp[ipISO][nelem][ipLo][ipTe] == 0. )
				rate = iso.DielecRecombVsTemp[ipISO][nelem][ipLo][ipTe+1];
			else
			{
				/* interpolate between tabulated points */
				rate = log10(iso.DielecRecombVsTemp[ipISO][nelem][ipLo][ipTe]) + 
					(phycon.alogte-TeDRCoef[ipTe])*
					(log10(iso.DielecRecombVsTemp[ipISO][nelem][ipLo][ipTe+1])-log10(iso.DielecRecombVsTemp[ipISO][nelem][ipLo][ipTe]))/
					(TeDRCoef[ipTe+1]-TeDRCoef[ipTe]);

				rate = pow( 10., rate );
			}
		}
	}
	else 
	{
		rate = 0.;
	}

	ASSERT( rate >= 0. && rate < 1.0e-12 );

	return rate*iso.lgDielRecom[ipISO];
}

/* TempInterp - interpolate on an array */
/** \todo	2	use a canned interpolation routine, no need for special one here */
STATIC double TempInterp( double* TempArray , double* ValueArray, long NumElements )
{
	static long int ipTe=-1;
	double rate = 0.;
	long i0;

	DEBUG_ENTRY( "TempInterp()" );

	ASSERT( fabs( 1. - (double)phycon.alogte/log10((double)phycon.te) ) < 0.0001 );

	if( ipTe < 0 )
	{
		/* te totally unknown */
		if( ( phycon.alogte < TempArray[0] ) || 
			( phycon.alogte > TempArray[NumElements-1] ) )
		{
			fprintf(ioQQQ," TempInterp called with te out of bounds \n");
			cdEXIT(EXIT_FAILURE);
		}
		ipTe = hunt_bisect( TempArray, NumElements, phycon.alogte );			
	}
	else if( phycon.alogte < TempArray[ipTe] )
	{
		/* temp is too low, must also lower ipTe */
		ASSERT( phycon.alogte > TempArray[0] );
		/* decrement ipTe until it is correct */
		while( ( phycon.alogte < TempArray[ipTe] ) && ipTe > 0)
			--ipTe;
	}
	else if( phycon.alogte > TempArray[ipTe+1] )
	{
		/* temp is too high */
		ASSERT( phycon.alogte <= TempArray[NumElements-1] );
		/* increment ipTe until it is correct */
		while( ( phycon.alogte > TempArray[ipTe+1] ) && ipTe < NumElements-1)
			++ipTe;
	}

	ASSERT( (ipTe >=0) && (ipTe < NumElements-1) );

	/* ipTe should now be valid */
	ASSERT( ( phycon.alogte >= TempArray[ipTe] )
		&& ( phycon.alogte <= TempArray[ipTe+1] ) && ( ipTe < NumElements-1 ) );

	if( ValueArray[ipTe+1] == 0. && ValueArray[ipTe] == 0. )
	{
		rate = 0.;
	}
	else
	{
		/* Do a four-point interpolation */
		const int ORDER = 3; /* order of the fitting polynomial */
		i0 = max(min(ipTe-ORDER/2,NumElements-ORDER-1),0);
		rate = lagrange( &TempArray[i0], &ValueArray[i0], ORDER+1, phycon.alogte );
	}

	return rate;
}
