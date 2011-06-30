/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeCollid evaluate collisional rates */
/*HeCSInterp interpolate on He1 collision strengths */
/*AtomCSInterp do the atom	*/
/*IonCSInterp do the ions	*/
/*CS_l_mixing_PS64 - find rate for l-mixing collisions by protons, for neutrals */
#include "cddefines.h"
#include "atmdat.h"
#include "conv.h"
#include "dense.h"
#include "helike.h"
#include "helike_cs.h"
#include "hydro_vs_rates.h"
#include "iso.h"
#include "lines_service.h"
#include "opacity.h"
#include "phycon.h"
#include "physconst.h"
#include "rfield.h"
#include "taulines.h"
#include "thirdparty.h"
#include "trace.h"

/* returns thermally-averaged Seaton 62 collision strength. */
STATIC double S62_Therm_ave_coll_str( double proj_energy_overKT, long nelem, long Collider, double deltaE, double osc_strength, double temp,
	double stat_weight, double I_energy_eV );

/* all of these are used in the calculation of Stark collision strengths
 * following the algorithms in Vrinceanu & Flannery (2001). */
STATIC double collision_strength_VF01( long ipISO, long nelem, long n, long l, long lp, long s, long Collider, 
	double ColliderCharge, double temp, double velOrEner, bool lgParamIsRedVel );
STATIC double L_mix_integrand_VF01( long n, long l, long lp, double bmax, double red_vel, double an, double ColliderCharge, double alpha );
STATIC double StarkCollTransProb_VF01( long int n, long int l, long int lp, double alpha, double deltaPhi);


class my_Integrand_S62
{
public:
	long nelem, Collider;
	double deltaE, osc_strength, temp, stat_weight, I_energy_eV;

	double operator() (double proj_energy_overKT)
	{
		double col_str = S62_Therm_ave_coll_str( proj_energy_overKT, nelem, Collider, deltaE, osc_strength,
			temp, stat_weight, I_energy_eV );
		return col_str;
	}
};

class my_Integrand_VF01_E
{
public:
	long ipISO, nelem, n, l, lp, s, Collider;
	double ColliderCharge, temp, velOrEner;
	bool lgParamIsRedVel;

	double operator() (double EOverKT)
	{
		double col_str = collision_strength_VF01( ipISO, nelem, n, l, lp, s, Collider,
			ColliderCharge, temp, EOverKT * temp / TE1RYD, lgParamIsRedVel );
		return exp( -1.*EOverKT ) * col_str;
	}
};

class my_Integrand_VF01_alpha
{
public:
	long n, l, lp;
	double bmax, red_vel, an, ColliderCharge;

	double operator() (double alpha)
	{
		double integrand = L_mix_integrand_VF01( n, l, lp,
			bmax, red_vel, an, ColliderCharge, alpha );
		return integrand;
	}
};


/* These are masses relative to the proton mass of the electron, proton, he+, and alpha particle. */
static const double ColliderMass[4] = {ELECTRON_MASS/PROTON_MASS, 1.0, 4.0, 4.0};
static const double ColliderCharge[4] = {1.0, 1.0, 1.0, 2.0};

void HeCollidSetup( void )
{
	/* this must be longer than data path string, set in path.h/cpu.cpp */
	long i, i1, j, nelem, ipHi, ipLo;
	bool lgEOL, lgHIT;
	FILE *ioDATA;
	long ipISO = ipHE_LIKE;

#	define chLine_LENGTH 1000
	char chLine[chLine_LENGTH];

	DEBUG_ENTRY( "HeCollidSetup()" );

	/* get the collision strength data for the He 1 lines */
	if( trace.lgTrace )
		fprintf( ioQQQ," HeCollidSetup opening he1_cs.dat:");

	ioDATA = open_data( "he1_cs.dat", "r" );

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " HeCollidSetup could not read first line of he1_cs.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	/*i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	i3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);*/

	/* the following is to check the numbers that appear at the start of he1_cs.dat */
	if( i1 !=COLLISMAGIC )
	{
		fprintf( ioQQQ, 
			" HeCollidSetup: the version of he1_cs.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			" HeCollidSetup: I expected to find the number %i and got %li instead.\n" ,
			COLLISMAGIC, i1 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* get the array of temperatures */
	lgHIT = false;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		if( chLine[0] != '#')
		{
			lgHIT = true;
			iso.nCS[ipISO] = 0;
			lgEOL = false;
			i = 1;
			while( !lgEOL && iso.nCS[ipISO] < HE1CSARRAY)
			{
				iso.CSTemp[ipISO][iso.nCS[ipISO]] = FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
				++iso.nCS[ipISO];
			}
			break;
		}
	}
	--iso.nCS[ipISO];
	if( !lgHIT )
	{
		fprintf( ioQQQ, " HeCollidSetup could not find line in CS temperatures in he1_cs.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* create space for array of CS values, if not already done */
	iso.HeCS.reserve( NISO );
	iso.HeCS.reserve( ipISO, LIMELM );

	for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
	{
		/** \todo	2	- this structure is currently only used for helium itself...
		 * stuff numbers in for other elements, or drop the [nelem] dimension off
		 * of iso.HeCS[ipISO]	*/
		if( nelem != ipHELIUM )
			continue;

		/* only grab core for elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem])
		{
			iso.HeCS.reserve( ipISO, nelem, iso.numLevels_max[ipHE_LIKE][nelem]-iso.nCollapsed_max[ipHE_LIKE][nelem] );

			for( ipHi=ipHe2s3S; ipHi < iso.numLevels_max[ipHE_LIKE][nelem] - iso.nCollapsed_max[ipHE_LIKE][nelem];++ipHi )
			{
				iso.HeCS.reserve( ipISO, nelem, ipHi, ipHi );

				for( ipLo=ipHe1s1S; ipLo<ipHi; ++ipLo )
					iso.HeCS.reserve( ipISO, nelem, ipHi, ipLo, iso.nCS[ipISO] );
			}
		}
	}

	iso.HeCS.alloc();
	iso.HeCS.zero();

	/* now read in the CS data */
	lgHIT = false;
	nelem = ipHELIUM;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		char *chTemp;
		/* only look at lines without '#' in first col */
		if( (chLine[0] == ' ') || (chLine[0]=='\n') )
			break;
		if( chLine[0] != '#')
		{
			lgHIT = true;

			/* get lower and upper level index */
			j = 1;
			ipLo = (long)FFmtRead(chLine,&j,INPUT_LINE_LENGTH,&lgEOL);
			ipHi = (long)FFmtRead(chLine,&j,INPUT_LINE_LENGTH,&lgEOL);
			ASSERT( ipLo < ipHi );
			if( ipHi >= iso.numLevels_max[ipHE_LIKE][nelem] - iso.nCollapsed_max[ipHE_LIKE][nelem] )
				continue;
			else
			{
				chTemp = chLine;
				/* skip over 4 tabs to start of cs data */
				for( i=0; i<3; ++i )
				{
					if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
					{
						fprintf( ioQQQ, " HeCollidSetup no init cs\n" );
						cdEXIT(EXIT_FAILURE);
					}
					++chTemp;
				}

				for( i=0; i<iso.nCS[ipISO]; ++i )
				{
					double a;
					if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
					{
						fprintf( ioQQQ, " HeCollidSetup not scan cs, current indices: %li %li\n", ipHi, ipLo );
						cdEXIT(EXIT_FAILURE);
					}
					++chTemp;
					sscanf( chTemp , "%le" , &a );
					iso.HeCS[ipISO][nelem][ipHi][ipLo][i] = (realnum)a;
				}
			}
		}
	}

	/* close the data file */
	fclose( ioDATA );

	return;
}

/* Choose either AtomCSInterp or IonCSInterp */
realnum HeCSInterp(long int nelem,
				 long int ipHi,
				 long int ipLo,
				 long int Collider )
{
	realnum cs, factor1;

	/* This variable is for diagnostic purposes:
	 * a string used in the output to designate where each cs comes from.	*/	
	const char *where = "      ";

	DEBUG_ENTRY( "HeCSInterp()" );

	if( !iso.lgColl_excite[ipHE_LIKE] )
	{
		return (realnum)1E-10;
	}

	if( nelem == ipHELIUM )
	{
		/* do for helium */
		cs = AtomCSInterp( nelem, ipHi , ipLo, &factor1, &where, Collider );
	}
	else
	{
		/* get collision strengths for an ion */
		cs = IonCSInterp( nelem, ipHi , ipLo, &factor1, &where, Collider );
	}

	ASSERT( cs >= 0.f );

	/* in many cases the correction factor for split states has already been made,
	 * if not then factor is still negative */
	/* Remove the second test here when IonCSInterp is up to par with AtomCSInterp */
	ASSERT( factor1 >= 0.f || nelem!=ipHELIUM );
	if( factor1 < 0.f )
	{
		ASSERT( iso.lgCS_Vriens[ipHE_LIKE] );

		factor1 = 1.f;
	}

	/* take factor into account, usually 1, ratio of stat weights if within 2 3P 
	 * and with collisions from collapsed to resolved levels */
	cs *= factor1;

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/

		if( DEBUG_LOC && ( nelem==ipOXYGEN ) && (cs > 0.f) && (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n == 2) 
			&& ( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n <= 2 ) )
			fprintf(ioQQQ,"%li\t%li\t%li\t-\t%li\t%li\t%li\t%.2e\t%s\n", 
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n , StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S ,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l ,	StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n ,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S , StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l , cs,where);
	}

	return MAX2(cs, 1.e-10f);
}

realnum AtomCSInterp(long int nelem,
				   long int ipHi,
				   long int ipLo,
				   realnum *factor1,
				   const char **where,
				   long int Collider )
{
	long ipArray;
	realnum cs;

	DEBUG_ENTRY( "AtomCSInterp()" );

	ASSERT( nelem == ipHELIUM );

	/* init values, better be positive when we exit */
	cs = -1.f; 

	/* this may be used for splitting up the collision strength within 2 3P when
	 * the lower level is withint 2 3P, and for collisions between resolved and collapsed levels.
	 * It may be set somewhere in routine, so set to negative value here as flag saying not set */
	*factor1 = -1.f;

	/* for most of the helium iso sequence, the order of the J levels within 2 3P 
	 * in increasing energy, is 0 1 2 - the exception is atomic helium itself,
	 * which is swapped, 2 1 0 */

	/* this branch is for upper and lower levels within 2p3P */
	if( ipLo >= ipHe2p3P0 && ipHi <= ipHe2p3P2 && Collider==ipELECTRON )
	{
		*factor1 = 1.f;
		/* atomic helium, use Berrington private comm cs */

		/* >>refer	he1	cs	Berrington, Keith, 2001, private communication - email follows
		> Dear Gary,
		> I could not find any literature on the He fine-structure
		> problem (but I didn't look very hard, so there may be 
		> something somewhere). However, I did a quick R-matrix run to 
		> see what the magnitude of the collision strengths are... At 
		> 1000K, I get the effective collision strength for 2^3P J=0-1, 
		>  0-2, 1-2 as 0.8,0.7,2.7; for 10000K I get 1.2, 2.1, 6.0
		*/
		/* indices are the same and correct, only thing special is that energies are in inverted order...was right first time.	*/
		if( ipLo == ipHe2p3P0 && ipHi == ipHe2p3P1 )
		{
			cs = 1.2f;
		}
		else if( ipLo == ipHe2p3P0 && ipHi == ipHe2p3P2 )
		{
			cs = 2.1f;
		}
		else if( ipLo == ipHe2p3P1 && ipHi == ipHe2p3P2 )
		{
			cs = 6.0f;
		}
		else
		{
			cs = 1.0f;
			TotalInsanity();
		}

		*where = "Berr  ";
		/* statistical weights included */
	}
	/* >>chng 02 feb 25, Bray data should come first since it is the best we have.	*/
	/* this branch is the Bray et al data, for n <= 5, where quantal calcs exist 
	 * must exclude ipLo >= ipHe2p1P because they give no numbers for those	*/
	else if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n <= 5 && 
		( ipHi < iso.numLevels_max[ipHE_LIKE][nelem] - iso.nCollapsed_max[ipHE_LIKE][nelem] ) &&
		iso.HeCS[ipHE_LIKE][nelem][ipHi][ipLo][0] >= 0.f && Collider== ipELECTRON )
	{
		ASSERT( *factor1 == -1.f );
		ASSERT( ipLo < ipHi );
		ASSERT( ipHe2p3P0 == 3 );

		/* ipLo is within 2^3P	*/
		if( ipLo >= ipHe2p3P0 && ipLo <= ipHe2p3P2 )
		{
			/* *factor1 is ratio of statistical weights of level to term */

			/* ipHe2p3P0, ipHe2p3P1, ipHe2p3P2 have indices 3,4,and 5, but j=0,1,and 2.	*/
			*factor1 = (2.f*((realnum)ipLo-3.f)+1.f) / 9.f;

			/* ipHi must be above ipHe2p3P2 since 2p3Pj->2p3Pk taken care of above	*/
			ASSERT( ipHi > ipHe2p3P2 );
		}
		/* ipHi is within 2^3P	*/
		else if( ipHi >= ipHe2p3P0 && ipHi <= ipHe2p3P2 )
		{
			ASSERT( ipLo < ipHe2p3P0 );

			*factor1 = (2.f*((realnum)ipHi-3.f)+1.f) / 9.f;
		}
		/* neither are within 2^3P...no splitting necessary	*/
		else 
		{
			*factor1 = 1.f;
		}

		/* SOME OF THESE ARE NOT N-CHANGING!	*/
		/* Must be careful about turning each one on or off.	*/

		/* this is the case where we have quantal calculations */
		/* >>refer	He1	cs	Bray, I., Burgess, A., Fursa, D.V., & Tully, J.A., 2000, A&AS, 146, 481-498 */
		/* check whether we are outside temperature array bounds,
		 * and use extreme value if we are */
		if( phycon.alogte <= iso.CSTemp[ipHE_LIKE][0] )
		{
			cs = iso.HeCS[ipHE_LIKE][nelem][ipHi][ipLo][0];
		}
		else if( phycon.alogte >= iso.CSTemp[ipHE_LIKE][iso.nCS[ipHE_LIKE]-1] )
		{
			cs = iso.HeCS[ipHE_LIKE][nelem][ipHi][ipLo][iso.nCS[ipHE_LIKE]-1];
		}
		else
		{
			realnum flow; 

			/* find which array element within the cs vs temp array */
			ipArray = (long)((phycon.alogte - iso.CSTemp[ipHE_LIKE][0])/(iso.CSTemp[ipHE_LIKE][1]-iso.CSTemp[ipHE_LIKE][0]));
			ASSERT( ipArray < iso.nCS[ipHE_LIKE] );
			ASSERT( ipArray >= 0 );
			/* when taking the average, this is the fraction from the lower temperature value */
			flow = (realnum)( (phycon.alogte - iso.CSTemp[ipHE_LIKE][ipArray])/
				(iso.CSTemp[ipHE_LIKE][ipArray+1]-iso.CSTemp[ipHE_LIKE][ipArray]));
			ASSERT( (flow >= 0.f) && (flow <= 1.f) );

			cs = iso.HeCS[ipHE_LIKE][nelem][ipHi][ipLo][ipArray] * (1.f-flow) +
				iso.HeCS[ipHE_LIKE][nelem][ipHi][ipLo][ipArray+1] * flow;
		}

		*where = "Bray ";

		/* options to kill collisional excitation and/or l-mixing	*/
		if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n == StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n )
			/* iso.lgColl_l_mixing turned off with atom he-like l-mixing collisions off command */
			cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
		else
		{
			/* iso.lgColl_excite turned off with atom he-like collisional excitation off command */
			cs *= (realnum)iso.lgColl_excite[ipHE_LIKE];
		}

		ASSERT( cs >= 0.f );
		/* statistical weights included */
	}
	/* this branch, n-same, l-changing collision, and not case of He with quantal data */
	else if( (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n == StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n ) &&
		(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S == StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S ) )
	{
		ASSERT( *factor1 == -1.f );
		*factor1 = 1.f;

		/* ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n >= 3 ); */
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n <= iso.n_HighestResolved_max[ipHE_LIKE][nelem] );

		if( (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l <=2) &&
			abs(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l - StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l)== 1 )
		{
			/* Use the method given in 
			 * >>refer He	CS	Seaton, M. J. 1962, Proc. Phys. Soc. 79, 1105 
			 * statistical weights included */
			cs = (realnum)CS_l_mixing_S62(ipHE_LIKE, nelem, ipLo, ipHi, (double)phycon.te, Collider); 
		}
		else if( iso.lgCS_Vrinceanu[ipHE_LIKE] )
		{
			if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l >=3 &&
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l >=3 )
			{
				/* Use the method given in 
				 * >>refer He	CS	Vrinceanu, D. \& Flannery, M. R. 2001, PhysRevA 63, 032701 
				 * statistical weights included */
				cs = (realnum)CS_l_mixing_VF01( ipHE_LIKE,
					nelem,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S,
					(double)phycon.te,
					Collider );
			}
			else
			{
				cs = 0.f;
			}
		}
		/* this branch, l changing by one */
		else if( abs(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l - StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l)== 1)
		{
			/* >>refer	He	cs	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
			/* statistical weights included */
			cs = (realnum)CS_l_mixing_PS64( 
				nelem,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].lifetime,
				nelem+1.-ipHE_LIKE,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].g,
				Collider);
		}
		else
		{
			/* l changes by more than 1, but same-n collision */
			cs = 0.f;
		}

		/* ipLo is within 2^3P	*/
		if( ipLo >= ipHe2p3P0 && ipLo <= ipHe2p3P2 )
		{
			*factor1 = (2.f*((realnum)ipLo-3.f)+1.f) / 9.f;
		}

		/* ipHi is within 2^3P	*/
		if( ipHi >= ipHe2p3P0 && ipHi <= ipHe2p3P2 )
		{
			*factor1 = (2.f*((realnum)ipHi-3.f)+1.f) / 9.f;
		}

		*where = "lmix  ";
		cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
	}

	/* this is an atomic n-changing collision with no quantal calculations */
	else if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n != StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n )
	{
		ASSERT( *factor1 == -1.f );
		/* this is an atomic n-changing collision with no quantal calculations */
		/* gbar g-bar goes here */

		/* >>chng 02 jul 25, add option for fits to quantal cs data */
		if( iso.lgCS_Vriens[ipHE_LIKE] )
		{			
			/* >>refer He CS	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
			 * statistical weight IS included in the routine */
			cs = (realnum)CS_VS80(
							ipHE_LIKE,
							nelem,
							ipHi,
							ipLo,
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul,
							phycon.te,
							Collider );

			*factor1 = 1.f;
			*where = "Vriens";
		}
		else if( iso.lgCS_None[ipHE_LIKE] )
		{
			cs = 0.f;
			*factor1 = 1.f;
			*where = "no gb";
		}
		else if( iso.nCS_new[ipHE_LIKE] )
		{
			*factor1 = 1.f;
			/* Don't know if stat weights are included in this, but they're probably
			 * wrong anyway since they are based in part on the former (incorrect)
			 * implementation of Vriens and Smeets rates */

			/* two different fits, allowed and forbidden */
			if( Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul > 1. )
			{
				/* permitted lines - large A */
				double x = 
					log10(MAX2(34.7,Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyWN));

				if( iso.nCS_new[ipHE_LIKE] == 1 )
				{
					/* this is the broken power law fit, passing through both quantal
					 * calcs at high energy and asymptotically goes to VS at low energies */
					if( x < 4.5 )
					{
						/* low energy fit for permitted transitions */
						cs = (realnum)pow( 10. , -1.45*x + 6.75);
					}
					else
					{
						/* higher energy fit for permitted transitions */
						cs = (realnum)pow( 10. , -3.33*x+15.15);
					}
				}
				else if( iso.nCS_new[ipHE_LIKE] == 2 )
				{
					/* single parallel fit for permitted transitions, runs parallel to VS */
					cs = (realnum)pow( 10. , -2.3*x+10.3);
				}
				else
					TotalInsanity();
			}
			else
			{
				/* fit for forbidden transitions */
				if( Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyWN < 25119.f )
				{
					cs = 0.631f; 
				}
				else
				{
					cs = (realnum)pow(10., 
						-3.*log10(Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyWN)+12.8);
				}
			}

			*where = "newgb";
		}
		else
			TotalInsanity();

		/* ipLo is within 2^3P	*/
		if( ipLo >= ipHe2p3P0 && ipLo <= ipHe2p3P2 )
		{
			*factor1 = (2.f*((realnum)ipLo-3.f)+1.f) / 9.f;
		}

		/* ipHi is within 2^3P	*/
		if( ipHi >= ipHe2p3P0 && ipHi <= ipHe2p3P2 )
		{
			*factor1 = (2.f*((realnum)ipHi-3.f)+1.f) / 9.f;
		}

		/* options to turn off collisions */
		cs *= (realnum)iso.lgColl_excite[ipHE_LIKE];

	}
	else
	{
		/* If spin changing collisions are prohibited in the l-mixing routine,
		 * they will fall here, and will have been assigned no collision strength.	
		 * Assign zero for now.	*/
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S != StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S );
		cs = 0.f;
		*factor1 = 1.f;
	}

	ASSERT( cs >= 0.f );

	return(cs);

}

/* IonCSInterp interpolate on collision strengths for element nelem */
realnum IonCSInterp( long nelem , long ipHi , long ipLo, realnum *factor1, const char **where, long Collider  )
{
	realnum cs;

	DEBUG_ENTRY( "IonCSInterp()" );

	ASSERT( nelem > ipHELIUM );
	ASSERT( nelem < LIMELM );

	/* init values, better be positive when we exit */
	cs = -1.f; 

	/* this may be used for splitting up the collision strength for collisions between
	 * resolved and collapsed levels.  It may be set somewhere in routine, so set to 
	 * negative value here as flag saying not set */
	*factor1 = -1.f;


	/* >>chng 02 mar 04,  the approximation here for transitions within 2p3P was not needed,
	 * because the Zhang data give these transitions.  They are of the same order, but are 
	 * specific to the three transitions	*/

	/* this branch is ground to n=2 or from n=2 to n=2, for ions only	*/
	/*>>refer Helike	CS	Zhang, Honglin, & Sampson, Douglas H. 1987, ApJS 63, 487	*/
	if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n==2 
		&& StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n<=2 && Collider==ipELECTRON)
	{
		*where = "Zhang";
		*factor1 = 1.;

		/* Collisions from gound	*/
		if( ipLo == ipHe1s1S )
		{
			switch( ipHi )
			{
			case 1:	/* to 2tripS	*/
				cs = 0.25f/POW2(nelem+1.f);
				break;
			case 2: /* to 2singS	*/
				cs = 0.4f/POW2(nelem+1.f);
				break;
			case 3: /* to 2tripP0	*/
				cs = 0.15f/POW2(nelem+1.f);
				break;
			case 4: /* to 2tripP1	*/
				cs = 0.45f/POW2(nelem+1.f);
				break;
			case 5: /* to 2tripP2	*/
				cs = 0.75f/POW2(nelem+1.f);
				break;
			case 6: /* to 2singP	*/
				cs = 1.3f/POW2(nelem+1.f);
				break;
			default:
				TotalInsanity();
				break;
			}
			cs *= (realnum)iso.lgColl_excite[ipHE_LIKE];
		}
		/* collisions from 2tripS to n=2	*/
		else if( ipLo == ipHe2s3S )
		{
			switch( ipHi )
			{
			case 2: /* to 2singS	*/
				cs = 2.75f/POW2(nelem+1.f);
				break;
			case 3: /* to 2tripP0	*/
				cs = 60.f/POW2(nelem+1.f);
				break;
			case 4: /* to 2tripP1	*/
				cs = 180.f/POW2(nelem+1.f);
				break;
			case 5: /* to 2tripP2	*/
				cs = 300.f/POW2(nelem+1.f);
				break;
			case 6: /* to 2singP	*/
				cs = 5.8f/POW2(nelem+1.f);
				break;
			default:
				TotalInsanity();
				break;
			}
			cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
		}
		/* collisions from 2singS to n=2	*/
		else if( ipLo == ipHe2s1S )
		{
			switch( ipHi )
			{
			case 3: /* to 2tripP0	*/
				cs = 0.56f/POW2(nelem+1.f);
				break;
			case 4: /* to 2tripP1	*/
				cs = 1.74f/POW2(nelem+1.f);
				break;
			case 5: /* to 2tripP2	*/
				cs = 2.81f/POW2(nelem+1.f);
				break;
			case 6: /* to 2singP	*/
				cs = 190.f/POW2(nelem+1.f);
				break;
			default:
				TotalInsanity();
				break;
			}
			cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
		}
		/* collisions from 2tripP0 to n=2	*/
		else if( ipLo == ipHe2p3P0 )
		{
			switch( ipHi )
			{
			case 4: /* to 2tripP1	*/
				cs = 8.1f/POW2(nelem+1.f);
				break;
			case 5: /* to 2tripP2	*/
				cs = 8.2f/POW2(nelem+1.f);
				break;
			case 6: /* to 2singP	*/
				cs = 3.9f/POW2(nelem+1.f);
				break;
			default:
				TotalInsanity();
				break;
			}
			cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
		}
		/* collisions from 2tripP1 to n=2	*/
		else if( ipLo == ipHe2p3P1 )
		{
			switch( ipHi )
			{
			case 5: /* to 2tripP2	*/
				cs = 30.f/POW2(nelem+1.f);
				break;
			case 6: /* to 2singP	*/
				cs = 11.7f/POW2(nelem+1.f);
				break;
			default:
				TotalInsanity();
				break;
			}
			cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
		}
		/* collisions from 2tripP2 to n=2	*/
		else if( ipLo == ipHe2p3P2 )
		{
			/* to 2singP	*/
			cs = 19.5f/POW2(nelem+1.f) * (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
		}
		else
			TotalInsanity();

		/* statistical weights included */
	}

	/* this branch, n-same, l-changing collisions */
	else if( (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n == StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n ) &&
		(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S == StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S ) )
	{
		ASSERT( *factor1 == -1.f );
		*factor1 = 1.f;

		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n <= iso.n_HighestResolved_max[ipHE_LIKE][nelem] );

		if( (StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l <=2) &&
			abs(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l - StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l)== 1 )
		{
			/* Use the method given in 
			 * >>refer He	CS	Seaton, M. J. 1962, Proc. Phys. Soc. 79, 1105 
			 * statistical weights included */
			cs = (realnum)CS_l_mixing_S62(ipHE_LIKE, nelem, ipLo, ipHi, (double)phycon.te, Collider); 
		}
		else if( iso.lgCS_Vrinceanu[ipHE_LIKE] )
		{
			if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l >=3 &&
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l >=3 )
			{
				/* Use the method given in 
				 * >>refer He	CS	Vrinceanu, D. \& Flannery, M. R. 2001, PhysRevA 63, 032701 
				 * statistical weights included */
				cs = (realnum)CS_l_mixing_VF01( ipHE_LIKE,
					nelem,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l,
					StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S,
					(double)phycon.te,
					Collider );
			}
			else
			{
				cs = 0.f;
			}
		}
		/* this branch, l changing by one */
		else if( abs(StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l - StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l)== 1)
		{
			/* >>refer	He	cs	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
			/* statistical weights included */
			cs = (realnum)CS_l_mixing_PS64( 
				nelem,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].lifetime,
				nelem+1.-ipHE_LIKE,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l,
				StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].g,
				Collider);
		}
		else
		{
			/* l changes by more than 1, but same-n collision */
			cs = 0.f;
		}

		/* ipHi is within 2^3P, do not need to split on ipLo.	*/
		if( ipHi >= ipHe2p3P0 && ipHi <= ipHe2p3P2 )
		{
			*factor1 = (2.f*((realnum)ipHi-3.f)+1.f) / 9.f;
		}

		*where = "lmix  ";
		cs *= (realnum)iso.lgColl_l_mixing[ipHE_LIKE];
	}

	/* this branch, n changing collisions for ions */
	else if( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n != StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n )
	{
		if( iso.lgCS_Vriens[ipHE_LIKE] )
		{
			/* this is the default branch */
			/* >>refer He CS	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
			 * statistical weight is NOT included in the routine */
			cs = (realnum)CS_VS80(
							ipHE_LIKE,
							nelem,
							ipHi,
							ipLo,
							Transitions[ipHE_LIKE][nelem][ipHi][ipLo].Emis->Aul,
							phycon.te,
							Collider );

			*factor1 = 1.f;
			*where = "Vriens";
		}
		else
		{
			/* \todo 2 this branch and above now do the same thing.  Change something. */
			/* this branch is for testing and reached with command ATOM HE-LIKE COLLISIONS VRIENS OFF */
			fixit(); /* use Percival and Richards here */

			cs = 0.f;
			*where = "hydro";
		}
	}
	else
	{
		/* what's left are deltaN=0, spin changing collisions.
		 * These have not been accounted for.	*/
		/* Make sure what comes here is what we think it is.	*/
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n == StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n );
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S != StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S );
		cs = 0.f;
		*where = "spin ";
	}

	ASSERT( cs >= 0.f );

	return(cs);
}


/*CS_l_mixing_S62 - find rate for l-mixing collisions by protons, for neutrals */
/* The S62 stands for Seaton 1962 */
double CS_l_mixing_S62(
	long ipISO,
	long nelem /* the chemical element, 1 for He */,
	long ipLo /* lower level, 0 for ground */,
	long ipHi /* upper level, 0 for ground */,
	double temp,
	long Collider)
{
	/* >>refer	He	l-mixing	Seaton, M.J., 1962, Proc. Phys. Soc. */
	double coll_str;

	DEBUG_ENTRY( "CS_l_mixing_S62()" );

	if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont <= 0 )
	{
		return 0.;
	}

	my_Integrand_S62 func;

	func.temp = temp;
	func.stat_weight = StatesElemNEW[nelem][nelem-ipISO][ipLo].g;
	/* >>chng 05 sep 06, RP  - update energies of excited states */
	/* func.deltaE = EVRYD*(iso.xIsoLevNIonRyd[ipISO][nelem][ipLo] -
		iso.xIsoLevNIonRyd[ipISO][nelem][ipHi]); */
	func.deltaE = Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg/EN1EV;
	func.I_energy_eV = EVRYD*iso.xIsoLevNIonRyd[ipISO][nelem][0];
	func.Collider = Collider;
	func.nelem = nelem;

	ASSERT( TRANS_PROB_CONST*POW2(func.deltaE/WAVNRYD/EVRYD) > 0. );

	func.osc_strength = Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul/
		(TRANS_PROB_CONST*POW2(func.deltaE/WAVNRYD/EVRYD));

	Integrator<my_Integrand_S62,Gaussian32> S62;
	/* This returns a thermally averaged collision strength */
	coll_str = S62.sum( 0., 1., func );
	coll_str += S62.sum( 1., 10., func );
	ASSERT( coll_str > 0. );

	return coll_str;
}

/* The integrand for calculating the thermal average of collision strengths */
STATIC double S62_Therm_ave_coll_str( double proj_energy_overKT, long nelem, long Collider, double deltaE, double osc_strength, double temp,
	double stat_weight, double I_energy_eV )
{

	double integrand, cross_section, /*Rnot,*/ coll_str, zOverB2;
	double /* betanot, */ betaone, zeta_of_betaone, /* cs1, */ cs2;
	double Dubya, proj_energy;
	double reduced_mass;

	DEBUG_ENTRY( "S62_Therm_ave_coll_str()" );

	/* projectile energy in eV */
	proj_energy = proj_energy_overKT * phycon.te / EVDEGK;

	reduced_mass = dense.AtomicWeight[nelem]*ColliderMass[Collider]/
		(dense.AtomicWeight[nelem]+ColliderMass[Collider])*ATOMIC_MASS_UNIT;

	/* Rnot = 1.1229*H_BAR/sqrt(ELECTRON_MASS*deltaE*EN1EV)/Bohr_radius; in units of Bohr_radius */

	proj_energy *= ColliderMass[ipELECTRON]/ColliderMass[Collider];
	/* The deltaE here is to make sure that the collider has no less
	 * than the energy difference between the initial and final levels. */
	proj_energy += deltaE;
	Dubya = 0.5*(2.*proj_energy-deltaE);
	ASSERT( Dubya > 0. );

	/* betanot = sqrt(proj_energy/I_energy_eV)*(deltaE/2./Dubya)*Rnot; */

	ASSERT( I_energy_eV > 0. );
	ASSERT( osc_strength > 0. );

	/* from equation 33 */
	zOverB2 = 0.5*POW2(Dubya/deltaE)*deltaE/I_energy_eV/osc_strength;

	ASSERT( zOverB2 > 0. );

	if( zOverB2 > 100. )
	{
		betaone = sqrt( 1./zOverB2 );
	}
	else if( zOverB2 < 0.54 )
	{
		/* Low betaone approximation */
		betaone = (1./3.)*(log(PI)-log(zOverB2)+1.28);
		if( betaone > 2.38 )
		{
			/* average with this over approximation */
			betaone += 0.5*(log(PI)-log(zOverB2));
			betaone *= 0.5;
		}
	}
	else
	{
		long ip_zOverB2 = 0;
		double zetaOVERbeta2[101] = {
			1.030E+02,9.840E+01,9.402E+01,8.983E+01,8.583E+01,8.200E+01,7.835E+01,7.485E+01,
			7.151E+01,6.832E+01,6.527E+01,6.236E+01,5.957E+01,5.691E+01,5.436E+01,5.193E+01,
			4.961E+01,4.738E+01,4.526E+01,4.323E+01,4.129E+01,3.943E+01,3.766E+01,3.596E+01,
			3.434E+01,3.279E+01,3.131E+01,2.989E+01,2.854E+01,2.724E+01,2.601E+01,2.482E+01,
			2.369E+01,2.261E+01,2.158E+01,2.059E+01,1.964E+01,1.874E+01,1.787E+01,1.705E+01,
			1.626E+01,1.550E+01,1.478E+01,1.409E+01,1.343E+01,1.280E+01,1.219E+01,1.162E+01,
			1.107E+01,1.054E+01,1.004E+01,9.554E+00,9.094E+00,8.655E+00,8.234E+00,7.833E+00,
			7.449E+00,7.082E+00,6.732E+00,6.397E+00,6.078E+00,5.772E+00,5.481E+00,5.202E+00,
			4.937E+00,4.683E+00,4.441E+00,4.210E+00,3.989E+00,3.779E+00,3.578E+00,3.387E+00,
			3.204E+00,3.031E+00,2.865E+00,2.707E+00,2.557E+00,2.414E+00,2.277E+00,2.148E+00,
			2.024E+00,1.907E+00,1.795E+00,1.689E+00,1.589E+00,1.493E+00,1.402E+00,1.316E+00,
			1.235E+00,1.157E+00,1.084E+00,1.015E+00,9.491E-01,8.870E-01,8.283E-01,7.729E-01,
			7.206E-01,6.712E-01,6.247E-01,5.808E-01,5.396E-01};

		ASSERT( zOverB2 >= zetaOVERbeta2[100] );

		/* find beta in the table */
		for( long i=0; i< 100; ++i )
		{
			if( ( zOverB2 < zetaOVERbeta2[i] ) && ( zOverB2 >= zetaOVERbeta2[i+1] ) )
			{
				/* found the temperature - use it */
				ip_zOverB2 = i;
				break;
			}
		}

		ASSERT( (ip_zOverB2 >=0) && (ip_zOverB2 < 100) );

		betaone = (zOverB2 - zetaOVERbeta2[ip_zOverB2]) / 
			(zetaOVERbeta2[ip_zOverB2+1] - zetaOVERbeta2[ip_zOverB2]) *
			(pow(10., ((double)ip_zOverB2+1.)/100. - 1.) - pow(10., ((double)ip_zOverB2)/100. - 1.)) +
			pow(10., (double)ip_zOverB2/100. - 1.);

	}

	zeta_of_betaone = zOverB2 * POW2(betaone);

	/* cs1 = betanot * bessel_k0(betanot) * bessel_k1(betanot); */
	cs2 = 0.5*zeta_of_betaone + betaone * bessel_k0(betaone) * bessel_k1(betaone);

	/* cross_section = MIN2(cs1, cs2); */
	cross_section = cs2;

	/* cross section in units of PI * a_o^2 */
	cross_section *= 8. * (I_energy_eV/deltaE) * osc_strength * (I_energy_eV/proj_energy);

	/* convert to collision strength */
	coll_str = ConvCrossSect2CollStr( cross_section * PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM, stat_weight, proj_energy/EVRYD, reduced_mass );

	integrand = exp( -1.*(proj_energy-deltaE)*EVDEGK/temp ) * coll_str;
	return integrand;
}

/*CS_l_mixing_PS64 - find rate for l-mixing collisions by protons, for neutrals */
double CS_l_mixing_PS64(
	long nelem, /* the chemical element, 1 for He */
	double tau, /* the radiative lifetime of the initial level.	*/
	double target_charge,
	long n,
	long l,
	double gHi,
	long Collider)
{
	/* >>refer	H-like	l-mixing	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
	/* >>refer	He-like	l-mixing	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
	double cs, Dnl,
		TwoLogDebye, TwoLogRc1, 
		factor1, factor2, 
		bestfactor,	factorpart,
		reduced_mass, reduced_mass_2_emass,
		rate;
	const double ChargIncoming = ColliderCharge[Collider];

	DEBUG_ENTRY( "CS_l_mixing_PS64()" );

	/* In this routine, two different cutoff radii are calculated, and as per PS64,
	 * the least of these is selected.  We take the least positive result.	*/

	/* This reduced mass is in grams.	*/
	reduced_mass = dense.AtomicWeight[nelem]*ColliderMass[Collider]/
		(dense.AtomicWeight[nelem]+ColliderMass[Collider])*ATOMIC_MASS_UNIT;
	/* this mass always appears relative to the electron mass, so define it that way */
	reduced_mass_2_emass = reduced_mass / ELECTRON_MASS;

	/* equation 46 of PS64 */
	/* min on density added to prevent this from becoming large and negative
	 * at very high n_e - Pengelly & Seaton did not apply this above 1e12 cm^-3 */
	/* This is 2 times the log of the Debye radius.	*/
	//TwoLogDebye = 1.68 + log10( phycon.te / MIN2(1e11 , dense.eden ) );
	/* Brocklehurst (1971, equation 3.40) has 1.181 instead of 1.68.  This appears to be due to 
	 * Pengelly and Seaton neglecting one of the two factors of PI in their Equation 33 */
	TwoLogDebye = 1.181 + log10( phycon.te / MIN2(1e10 , dense.eden ) );

	/* This is 2 times the log of cutoff = 0.72v(tau), where tau is the lifetime of the level nl.	
	 * This is PS64 equation 45 (same as Brocklehurst 1971 equation 3.41)  */
	TwoLogRc1 = 10.95 + log10( phycon.te * tau * tau / reduced_mass_2_emass );

#if 0
	/* This is 2 times the log of cutoff = 1.12 * hbar v / deltaE, where v is reduced velocity = sqrt( 2kT/mu ). */	
	/* This is PS64 equation 29 */
	if( deltaE > 0. )
		TwoLogRc2 = 2. * log10( 1.12 * H_BAR * v / deltaE );
	else
		TwoLogRc2 = BIGDOUBLE;
#endif

	/* this is equation 44 of PS64 */
	Dnl = POW2( ChargIncoming / target_charge) * 6. * POW2( (double)n) *
		( POW2((double)n) - POW2((double)l) - l - 1);

	ASSERT( Dnl > 0. );
	ASSERT( phycon.te  / Dnl / reduced_mass_2_emass > 0. );

	factorpart = (11.54 + log10( phycon.te  / Dnl / reduced_mass_2_emass ) );

	if( (factor1 = factorpart + TwoLogDebye) <= 0.)
		factor1 = BIGDOUBLE;
	if( (factor2 = factorpart + TwoLogRc1) <= 0.)
		factor2 = BIGDOUBLE;

	/* Now we find the least positive result.	*/
	bestfactor = MIN2(factor1,factor2);

	ASSERT( bestfactor > 0. );

	/* If both factors are bigger than 100, toss out the result, and return SMALLFLOAT. */
	if( bestfactor > 100. )
	{
		return SMALLFLOAT;
	}

	/* This is the rate coefficient.   Units: cm^3 s-1	*/
	rate = 9.93e-6 * sqrt( reduced_mass_2_emass  ) * Dnl / phycon.sqrte * bestfactor;

	/***** NB NB NB NB 
	 * Brocklehurst (1971) has a factor of electron density in the rate coefficient (equation 3.38).  
	 * This is NOT a proper rate, as can be seen in his equations 2.2 and 2.4.  This differs
	 * from the formulism given by PS64. */
	//rate *= dense.eden;

	/* this is the TOTAL rate from nl to nl+/-1. So assume we can
	 * divide by two to get the rate nl to either nl+1 or nl-1. */
	if( l > 0 )
		rate /=2;

	/* convert rate to collision strength */
	/* NB - the term in parentheses corrects for the fact that COLL_CONST is only appropriate 
	 * for electron colliders and is off by reduced_mass_2_emass^-1.5 */
	cs = rate / ( COLL_CONST * pow(reduced_mass_2_emass, -1.5) ) *
		phycon.sqrte * gHi;

	ASSERT( cs > 0. );

	return cs;
}

/*CS_l_mixing - find collision strength for l-mixing collisions for neutrals */
/* The VF stands for Vrinceanu & Flannery 2001 */
/* >>refer	He	l-mixing	Vrinceanu, D. & Flannery, M. R. 2001, PhysRevA 63, 032701	*/
/* >>refer	He	l-mixing	Hezel, T. P., Burkhardt, C. E., Ciocca, M., He, L-W., */
/* >>refercon	Leventhal, J. J. 1992, Am. J. Phys. 60, 329 */
double CS_l_mixing_VF01(long int ipISO,
						long int nelem,
						long int n,
						long int l,
						long int lp,
						long int s,
						double temp,
						long int Collider )
{

	double coll_str;

	DEBUG_ENTRY( "CS_l_mixing_VF01()" );

	my_Integrand_VF01_E func;
	func.ipISO = ipISO;
	func.nelem = nelem;
	func.n = n;
	func.l = l;
	func.lp = lp;
	func.s = s;
	func.Collider = Collider;
	func.temp = temp;
	func.ColliderCharge = ColliderCharge[Collider];
	func.lgParamIsRedVel = false;
	ASSERT( func.ColliderCharge > 0. );

	Integrator<my_Integrand_VF01_E, Gaussian32> VF01_E;

	/* no need to do this for h-like */
	if( ipISO > ipH_LIKE )
	{
		ASSERT( l != 0 );
		ASSERT( lp != 0 );
	}

	if( !iso.lgCS_therm_ave[ipISO] )
	{
		/* Must do some thermal averaging for densities greater
		 * than about 10000 and less than about 1e10,
		 * because kT gives significantly different results.
		 * Still, do sparser integration than is done below */
		if( (dense.eden > 10000.) && (dense.eden < 1E10 ) )
		{
			coll_str = VF01_E.sum( 0.0, 6.0, func );
		}
		else
		{
			/* Do NOT average over Maxwellian */
			coll_str = collision_strength_VF01( ipISO, nelem, n, l, lp, s, Collider,
				ColliderCharge[Collider], temp, temp/TE1RYD, false );
		}
	}
	else
	{
		/* DO average over Maxwellian */
		coll_str = VF01_E.sum( 0.0, 1.0, func );
		coll_str += VF01_E.sum( 1.0, 10.0, func );
	}

	return coll_str;
}

STATIC double collision_strength_VF01( long ipISO, long nelem, long n, long l, long lp, long s, long Collider, double ColliderCharge, double temp, double velOrEner, bool lgParamIsRedVel )
{
	double cross_section, coll_str, RMSv, aveRadius, reduced_vel, E_Proj_Ryd;
	double ConstantFactors, reduced_mass, CSIntegral, stat_weight;
	double quantum_defect1, quantum_defect2, omega_qd1, omega_qd2, omega_qd;
	double reduced_b_max, reduced_b_min, alphamax, alphamin, step, alpha1 /*, alpha2*/;
	long ipLo, ipHi;

	DEBUG_ENTRY( "collision_strength_VF01()" );

	ASSERT( n > 0 );

	/* >>chng 06 may 30, move these up from below.  Also ipHi needs lp not l. */
	ipLo = iso.QuantumNumbers2Index[ipISO][nelem][n][l][s];
	ipHi = iso.QuantumNumbers2Index[ipISO][nelem][n][lp][s];
	stat_weight = StatesElemNEW[nelem][nelem-ipISO][ipLo].g;

	/* no need to do this for h-like */
	if( ipISO > ipH_LIKE )
	{
		/* these shut up the lint, already done above */
		ASSERT( l > 0 );
		ASSERT( lp > 0 );
	}

	/* This reduced mass is in grams.	*/
	reduced_mass = dense.AtomicWeight[nelem]*ColliderMass[Collider]/
		(dense.AtomicWeight[nelem]+ColliderMass[Collider])*ATOMIC_MASS_UNIT;
	ASSERT( reduced_mass > 0. );

	/* this is root mean squared velocity */
	/* use this as projectile velocity for thermally-averaged cross-section */
	aveRadius = (BOHR_RADIUS_CM/((double)nelem+1.-(double)ipISO))*POW2((double)n);
	ASSERT( aveRadius < 1.e-4 );
	/* >>chng 05 jul 14, as per exchange with Ryan Porter & Peter van Hoof, avoid
	 * roundoff error and give ability to go beyond zinc */
	/*ASSERT( aveRadius >=  BOHR_RADIUS_CM );*/
	ASSERT( aveRadius > 3.9/LIMELM * BOHR_RADIUS_CM );

	/* vn = n * H_BAR/ m / r = Z * e^2 / n / H_BAR 
	 * where Z is the effective charge. */
	RMSv = ((double)nelem+1.-(double)ipISO)*POW2(ELEM_CHARGE_ESU)/(double)n/H_BAR;
	ASSERT( RMSv > 0. );

	ASSERT( ColliderMass[Collider] > 0. );

	if( lgParamIsRedVel )
	{
		/* velOrEner is a reduced velocity */
		reduced_vel = velOrEner;
		/* The proton mass is needed here because the ColliderMass array is
		 * expressed in units of the proton mass, but here we need absolute mass. */
		E_Proj_Ryd = 0.5 * POW2( reduced_vel * RMSv ) * ColliderMass[Collider] *
			PROTON_MASS / EN1RYD;
	}
	else
	{	
		/* velOrEner is a projectile energy in Rydbergs */
		E_Proj_Ryd = velOrEner;
		reduced_vel = sqrt( 2.*E_Proj_Ryd*EN1RYD/ColliderMass[Collider]/PROTON_MASS )/RMSv;
	}

	/* put limits on the reduced velocity.   These limits should be more than fair. */
	ASSERT( reduced_vel > 1.e-10 );
	ASSERT( reduced_vel < 1.e10 );

	/* Factors outside integral	*/
	ConstantFactors = 4.5*PI*POW2(ColliderCharge*aveRadius/reduced_vel);

	/* Reduced here means in units of aveRadius: */
	reduced_b_min = 1.5 * ColliderCharge / reduced_vel;
	ASSERT( reduced_b_min > 1.e-10 );
	ASSERT( reduced_b_min < 1.e10 );

	if( ipISO == ipH_LIKE )
	{
		/* Debye radius: appears to be too large, results in 1/v^2 variation. */
		reduced_b_max = sqrt( BOLTZMANN*temp/ColliderCharge/dense.eden ) 
			/ (PI2*ELEM_CHARGE_ESU)/aveRadius;
	}
	else if( ipISO == ipHE_LIKE )
	{
		quantum_defect1  = (double)n- (double)nelem/sqrt(iso.xIsoLevNIonRyd[ipISO][nelem][ipLo]);
		quantum_defect2  = (double)n- (double)nelem/sqrt(iso.xIsoLevNIonRyd[ipISO][nelem][ipHi]);

		/* The magnitude of each quantum defect must be between zero and one. */
		ASSERT( fabs(quantum_defect1)  < 1.0 );
		ASSERT( fabs(quantum_defect1)  > 0.0 );
		ASSERT( fabs(quantum_defect2)  < 1.0 );
		ASSERT( fabs(quantum_defect2)  > 0.0 );

		/* The quantum defect precession frequencies */
		omega_qd1 = fabs( 5.* quantum_defect1 * (1.-0.6*POW2((double)l/(double)n)) / POW3( (double)n ) / (double)l );
		/* >>chng 06 may 30, this needs lp not l. */
		omega_qd2 = fabs( 5.* quantum_defect2 * (1.-0.6*POW2((double)lp/(double)n)) / POW3( (double)n ) / (double)lp );
		/* Take the average for the two levels, for reciprocity. */
		omega_qd = 0.5*( omega_qd1 + omega_qd2 );

		ASSERT( omega_qd > 0. );

		reduced_b_max = sqrt( 1.5 * ColliderCharge * n / omega_qd )/aveRadius;
	}
	else
		/* rethink this before using on other iso sequences. */
		TotalInsanity();

	reduced_b_max = MAX2( reduced_b_max, reduced_b_min );
	ASSERT( reduced_b_max > 0. );

	// set up the integrator.
	my_Integrand_VF01_alpha func;
	func.n = n;
	func.l = l;
	func.lp = lp;
	func.red_vel = reduced_vel;
	func.an = aveRadius;
	func.ColliderCharge = ColliderCharge;
	func.bmax = reduced_b_max*aveRadius;
	Integrator<my_Integrand_VF01_alpha, Gaussian32> VF01_alpha;

	alphamin = 1.5*ColliderCharge/(reduced_vel * reduced_b_max);
	alphamax = 1.5*ColliderCharge/(reduced_vel * reduced_b_min);

	ASSERT( alphamin > 0. );
	ASSERT( alphamax > 0. );

	alphamin = MAX2( alphamin, 1.e-30 );
	alphamax = MAX2( alphamax, 1.e-20 );

	CSIntegral = 0.;

	if( alphamax > alphamin )
	{

		step = (alphamax - alphamin)/5.;
		alpha1 = alphamin;
		CSIntegral += VF01_alpha.sum( alpha1, alpha1+step, func );
		CSIntegral += VF01_alpha.sum( alpha1+step, alpha1+4.*step, func );
	}

	/* Calculate cross section */
	cross_section = ConstantFactors * CSIntegral;

	/* convert to collision strength, cross section is already in cm^2 */
	coll_str = ConvCrossSect2CollStr( cross_section, stat_weight, E_Proj_Ryd, reduced_mass );

	coll_str = MAX2( SMALLFLOAT, coll_str);

	return coll_str;
}

STATIC double L_mix_integrand_VF01( long n, long l, long lp, double bmax, double red_vel, double an, double ColliderCharge, double alpha )
{
	double integrand, deltaPhi, b;

	DEBUG_ENTRY( "L_mix_integrand_VF01()" );

	ASSERT( alpha >= 1.e-30 );
	ASSERT( bmax > 0. );
	ASSERT( red_vel > 0. );

	/* >>refer He l-mixing Kazansky, A. K. & Ostrovsky, V. N. 1996, JPhysB: At. Mol. Opt. Phys. 29, 3651 */
	b = 1.5*ColliderCharge*an/(red_vel * alpha);
	/* deltaPhi is the effective angle swept out by the projector as viewed by the target.  */
	if( b < bmax )
	{
		deltaPhi = -1.*PI + 2.*asin(b/bmax);
	}
	else
	{
		deltaPhi = 0.;
	}
	integrand = 1./alpha/alpha/alpha;
	integrand *= StarkCollTransProb_VF01( n, l, lp, alpha, deltaPhi );

	return integrand;
}

STATIC double StarkCollTransProb_VF01( long n, long l, long lp, double alpha, double deltaPhi)
{
	double probability;
	double cosU1, cosU2, sinU1, sinU2, cosChiOver2, sinChiOver2, cosChi, A, B;

	DEBUG_ENTRY( "StarkCollTransProb_VF01()" );

	ASSERT( alpha > 0. );

	/* These are defined on page 11 of VF01 */ 
	cosU1 = 2.*POW2((double)l/(double)n) - 1.;
	cosU2 = 2.*POW2((double)lp/(double)n) - 1.;

	sinU1 = sqrt( 1. - cosU1*cosU1 );
	sinU2 = sqrt( 1. - cosU2*cosU2 );


	cosChiOver2 = (1. + alpha*alpha*cos( sqrt(1.+alpha*alpha) * deltaPhi ) )/(1.+alpha*alpha);
	sinChiOver2 = sqrt( 1. - cosChiOver2*cosChiOver2 );
	cosChi = 2. * POW2( cosChiOver2 ) - 1.;

	if( l == 0 )
	{
		if( -1.*cosU2 - cosChi < 0. )
		{
			probability = 0.;
		}
		else
		{
			/* Here is the initial state S case */
			ASSERT( sinChiOver2 > 0. );
			ASSERT( sinChiOver2*sinChiOver2 > POW2((double)lp/(double)n) );
			/* This is equation 35 of VF01.  There is a factor of hbar missing in the denominator of the
			 * paper, but it's okay if you use atomic units (hbar = 1). */
			probability = (double)lp/(POW2((double)n)*sinChiOver2*sqrt( POW2(sinChiOver2) - POW2((double)lp/(double)n) ) );
		}
	}
	else
	{
		double OneMinusCosChi = 1. - cosChi;

		if( OneMinusCosChi == 0. )
		{
			double hold = sin( deltaPhi / 2. );
			/* From approximation at bottom of page 10 of VF01. */
			OneMinusCosChi = 8. * alpha * alpha * POW2( hold );
		}

		if( OneMinusCosChi == 0. )
		{
			probability = 0.;
		}
		else
		{
			/* Here is the general case */
			A = (cosU1*cosU2 - sinU1*sinU2 - cosChi)/OneMinusCosChi;
			B = (cosU1*cosU2 + sinU1*sinU2 - cosChi)/OneMinusCosChi;

			ASSERT( B > A );

			/* These are the three cases of equation 34. */
			if( B <= 0. )
			{
				probability = 0.;
			}
			else
			{
				ASSERT( POW2( sinChiOver2 ) > 0. );

				probability = 2.*lp/(PI* /*H_BAR* */ n*n*POW2( sinChiOver2 ));

				if( A < 0. )
				{
					probability *= ellpk( -A/(B-A) );
					probability /= sqrt( B - A );
				}
				else
				{
					probability *= ellpk( A/B );
					probability /= sqrt( B );
				}
			}
		}

	}

	return probability;

}
