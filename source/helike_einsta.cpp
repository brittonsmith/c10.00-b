/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*he_1trans compute Aul for given line	*/
/*ritoa - converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
/*ForbiddenAuls calculates transition probabilities for forbidden transitions.	*/
/*scqdri - stands for Semi-Classical Quantum Defect Radial Integral	*/
/*Jint - used by scqdri	*/
/*AngerJ - used by scqdri */
/*DoFSMixing - applies a fine structure mixing approximation to A's.  To be replaced by 
 * method that treats the entire rate matrix.	*/

#include "cddefines.h" 
#include "physconst.h" 
#include "taulines.h"
#include "dense.h"
#include "trace.h"
#include "hydro_bauman.h"
#include "iso.h"
#include "helike.h"
#include "helike_einsta.h"
#include "hydroeinsta.h"

/* the array of transitions probabilities read from data file.  */
static double ***TransProbs;

/*ritoa converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
STATIC double ritoa( long li, long lf, long nelem, double k, double RI2 );

/* handles all forbidden transition probabilities. */
STATIC double ForbiddenAuls( long ipHi, long ipLo, long nelem );

/*
static long RI_ipHi, RI_ipLo, RI_ipLev;
static long globalZ;
*/

/* used as parameters in qg32 integration */
static double vJint , zJint;

/* the integrand used in the AngerJ function below. */
STATIC double Jint( double theta )
{
	/*	[ cos[vx - z sin[x]] ]  */
	double 
		d0 = ( 1.0 / PI ),
		d1 = vJint * theta,
		d2 = zJint * sin(theta),
		d3 = (d1 - d2),
		d4 = cos(d3),
		d5 = (d0 * d4);

	return( d5 );
}

/* AngerJ function. */
STATIC double AngerJ( double vv, double zz )
{
	long int rep = 0, ddiv, divsor;

	double y = 0.0;

	DEBUG_ENTRY( "AngerJ()" );

	/* Estimate number of peaks in integrand.  */                            
	/* Divide region of integration by number  */
	/*  peaks in region.                       */
	if( (fabs(vv)) - (int)(fabs(vv)) > 0.5 )
		ddiv = (int)(fabs(vv)) + 1;
	else 
		ddiv = (int)(fabs(vv));

	divsor  = ((ddiv == 0) ? 1 : ddiv);
	vJint = vv;
	zJint = zz;

	for( rep = 0; rep < divsor; rep++ )
	{
		double
		rl = (((double) rep)/((double) divsor)),
		ru = (((double) (rep+1))/((double) divsor)),
		x_low = (PI * rl),
		x_up  = (PI * ru);      

		y += qg32( x_low, x_up, Jint );
	}

	return( y );
}

/******************************************************************************/
/******************************************************************************/
/*                                                                            */
/*    Semi-Classical Quantum Defect Radial Integral                           */      
/*                                                                            */
/*   See for example                                                          */
/*     Atomic, Molecular & Optical Physics Handbook                           */
/*     Gordon W. F. Drake; Editor                                             */
/*     AIP Press                                                              */
/*     Woddbury, New York.                                                    */
/*     1996                                                                   */
/*                                                                            */
/* NOTE:: we do not include the Bohr Radius a_o in the                        */
/*           definition of of  R(n,L;n'L') as per Drake.                      */
/*                                                                            */
/*                                                                            */
/*                   1  (n_c)^2 | {      D_l max(L,L') }                      */
/*    R(n,L;n'L') = --- ------- | { 1 -  ------------- } J( D_n-1; -x )  -    */
/*                   Z   2 D_n  | {          n_c       }                      */
/*                                                                            */
/*                                                                            */
/*                    {      D_L max(L,L') }                                  */
/*                -   { 1 +  ------------- } J( D_n+1; -x )                   */
/*                    {          n_c       }                                  */
/*                                                                            */
/*                                                                            */
/*                     2                    |                                 */
/*                  + ---  sin(Pi D_n)(1-e) |                                 */
/*                     Pi                   |                                 */
/*                                                                            */
/*  where                                                                     */
/*        n_c = (2n*'n*)/(n*'+n*)                                             */
/*                                                                            */
/*       Here is the quantity Drake gives...                                  */
/*            n_c = ( 2.0 * nstar * npstar ) / ( nstar + npstar );            */
/*                                                                            */
/*       while V.A. Davidkin uses                                             */
/*            n_c = sqrt(  nstar * npstar  );                                 */
/*                                                                            */
/*        D_n = n*' - n*                                                      */
/*                                                                            */      
/*        D_L = L*' - L*                                                      */
/*                                                                            */
/*        x = e D_n                                                           */
/*                                                                            */
/*        Lmx  = max(L',L)                                                    */
/*                                                                            */
/*        e = sqrt( 1 - {Lmx/n_c}^2 )                                         */
/*                                                                            */
/*                                                                            */
/*       Here n* = n - qd where qd is the quantum defect                      */
/*                                                                            */
/******************************************************************************/
/******************************************************************************/
double scqdri(/* upper and lower quantum numbers...n's are effective	*/
			  double nstar, long int l,
			  double npstar, long int lp,
              double iz
              )
{
	double n_c = ((2.0 * nstar * npstar ) / ( nstar + npstar ));

	double D_n = (nstar - npstar);
	double D_l = (double) ( l - lp );
	double lg  = (double) ( (lp > l) ? lp : l);

	double h = (lg/n_c);
	double g = h*h;
	double f = ( 1.0 - g );
	double e = (( f >= 0.0) ? sqrt( f ) : 0.0 );

	double x  = (e * D_n);
	double z  = (-1.0 * x);
	double v1 = (D_n + 1.0);
	double v2 = (D_n - 1.0); 

	double d1,d2,d7,d8,d9,d34,d56,d6_1; 

	DEBUG_ENTRY( "scqdri()" );

	if( iz == 0.0 )
		iz += 1.0;

	if( D_n == 0.0 )
	{
		return( -1.0 );
	}

	if( D_n < 0.0 )
	{
		return( -1.0 );
	}

	if( f < 0.0 )
	{
		/* This can happen for certain  quantum defects   */
		/* in the lower n=1:l=0 state. In which case you  */
		/* probably should be using some other alogrithm  */
		/* or theory to calculate the dipole moment.      */
		return( -1.0 );
	}

	d1 = ( 1.0 / iz );

	d2 = (n_c * n_c)/(2.0 * D_n);

	d34 = (1.0 - ((D_l * lg)/n_c)) * AngerJ( v1, z );

	d56 = (1.0 + ((D_l * lg)/n_c)) * AngerJ( v2, z );

	d6_1 = PI * D_n;

	d7 = (2./PI) * sin( d6_1 ) * (1.0 - e);

	d8 = d1 * d2 * ( (d34) - (d56) + d7 );

	d9 = d8 * d8;

	ASSERT( D_n  > 0.0 );
	ASSERT( l  >= 0  );
	ASSERT( lp >= 0 );
	ASSERT( (l == lp + 1) || ( l == lp - 1) );
	ASSERT( n_c != 0.0 );
	ASSERT( f >= 0.0 );
	ASSERT( d9  > 0.0 );

	return( d9 );
}

STATIC double ForbiddenAuls( long ipHi, long ipLo, long nelem )
{
	double A;
	/* >>refer	Helike	2pho	Derevianko, A., & Johnson, W.R. 1997, Phys. Rev. A 56, 1288
	 * numbers are not explicitly given in this paper for Z=21-24,26,27,and 29.
	 * So numbers given here are interpolated.	*/
	double As2nuFrom1S[28] = {1940.,1.82E+04,9.21E+04,3.30E+05,9.44E+05,2.31E+06,5.03E+06,1.00E+07,
		1.86E+07,3.25E+07,5.42E+07,8.69E+07,1.34E+08,2.02E+08,2.96E+08,4.23E+08,5.93E+08,8.16E+08,
		1.08E+09,1.43E+09,1.88E+09,2.43E+09,3.25E+09,3.95E+09,4.96E+09,6.52E+09,7.62E+09,9.94E+09};
	/* Important clarification, according to Derevianko & Johnson (see ref above), 2^3S can decay
	 * to ground in one of two ways: through a two-photon process, or through a single-photon M1 decay,
	 * but the M1 rates are about 10^4 greater that the two-photon decays throughout the entire
	 * sequence.  Thus these numbers, are much weaker than the effective decay rate, but should probably
	 * be treated in as a two-photon decay at some point	*/
	double As2nuFrom3S[28] = {1.25E-06,5.53E-05,8.93E-04,8.05E-03,4.95E-02,2.33E-01,8.94E-01,2.95E+00,
		8.59E+00,2.26E+01,5.49E+01,1.24E+02,2.64E+02,5.33E+02,1.03E+03,1.91E+03,3.41E+03,5.91E+03,
		9.20E+03,1.50E+04,2.39E+04,3.72E+04,6.27E+04,8.57E+04,1.27E+05,2.04E+05,2.66E+05,4.17E+05};

	DEBUG_ENTRY( "ForbiddenAuls()" );

	int ipISO = ipHE_LIKE;

	if( (ipLo == ipHe1s1S) && (N_(ipHi) == 2) )
	{
		if( nelem == ipHELIUM )
		{
			/* All of these but the second and third one (values 51.02 and 1E-20) are from
			 * >>refer	HeI	As	Lach, G, & Pachucki, K, 2001, Phys. Rev. A 64, 042510
			,* 1E-20 is made up
			 * 51.3 is from the Derevianko & Johnson paper cited above.	*/
			double ForbiddenHe[5] = { 1.272E-4,	51.02,	1E-20,	177.58,	0.327 };

			fixit(); // adding the 2-photon decay 2^3S - 1^1S may be important in early universe
			A = ForbiddenHe[ipHi - 1];
			iso_put_error(ipHE_LIKE,nelem,ipHe2s3S ,ipHe1s1S,IPRAD, 0.01f, 0.20f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2s1S ,ipHe1s1S,IPRAD, 0.01f, 0.05f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2p3P0,ipHe1s1S,IPRAD, 0.00f, 0.00f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2p3P1,ipHe1s1S,IPRAD, 0.01f, 0.05f);
			iso_put_error(ipHE_LIKE,nelem,ipHe2p3P2,ipHe1s1S,IPRAD, 0.01f, 0.01f);
		}
		else
		{
			switch ( (int)ipHi )
			{
			case 1: /* Parameters for 2^3S to ground transition.	*/
				/* >>refer	Helike	As	Lin, C.D., Johnson, W.R., and Dalgarno, A. 1977, 
				 * >>refercon	Phys. Rev. A 15, 1, 010015	*/
				A = (3.9061E-7) * pow( (double)nelem+1., 10.419 ) + As2nuFrom3S[nelem-2];
				break;
			case 2: /* Parameters for 2^1S to ground transition.	*/
				A = As2nuFrom1S[nelem-2];
				break;
			case 3: /* Parameters for 2^3P0 to ground transition.	*/
				A = iso.SmallA;
				break;
			case 4: /* Parameters for 2^3P1 to ground transition.	*/
				/* >>chng 06 jul 25, only use the fit to Johnson et al. values up to and
				 * including argon, where there calculation stops.  For higher Z use below */
				if( nelem <= ipARGON )
				{
					A = ( 11.431 * pow((double)nelem, 9.091) );
				}
				else
				{
					/* a fit to the Lin et al. 1977. values, which go past zinc. */
					A = ( 383.42 * pow((double)nelem, 7.8901) );
				}
				break;
			case 5: /* Parameters for 2^3P2 to ground transition.	*/
				/* fit to Lin et al. 1977 values.  This fit is good
				 * to 7% for the range from carbon to iron. The Lin et al. values
				 * differ from the Hata and Grant (1984) values (only given from 
				 * calcium to copper) by less than 2%. */
				A = ( 0.11012 * pow((double)nelem, 7.6954) ); 
				break;
			default:
				TotalInsanity();
			}
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
		}

		return A;
	}

	/* The next two cases are fits to probabilities given in 
	 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
	 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
	/* originally astro.ph. 0201454 */
	/* The involve Triplet P and Singlet S.  Rates for Triplet S to Singlet P 
	 * do not seem to be available.	*/

	/* Triplet P to Singlet S...Delta n not equal to zero!	*/
	else if( nelem>ipHELIUM && L_(ipHi)==1 && S_(ipHi)==3 && 
		L_(ipLo)==0 && S_(ipLo)==1 && N_(ipLo) < N_(ipHi) )
	{
		A = 8.0E-3 * exp(9.283/sqrt((double)N_(ipLo))) * pow((double)nelem,9.091) /
			pow((double)N_(ipHi),2.877);
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	/* Singlet S to Triplet P...Delta n not equal to zero!	*/
	else if( nelem > ipHELIUM && L_(ipHi)==0 && S_(ipHi)==1 && 
		L_(ipLo)==1 && S_(ipLo)==3 && N_(ipLo) < N_(ipHi) )
	{
		A = 2.416 * exp(-0.256*N_(ipLo)) * pow((double)nelem,9.159) / pow((double)N_(ipHi),3.111);

		if( ( (ipLo == ipHe2p3P0) || (ipLo == ipHe2p3P2) ) )
		{
			/* This is divided by 3 instead of 9, because value calculated is specifically for 2^3P1.
			 * Here we assume statistical population of the other two.	*/
			A *= (2.*(ipLo-3)+1.0)/3.0;
		}
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else if( ( ipLo == ipHe2s3S ) && ( ipHi == ipHe3s3S ) )
	{
		double As_3TripS_to_2TripS[9] = {
			7.86E-9, 4.59E-6, 1.90E-4, 2.76E-3, 2.25E-2,
			1.27E-1, 5.56E-1, 2.01E0, 6.26E0 };

		/* These M1 transitions given by 
		 * >>refer He-like As Savukov, I.M., Labzowsky, and Johnson, W.R. 2005, PhysRevA, 72, 012504 */
		if( nelem <= ipNEON )
		{
			A = As_3TripS_to_2TripS[nelem-1];
			/* 20% error is based on difference between Savukov, Labzowsky, and Johnson (2005)
			 * and Lach and Pachucki (2001) for the helium transition. */
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.2f,0.2f);
		}
		else
		{
			/* This is an extrapolation to the data given above.  The fit reproduces
			 * the above values to 10% or better. */
			A= 7.22E-9*pow((double)nelem, 9.33);
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.3f,0.3f);
		}
	}

	else if( ( ipLo == ipHe2s3S ) && ( ipHi == ipHe2p1P ) )
	{
		/* This transition,1.549 , given by Lach and Pachucki, 2001 for the atom */
		if( nelem == ipHELIUM )
		{
			A = 1.549;
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
		}
		else
		{
			/* This is a fit to data given in
			 * >>refer	He-like	As	Savukov, I.M., Johnson, W.R., & Safronova, U.I. 
			 * >>refercon	astro-ph 0205163	*/
			A= 0.1834*pow((double)nelem, 6.5735);
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
		}
	}

	else if( nelem==ipHELIUM && ipHi==4 && ipLo==3 )
	{
		/* >>refer	He	As	Bulatov, N.N. 1976, Soviet Astronomy, 20, 315 */
		fixit();
		/* This is the 29.6 GHz line that can be seen in radio galaxies. */
		/** \todo	2	find a transition probability for this 2^3P0 - 2^3P1 transition.
		 * It will require a bit of trickery to insert into the rate matrix, 
		 * because of the fact that the lower level has a higher index.  
		 * See discussion "Energy order within 2 3P" near the top of helike.c */
		A = 5.78E-12;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);

	}

	else if( nelem==ipHELIUM && ipHi==5 && ipLo==4 )
	{
		fixit();
		/* This is the 3 GHz line that can be seen in radio galaxies. */
		/** \todo	2	find a transition probability for this 2^3P1 - 2^3P2 transition.
		 * It will require a bit of trickery to insert into the rate matrix, 
		 * because of the fact that the lower level has a higher index.  
		 * See discussion "Energy order within 2 3P" near the top of helike.c */
		A = 3.61E-15;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	else
	{
		/* Current transition is not supported.	*/
		A = iso.SmallA;
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
	}

	ASSERT( A > 0.);

	return A;   
}

/* Calculates Einstein A for a given transition.	*/
double he_1trans( 
			   /* charge on c scale, Energy is wavenumbers, Einstein A	*/
			   long nelem , double Enerwn ,
			   /* quantum numbers of upper level:	*/
			   double Eff_nupper, long lHi, long sHi, long jHi,
			   /* and of lower level: */
			   double Eff_nlower, long lLo, long sLo, long jLo )
			   /* Note j is only necessary for 2 triplet P...for all other n,l,s,
			    * j is completely ignored.	*/
{
	int ipISO = ipHE_LIKE;
	double RI2, Aul;
	long nHi, nLo, ipHi, ipLo;

	DEBUG_ENTRY( "he_1trans()" );

	ASSERT(nelem > ipHYDROGEN);

	/* Since 0.4 is bigger than any defect, adding that to the effective principle quantum number,
	 * and truncating to an integer will produce the principal quantum number.	*/
	nHi = (int)(Eff_nupper + 0.4);
	nLo = (int)(Eff_nlower + 0.4);

	/* Make sure this worked correctly.	*/
	ASSERT( fabs(Eff_nupper-(double)nHi) < 0.4 );
	ASSERT( fabs(Eff_nlower-(double)nLo) < 0.4 );

	ipHi = iso.QuantumNumbers2Index[ipHE_LIKE][nelem][nHi][lHi][sHi];
	if( (nHi==2) && (lHi==1) && (sHi==3) )
	{
		ASSERT( (jHi>=0) && (jHi<=2) );
		ipHi -= (2 - jHi);
	}

	ipLo = iso.QuantumNumbers2Index[ipHE_LIKE][nelem][nLo][lLo][sLo];
	if( (nLo==2) && (lLo==1) && (sLo==3) )
	{
		ASSERT( (jLo>=0) && (jLo<=2) );
		ipLo -= (2 - jLo);
	}

	ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].n == nHi );
	if( nHi <= iso.n_HighestResolved_max[ipHE_LIKE][nelem] )
	{
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].l == lHi );
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHi].S == sHi );
	}
	ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].n == nLo );
	if( nLo <= iso.n_HighestResolved_max[ipHE_LIKE][nelem] )
	{
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].l == lLo );
		ASSERT( StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].S == sLo );
	}

	/* First do allowed transitions	*/
	if( (sHi == sLo) && (abs((int)(lHi - lLo)) == 1) )
	{
		Aul = -2.;

		/* For clarity, let's do this in two separate chunks...one for helium, one for everything else.	*/
		if( nelem == ipHELIUM )
		{
			/* Retrieve transition probabilities for Helium.	*/
			/* >>refer He	As	Drake, G.W.F., Atomic, Molecular, and Optical Physics Handbook */
			if( ipHi <= MAX_TP_INDEX && N_(ipHi) <= iso.n_HighestResolved_max[ipHE_LIKE][ipHELIUM] )
			{
				/*Must not be accessed by collapsed levels!	*/
				ASSERT( ipHi < ( iso.numLevels_max[ipHE_LIKE][ipHELIUM] - iso.nCollapsed_max[ipHE_LIKE][ipHELIUM] ) );
				ASSERT( ipLo < ( iso.numLevels_max[ipHE_LIKE][ipHELIUM] - iso.nCollapsed_max[ipHE_LIKE][ipHELIUM] ) );
				ASSERT( ipHi > 2 );

				Aul = TransProbs[nelem][ipHi][ipLo];

				iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.0001f,0.002f);
			}

			if( Aul < 0. )
			{
				/* Here are the Lyman transitions.	*/
				if( ipLo == ipHe1s1S )
				{
					ASSERT( (lHi == 1) && (sHi == 1) );

					/* these fits calculated from Drake A's (1996) */
					if( nLo == 1 )
						Aul = (1.59208e10) / pow(Eff_nupper,3.0);
					ASSERT( Aul > 0.);
					iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.0002f,0.002f);
				}

				/* last resort for transitions involving significant defects, 
				 * except that highest lLo are excluded */
				else if( lHi>=2 && lLo>=2 && nHi>nLo )
				{
					Aul = H_Einstein_A(nHi ,lHi , nLo , lLo , nelem);
					ASSERT( Aul > 0.);
					iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.006f,0.04f);
				}
				/* These fits come from extrapolations of Drake's oscillator strengths
				 * out to the series limit.  We also use this method to obtain threshold
				 * photoionization cross-sections for the lower level of each transition here. */
				/* See these two references :
				 * >>refer	He	As	Hummer, D. G. \& Storey, P. J. 1998, MNRAS 297, 1073 
				 * >>refer	Seaton's Theorem	Seaton, M. J. 1958, MNRAS 118, 504 */
				else if( N_(ipHi)>10 && N_(ipLo)<=5 && lHi<=2 && lLo<=2 )
				{
					int paramSet=0;
					double emisOscStr, x, a, b, c;
					double extrapol_Params[2][4][4][3] = {
						/* these are for singlets */
						{	
							{	/* these are P to S */
								{	0.8267396	,	1.4837624	,	-0.4615955	},
								{	1.2738405	,	1.5841806	,	-0.3022984	},
								{	1.6128996	,	1.6842538	,	-0.2393057	},
								{	1.8855491	,	1.7709125	,	-0.2115213	}},
							{	/* these are S to P */
								{	-1.4293664	,	2.3294080	,	-0.0890470	},
								{	-0.3608082	,	2.3337636	,	-0.0712380	},
								{	0.3027974	,	2.3326252	,	-0.0579008	},
								{	0.7841193	,	2.3320138	,	-0.0497094	}},
							{	/* these are D to P */
								{	1.1341403	,	3.1702435	,	-0.2085843	},
								{	1.7915926	,	2.4942946	,	-0.2266493	},
								{	2.1979400	,	2.2785377	,	-0.1518743	},
								{	2.5018229	,	2.1925720	,	-0.1081966	}},
							{	/* these are P to D */
								{	0.0000000	,	0.0000000	,	0.0000000	},
								{	-2.6737396	,	2.9379143	,	-0.3805367	},
								{	-1.4380124	,	2.7756396	,	-0.2754625	},
								{	-0.6630196	,	2.6887253	,	-0.2216493	}},
						},
						/* these are for triplets */
						{	
							{	/* these are P to S */
								{	0.3075287	,	0.9087130	,	-1.0387207	},
								{	0.687069	,	1.1485864	,	-0.6627317	},
								{	0.9776064	,	1.3382024	,	-0.5331906	},
								{	1.2107725	,	1.4943721	,	-0.4779232	}},
							{	/* these are S to P */
								{	-1.3659605	,	2.3262253	,	-0.0306439	},
								{	-0.2899490	,	2.3279391	,	-0.0298695	},
								{	0.3678878	,	2.3266603	,	-0.0240021	},
								{	0.8427457	,	2.3249540	,	-0.0194091	}},
							{	/* these are D to P */
								{	1.3108281	,	2.8446367	,	-0.1649923	},
								{	1.8437692	,	2.2399326	,	-0.2583398	},
								{	2.1820792	,	2.0693762	,	-0.1864091	},
								{	2.4414052	,	2.0168255	,	-0.1426083	}},
							{	/* these are P to D */
								{	0.0000000	,	0.0000000	,	0.0000000	},
								{	-1.9219877	,	2.7689624	,	-0.2536072	},
								{	-0.7818065	,	2.6595150	,	-0.1895313	},
								{	-0.0665624	,	2.5955623	,	-0.1522616	}},
						}
					};

					if( lLo==0 )
					{
						paramSet = 0;
					}
					else if( lLo==1 && lHi==0 )
					{
						paramSet = 1;
					}
					else if( lLo==1 && lHi==2 )
					{
						paramSet = 2;
					}
					else if( lLo==2 )
					{
						paramSet = 3;
						ASSERT( lHi==1 );
					}

					ASSERT( (int)((sHi-1)/2) == 0 || (int)((sHi-1)/2) == 1 );
					a = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][0];
					b = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][1];
					c = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][2];
					ASSERT( Enerwn > 0. );
					x = log( iso.xIsoLevNIonRyd[ipHE_LIKE][nelem][ipLo]*RYD_INF/Enerwn );

					emisOscStr = exp(a+b*x+c*x*x)/pow(Eff_nupper,3.)*
						(2.*lLo+1)/(2.*lHi+1);

					Aul = TRANS_PROB_CONST*Enerwn*Enerwn*emisOscStr;

					if( (ipLo == ipHe2p3P0) || (ipLo == ipHe2p3P1) || (ipLo == ipHe2p3P2) )
					{
						Aul *= (2.*(ipLo-3)+1.0)/9.0;
					}

					ASSERT( Aul > 0. );
					iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.0002f,0.002f);
				}
				else
				{
					/* Calculate the radial integral from the quantum defects.	*/
					RI2 = scqdri(Eff_nupper,lHi,Eff_nlower,lLo,(double)(ipHELIUM));
					ASSERT( RI2 > 0. );
					/* Convert radial integral to Aul.	*/
					Aul = ritoa(lHi,lLo,ipHELIUM,Enerwn,RI2);
					/* radial integral routine does not recognize fine structure.
					 * Here we split 2^3P.	*/
					if( (ipLo == ipHe2p3P0) || (ipLo == ipHe2p3P1) || (ipLo == ipHe2p3P2) )
					{
						Aul *= (2.*(ipLo-3)+1.0)/9.0;
					}

					ASSERT( Aul > 0. );
					iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.01f,0.07f);
				}
			}
		}

		/* Heavier species	*/
		else
		{
			/* Retrieve transition probabilities for Helike ions.	*/
			/* >>refer He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
			 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J, originally astro.ph. 0201454 */
			if( ipHi <= MAX_TP_INDEX && N_(ipHi) <= iso.n_HighestResolved_max[ipHE_LIKE][nelem] )
			{
				/*Must not be accessed by collapsed levels!	*/
				ASSERT( ipHi < ( iso.numLevels_max[ipHE_LIKE][nelem] - iso.nCollapsed_max[ipHE_LIKE][nelem] ) );
				ASSERT( ipLo < ( iso.numLevels_max[ipHE_LIKE][nelem] - iso.nCollapsed_max[ipHE_LIKE][nelem] ) );
				ASSERT( ipHi > 2 );

				Aul = TransProbs[nelem][ipHi][ipLo];
				iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
			}

			if( Aul < 0. )
			{
				/* Do same-n transitions. */
				if( nLo==nHi && lHi<=2 && lLo<=2 )
				{
					/* These are 2p3Pj to 2s3S fits to (low Z) Porquet & Dubau (2000) & 
					 * (high Z) NIST Atomic Spectra Database.	*/
					if( ipLo == ipHe2s3S )
					{
						if(ipHi == ipHe2p3P0)
							Aul = 3.31E7 + 1.13E6 * pow((double)nelem+1.,1.76);
						else if(ipHi == ipHe2p3P1)
							Aul = 2.73E7 + 1.31E6 * pow((double)nelem+1.,1.76);
						else if(ipHi == ipHe2p3P2)
							Aul = 3.68E7 + 1.04E7 * exp(((double)nelem+1.)/5.29);
						else
						{
							fprintf(ioQQQ," PROBLEM Impossible value for ipHi in he_1trans.\n");
							TotalInsanity();
						}
					}

					/* These are 2p1P to 2s1S fits to data from TOPbase.	*/
					else if( ( ipLo == ipHe2s1S ) && ( ipHi == ipHe2p1P) )
					{
						Aul = 5.53E6 * exp( 0.171*(nelem+1.) );
					}

					else 
					{
						/* This case should only be entered if n > 2.  Those cases were done above.	*/
						ASSERT( nLo > 2 );

						/* Triplet P to triplet S, delta n = 0	*/
						if( (lHi == 1) && (sHi == 3) && (lLo == 0) && (sLo == 3))
						{
							Aul = 0.4 * 3.85E8 * pow((double)nelem,1.6)/pow((double)nHi,5.328);
						}
						/* Singlet P to singlet D, delta n = 0	*/
						else if( (lHi == 1) && (sHi == 1) && (lLo == 2) && (sLo == 1))
						{
							Aul = 1.95E4 * pow((double)nelem,1.6) / pow((double)nHi, 4.269);
						}
						/* Singlet P to singlet S, delta n = 0	*/
						else if( (lHi == 1) && (sHi == 1) && (lLo == 0) )
						{
							Aul = 6.646E7 * pow((double)nelem,1.5) / pow((double)nHi, 5.077);
						}
						else 
						{
							ASSERT( (lHi == 2) && (sHi == 3)  && (lLo == 1) );
							Aul = 3.9E6 * pow((double)nelem,1.6) / pow((double)nHi, 4.9);
							if( (lHi >2) || (lLo > 2) )
								Aul *= (lHi/2.);
							if(lLo > 2)
								Aul *= (1./9.);
						}
					}
					ASSERT( Aul > 0.);
				}

				/* assume transitions involving F and higher orbitals are hydrogenic.	*/
				else if( (nHi > nLo) && ((lHi > 2) || (lLo > 2)) )
				{
					Aul = H_Einstein_A(nHi ,lHi , nLo , lLo , nelem);
					ASSERT( Aul > 0.);
				}

				/* These transitions are of great importance, but the below radial integral 
				 * routine fails to achieve desirable accuracy, so these are fits as produced 
				 * from He A's for nupper through 9.  They are transitions to ground and 
				 * 2, 3, and 4 triplet S.	*/
				else if( ( ipLo == 0 ) || ( ipLo == ipHe2s1S ) || ( ipLo == ipHe2s3S ) 
					|| ( ipLo == ipHe3s3S ) || ( ipLo == ipHe4s3S ) )
				{
					/* Here are the Lyman transitions.	*/
					if( ipLo == 0 )
					{
						ASSERT( (lHi == 1) && (sHi == 1) );

						/* In theory, this Z dependence should be Z^4, but values from TOPbase 
						 * suggest 3.9 is a more accurate exponent.	Values from 
						 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
						 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
						/* originally astro.ph. 0201454  */
						Aul = 1.375E10 * pow((double)nelem, 3.9) / pow((double)nHi,3.1);
					}

					/* Here are the Balmer transitions.	*/
					else if( ipLo == ipHe2s1S )
					{
						ASSERT( (lHi == 1) && (sHi == 1) );

						/* More fits, as above. */
						Aul = 5.0e8 * pow((double)nelem,4.) / pow((double)nHi, 2.889);
					}

					/* Here are transitions down to triplet S	*/
					else
					{
						ASSERT( (lHi == 1) && (sHi == 3) );

						/* These fits also as above. */
						if( nLo == 2 )
							Aul = 1.5 * 3.405E8 * pow((double)nelem,4.) / pow((double)nHi, 2.883);
						else if( nLo == 3 )
							Aul = 2.5 * 4.613E7 * pow((double)nelem,4.) / pow((double)nHi, 2.672);
						else 
							Aul = 3.0 * 1.436E7 * pow((double)nelem,4.) / pow((double)nHi, 2.617);
					}

					ASSERT( Aul > 0.);
				}

				/* Every other allowed transition is calculated as follows.	*/
				else
				{
					/* Calculate the radial integral from the quantum defects.	*/
					RI2 = scqdri(Eff_nupper,lHi,Eff_nlower,lLo,(double)(nelem));
					/* Convert radial integral to Aul.	*/
					Aul = ritoa(lHi,lLo,nelem,Enerwn,RI2);
					/* radial integral routine does not recognize fine structure.
					 * Here we split 2^3P.	*/
					if( ( (ipLo == ipHe2p3P0) || (ipLo == ipHe2p3P1) || (ipLo == ipHe2p3P2) ) && (Aul > iso.SmallA) )
					{
						Aul *= (2.*(ipLo-3)+1.0)/9.0;
					}

				}
				iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.1f,0.1f);
			}
		}
	}

	/* Now do forbidden transitions from 2-1 ... */
	/* and those given by  
	 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
	 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
	/* originally astro.ph. 0201454  
	 * for heavy elements. These are triplet P to singlet S, 
	 * going either up or down...Triplet S to Singlet P are not included, as they are far weaker.	*/
	else 
	{
		ASSERT( (sHi != sLo) || (abs((int)(lHi - lLo)) != 1) );
		Aul = ForbiddenAuls(ipHi, ipLo, nelem);
		ASSERT( Aul > 0. );
	}

	Aul = MAX2( Aul, iso.SmallA );
	ASSERT( Aul >= iso.SmallA );

	/* negative energy for a transition with substantial transition probability
	 * would be major logical error - but is ok for same-n l transitions */
	if( Enerwn < 0. && Aul > iso.SmallA )
	{
		fprintf( ioQQQ," he_1trans hit negative energy, nelem=%li, val was %f \n", nelem ,Enerwn );
	}

	return Aul;
}

void DoFSMixing( long nelem, long ipLoSing, long ipHiSing )
{
	/* Every bit of this routine is based upon the singlet-triplet mixing formalism given in
	 * >>refer He	FSM	Drake, G. W. F. 1996, in Atomic, Molecular, \& Optical Physics Handbook,
	 * >>refercon	ed. G. W. F. Drake (New York: AIP Press).	
	 * That formalism mixes the levels themselves, but since this code is not J-resolved, we simulate
	 * that by mixing only the transition probabilities.  We find results comparable to those calculated
	 * in the fully J-resolved model spearheaded by Rob Bauman, and described in
	 * >>refer	He	FSM	Bauman, R. P., Porter, R. L., Ferland, G. J., \& MacAdam, K. B. 2005, ApJ, accepted */
	long int nHi, lHi, sHi, nLo, lLo, sLo, ipHiTrip, ipLoTrip;
	double Ass, Att, Ast, Ats;
	double SinHi, SinLo, CosHi, CosLo;
	double HiMixingAngle, LoMixingAngle , error;
	double Kss, Ktt, Kts, Kst, fss, ftt, fssNew, fttNew, ftsNew, fstNew, temp;

	DEBUG_ENTRY( "DoFSMixing()" );

	nHi = StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHiSing].n;
	lHi = StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHiSing].l;
	sHi = StatesElemNEW[nelem][nelem-ipHE_LIKE][ipHiSing].S;
	nLo = StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLoSing].n;
	lLo = StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLoSing].l;
	sLo = StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLoSing].S;

	if( ( sHi == 3 || sLo == 3 ) ||
	    ( abs(lHi - lLo) != 1 ) ||
	    ( nLo < 2 ) ||
	    ( lHi <= 1 || lLo <= 1 ) ||
	    ( nHi == nLo && lHi == 1 && lLo == 2 ) ||
	    ( nHi > nLo && lHi != 1 && lLo != 1 ) )
	{
		return;
	}

	ASSERT( lHi > 0 );
	/*ASSERT( (lHi > 1) && (lLo > 1) );*/

	ipHiTrip = iso.QuantumNumbers2Index[ipHE_LIKE][nelem][nHi][lHi][3];
	ipLoTrip = iso.QuantumNumbers2Index[ipHE_LIKE][nelem][nLo][lLo][3];

	if( lHi == 2 )
	{
		HiMixingAngle = 0.01;
	}
	else if( lHi == 3 )
	{
		HiMixingAngle = 0.5;
	}
	else
	{
		HiMixingAngle = PI/4.;
	}

	if( lLo == 2 )
	{
		LoMixingAngle = 0.01;
	}
	else if( lLo == 3 )
	{
		LoMixingAngle = 0.5;
	}
	else
	{
		LoMixingAngle = PI/4.;
	}

	/* These would not work correctly if l<=1 were included in this treatment!	*/
	ASSERT( ipHiTrip > ipLoTrip );
	ASSERT( ipHiTrip > ipLoSing );
	ASSERT( ipHiSing > ipLoTrip );
	ASSERT( ipHiSing > ipLoSing );

	SinHi = sin( HiMixingAngle );
	SinLo = sin( LoMixingAngle );
	CosHi = cos( HiMixingAngle );
	CosLo = cos( LoMixingAngle );

	Kss = Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].EnergyWN;
	Ktt = Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].EnergyWN;
	Kst = Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoTrip].EnergyWN;
	Kts = Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoSing].EnergyWN;

	fss = Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Emis->Aul/TRANS_PROB_CONST/POW2(Kss)/
		Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Lo->g*Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Hi->g;

	ftt = Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Emis->Aul/TRANS_PROB_CONST/POW2(Ktt)/
		Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Lo->g*Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Hi->g;

	temp = sqrt(fss/Kss)*CosHi*CosLo + sqrt(ftt/Ktt)*SinHi*SinLo;
	fssNew = Kss*POW2( temp );
	temp = sqrt(fss/Kss)*SinHi*SinLo + sqrt(ftt/Ktt)*CosHi*CosLo;
	fttNew = Ktt*POW2( temp );
	temp = sqrt(fss/Kss)*CosHi*SinLo - sqrt(ftt/Ktt)*SinHi*CosLo;
	fstNew = Kst*POW2( temp );
	temp = sqrt(fss/Kss)*SinHi*CosLo - sqrt(ftt/Ktt)*CosHi*SinLo;
	ftsNew = Kts*POW2( temp );

	Ass = TRANS_PROB_CONST*POW2(Kss)*fssNew*Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Lo->g/
		Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Hi->g;

	Att = TRANS_PROB_CONST*POW2(Ktt)*fttNew*Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Lo->g/
		Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Hi->g;

	Ast = TRANS_PROB_CONST*POW2(Kst)*fstNew*Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoTrip].Lo->g/
		Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoTrip].Hi->g;

	Ats = TRANS_PROB_CONST*POW2(Kts)*ftsNew*Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoSing].Lo->g/
		Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoSing].Hi->g;

	error = fabs( ( Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Emis->Aul+
		Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Emis->Aul )/
		(Ass+Ast+Ats+Att) - 1.f );

	if( error > 0.001 )
	{
		fprintf( ioQQQ, "FSM error %e LS %li HS %li LT %li HT %li Ratios Ass %e Att %e Ast %e Ats %e\n", error,
			ipLoSing, ipHiSing, ipLoTrip, ipHiTrip,
			Ass/Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Emis->Aul,
			Att/Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Emis->Aul,
			Ast/Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoTrip].Emis->Aul,
			Ats/Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoSing].Emis->Aul );
	}

	Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoSing].Emis->Aul = (realnum)Ass;
	Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoTrip].Emis->Aul = (realnum)Att;
	Transitions[ipHE_LIKE][nelem][ipHiSing][ipLoTrip].Emis->Aul = (realnum)Ast;
	Transitions[ipHE_LIKE][nelem][ipHiTrip][ipLoSing].Emis->Aul = (realnum)Ats;

	return;
}

/*ritoa converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
STATIC double ritoa(long li, long lf, long nelem, double k, double RI2)
{
   	/*	Variables are as follows:				*/
	/*	lg = larger of li and lf				*/
	/*	fmean = mean oscillator strength		*/
	/*		for a given level.					*/
	/*	mu = reduced mass of optical electron.	*/
	/*	EinsteinA = Einstein emission coef.		*/
	/*	w = angular frequency of transition.	*/
	/*	RI2_cm = square of rad. int. in cm^2.	*/
	long lg;
	double fmean,mu,EinsteinA,w,RI2_cm;

	DEBUG_ENTRY( "ritoa()" );

	mu = ELECTRON_MASS/(1+ELECTRON_MASS/(dense.AtomicWeight[nelem]*ATOMIC_MASS_UNIT));

	w = 2. * PI * k * SPEEDLIGHT;

	RI2_cm = RI2 * BOHR_RADIUS_CM * BOHR_RADIUS_CM;

	lg = (lf > li) ? lf : li;

	fmean = 2.0*mu*w*lg*RI2_cm/((3.0*H_BAR) * (2.0*li+1.0));

	EinsteinA = TRANS_PROB_CONST*k*k*fmean;

	/* ASSERT( EinsteinA > SMALLFLOAT ); */

	return EinsteinA;
}

realnum helike_transprob( long nelem, long ipHi, long ipLo )
{
	double Aul, Aul1;
	double Enerwn, n_eff_hi, n_eff_lo;
	long ipISO = ipHE_LIKE;
	/* charge to 4th power, needed for scaling laws for As*/
	double z4 = POW2((double)nelem);
	z4 *= z4;

	DEBUG_ENTRY( "helike_transprob()" );

	Enerwn = Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN;
	n_eff_hi = N_(ipHi) - helike_quantum_defect(nelem,ipHi);
	n_eff_lo = N_(ipLo) - helike_quantum_defect(nelem,ipLo);

	if( ipHi >= iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem] )
	{
		if( ipLo >= iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem] )
		{
			/* Neither upper nor lower is resolved...Aul's are easy.	*/
			Aul = HydroEinstA( N_(ipLo), N_(ipHi) )*z4;
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);

			ASSERT( Aul > 0.);
		}
		else 
		{
			/* Lower level resolved, upper not. First calculate Aul
			 * from upper level with ang mom one higher.	*/
			Aul = he_1trans( nelem, Enerwn, n_eff_hi, L_(ipLo)+1,
				S_(ipLo), -1, n_eff_lo, L_(ipLo), S_(ipLo), ipLo-3);

			iso.CachedAs[ipISO][nelem][ N_(ipHi)-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][0] = (realnum)Aul;

			Aul *= (2.*L_(ipLo)+3.) * S_(ipLo) / (4.*(double)N_(ipHi)*(double)N_(ipHi));

			if( L_(ipLo) != 0 )
			{
				/* For all l>0, add in transitions from upper level with ang mom one lower.	*/
				Aul1 = he_1trans( nelem, Enerwn, n_eff_hi, L_(ipLo)-1,
					S_(ipLo), -1, n_eff_lo, L_(ipLo), S_(ipLo), ipLo-3);

				iso.CachedAs[ipISO][nelem][ N_(ipHi)-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][1] = (realnum)Aul1;

				/* now add in this part after multiplying by stat weight for lHi = lLo-1.	*/
				Aul += Aul1*(2.*L_(ipLo)-1.) * S_(ipLo) / (4.*(double)N_(ipHi)*(double)N_(ipHi));
			}
			else
				iso.CachedAs[ipISO][nelem][ N_(ipHi)-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][1] = 0.f;

			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
			ASSERT( Aul > 0.);
		}
	}
	else
	{
		/* Both levels are resolved...do the normal bit.	*/
		if( Enerwn < 0. )
		{
			Aul = he_1trans( nelem, -1.*Enerwn,  n_eff_lo,
				L_(ipLo), S_(ipLo), ipLo-3, n_eff_hi, L_(ipHi), S_(ipHi), ipHi-3);
		}
		else
		{
			Aul = he_1trans( nelem, Enerwn,  n_eff_hi,
				L_(ipHi), S_(ipHi), ipHi-3, n_eff_lo, L_(ipLo), S_(ipLo), ipLo-3);
		}
	}

	return (realnum)Aul;
}

/* This routine is pure infrastructure and bookkeeping, no physics. */
void HelikeTransProbSetup( void )
{

	const int chLine_LENGTH = 1000;
	char chLine[chLine_LENGTH];

	FILE *ioDATA;
	bool lgEOL;

	long nelem, ipLo, ipHi, i, i1, i2, i3;

	DEBUG_ENTRY( "HelikeTransProbSetup()" );

	TransProbs = (double ***)MALLOC(sizeof(double **)*(unsigned)LIMELM );

	for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
	{

		TransProbs[nelem] = (double**)MALLOC(sizeof(double*)*(unsigned)(MAX_TP_INDEX+1) );

		for( ipLo=ipHe1s1S; ipLo <= MAX_TP_INDEX;++ipLo )
		{
			TransProbs[nelem][ipLo] = (double*)MALLOC(sizeof(double)*(unsigned)MAX_TP_INDEX );
		}
	}

	/********************************************************************/
	/*************** Read in data from he_transprob.dat	*****************/
	if( trace.lgTrace )
		fprintf( ioQQQ," HelikeTransProbSetup opening he_transprob.dat:");

	ioDATA = open_data( "he_transprob.dat", "r" );

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " HelikeTransProbSetup could not read first line of he_transprob.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	if( i1 !=TRANSPROBMAGIC || i2 != N_HE1_TRANS_PROB )
	{
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: the version of he_transprob.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: I expected to find the number %i %i and got %li %li instead.\n" ,
			 TRANSPROBMAGIC, N_HE1_TRANS_PROB, i1, i2 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* Initialize TransProbs[nelem][ipLo][ipHi] before filling it in.	*/
	for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{	
		for( ipHi=0; ipHi <= MAX_TP_INDEX; ipHi++ )
		{
			for( ipLo=0; ipLo < MAX_TP_INDEX; ipLo++ )
			{
				TransProbs[nelem][ipHi][ipLo] = -1.;
			}
		}
	}

	for( ipLo=1; ipLo <= N_HE1_TRANS_PROB; ipLo++ )
	{
		char *chTemp;

		/* get next line image */
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			BadRead();

		while( chLine[0]=='#' )
		{
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				BadRead();
		}

		i3 = 1;
		i1 = (long)FFmtRead(chLine,&i3,INPUT_LINE_LENGTH,&lgEOL);
		i2 = (long)FFmtRead(chLine,&i3,INPUT_LINE_LENGTH,&lgEOL);
		/* check that these numbers are correct */
		if( i1<0 || i2<=i1 )
		{
			fprintf( ioQQQ, " HelikeTransProbSetup detected insanity in he_transprob.dat.\n");
			cdEXIT(EXIT_FAILURE);
		}

		chTemp = chLine;

		/* skip over 2 tabs to start of data */
		for( i=0; i<1; ++i )
		{
			if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
			{
				fprintf( ioQQQ, " HelikeTransProbSetup could not init he_transprob\n" );
				cdEXIT(EXIT_FAILURE);
			}
			++chTemp;
		}

		/* now read in the data */
		for( nelem = ipHELIUM; nelem < LIMELM; nelem++ )
		{
			if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
			{
				fprintf( ioQQQ, " HelikeTransProbSetup could not scan he_transprob\n" );
				cdEXIT(EXIT_FAILURE);
			}
			++chTemp;

			sscanf( chTemp , "%le" , &TransProbs[nelem][i2][i1] );

			/** \todo	2	this test is out of place, where should it go? */
			if( lgEOL )
			{
				fprintf( ioQQQ, " HelikeTransProbSetup detected insanity in he_transprob.dat.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* check that ending magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " HelikeTransProbSetup could not read last line of he_transprob.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	i1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
	if( i1 !=TRANSPROBMAGIC || i2 != N_HE1_TRANS_PROB )
	{
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: the version of he_transprob.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			 " HelikeTransProbSetup: I expected to find the number %i %i and got %li %li instead.\n" ,
			 TRANSPROBMAGIC, N_HE1_TRANS_PROB, i1, i2 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* close the data file */
	fclose( ioDATA );
	return;
}
