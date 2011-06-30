/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroEinstA calculates Einstein A's from  osillator strengths*/
#include "cddefines.h"
#include "hydro_bauman.h"
#include "hydrooscilstr.h"
#include "hydroeinsta.h"
#include "iso.h"
#include "physconst.h"
#include "taulines.h"

double HydroEinstA(long int n1, 
	  long int n2)
{
	long int lower, iupper;
	double EinstA_v, 
	  ryd, 
	  xl, 
	  xmicron, 
	  xu;

	DEBUG_ENTRY( "HydroEinstA()" );
	/* (lower,upper) of Johnson 1972.  */

	/* strictly n -> n' transition probabilities
	 * no attempt to distribute according to l,l' */

	/* sort out the order of upper and lower, so can be called either way */
	lower = MIN2( n1 , n2 );
	iupper = MAX2( n1, n2 );
	if( lower < 1 || lower == iupper )
	{
		fprintf(ioQQQ," HydroEinstA called with impossible ns, =%li %li\n", lower, iupper);
		cdEXIT(EXIT_FAILURE);
	}

	xl = (double)lower;
	xu = (double)iupper;
	ryd = 1./POW2(xl) - 1./POW2(xu);
	xmicron = 1.E4/(ryd*RYD_INF);
	EinstA_v = HydroOscilStr(xl,xu)*TRANS_PROB_CONST*1e8f/(POW2(xmicron))*xl*xl/xu/xu;
	return( EinstA_v );
}

realnum hydro_transprob( long nelem, long ipHi, long ipLo )
{
	double Aul, Aul1;
	long ipISO = ipH_LIKE;
	/* charge to 4th power, needed for scaling laws for As*/
	double z4 = POW4((double)nelem+1.);
	DEBUG_ENTRY( "hydro_transprob()" );

	if( ipHi >= iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem] )
	{
		if( ipLo >= iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem] )
		{
			/* Neither upper nor lower is resolved...Aul's are easy.	*/
			Aul = HydroEinstA( N_(ipLo), N_(ipHi) )*z4;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);

			ASSERT( Aul > 0.);
		}
		else 
		{
			/* Lower level resolved, upper not. First calculate Aul
			 * from upper level with ang mom one higher.	*/
			Aul = H_Einstein_A( N_(ipHi), L_(ipLo)+1, N_(ipLo), L_(ipLo), nelem+1 );

			iso.CachedAs[ipISO][nelem][ N_(ipHi)-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][0] = (realnum)Aul;

			Aul *= (2.*L_(ipLo)+3.) * 2. / (2.*(double)N_(ipHi)*(double)N_(ipHi));

			if( L_(ipLo) != 0 )
			{
				/* For all l>0, add in transitions from upper level with ang mom one lower.	*/
				Aul1 = H_Einstein_A( N_(ipHi), L_(ipLo)-1, N_(ipLo), L_(ipLo), nelem+1 );

				iso.CachedAs[ipISO][nelem][ N_(ipHi)-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][1] = (realnum)Aul1;

				/* now add in this part after multiplying by stat weight for lHi = lLo-1.	*/
				Aul += Aul1*(2.*L_(ipLo)-1.) * 2. / (2.*(double)N_(ipHi)*(double)N_(ipHi));
			}
			else
				iso.CachedAs[ipISO][nelem][ N_(ipHi)-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][1] = 0.f;

			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
			ASSERT( Aul > 0.);
		}
	}
	else
	{
		if(  N_(ipHi) == N_(ipLo)  )
		{	
			/** \todo 1 define quantum defects and use scqdri to 
			 * calculate A's if levels are not exactly degenerate. */
			Aul = SMALLFLOAT;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else if( ipLo == 0 && ipHi == 1 )
		{
			/* 2s two photon */
			Aul =  8.226*pow((double)(nelem+1.),6.);
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else if( ipLo == 0 && ipHi == 2 )
		{
			Aul = 6.265e8*z4;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else if( abs( L_(ipLo) - L_(ipHi) )== 1 )
		{
			Aul = H_Einstein_A( N_(ipHi), L_(ipHi), N_(ipLo), L_(ipLo), nelem+1 );
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else
		{
			ASSERT( N_(ipHi) > N_(ipLo) );
			ASSERT( (L_(ipHi) == L_(ipLo)) || 
				( abs(L_(ipHi)-L_(ipLo)) > 1) );
			Aul = SMALLFLOAT;
			iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
	}

	return (realnum)Aul;
}
