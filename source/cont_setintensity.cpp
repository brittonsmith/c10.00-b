/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ContSetIntensity derive intensity of incident continuum */
/*extin do extinction of incident continuum as set by extinguish command */
/*sumcon sums L and Q for net incident continuum */
/*ptrcer show continuum pointers in real time following drive pointers command */
/*conorm normalize continuum to proper intensity */
/*qintr integrates Q for any continuum between two limits, used for normalization */
/*pintr integrates L for any continuum between two limits, used for normalization */
#include "cddefines.h"
#include "physconst.h"
#include "iso.h"
#include "noexec.h"
#include "ionbal.h"
#include "hextra.h"
#include "trace.h"
#include "dense.h"
#include "oxy.h"
#include "prt.h"
#include "heavy.h"
#include "rfield.h"
#include "phycon.h"
#include "called.h"
#include "hydrogenic.h"
#include "timesc.h"
#include "secondaries.h"
#include "opacity.h"
#include "thermal.h"
#include "ipoint.h"
#include "atmdat.h"
#include "rt.h"
#include "radius.h"
#include "geometry.h"
#include "grainvar.h"
#include "continuum.h"
#include "taulines.h"

/* these are weights used for continuum integration */
static const double aweigh[4]={-0.4305682,-0.1699905, 0.1699905, 0.4305682}; 
static const double fweigh[4]={ 0.1739274, 0.3260726, 0.3260726, 0.1739274};

/*conorm normalize continuum to proper intensity */
STATIC void conorm(void);

/*pintr integrates L for any continuum between two limits, used for normalization */
STATIC double pintr(double penlo, 
  double penhi);

/*qintr integrates Q for any continuum between two limits, used for normalization */
STATIC double qintr(double *qenlo, 
  double *qenhi);


/*sumcon sums L and Q for net incident continuum */
STATIC void sumcon(long int il, 
  long int ih, 
  realnum *q, 
  realnum *p, 
  realnum *panu);

/*extin do extinction of incident continuum as set by extinguish command */
STATIC void extin(realnum *ex1ryd);

/*ptrcer show continuum pointers in real time following drive pointers command */
STATIC void ptrcer(void);

void IncidentContinuumHere( void )
{

	DEBUG_ENTRY( "IncidentContinuumHere()" );
	double frac_beam_time;
	/* fraction of beamed continuum that is constant */
	double frac_beam_const;
	/* fraction of continuum that is isotropic */
	double frac_isotropic;

	double BigLog = 0.;
	for( long i = 0; i<rfield.nflux; ++i )
	{
		double flux_org = rfield.flux[0][i];
		double flux_now = ffun( rfield.anu[i] , 
			&frac_beam_time ,
			&frac_beam_const ,
			&frac_isotropic )*rfield.widflx[i]*rfield.ExtinguishFactor[i];

		double ratio = 1.;
		if( flux_org>SMALLFLOAT && flux_now>SMALLFLOAT )
		{
			ratio = fabs( log10( flux_now / flux_org ) );
			BigLog = max( ratio , BigLog );
		}
	}

	if( BigLog > 0.01 )
		fprintf(ioQQQ , "DEBUG diff continua %.2e\n", BigLog );
	return;
}
/* called by Cloudy to set up continuum 
 * return value is zero if  ok, 1 if no radiation - main would then stop */
int ContSetIntensity(void)
{
	bool lgCheckOK;

	long int i, 
	  ip, 
	  j, 
	  n;

	realnum EdenHeav, 
	  ex1ryd, 
	  factor, 
	  occ1, 
	  p, 
	  p1, 
	  p2, 
	  p3, 
	  p4, 
	  p5, 
	  p6, 
	  p7, 
	  p8,
	  pgn, 
	  phe, 
	  pheii, 
	  qgn;

	realnum xIoniz;

	double HCaseBRecCoeff,
	  wanu[4],
	  alf, 
	  bet, 
	  fntest, 
	  fsum, 
	  ecrit, 
	  tcompr, 
	  tcomp, 
	  RatioIonizToRecomb, 
	  r3ov2;

	double amean, 
	  amean2, 
	  amean3, 
	  peak, 
	  wfun[4];

	/* fraction of beamed continuum that is varies with time */
	double frac_beam_time , frac_beam_time1;
	/* fraction of beamed continuum that is constant */
	double frac_beam_const , frac_beam_const1;
	/* fraction of continuum that is isotropic */
	double frac_isotropic , frac_isotropic1;

	long int nelem , ion;

	DEBUG_ENTRY( "ContSetIntensity()" );

	/* set continuum */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " ContSetIntensity called.\n" );
	}

	/* find normalization factors for the continua - this decides whether continuum is
	 * per unit area of luminosity, and which is desired final product */
	conorm();

	/* define factors to convert rfeld.flux array into photon occupation array OCCNUM
	 * by multiplication */
	factor = (realnum)(EN1RYD/PI4/FR1RYD/HNU3C2);

	/*------------------------------------------------------------- */
	lgCheckOK = true;

	for( i=0; i < rfield.nupper; i++ )
	{
		/* this was original anu array with no continuum averaging */
		rfield.anu[i] = rfield.AnuOrg[i];
		rfield.ContBoltz[i] = 0.;
		fsum = 0.;
		amean = 0.;
		amean2 = 0.;
		amean3 = 0.;
		frac_beam_time = 0.;
		frac_beam_const = 0.;
		frac_isotropic = 0.;

		for( j=0; j < 4; j++ )
		{
			/* aweigh is symmetric about 0.5 widflx */
			wanu[j] = rfield.anu[i] + rfield.widflx[i]*aweigh[j];
			/* >>chng 02 jul 16, add test on continuum limits -
			 * this was exceeded when resolution set very large */
			wanu[j] = MAX2( wanu[j] , rfield.emm );
			wanu[j] = MIN2( wanu[j] , rfield.egamry );
			/* >>chng 06 feb 03, the continuum binning can change dramatically
			 * at some energies - make sure that this cell does not overextend the
			 * boundaries of its neighbors */
			if( i > 0 && i < rfield.nupper-1 )
			{
				wanu[j] = MAX2( wanu[j] , rfield.anu[i-1] + 0.5*rfield.widflx[i-1] );
				wanu[j] = MIN2( wanu[j] , rfield.anu[i+1] - 0.5*rfield.widflx[i+1]  );
			}

			wfun[j] = fweigh[j]*ffun(wanu[j] , 
				&frac_beam_time1 ,
				&frac_beam_const1 ,
				&frac_isotropic1 );
			/*if( i==76 )
				fprintf(ioQQQ,"DEBUG ffun %li %.4e at %.4e R\n", j , 
				ffun(wanu[j]) , wanu[j] );*/
			fsum += wfun[j];
			amean += wanu[j]*wfun[j];
			amean2 += wanu[j]*wanu[j]*wfun[j];
			amean3 += wanu[j]*wanu[j]*wanu[j]*wfun[j];
			frac_beam_time += fweigh[j]*frac_beam_time1;
			frac_beam_const += fweigh[j]*frac_beam_const1;
			frac_isotropic += fweigh[j]*frac_isotropic1;
		}

		ASSERT( fabs( 1.-frac_beam_time-frac_beam_const-frac_isotropic)<
			10.*FLT_EPSILON);
		/* This is a fix for ticket #1 */
		if( fsum*rfield.widflx[i] > BIGFLOAT )
		{
			fprintf( ioQQQ, "\n Cannot continue.  The continuum is far too intense.\n" );
			for( j=0; j < rfield.nShape; j++ )
			{
				if( (wfun[j]*rfield.widflx[i] > BIGFLOAT) && ( rfield.nShape > 1 ) )
				{
					fprintf( ioQQQ, " Problem is with source number %li\n", j );
					break;
				}
			}
			cdEXIT(EXIT_FAILURE);
		}

		rfield.flux[0][i] = (realnum)(fsum*rfield.widflx[i]);

		/* save separation into isotropic constant and beamed, and possibly 
		 * time-variable beamed continua */
		rfield.flux_beam_const[i] = rfield.flux[0][i] * (realnum)frac_beam_const;
		rfield.flux_beam_time[i] = rfield.flux[0][i] * (realnum)frac_beam_time;
		rfield.flux_isotropic[i] = rfield.flux[0][i] * (realnum)frac_isotropic;

		if( rfield.flux[0][i] > 0. )
		{
			/*fprintf(ioQQQ,"DEBUG %4li %e %e %.3e %.3e\n", 
				i , rfield.anu[i] , (amean2/amean) , amean2 , amean );*/
			rfield.anu[i] = (realnum)(amean2/amean);
			rfield.anu2[i] = (realnum)(amean3/amean);
			/* mesh must be strictly monotonically increasing - make it so */
			if( i > 0 && rfield.anu[i] <= rfield.anu[i-1] )
			{
				/* prevent roundoff from allowing i cell to lie below i-1
				 * cell when continuum mesh is very fine. */
				/* use 2*epsilon to protect against unusual rounding modes */
				rfield.anu[i] = rfield.anu[i-1]*(1.f+2.f*FLT_EPSILON);
				rfield.anu2[i] = pow2(rfield.anu[i]);
			}
			ASSERT( i==0 || rfield.anu[i] > rfield.anu[i-1] );
			/* define array of LOG10( nu(ryd) ) these are not guaranteed to
			 * be monotonically increasing due to loss of precision in log
			 * of float */
			rfield.anulog[i] = (realnum)log10(rfield.anu[i]);
		}

		else if( rfield.flux[0][i] == 0. )
		{
			rfield.anu2[i] = rfield.anu[i]*rfield.anu[i];
			rfield.anulog[i] = (realnum)log10(rfield.anu[i]);
		}

		else
		{
			rfield.anu2[i] = rfield.anu[i]*rfield.anu[i];
			fprintf( ioQQQ, " negative continuum returned at%6ld%10.2e%10.2e\n", 
			  i, rfield.anu[i], rfield.flux[0][i] );
			lgCheckOK = false;
		}
		rfield.anu3[i] = rfield.anu2[i]*rfield.anu[i];

		rfield.ConEmitReflec[0][i] = 0.;
		rfield.ConEmitOut[0][i] = 0.;
		rfield.convoc[i] = factor/rfield.widflx[i]/rfield.anu2[i];

		/* following are Compton exchange factors from Tarter */
		alf = 1./(1. + rfield.anu[i]*(1.1792e-4 + 7.084e-10*rfield.anu[i]));
		bet = 1. - alf*rfield.anu[i]*(1.1792e-4 + 2.*7.084e-10*rfield.anu[i])/
		  4.;
		rfield.csigh[i] = (realnum)(alf*rfield.anu2[i]*3.858e-25);
		rfield.csigc[i] = (realnum)(alf*bet*rfield.anu[i]*3.858e-25);
	}

	/*i = 76;
	fprintf(ioQQQ,"\nDEBUG %4li %e \n", 
		i , rfield.anu[i] );*/

	/* >>chng 05 jul 01 add this
	 * finished with stored continua - return these vectors */
#if 0
	/* commented out since we must conserve energy, and continuum was set with old widflx */
	/* now fix widflx array so that it is correct */
	for( i=1; i<rfield.nupper-1; ++i )
	{
		/*rfield.widflx[i] = rfield.anu[i+1] - rfield.anu[i];*/
		rfield.widflx[i] = ((rfield.anu[i+1] - rfield.anu[i]) + (rfield.anu[i] - rfield.anu[i-1]))/2.f;
	}
#endif

	if( !lgCheckOK )
	{
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace && trace.lgComBug )
	{
		fprintf( ioQQQ, "\n\n Compton heating, cooling coefficients \n" );
		for( i=0; i < rfield.nupper; i += 2 )
		{
			fprintf( ioQQQ, "%6ld%10.2e%10.2e%10.2e", i, rfield.anu[i], 
			  rfield.csigh[i], rfield.csigc[i] );
		}
		fprintf( ioQQQ, "\n" );
	}

	/* option to check frequencies in real time, drive pointers command,
	 * routine is below, is file static */
	if( trace.lgPtrace )
		ptrcer();

	/* extinguish continuum if set on */
	extin(&ex1ryd);

	/* now find peak of hydrogen ionizing continuum - for PDR calculations
	 * this will remain equal to 1 since the loop will not execute */
	prt.ipeak = 1;
	peak = 0.;

	for( i=iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1; i < rfield.nupper; i++ )
	{
		if( rfield.flux[0][i]*rfield.anu[i]/rfield.widflx[i] > (realnum)peak )
		{
			/* prt.ipeak points to largest f_nu at H-ionizing energies
			 * and is passed to other parts of code */
			/* i+1 to keep ipeak on fortran version of energy array */
			prt.ipeak = i+1;
			peak = rfield.flux[0][i]*rfield.anu[i]/rfield.widflx[i];
		}
	}

	/* find highest energy to consider in continuum flux array
	 * peak is the peak product nu*flux */
	peak = rfield.flux[0][prt.ipeak-1]/rfield.widflx[prt.ipeak-1]*
	  rfield.anu2[prt.ipeak-1];

	/* say what type of cpu this is, if desired */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " ContSetIntensity: The peak of the H-ion continuum is at %.3e Ryd - its value is %.2e\n", 
		  rfield.anu[prt.ipeak-1] , peak);
	}

	if( peak > 1e38 )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER The continuum is too intense to compute. Use a fainter continuum. (This is the nu*f_nu test)\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* FluxFaint set in zero.c; normally 1e-10 */
	/* this will be faintest level of continuum we want to consider.
	 * peak was set above, is peak of hydrogen ionizing radiation field, 
	 * and is zero if no H-ionizing radiation */
	fntest = peak*rfield.FluxFaint;
	{
		enum {DEBUG_LOC=false};
		/* print flux array then quit */
		if( DEBUG_LOC )
		{
			for( i=0; i<rfield.nupper; ++i )
			{
				fprintf(ioQQQ," consetintensityBUGGG\t%.2e\t%.2e\n" , 
					rfield.anu[i] , rfield.flux[0][i]/rfield.widflx[i] );
			}
			cdEXIT(EXIT_SUCCESS);
		}
	}

	if( fntest > 0. )
	{
		/* this test is not done in pdr conditions where NO H-ionizing radiation,
		 * since fntest is zero*/
		i = rfield.nupper;
		/* >>chng 05 aug 16, change ipeak to ipeak+3, to avoid possible one-off bugs
		 * where continuum barely goes into H-ionizing radiation */
		while( i > prt.ipeak+3 && 
			rfield.flux[0][i-1]*rfield.anu2[i-1]/rfield.widflx[i-1] < (realnum)fntest )
		{
			--i;
		}
	}
	else
	{
		/* when no H-ionizing radiation set to Lyman edge */
		/* >>chng 05 aug 16, change ipeak to ipeak+3, to avoid possible one-off bugs
		 * where continuum barely goes into H-ionizing radiation */
		i = iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]+3;
	}

	/* 
	 * this line of code dates from 1979 and IOA Cambridge.  removed july 19 95
	 * I think it was the last line of the original Cambridge Fortran source
	   nflux = MAX( ineon(1)+4 , i )
	 */

	/* >>chng 99 apr 28, reinstate the rfield.FluxFaint limit with nflux */
	rfield.nflux = i;

	/* trim down nflux, was set to rfield.nupper, the dimension of all vectors, in zero.c,
	 * in ContCreatePointers was set to nupper, the number of cells needed to get up to the
	 * high-energy limit of the code */
	while( rfield.flux[0][rfield.nflux-1] < SMALLFLOAT && rfield.nflux > 1 )
	{
		--rfield.nflux;
	}
	/* make sure we go high enough, avoid 1-off bugs as in draine field */
	++rfield.nflux;

	if( rfield.nflux == 1 )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER This incident continuum appears to have no radiation.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* >>chng 04 oct 10, add this limit - arrays will malloc to nupper, but will add unit
	 * continuum to [nflux] - this must be within array */
	rfield.nflux = MIN2( rfield.nupper-1 , rfield.nflux );

	/* >>chng 05 mar 01, nfine should not extend beyond continuum 
	 * make sure fine opacity scale does not extend beyond continuum we will use */
	rfield.nfine = rfield.nfine_malloc;
	while( rfield.nfine > 0 && rfield.fine_anu[rfield.nfine-1] > rfield.anu[rfield.nflux-1] )
	{
		--rfield.nfine;
	}

	/* check that continuum defined everywhere - look for zero's and comment if present */
	continuum.lgCon0 = false;
	ip = 0;
	for( i=0; i < rfield.nflux; i++ )
	{
		if( rfield.flux[0][i] == 0. )
		{
			if( called.lgTalk && !continuum.lgCon0 )
			{
				fprintf( ioQQQ, " NOTE Setcon: continuum has zero intensity starting at %11.4e Ryd.\n", 
				  rfield.anu[i] );
				continuum.lgCon0 = true;
			}
			++ip;
		}
	}

	if( continuum.lgCon0 && called.lgTalk )
	{
		fprintf( ioQQQ, 
			"%6ld cells in the incident continuum have zero intensity.  Problems???\n\n", 
		  ip );
	}

	/* check for devastating error in the continuum mesh or intensity */
	lgCheckOK = true;
	for( i=1; i < rfield.nflux; i++ )
	{
		if( rfield.flux[0][i] < 0. )
		{
			fprintf( ioQQQ, 
				" PROBLEM DISASTER Continuum has negative intensity at %.4e Ryd=%.2e %4.4s %4.4s\n", 
			  rfield.anu[i], rfield.flux[0][i], rfield.chLineLabel[i], rfield.chContLabel[i] );
			lgCheckOK = false;
		}
		else if( rfield.anu[i] <= rfield.anu[i-1] )
		{
			fprintf( ioQQQ, 
				" PROBLEM DISASTER cont_setintensity - internal error - continuum energies not in increasing order: energies follow\n" );
			fprintf( ioQQQ, 
				"%ld %e %ld %e %ld %e\n", 
			  i -1 , rfield.anu[i-1], i, rfield.anu[i], i +1, rfield.anu[i+1] );
			lgCheckOK = false;
		}
	}

	/* either of the ways lgCheckOK would be set true would be a major internal error */
	if( !lgCheckOK )
	{
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* turn off recoil ionization if high energy < 190R */
	if( rfield.anu[rfield.nflux-1] <= 190 )
	{
		ionbal.lgCompRecoil = false;
	}

	/* sum photons and energy, save mean */

	/* sum from low energy to Balmer edge */
	sumcon(1,iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH2p]-1,&rfield.qrad,&prt.pradio,&p1);

	/* sum over Balmer continuum */
	sumcon(iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH2p],iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1,&rfield.qbal,&prt.pbal,&p2);

	/* sum from Lyman edge to HeI edge */
	sumcon(iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s],
		iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0]-1,&prt.q,&p,&p3);

	/* sum from HeI to HeII edges */
	sumcon(iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0],
		iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s]-1,&rfield.qhe,&phe,&p4);

	/* sum from Lyman edge to carbon k-shell */
	sumcon(iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s],opac.ipCKshell-1,&rfield.qheii,&pheii,&p5);

	/* sum from c k-shell to gamma ray - where pairs start */
	sumcon( opac.ipCKshell , rfield.ipEnerGammaRay-1 , &prt.qx,
		&prt.xpow ,  &p6);

	/* complete sum up to high-energy limit */
	sumcon(rfield.ipEnerGammaRay,rfield.nflux,&prt.qgam,&prt.GammaLumin, &p7);

	/* find to estimate photoerosion timescale */
	n = ipoint(7.35e5);
	sumcon(n,rfield.nflux,&qgn,&pgn,&p8);
	timesc.TimeErode = qgn;

	/* find Compton temp */
	tcompr = (p1 + p2 + p3 + p4 + p5 + p6 + p7)/(prt.pradio + prt.pbal + 
	  p + phe + pheii + prt.xpow + prt.GammaLumin);

	tcomp = tcompr/(4.*6.272e-6);

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" mean photon energy=%10.3eR =%10.3eK low, high nu=%12.4e%12.4e\n", 
		  tcompr, tcomp, rfield.anu[0] - rfield.widflx[0]/2., rfield.anu[rfield.nflux-1] + 
		  rfield.widflx[rfield.nflux-1]/2. );
	}

	/* this is total power in ionizing radiation, per unit area */
	prt.powion = p + phe + pheii + prt.xpow + prt.GammaLumin;

	/* this is the total radiation intensity, erg cm-2 s-1 */
	continuum.TotalLumin = prt.pradio + prt.powion + prt.pbal;

	/* this is placed into the line stack on the first zone, then
	 * reset to zero, to end up with luminosity in the emission lines array.
	 * at end of iteration it is reset to TotalLumin */
	continuum.totlsv = continuum.TotalLumin;

	/* total H-ionizing photon number, */
	rfield.qhtot = prt.q + rfield.qhe + rfield.qheii + prt.qx + prt.qgam;

	/* ftotal photon number, all energies  */
	rfield.qtot = rfield.qhtot + rfield.qbal + rfield.qrad;

	if( prt.powion <= 0. && called.lgTalk )
	{
		rfield.lgHionRad = true;
		fprintf( ioQQQ, " NOTE There is no hydrogen-ionizing radiation.\n" );
		fprintf( ioQQQ, " Was this intended?\n\n" );
		/* check if any Balmer ionizing radiation since even metals will be 
		 * totally neutral */
		if( prt.pbal <=0. && called.lgTalk  )
		{
			fprintf( ioQQQ, " NOTE There is no Balmer continuum radiation.<<<<\n" );
			fprintf( ioQQQ, " Was this intended?\n\n" );
		}
	}

	else
	{
		rfield.lgHionRad = false;
	}

	/* option to add energy deposition due to fast neutrons, 
	 * entered as fraction of total photon luminosity */
	if( hextra.lgNeutrnHeatOn )
	{
		/* hextra.totneu is erg cm-2 s-1 in neutrons 
		 * hextra.effneu - efficiency default is unity */
		hextra.totneu = (realnum)(hextra.effneu*continuum.TotalLumin*
			pow((realnum)10.f,hextra.frcneu));
	}
	else
	{
		hextra.totneu = (realnum)0.;
	}

	/* temp correspond to energy density, printed in STARTR */
	phycon.TEnerDen = pow(continuum.TotalLumin/SPEEDLIGHT/7.56464e-15,0.25);

	/* sanity check for single blackbody, that energy density temperature
	 * is smaller than black body temperature */
	if( rfield.nShape==1 && 
		strcmp( rfield.chSpType[rfield.nShape-1], "BLACK" )==0 )
	{
		/* single black body, now confirm that TEnerDen is less than this temperature,
		 * >>>chng 99 may 02,
		 * in lte these are very very close, factor of 1.00001 allows for numerical
		 * errors, and apparently slightly different atomic coef in different parts
		 * of code.  eventaully all mustuse physonst.h and agree exactly */
		if( phycon.TEnerDen > 1.0001f*rfield.slope[rfield.nShape-1] )
		{
			fprintf( ioQQQ,
				"\n WARNING:  The energy density temperature (%g) is greater than the"
				" black body temperature (%g).  This is unphysical.\n\n",
				phycon.TEnerDen , rfield.slope[rfield.nShape-1]);
		}
	}

	/* incident continuum nu*f_nu at Hbeta and Ly-alpha */
	continuum.cn4861 = (realnum)(ffun(0.1875)*HPLANCK*FR1RYD*0.1875*0.1875);
	continuum.cn1216 = (realnum)(ffun(0.75)*HPLANCK*FR1RYD*0.75*0.75);
	continuum.sv4861 = continuum.cn4861;
	continuum.sv1216 = continuum.cn1216;
	/* flux density in V, erg / s / cm2 / hz */
	continuum.fluxv = (realnum)(ffun(0.1643)*HPLANCK*0.1643);
	continuum.fbeta = (realnum)(ffun(0.1875)*HPLANCK*0.1875*6.167e14);

	/* flux density nu*Fnu = erg / s / cm2
	 * EX1RYD is optional extinction factor at 1 ryd */
	prt.fx1ryd = (realnum)(ffun(1.000)*HPLANCK*ex1ryd*FR1RYD);

	realnum plsFrqConstant = (realnum)(ELEM_CHARGE_ESU/sqrt(PI*ELECTRON_MASS)/FR1RYD);
	ASSERT( plsFrqConstant > 2.7e-12 && plsFrqConstant < 2.8e-12 );
	
	/* check for plasma frequency - then zero out incident continuum
	 * for energies below this
	 * this is critical electron density, shielding of incident continuum
	 * if electron density is greater than this */
	ecrit = POW2(rfield.anu[0]/plsFrqConstant);

	if( dense.gas_phase[ipHYDROGEN]*1.2 > ecrit )
	{
		rfield.lgPlasNu = true;
                // This should be electron density, but we don't know that yet.
                // We use n_H (with 1.2 factor for He) as an ansatz.
		rfield.plsfrq = plsFrqConstant*sqrt(dense.gas_phase[ipHYDROGEN]*1.2f);
		rfield.plsfrqmax = rfield.plsfrq;
		rfield.ipPlasma = ipoint(rfield.plsfrq);

		/* save max pointer too */
		rfield.ipPlasmax = rfield.ipPlasma;

		/* now loop over all affected energies, setting incident continuum
		 * to zero there, and counting all as reflected */
		/* >>chng 01 jul 14, from i < ipPlasma to ipPlasma-1 - 
		 * when ipPlasma is 1 plasma freq is not on energy scale */
		for( i=0; i < rfield.ipPlasma-1; i++ )
		{
			/* count as reflected incident continuum */
			rfield.ConRefIncid[0][i] = rfield.flux[0][i];
			/* set continuum to zero there */
			rfield.flux_beam_const[i] = 0.;
			rfield.flux_beam_time[i] = 0.;
			rfield.flux_isotropic[i] = 0.;
			rfield.flux[0][i] = rfield.flux_beam_const[i] + rfield.flux_beam_time[i] +
				rfield.flux_isotropic[i];
		}
	}
	else
	{
		rfield.lgPlasNu = false;
		/* >>chng 01 jul 14, from 0 to 1 - 1 is the first array element on the F scale,
		 * ipoint would return this, so rest of code assumes ipPlasma is 1 plus correct index */
		rfield.ipPlasma = 1;
		rfield.plsfrqmax = 0.;
		rfield.plsfrq = 0.;
	}

	if( rfield.ipPlasma > 1 && called.lgTalk )
	{
		fprintf( ioQQQ, 
			"           !The plasma frequency is %.2e Ryd.  The incident continuum is set to 0 below this.\n", 
		  rfield.plsfrq );
	}

	rfield.occmax = 0.;
	rfield.tbrmax = 0.;
	for( i=0; i < rfield.nupper; i++ )
	{
		/* set up occupation number array */
		rfield.OccNumbIncidCont[i] = rfield.flux[0][i]*rfield.convoc[i];
		if( rfield.OccNumbIncidCont[i] > rfield.occmax )
		{
			rfield.occmax = rfield.OccNumbIncidCont[i];
			rfield.occmnu = rfield.anu[i];
		}
		/* following product is continuum brightness temperature */
		if( rfield.OccNumbIncidCont[i]*TE1RYD*rfield.anu[i] > rfield.tbrmax )
		{
			rfield.tbrmax = (realnum)(rfield.OccNumbIncidCont[i]*TE1RYD*rfield.anu[i]);
			rfield.tbrmnu = rfield.anu[i];
		}
		/* save continuum for next iteration */
		rfield.flux_total_incident[0][i] = rfield.flux[0][i];
		rfield.flux_beam_const_save[i] = rfield.flux_beam_const[i];
		rfield.flux_time_beam_save[i] = rfield.flux_beam_time[i];
		rfield.flux_isotropic_save[i] = rfield.flux_isotropic[i];
		/*fprintf(ioQQQ,"DEBUG type cont %li\t%.3e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			i, rfield.anu[i],
			rfield.flux[0][i],rfield.flux_beam_const[i],rfield.flux_beam_time[i],
			rfield.flux_isotropic[i]);
		fflush(ioQQQ);*/
	}

	/* if continuum brightness temp is large, where does it fall below
	 * 1e4k??? */
	if( rfield.tbrmax > 1e4 )
	{
		i = ipoint(rfield.tbrmnu)-1;
		while( i < rfield.nupper-1 && (rfield.OccNumbIncidCont[i]*TE1RYD*
		  rfield.anu[i] > 1e4) )
		{
			++i;
		}
		rfield.tbr4nu = rfield.anu[i];
	}
	else
	{
		rfield.tbr4nu = 0.;
	}

	/* if continuum occ num is large, where does it fall below 1? */
	if( rfield.occmax > 1. )
	{
		i = ipoint(rfield.occmnu)-1;
		while( i < rfield.nupper && (rfield.OccNumbIncidCont[i] > 1.) )
		{
			++i;
		}
		rfield.occ1nu = rfield.anu[i];
	}
	else
	{
		rfield.occ1nu = 0.;
	}

	/* remember if incident radiation field is less than 10*Habing ISM */
	/* >>chng 06 aug 01, change this test from continuum.TotalLumin to 
	 * energy in balmer and ionizing continuum, since this is the true habing field
	 * and is the continuum that interacts with gas.  When CMB set this
	 * tests on total did not trigger due to cold blackbody, which has little
	 * effect on gas, other than compton */
	if( (prt.powion + prt.pbal) < 1.8e-2 )
	{
		/* thermal.ConstTemp def is zero, set pos when constant temperature is set */
		rfield.lgHabing = true;
		/* >>chng 06 aug 01 also print warning if substantially below Habing, this may be a mistake */
		if( ((prt.powion + prt.pbal) < 1.8e-12) &&
			/* this is test for not constant temperature */
			(!thermal.lgTemperatureConstant) )
		{
			fprintf( ioQQQ, "\n >>>\n"
							" >>> NOTE The incident continuum is surprisingly faint.\n" );
			fprintf( ioQQQ, 
				" >>> The total energy in the Balmer and Lyman continua is %.2e erg cm-2 s-1.\n"
				,(prt.powion + prt.pbal));
			fprintf( ioQQQ, " >>> This is many orders of magnitude fainter than the ISM galactic background.\n" );
			fprintf( ioQQQ, " >>> This seems unphysical - please check that the continuum intensity has been properly set.\n" );
			fprintf( ioQQQ, " >>> YOU MAY BE MAKING A BIG MISTAKE!!\n >>>\n\n\n\n" );
		}
	}

	/* fix ionization parameter (per hydrogen) at inner edge */
	rfield.uh = (realnum)(rfield.qhtot/dense.gas_phase[ipHYDROGEN]/SPEEDLIGHT);
	rfield.uheii = (realnum)((rfield.qheii + prt.qx)/dense.gas_phase[ipHYDROGEN]/SPEEDLIGHT);
	if( rfield.uh > 1e10 )
	{
		fprintf( ioQQQ, "\n\n"
						" CAUTION The incident radiation field is surprisingly intense.\n" );
		fprintf( ioQQQ, 
			" The dimensionless hydrogen ionization parameter is %.2e.\n"
			, rfield.uh );
		fprintf( ioQQQ, " This is many orders of magnitude brighter than commonly seen.\n" );
		fprintf( ioQQQ, " This seems unphysical - please check that the radiation field intensity has been properly set.\n" );
		fprintf( ioQQQ, " YOU MAY BE MAKING A BIG MISTAKE!!\n\n\n\n\n" );
	}

	/* guess first temperature and neutral h density */
	if( thermal.ConstTemp > 0. )
	{
		phycon.te = thermal.ConstTemp;
	}
	else
	{
		if( rfield.uh > 0. )
		{
			phycon.te = (20000.+log10(rfield.uh)*5000.);
			phycon.te = MAX2(8000. , phycon.te );
		}
		else
		{
			phycon.te = 1000.;
		}
	}

	/* this is an option to stop after printing header only */
	if( noexec.lgNoExec )
		return(0);

	/* estimate secondary ionization rate - probably 0, but possible extra
	 * SetCsupra set with "set csupra" */
	/* >>>chng 99 apr 29, added cosmic ray ionization since this is used to get
	 * helium ionization fraction, and was zero in pdr, so He turned off at start,
	 * and never turned back on */
	/* coef on cryden is from highen.c */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		for( ion=0; ion<nelem+1; ++ion )
		{
			secondaries.csupra[nelem][ion] = 
				secondaries.SetCsupra + hextra.cryden*2e-9f;
		}
	}

	/*********************************************************************
	 *                                                                   *
	 * estimate hydrogen's level of ionization                           *
	 *                                                                   *
	 *********************************************************************/

	/* create fake ionization balance, but will conserve number of hydrogens */
	dense.xIonDense[ipHYDROGEN][0] = 0.;
	dense.xIonDense[ipHYDROGEN][1] = dense.gas_phase[ipHYDROGEN];
	/* this must be zero since PresTotCurrent will do radiation pressure due to H */
	StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop = 0.;

	/* "extra" electrons from command line, or assumed residual electrons */
	double EdenExtraLocal = dense.EdenExtra + 
		/* if we are in a molecular cloud the current logic could badly fail
		* do not let electron density fall below 1e-7 of H density */
		1e-7*dense.gas_phase[ipHYDROGEN];
	dense.eden = dense.xIonDense[ipHYDROGEN][1] + EdenExtraLocal;

	/* hydrogen case B recombination coefficient */
	HCaseBRecCoeff = (-9.9765209 + 0.158607055*phycon.telogn[0] + 0.30112749*
	  phycon.telogn[1] - 0.063969007*phycon.telogn[2] + 0.0012691546*
	  phycon.telogn[3])/(1. + 0.035055422*phycon.telogn[0] - 
	  0.037621619*phycon.telogn[1] + 0.0076175048*phycon.telogn[2] - 
	  0.00023432613*phycon.telogn[3]);
	HCaseBRecCoeff = pow(10.,HCaseBRecCoeff)/phycon.te;

	double CollIoniz = t_ADfA::Inst().coll_ion(1,1,phycon.te);
	double OtherIonization = rfield.qhtot*2e-18 + secondaries.csupra[ipHYDROGEN][0];
	double newEden = dense.eden;
	long loopCount = 0;

	do
	{
		/* update electron density */
		dense.eden = newEden;
		double RatioIoniz = 
			(CollIoniz*dense.eden+OtherIonization)/(HCaseBRecCoeff*dense.eden);
		if( RatioIoniz<1e-3 )
		{
			/* very low ionization solution */
			dense.xIonDense[ipHYDROGEN][1] = (realnum)(
				dense.gas_phase[ipHYDROGEN]*RatioIoniz);
			dense.xIonDense[ipHYDROGEN][0] = dense.gas_phase[ipHYDROGEN] - 
				dense.xIonDense[ipHYDROGEN][1];
			ASSERT( dense.xIonDense[ipHYDROGEN][0]>=0. &&
				dense.xIonDense[ipHYDROGEN][0]<=dense.gas_phase[ipHYDROGEN]);
			ASSERT( dense.xIonDense[ipHYDROGEN][1]>=0. &&
				dense.xIonDense[ipHYDROGEN][1]<dense.gas_phase[ipHYDROGEN]);
		}
		else if( RatioIoniz>1e3 )
		{
			/* very high ionization solution */
			dense.xIonDense[ipHYDROGEN][0] = (realnum)(
				dense.gas_phase[ipHYDROGEN]/RatioIoniz);
			dense.xIonDense[ipHYDROGEN][1] = dense.gas_phase[ipHYDROGEN] - 
				dense.xIonDense[ipHYDROGEN][0];
			ASSERT( dense.xIonDense[ipHYDROGEN][0]>=0. &&
				dense.xIonDense[ipHYDROGEN][0]<dense.gas_phase[ipHYDROGEN]);
			ASSERT( dense.xIonDense[ipHYDROGEN][1]>=0. &&
				dense.xIonDense[ipHYDROGEN][1]<=dense.gas_phase[ipHYDROGEN]);
		}
		else
		{
			/* intermediate ionization - solve quadratic */
			double alpha = HCaseBRecCoeff + CollIoniz ;
			double beta = HCaseBRecCoeff*EdenExtraLocal + OtherIonization +
				(EdenExtraLocal - dense.gas_phase[ipHYDROGEN])*CollIoniz;
			double gamma = -dense.gas_phase[ipHYDROGEN]*(OtherIonization+EdenExtraLocal*CollIoniz);

			double discriminant = POW2(beta) - 4.*alpha*gamma;
			if( discriminant <0 )
			{
				/* oops */
				fprintf(ioQQQ," DISASTER PROBLEM cont_initensity found negative discriminant.\n");
				TotalInsanity();
			}

			dense.xIonDense[ipHYDROGEN][1] = (realnum)((-beta + sqrt(discriminant))/(2.*alpha));
			if( dense.xIonDense[ipHYDROGEN][1]> dense.gas_phase[ipHYDROGEN] )
			{
				/* oops */
				fprintf(ioQQQ," DISASTER PROBLEM cont_initensity found n(H+)>n(H).\n");
				TotalInsanity();
			}
			dense.xIonDense[ipHYDROGEN][0] = dense.gas_phase[ipHYDROGEN] -
				dense.xIonDense[ipHYDROGEN][1];
			if( dense.xIonDense[ipHYDROGEN][0]<= 0 )
			{
				/* oops */
				fprintf(ioQQQ," DISASTER PROBLEM cont_initensity found n(H0)<0.\n");
				TotalInsanity();
			}
		}


		if( dense.xIonDense[ipHYDROGEN][1] > 1e-30 )
		{
			StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop = dense.xIonDense[ipHYDROGEN][0];
		}
		else
		{
			StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop = 0.;
		}

		/* now save estimates of whether induced recombination is going
		 * to be important -this is a code pacesetter since GammaBn is slower
		 * than GammaK */
		hydro.lgHInducImp = false;
		for( i=ipH1s; i < iso.numLevels_max[ipH_LIKE][ipHYDROGEN]; i++ )
		{
			if( rfield.OccNumbIncidCont[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][i]-1] > 0.01 )
				hydro.lgHInducImp = true;
		}

		/*******************************************************************
		 *                                                                 *
		 * estimate helium's level of ionization                           *
		 *                                                                 *
		 *******************************************************************/

		/* only if helium is turned on */
		if( dense.lgElmtOn[ipHELIUM] )
		{
			/* next estimate level of helium singly ionized */
			xIoniz = (realnum)t_ADfA::Inst().coll_ion(2,2,phycon.te);
			/* >>chng 05 jan 05, add cosmic rays */
			xIoniz = (realnum)(xIoniz*dense.eden + rfield.qhe*1e-18 +  secondaries.csupra[ipHELIUM][1] );
			double RecTot = HCaseBRecCoeff*dense.eden;
			RatioIonizToRecomb = xIoniz/RecTot;

			/* now estimate level of helium doubly ionized */
			xIoniz = (realnum)t_ADfA::Inst().coll_ion(2,1,phycon.te);
			/* >>chng 05 jan 05, add cosmic rays */
			xIoniz = (realnum)(xIoniz*dense.eden + rfield.qheii*1e-18 +  ionbal.CosRayIonRate );

			/* rough charge dependence */
			RecTot *= 4.;
			r3ov2 = xIoniz/RecTot;

			/* now set level of helium ionization */
			if( RatioIonizToRecomb > 0. )
			{
				dense.xIonDense[ipHELIUM][1] = (realnum)(dense.gas_phase[ipHELIUM]/(1./RatioIonizToRecomb + 1. + r3ov2));
				dense.xIonDense[ipHELIUM][0] = (realnum)(dense.xIonDense[ipHELIUM][1]/RatioIonizToRecomb);
			}
			else
			{
				/* no He ionizing radiation */
				dense.xIonDense[ipHELIUM][1] = 0.;
				dense.xIonDense[ipHELIUM][0] = dense.gas_phase[ipHELIUM];
			}

			dense.xIonDense[ipHELIUM][2] = (realnum)(dense.xIonDense[ipHELIUM][1]*r3ov2);

			if( dense.xIonDense[ipHELIUM][2] > 1e-30 )
			{
				StatesElemNEW[ipHELIUM][ipHELIUM-ipH_LIKE][ipH1s].Pop = dense.xIonDense[ipHELIUM][1];
			}
			else
			{
				StatesElemNEW[ipHELIUM][ipHELIUM-ipH_LIKE][ipH1s].Pop = 0.;
			}
		}
		else
		{
			/* case where helium is turned off */
			dense.xIonDense[ipHELIUM][1] = 0.;
			dense.xIonDense[ipHELIUM][0] = 0.;
			dense.xIonDense[ipHELIUM][2] = 0.;
			StatesElemNEW[ipHELIUM][ipHELIUM-ipH_LIKE][ipH1s].Pop = 0.;
		}
		
		/* update electron density */
		newEden = dense.xIonDense[ipHYDROGEN][1] + EdenExtraLocal + dense.xIonDense[ipHELIUM][1] + 2.*dense.xIonDense[ipHELIUM][2];

		loopCount++;
	}
	/* repeat above until guessed and calculated eden agree to at least 0.1%. */
	while( loopCount < 10 && fabs(newEden/dense.eden - 1.) > 0.001 );

	if( dense.xIonDense[ipHYDROGEN][0]<0.)
		TotalInsanity();
	else if( dense.xIonDense[ipHYDROGEN][0] == 0. )
	{
		fprintf(ioQQQ,"PROBLEM the derived atomic hydrogen density is zero.\n");
		if( dense.gas_phase[ipHYDROGEN]<1e-5 && rfield.uh > 1e10)
		{
			fprintf(ioQQQ,"This is almost certainly due to floating point "
				"limits on this computer.\nThe ionization parameter is very large,"
				" the density is very small,\nand the H^0 density cannot be"
				" stored in a float.\n");
			cdEXIT( EXIT_FAILURE );
		}
	}
	ASSERT( dense.xIonDense[ipHYDROGEN][0] >0 && dense.xIonDense[ipHYDROGEN][1]>= 0.);

	/* update electron density */
	dense.eden = newEden;

	if( dense.eden <= SMALLFLOAT )
	{
		/* no electrons! */
		fprintf(ioQQQ,"\n PROBLEM DISASTER - this simulation has no source"
			" of ionization.  The electron density is zero.  Consider "
			"adding a source of ionization such as cosmic rays.\n\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* fix range of stages of ionization */
	for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			// IonHigh[n] is the highest stage of ionization present
			// the IonHigh array index is on the C scale, so [0] is hydrogen
			// the value is also on the C scale, so element [nelem] can range
			// from 0 to nelem+1 
			dense.IonHigh[nelem] = nelem + 1;

			dense.IonLow[nelem] = 0;
			// for very intense radiation fields very heavy elements, N>Fe, will fail
			// in ion_solver due to ill conditioned matrix.  all populations are in
			// fully stripped state.  Start will fully stripped ion distribution in this case.
			if( rfield.uh > 1e15 )
			{
				//trim down highest stage to be within incident radiation field
				while ( rfield.anu[Heavy.ipHeavy[nelem][dense.IonHigh[nelem]-1]] > 
					rfield.anu[rfield.nflux] && dense.IonHigh[nelem]>1 )
					--dense.IonHigh[nelem];
				dense.IonLow[nelem] = max( 0 , dense.IonHigh[nelem]-1 );
			}

			// make low-stage populations zero
			for( ion=0; ion<dense.IonLow[nelem]; ++ion )
				dense.xIonDense[nelem][ion] = 0.;

			/* >>chng 04 jan 13, add this test, caught by Orly Gnat */
			/* check on actual zero abundances of lower stages - this will only 
			 * happen when ionization is set with element ionization command */
			if( dense.lgSetIoniz[nelem] )
			{
				while( dense.SetIoniz[nelem][dense.IonLow[nelem]] < dense.density_low_limit )
					++dense.IonLow[nelem];
				while( dense.SetIoniz[nelem][dense.IonHigh[nelem]] < dense.density_low_limit )
					--dense.IonHigh[nelem];
			}
		}
		else
		{
			/* this element is turned off, make stages impossible */
			dense.IonLow[nelem] = -1;
			dense.IonHigh[nelem] = -1;
		}
	}

	/* estimate electrons from heavies, assuming each at least
	 * 1 times ionized */
	EdenHeav = 0.;
	realnum atomFrac = dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN];
	realnum firstIonFrac = dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN];
	for( nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
	{
		// for very high ionization parameters, logic above may have raise
		// lowest ion stage
		if( dense.lgElmtOn[nelem] )
		{
			if( dense.IonLow[nelem] == 0)
				dense.xIonDense[nelem][0] = dense.gas_phase[nelem] * atomFrac;
			if( dense.IonLow[nelem] <= 1)
			{
				dense.xIonDense[nelem][1] = dense.gas_phase[nelem] * firstIonFrac;
				EdenHeav += dense.xIonDense[nelem][1];
			}
		}
	}

	/* modify estimate of electron density */
	dense.eden += EdenHeav;

	/* >>chng 05 jan 05, insure positive eden */
	dense.eden = MAX2( SMALLFLOAT , dense.eden );

	if( dense.EdenSet > 0. )
	{
		dense.eden = dense.EdenSet;
	}

	dense.EdenHCorr = dense.eden;

	if( dense.eden < 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Negative electron density results in ContSetIntensity.\n" );
		fprintf( ioQQQ, "%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e\n", 
		  dense.eden, dense.xIonDense[ipHYDROGEN][1], dense.xIonDense[ipHELIUM][1], 
		  dense.xIonDense[ipHELIUM][2], dense.gas_phase[ipCARBON], dense.EdenExtra );
		TotalInsanity();
	}

	if( dense.EdenSet > 0. )
	{
		dense.EdenTrue = dense.EdenSet;
	}
	else
		dense.EdenTrue = dense.eden;

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" ContSetIntensity sets initial EDEN to %.4e, contributors H+=%.2e He+, ++= %.2e %.2e Heav %.2e extra %.2e\n", 
		  dense.eden ,
		  dense.xIonDense[ipHYDROGEN][1],
		  dense.xIonDense[ipHELIUM][1],
          2.*dense.xIonDense[ipHELIUM][2],
		  EdenHeav,
		  dense.EdenExtra);
	}

	// photon occupation number at 1 Ryd - used for printout
	occ1 = (realnum)(prt.fx1ryd/HNU3C2/PI4/FR1RYD);
	if( occ1 > 1. )
		rfield.lgOcc1Hi = true;
	else
		rfield.lgOcc1Hi = false;

	if( trace.lgTrace && trace.lgConBug )
	{
		fprintf(ioQQQ,"\ntrace continuum print of %li incident spectral "
			"components\n", rfield.nShape);
		fprintf(ioQQQ,"  #   type Illum Beam? 1/cos TimeVary?\n");
		for( rfield.ipSpec=0; rfield.ipSpec<rfield.nShape; ++rfield.ipSpec )
		{
			fprintf(ioQQQ,"%3li %6s %5i     %c %.3f        %c\n", 
				rfield.ipSpec ,
				rfield.chSpType[rfield.ipSpec],
				rfield.Illumination[rfield.ipSpec],
				TorF(rfield.lgBeamed[rfield.ipSpec]),
				rfield.OpticalDepthScaleFactor[rfield.ipSpec],
				TorF(rfield.lgTimeVary[rfield.ipSpec]));
		}
		fprintf(ioQQQ,"\n");

		/*  print some useful pointers to ionization edges */
		fprintf( ioQQQ, " H2,1=%5ld%5ld NX=%5ld IRC=%5ld\n", 
		  iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][2], 
		  iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s],
		  opac.ipCKshell, 
		  ionbal.ipCompRecoil[ipHYDROGEN][0] );

		fprintf( ioQQQ, " CARBON" );
		for( i=0; i < 6; i++ )
			fprintf( ioQQQ, "%5ld", Heavy.ipHeavy[ipCARBON][i] );
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " OXY" );
		for( i=0; i < 8; i++ )
			fprintf( ioQQQ, "%5ld", Heavy.ipHeavy[ipOXYGEN][i] );
		fprintf( ioQQQ, "%5ld%5ld%5ld\n", opac.ipo3exc[0], 
		  oxy.i2d, oxy.i2p );

		fprintf( ioQQQ, 
			"\n\n                                                  PHOTONS PER CELL (NOT RYD)\n" );
		fprintf( ioQQQ, 
			"\n\n                                        nu, flux, wid, occ \n" );
		fprintf( ioQQQ, 
			" " );

		for( i=0; i < rfield.nflux; i++ )
		{
			fprintf( ioQQQ, "%4ld%10.2e%10.2e%10.2e%10.2e", i, 
			  rfield.anu[i], rfield.flux[0][i], rfield.widflx[i], 
			  rfield.OccNumbIncidCont[i] );
		}
		fprintf( ioQQQ, " \n" );
	}

	/* zero out some continua related to the ots rates,
	 * prototype and routine in RT_OTS_Update.  This is done here since summed cont will
	 * be set to rfield */
	RT_OTS_Zero();

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " ContSetIntensity returns, nflux=%5ld anu(nflux)=%11.4e eden=%10.2e\n", 
		  rfield.nflux, rfield.anu[rfield.nflux-1], dense.eden );
	}

	return(0);
}

/*sumcon sums L and Q for net incident continuum */
STATIC void sumcon(long int il, 
  long int ih, 
  realnum *q, 
  realnum *p, 
  realnum *panu)
{
	long int i, 
	  iupper; /* used as upper limit to the sum */

	DEBUG_ENTRY( "sumcon()" );

	*q = 0.;
	*p = 0.;
	*panu = 0.;

	/* soft continua may not go as high as the requested bin */
	iupper = MIN2(rfield.nflux,ih);

	/* n.b. - in F77 loop IS NOT executed when IUPPER < IL */
	for( i=il-1; i < iupper; i++ )
	{
		/* sum photon number */
		*q += rfield.flux[0][i];
		/* en1ryd is needed to stop overflow */
		/* sum flux */
		*p += (realnum)(rfield.flux[0][i]*(rfield.anu[i]*EN1RYD));
		/* this sum needed for means */
		*panu += (realnum)(rfield.flux[0][i]*(rfield.anu2[i]*EN1RYD));
	}

	return;
}

/*ptrcer show continuum pointers in real time following drive pointers command */
STATIC void ptrcer(void)
{
	char chCard[INPUT_LINE_LENGTH];
	/* in case of checking everything, will write errors to this file */
	FILE * ioERRORS=NULL;
	bool lgEOL;
	char chKey;
	long int i, 
	  ipnt, 
	  j;
	double pnt, 
	  t1, 
	  t2;

	DEBUG_ENTRY( "ptrcer()" );

	fprintf( ioQQQ, " There are two ways to do this:\n");
	fprintf( ioQQQ, " do you want me to test all the pointers (enter y)\n");
	fprintf( ioQQQ, " or do you want to enter energies yourself? (enter n)\n" );

	if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
	{
		fprintf( ioQQQ, " error getting input \n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* this must be either y or n */
	chKey = chCard[0];

	if( chKey == 'n' )
	{
		/* this branch, enter energies by hand, and see what happens */
		fprintf( ioQQQ, " Enter energy (Ryd); 0 to stop; negative is log.\n" );
		pnt = 1.;
		while( pnt!=0. )
		{
			if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " error getting input2 \n" );
				cdEXIT(EXIT_FAILURE);
			}
			/* now get the number off the line */
			i = 1;
			pnt = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);

			/* bail if no number at all, or it is zero*/
			if( lgEOL || pnt==0. )
			{
				break;
			}

			/* if number negative then interpret as log */
			if( pnt < 0. )
			{
				pnt = pow(10.,pnt);
			}

			/* get pointer to call */
			ipnt = ipoint(pnt);
			fprintf( ioQQQ, " Cell num%4ld center:%10.2e width:%10.2e low:%10.2e hi:%10.2e convoc:%10.2e\n", 
			  ipnt, rfield.anu[ipnt-1], rfield.widflx[ipnt-1], 
			  rfield.anu[ipnt-1] - rfield.widflx[ipnt-1]/2., 
			  rfield.anu[ipnt-1] + rfield.widflx[ipnt-1]/2., 
			  rfield.convoc[ipnt-1] );
		}
	}

	else if( chKey == 'y' )
	{
		/* first check that ipoint will not crash due to out of range call*/
		if( rfield.anu[0] - rfield.widflx[0]/2.*0.9 < continuum.filbnd[0] )
		{
			fprintf( ioQQQ," ipoint would crash since lowest desired energy of %e ryd is below limit of %e\n",
				rfield.anu[0] - rfield.widflx[0]/2.*0.9 , continuum.filbnd[0] );
			fprintf( ioQQQ," width of cell is %e\n",rfield.widflx[0]);
			cdEXIT(EXIT_FAILURE);
		}

		else if( rfield.anu[rfield.nflux-1] + rfield.widflx[rfield.nflux-1]/2.*0.9 > 
			continuum.filbnd[continuum.nrange] )
		{
			fprintf( ioQQQ," ipoint would crash since highest desired energy of %e ryd is above limit of %e\n",
				rfield.anu[rfield.nflux-1] + rfield.widflx[rfield.nflux-1]/2.*0.9 , 
				continuum.filbnd[continuum.nrange-1] );
			fprintf( ioQQQ," width of cell is %e\n",rfield.widflx[rfield.nflux]);
			fprintf( ioQQQ," this, previous cells are %e %e\n",
				rfield.anu[rfield.nflux-1],rfield.anu[rfield.nflux-2]);
			cdEXIT(EXIT_FAILURE);
		}

		/* this branch check everything, write errors to error file */
		fprintf( ioQQQ, " errors output on errors.txt\n");
		fprintf( ioQQQ, " IP(cor),IP(fount),nu lower, upper of found, desired cell.\n" );

		/* error file not open, set to null so we can check later */
		ioERRORS = NULL;
		for( i=0; i < rfield.nflux-1; i++ )
		{
			t1 = rfield.anu[i] - rfield.widflx[i]/2.*0.9;
			t2 = rfield.anu[i] + rfield.widflx[i]/2.*0.9;

			j = ipoint(t1);
			if( j != i+1 )
			{
				/* open file for errors if not already open */
				if( ioERRORS == NULL )
				{
					ioERRORS = open_data( "errors.txt", "w", AS_LOCAL_ONLY );
					fprintf( ioQQQ," created errors.txt file with error summary\n");
				}

				fprintf( ioQQQ, " Pointers do not agree for lower bound of cell%4ld, %e\n",
					i, rfield.anu[i]);
				fprintf( ioERRORS, " Pointers do not agree for lower bound of cell%4ld, %e\n",
					i, rfield.anu[i] );
			}

			j = ipoint(t2);
			if( j != i+1 )
			{
				/* open file for errors if not already open */
				if( ioERRORS == NULL )
				{
					ioERRORS = open_data( "errors.txt", "w", AS_LOCAL_ONLY );
					fprintf( ioQQQ," created errors.txt file with error summary\n");
				}
				fprintf( ioQQQ, " Pointers do not agree for upper bound of cell%4ld, %e\n", 
					i , rfield.anu[i]);
				fprintf( ioERRORS, " Pointers do not agree for upper bound of cell%4ld, %e\n", 
					i , rfield.anu[i]);
			}

		}
	}

	else
	{
		fprintf( ioQQQ, "I do not understand this key, sorry.  %c\n", chKey );
		cdEXIT(EXIT_FAILURE);
	}

	if( ioERRORS!=NULL )
		fclose( ioERRORS );
	cdEXIT(EXIT_FAILURE);
}

/*extin do extinction of incident continuum as set by extinguish command */
STATIC void extin(realnum *ex1ryd)
{

	DEBUG_ENTRY( "extin()" );

	/* modify input continuum by leaky absorber
	 * power law fit to 
	 * >>refer	XUV	extinction	Cruddace et al. 1974 ApJ 187, 497. */
	if( rfield.ExtinguishColumnDensity == 0. )
	{
		*ex1ryd = 1.;

		for( long i=0; i<rfield.nupper; ++i )
			rfield.ExtinguishFactor[i] = 1.;
	}
	else
	{
		double absorb = 1. - rfield.ExtinguishLeakage;
		double factor = rfield.ExtinguishColumnDensity*
			rfield.ExtinguishConvertColDen2OptDepth;
		/* extinction at 1 4 Ryd */
		*ex1ryd = (realnum)(rfield.ExtinguishLeakage + absorb*sexp(factor));

		// low energy limit of extinction
		long low = ipoint(rfield.ExtinguishLowEnergyLimit);

		for( long i=0; i<low-1; ++i )
			rfield.ExtinguishFactor[i] = 1.;

		for( long i=low-1; i < rfield.nupper; i++ )
		{
			realnum extfactor = (realnum)(rfield.ExtinguishLeakage + absorb*
				sexp(factor*(pow(rfield.anu[i],rfield.ExtinguishEnergyPowerLow))));

			rfield.ExtinguishFactor[i] = extfactor;
			rfield.flux_beam_const[i] *= extfactor;
			rfield.flux_beam_time[i] *= extfactor;
			rfield.flux_isotropic[i] *= extfactor;

			rfield.flux[0][i] = rfield.flux_beam_const[i] + rfield.flux_beam_time[i] +
				rfield.flux_isotropic[i];
		}
	}
	return;
}

/*conorm normalize continuum to proper intensity */
STATIC void conorm(void)
{
	long int i;
	double xLog_radius_inner, 
	  diff, 
	  f, 
	  flx1, 
	  flx2, 
	  pentrd, 
	  qentrd;

	DEBUG_ENTRY( "conorm()" );

	xLog_radius_inner = log10(radius.Radius);

	/* this is a sanity check, it can't happen */
	for( i=0; i < rfield.nShape; i++ )
	{
		if( strcmp(rfield.chRSpec[i],"UNKN") == 0 )
		{
			fprintf( ioQQQ, " UNKN spectral normalization cannot happen.\n" );
			fprintf( ioQQQ, " conorm punts.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		else if( strcmp(rfield.chRSpec[i],"SQCM") != 0 && 
			      strcmp(rfield.chRSpec[i],"4 PI") != 0 )
		{
			fprintf( ioQQQ, " chRSpec must be SQCM or 4 PI, and it was %4.4s.  This cannot happen.\n", 
			  rfield.chRSpec[i] );
			fprintf( ioQQQ, " conorm punts.\n" );
			cdEXIT(EXIT_FAILURE);
		}


		/* this sanity check makes sure that atlas.mod or werner.mod grids
		 * are for the current version of the code */
		if( strcmp(rfield.chSpType[i],"VOLK ") == 0 )
		{
			if( !fp_equal( rfield.RSFCheck[i], continuum.ResolutionScaleFactor ) )
			{
				fprintf( ioQQQ,"\n\n PROBLEM DISASTER At least one of the compiled stellar atmosphere"
					 " grids has been compiled with a different energy grid resolution factor.\n" );
				fprintf( ioQQQ, " Please recompile this file using the COMPILE STARS command "
					 "and make sure that you use the correct SET CONTINUUM RESOLUTION factor.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( strcmp(rfield.chSpType[i],"READ ") == 0 )
		{
			if( !fp_equal( rfield.RSFCheck[i], continuum.ResolutionScaleFactor ) )
			{
				fprintf( ioQQQ,"\n\n PROBLEM DISASTER The file read by the TABLE READ command "
					 "has been compiled with a different energy grid resolution factor.\n" );
				fprintf( ioQQQ, " Please recompile this file using the SAVE TRANSMITTED CONTINUUM "
					 "command and use the correct SET CONTINUUM RESOLUTION factor.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* this sanity check is that the grains we have read in from opacity files agree
	 * with the energy grid in this version of cloudy */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( !fp_equal( gv.bin[nd]->RSFCheck, continuum.ResolutionScaleFactor ) )
		{
			fprintf( ioQQQ,"\n\n PROBLEM DISASTER At least one of the grain opacity files "
				 "has been compiled with a different energy grid resolution factor.\n" );
			fprintf( ioQQQ, " Please recompile this file using the COMPILE GRAINS command "
				 "and make sure that you use the correct SET CONTINUUM RESOLUTION factor.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* default is is to predict line intensities, 
	 * but if any continuum specified as luminosity then would override this -
	 * following two values are correct for intensities */
	radius.pirsq = 0.;
	radius.lgPredLumin = false;

	/* check whether ANY luminosities are present */
	for( i=0; i < rfield.nShape; i++ )
	{
		if( strcmp(rfield.chRSpec[i],"4 PI") == 0 )
		{
			radius.pirsq = (realnum)(1.0992099 + 2.*xLog_radius_inner);
			radius.lgPredLumin = true;
			/* convert down to intensity */
			rfield.totpow[i] -= radius.pirsq;

			if( trace.lgTrace )
			{
				fprintf( ioQQQ, 
					" conorm converts continuum %ld from luminosity to intensity.\n", 
					i );
			}
		}
	}

	/* if total luminosities are present, must have specified a starting radius */
	if( radius.lgPredLumin && !radius.lgRadiusKnown )
	{
		fprintf(ioQQQ,"PROBLEM DISASTER conorm: - A continuum source was specified as a luminosity, but the inner radius of the cloud was not set.\n");
		fprintf(ioQQQ,"Please set an inner radius.\nSorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* convert ionization parameters to number of photons, called "q(h)"
	 * at this stage q(h) and "PHI " are the same */
	for( i=0; i < rfield.nShape; i++ )
	{
		if( strcmp(rfield.chSpNorm[i],"IONI") == 0 )
		{
			/* the log of the ionization parameter was stored here, this converts
			 * it to the log of the number of photons per sq cm */
			rfield.totpow[i] += log10(dense.gas_phase[ipHYDROGEN]) + log10(SPEEDLIGHT);
			strcpy( rfield.chSpNorm[i], "Q(H)" );
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, 
					" conorm converts continuum %ld from ionizat par to q(h).\n", 
					i );
			}
		}
	}

	/* convert x-ray ionization parameter xi to intensity */
	for( i=0; i < rfield.nShape; i++ )
	{
		if( strcmp(rfield.chSpNorm[i],"IONX") == 0 )
		{
			/* this converts it to an intensity */
			rfield.totpow[i] += log10(dense.gas_phase[ipHYDROGEN]) - log10(PI4);
			strcpy( rfield.chSpNorm[i], "LUMI" );
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm converts continuum%3ld from x-ray ionizat par to I.\n", 
				  i );
			}
		}
	}

	/* indicate whether we ended up with luminosity or intensity */
	if( trace.lgTrace )
	{
		if( radius.lgPredLumin )
		{
			fprintf( ioQQQ, " Cloudy will predict lumin into 4pi\n" );
		}
		else
		{
			fprintf( ioQQQ, " Cloudy will do surface flux for lumin\n" );
		}
	}

	/* if intensity per unit area is predicted then geometric
	 * covering factor must be unity
	 * variable can also be set elsewhere */
	if( !radius.lgPredLumin )
	{
		geometry.covgeo = 1.;
	}

	/* main loop over all continuum shapes to find continuum normalization 
	 * for each one */
	for( i=0; i < rfield.nShape; i++ )
	{
		rfield.ipSpec = i;

		/* check that, if laser, bounds include laser energy */
		if( strcmp(rfield.chSpType[rfield.ipSpec],"LASER") == 0 )
		{
			if( !( rfield.range[rfield.ipSpec][0] < rfield.slope[rfield.ipSpec] &&
				rfield.range[rfield.ipSpec][1] > rfield.slope[rfield.ipSpec]) )
			{
				fprintf(ioQQQ,"PROBLEM DISASTER, a continuum source is a laser at %f Ryd, but the intensity was specified over a range from %f to %f Ryd.\n",
					rfield.slope[rfield.ipSpec],
					rfield.range[rfield.ipSpec][0],
					rfield.range[rfield.ipSpec][1]);
				fprintf(ioQQQ,"Please specify the continuum flux where the laser is active.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		if( trace.lgTrace )
		{
			long int jj;
			fprintf( ioQQQ, " conorm continuum number %ld is shape %s range is %.2e %.2e\n", 
			  i, 
			  rfield.chSpType[i], 
			  rfield.range[i][0], 
			  rfield.range[i][1] );
			fprintf( ioQQQ, "the continuum points follow\n");
			jj = 0;
			if( rfield.lgContMalloc[rfield.ipSpec] )
			{
				while( rfield.tNu[rfield.ipSpec][jj].Ryd() != 0. && jj < 100 )
				{
					fprintf( ioQQQ, "%li %e %e\n",
						jj ,
						rfield.tNu[rfield.ipSpec][jj].Ryd(),
						rfield.tslop[rfield.ipSpec][jj] );
					++jj;
				}
			}
		}

		if( strcmp(rfield.chSpNorm[i],"RATI") == 0 )
		{
			/* option to scale relative to previous continua
			 * this must come first since otherwise may be too late
			 * BUT ratio cannot be the first continuum source */
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm this is ratio to 1st con\n" );
			}

			/* check that this is not the first continuum source, we must ratio */
			if( i == 0 )
			{
				fprintf( ioQQQ, " I cant form a ratio if continuum is first source.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* first find photon flux and Q of previous continuum */
			rfield.ipSpec -= 1;
			flx1 = ffun1(rfield.range[i][0])*rfield.spfac[rfield.ipSpec]*
			  rfield.range[i][0];

			/* check that previous continua were not zero where ratio is formed */
			if( flx1 <= 0. )
			{
				fprintf( ioQQQ, " Previous continua were zero where ratio is desired.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* return pointer to previous (correct) value, find F, Q */
			rfield.ipSpec += 1;

			/* we want a continuum totpow as powerful, flx is now desired flx */
			flx1 *= rfield.totpow[i];

			/*.        find flux of this new continuum at that point */
			flx2 = ffun1(rfield.range[i][1])*rfield.range[i][1];

			/* this is ratio of desired to actual */
			rfield.spfac[i] = flx1/flx2;
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm ratio will set scale fac to%10.3e at%10.2e Ryd.\n", 
				  rfield.totpow[i], rfield.range[i][0] );
			}
		}

		else if( strcmp(rfield.chSpNorm[i],"FLUX") == 0 )
		{
			/* specify flux density
			 * option to use arbitrary frequency or range */
			f = ffun1(rfield.range[i][0]); 

			/* make sure this is positive, could be zero if out of range of table,
			 * or for various forms of insanity */
			if( f<=0. )
			{
				fprintf( ioQQQ, "\n\n PROBLEM DISASTER\n The intensity of continuum source %ld is non-positive at the energy used to normalize it (%.3e Ryd).  Something is seriously wrong.\n", 
					i , rfield.range[i][0]);
				/* is this a table?  if so, is ffun within its bounds? */
				if( strcmp(rfield.chSpType[i],"INTER") == 0 )
					fprintf( ioQQQ, " This continuum shape given by a table of points - check that the intensity is specified at an energy within the range of that table.\n");
				fprintf( ioQQQ, " Also check that the numbers used to specify the shape and intensity do not under or overflow on this cpu.\n\n");

				cdEXIT(EXIT_FAILURE);
			}

			/* now convert to log in continuum units we shall use */
			f = MAX2(1e-37,f); 
			f = log10(f) + log10(rfield.range[i][0]*EN1RYD/FR1RYD);

			f = rfield.totpow[i] - f;
			/* >>chng 96 dec 31, added following test */
			if( f > 35. )
			{
				fprintf( ioQQQ, " PROBLEM DISASTER Continuum source %ld is too intense for this cpu - is it normalized correctly?\n", 
					 i );
				cdEXIT(EXIT_FAILURE);
			}

			rfield.spfac[i] = pow(10.,f);
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm will set log fnu to%10.3e at%10.2e Ryd.  Factor=%11.4e\n", 
				  rfield.totpow[i], rfield.range[i][0], rfield.spfac[i] );
			}
		}

		else if( strcmp(rfield.chSpNorm[i],"Q(H)") == 0 || 
			strcmp(rfield.chSpNorm[i],"PHI ") == 0 )
		{
			/* some type of photon density entered */
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm calling qintr range=%11.3e %11.3e desired val is %11.3e\n", 
				  rfield.range[i][0], 
				  rfield.range[i][1] ,
				  rfield.totpow[i]);
			}

			/* the total number of photons over the specified range in
			 * the arbitrary system of units that the code save the continuum shape */
			ASSERT( rfield.range[i][0] < rfield.range[i][1] );
			qentrd = qintr(&rfield.range[i][0],&rfield.range[i][1]);
			/* this is the log of the scale factor that must multiply the 
			 * continuum shape to get the final set of numbers */
			diff = rfield.totpow[i] - qentrd;

			/* >>chng 03 mar 13, from diff < -25 to <-35,
			 * tripped for very low U models used for H2 simulations */
			/*if( diff < -25. || diff > 35. )*/
			if( diff < -35. || diff > 35. )
			{
				fprintf( ioQQQ, " PROBLEM DISASTER Continuum source specified is too extreme.\n" );
				fprintf( ioQQQ, 
					" The integral over the continuum shape gave (log) %.3e photons, and the command requested (log) %.3e.\n" ,
					qentrd , rfield.totpow[i]);
				fprintf( ioQQQ, 
					" The difference in the log is %.3e.\n" ,
					diff );
				if( diff>0. )
				{
					fprintf( ioQQQ, " The continuum source is too bright.\n" );
				}
				else
				{
					fprintf( ioQQQ, " The continuum source is too faint.\n" );
				}
				/* explain what happened */
				fprintf( ioQQQ, " The usual cause for this problem is an incorrect continuum intensity/luminosity or radius command.\n" );
				fprintf( ioQQQ, " There were a total of %li continuum shape commands entered - the problem is with number %li.\n",
					rfield.nShape , i+1 );
				cdEXIT(EXIT_FAILURE);
			}

			else
			{
				rfield.spfac[i] = pow(10.,diff);
			}

			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm finds Q over range from%11.4e-%11.4e Ryd, integral= %10.4e Factor=%11.4e\n", 
				  rfield.range[i][0], 
				  rfield.range[i][1], 
				  qentrd ,
				  rfield.spfac[i] );
			}
		}

		else if( strcmp(rfield.chSpNorm[i],"LUMI") == 0 )
		{
			/* luminosity entered, special since default is TOTAL lumin */
			/*pintr integrates L for any continuum between two limits, used for normalization,
			 * return units are log of ryd cm-2 s-1, last log conv to ergs */
			pentrd = pintr(rfield.range[i][0],rfield.range[i][1]) + 
			  log10(EN1RYD);
			f = rfield.totpow[i] - pentrd;

			/* >>chng 96 dec 31, added following test */
			if( f > 35. )
			{
				fprintf( ioQQQ, " PROBLEM DISASTER Continuum source %ld is too intense for this cpu - is it normalized correctly?\n", 
					 i );
				cdEXIT(EXIT_FAILURE);
			}

			rfield.spfac[i] = pow(10.,f);
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, " conorm finds luminosity range is %10.3e to %9.3e Ryd, factor is %11.4e\n", 
				  rfield.range[i][0], rfield.range[i][1], 
				  rfield.spfac[i] );
			}
		}

		else
		{
			fprintf( ioQQQ, "PROBLEM DISASTER What chSpNorm label is this? =%s=\n", rfield.chSpNorm[i]);
			TotalInsanity();
		}

		/* sec part after .or. added June 93 because sometimes spfac=0
		 * got past first test */
		if( 1./rfield.spfac[i] == 0. || rfield.spfac[i] == 0. )
		{
			fprintf( ioQQQ, "PROBLEM DISASTER conorm finds infinite continuum scale factor.\n" );
			fprintf( ioQQQ, "The continuum is too intense to compute with this cpu.\n" );
			fprintf( ioQQQ, "Were the intensity and luminosity commands switched?\n" );
			fprintf( ioQQQ, "Sorry, but I cannot go on.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* this is conversion factor for final units of line intensities or luminosities in printout,
	 * will be intensities (==0) unless luminosity is to be printed, or flux at Earth 
	 * pirsq is the log of 4 pi r_in^2 */
	radius.Conv2PrtInten = radius.pirsq;

	/* >>chng 02 apr 25, add option for slit on aperture command */
	if( geometry.iEmissPower == 1  )
	{
		if( radius.lgPredLumin )
		{
			/* factor should be divided by 2 r_in (so that 2pi*r_in remains) */
			radius.Conv2PrtInten -= (log10(2.) + xLog_radius_inner);
		}
		else if( !radius.lgPredLumin )
		{
			/* this is an error - slit requested but radius is not known */
			fprintf( ioQQQ, "PROBLEM DISASTER conorm: Aperture slit specified, but not predicting luminosity.\n" );
			fprintf( ioQQQ, "conorm: Please specify an inner radius to determine L.\nSorry\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	if( geometry.iEmissPower == 0 && radius.lgPredLumin )
	{
		/* leave Conv2PrtInten at zero if not predicting luminosity */
		radius.Conv2PrtInten = log10(2.);
	}

	/* this is option to give final absolute results as flux observed at Earth */
	if( radius.distance > 0. && radius.lgRadiusKnown && prt.lgPrintFluxEarth )
	{
		/* this implements the conversion from Q_alpha -> F_alpha
		 * described in the section on the APERTURE command in Hazy 1 */
		if( geometry.iEmissPower == 0 )
			radius.Conv2PrtInten += log10(double(geometry.size)/SQAS_SKY);
		else if( geometry.iEmissPower == 1 )
			radius.Conv2PrtInten += log10(double(geometry.size)/(PI4*AS1RAD*radius.distance));
		else if( geometry.iEmissPower == 2 )
			radius.Conv2PrtInten -= log10( PI4*pow2(radius.distance) );
		else
			TotalInsanity();
	}

	/* normally lines are into 4pi, this is option to do per sr or arcsec^2 */
	if( prt.lgSurfaceBrightness )
	{
		if( radius.pirsq != 0. )
		{
			/* make sure we are predicting line intensities, not luminosity */
			fprintf( ioQQQ, " PROBLEM DISASTER Sorry, but both luminosity and surface brightness have been requested for lines.\n" );
			fprintf( ioQQQ, " the PRINT LINE SURFACE BRIGHTNESS command can only be used when lines are predicted per unit cloud area.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if( prt.lgSurfaceBrightness_SR )
		{
			/* we want final units to be per sr */
			radius.Conv2PrtInten -= log10( PI4 );
		}
		else
		{
			/* we want final units to be per square arcsec */
			radius.Conv2PrtInten -= log10(SQAS_SKY);
		}
	}
	return;
}

/*qintr integrates Q for any continuum between two limits, used for normalization */
STATIC double qintr(double *qenlo, 
  double *qenhi)
{
	long int i, 
	  ipHi, 
	  ipLo, 
	  j;
	double qintr_v, 
	  sum, 
	  wanu;

	DEBUG_ENTRY( "qintr()" );

	/* returns LOG of number of photons over energy interval */

	/* this is copy of logic that occurs three times across code */
	ASSERT(*qenhi > *qenlo);
	ipLo = ipoint(*qenlo);
	ipHi = ipoint(*qenhi);
	/* this is actual sum of photons within band */
	sum = 0.;
	for( i=ipLo-1; i < (ipHi - 1); i++ )
	{
		/*sum += ffun1(rfield.anu[i])*rfield.widflx[i];*/
		for( j=0; j < 4; j++ )
		{
			wanu = rfield.anu[i] + rfield.widflx[i]*aweigh[j];
			/* >>chng 02 jul 16, add test on continuum limits -
			* this was exceeded when resolution set very large */
			wanu = MAX2( wanu , rfield.emm );
			wanu = MIN2( wanu , rfield.egamry );
			sum += fweigh[j]*ffun1(wanu)*rfield.widflx[i];
		}
	}

	if( sum <= 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Photon number sum in QINTR is %.3e\n", 
		  sum );
		fprintf( ioQQQ, " This source has no ionizing radiation, and the number of ionizing photons was specified.\n" );
		fprintf( ioQQQ, " This was continuum source number%3ld\n", 
		  rfield.ipSpec );
		fprintf( ioQQQ, " Sorry, but I cannot go on.  ANU and FLUX arrays follow.  Enjoy.\n" );
		fprintf( ioQQQ, "\n\n This error is also caused by an old table read file whose energy mesh does not agree with the code.\n" );
		for( i=0; i < rfield.nupper; i++ )
		{
			fprintf( ioQQQ, "%.2e\t%.2e\n", 
				rfield.anu[i], 
				rfield.flux[0][i] );
		}
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		qintr_v = log10(sum);
	}
	return( qintr_v );
}

/*pintr integrates L for any continuum between two limits, used for normalization,
 * return units are log of ryd cm-2 s-1 */
STATIC double pintr(double penlo, 
  double penhi)
{
	long int i, 
	  j;
	double fsum, 
	  pintr_v, 
	  sum, 
	  wanu, 
	  wfun;
	long int ip1 , ip2;

	DEBUG_ENTRY( "pintr()" );

	/* computes log of luminosity in radiation over some integral
	 * answer is in Ryd per sec */

	sum = 0.;
	/* >>chng 02 oct 27, do not call qg32, do same type sum as
		* final integration */
	/* laser is special since delta function, this is center of laser */
	/* >>chng 01 jul 01, was +-21 cells, change to call to ipoint */
	ip1 = ipoint( penlo );

	ip2 = ipoint( penhi );

	for( i=ip1-1; i < ip2-1; i++ )
	{
		fsum = 0.;
		for( j=0; j < 4; j++ )
		{
			wanu = rfield.anu[i] + rfield.widflx[i]*aweigh[j];
			/*++iiii;fprintf(ioQQQ,"DEBUG iii %li %e \n",iiii, wanu );*/
			wfun = fweigh[j]*ffun1(wanu)*wanu;
			fsum += wfun;
		}
		sum += fsum*rfield.widflx[i];
	}

	if( sum > 0. )
	{
		pintr_v = log10(sum);
	}
	else
	{
		pintr_v = -38.;
	}

	return( pintr_v );
}
