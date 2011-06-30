/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*esc_CRDwing_1side fundamental escape probability radiative transfer routine, for complete redistribution */
/*esc_PRD_1side fundamental escape probability radiative transfer routine for incomplete redistribution */
/*RTesc_lya escape prob for hydrogen atom Lya, using Hummer and Kunasz results,
 * called by hydropesc */
/*esc_PRD escape probability radiative transfer for incomplete redistribution */
/*esc_CRDwing escape probability for CRD with wings */
/*esc_CRDcore escape probability for CRD with no wings */
/*esca0k2 derive Hummer's K2 escape probability for Doppler core only */
/*RT_DestProb returns line destruction probability due to continuum opacity */
/*RT_LineWidth determine half width of any line with known optical depths */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "conv.h"
#include "rfield.h"
#include "opacity.h"
#include "lines_service.h"
#include "taulines.h"
#include "doppvel.h"
#include "pressure.h"
#include "wind.h"
#include "rt.h"

static const double SCALE = 2.;

/*escConE2 one of the forms of the continuum escape probability */
class my_Integrand_escConE2
{
public:
	double chnukt_ContTkt, chnukt_ctau;

	double operator() (double x)
	{
		return exp(-chnukt_ContTkt*(x-1.))/x*e2(chnukt_ctau/POW3(x));
	}
};

/* a continuum escape probability */
class my_Integrand_conrec
{
public:
	double chnukt_ContTkt;

	double operator() (double x)
	{
		return exp(-chnukt_ContTkt*(x-1.))/x;
	}
};


/*escmase escape probability for negative (masing) optical depths,*/
STATIC double escmase(double tau);
/*RTesc_lya_1side fit Hummer and Kunasz escape probability for hydrogen atom Lya */
STATIC void RTesc_lya_1side(double taume, 
  double beta, 
  realnum *esc, 
  realnum *dest,
  /* position of line in frequency array on c scale */
  long ipLine );

double esc_PRD_1side(double tau, 
  double a)
{
	double atau, 
	  b, 
	  escinc_v;

	DEBUG_ENTRY( "esc_PRD_1side()" );
	ASSERT( a>0.0 );

	/* this is one of the three fundamental escape probability routines
	 * the three are esc_CRDwing_1side, esc_PRD_1side, and RTesc_lya
	 * it computes esc prob for incomplete redistribution
	 * */
#	if 0
	if( strcmp(rt.chEscFunSubord,"SIMP") == 0 )
	{
		/* this set with "escape simple" command, used for debugging */
		escinc_v = 1./(1. + tau);
		return escinc_v;
	}
#	endif

	if( tau < 0. )
	{
		/* line mased */
		escinc_v = escmase(tau);
	}
	else if( tau < 10. )
	{
		/* linear part of doppler core */
		escinc_v = 1./(1. + 1.6*tau);
	}
	else
	{
		/* first find coefficient b(tau) */
		atau = a*tau;
		if( atau > 1. )
		{
			b = 1.6 + (3.*pow(2.*a,-0.12))/(1. + atau);
		}
		else
		{
			b = 1.6 + (3.*pow(2.*a,-0.12))/(1. + 1./sqrt(atau));
		}
		b = MIN2(6.,b);

		escinc_v = 1./(1. + b*tau);
	}
	return escinc_v;
}

/*esc_CRDwing_1side fundamental escape probability radiative transfer routine, for complete redistribution */
double esc_CRDwing_1side(double tau, 
  double a )
{
	double esccom_v;

	DEBUG_ENTRY( "esc_CRDwing_1side()" );

	/* this is one of the three fundamental escape probability routines
	 * the three are esc_CRDwing_1side, esc_PRD_1side, and RTesc_lya
	 * it computes esc prob for complete redistribution with wings
	 * computes escape prob for complete redistribution in one direction
	 * */

	/* this is the only case that this routine computes,
	 * and is the usual case for subordinate lines, 
	 * complete redistribution with damping wings */
	esccom_v = esca0k2(tau);
	if( tau > 1e3 )
	{
		esccom_v += 0.333*sqrt(a/(SQRTPI*tau));
	}
	return esccom_v;
}

/*RTesc_lya escape prob for hydrogen atom Lya, using 
 >>refer	La	escp	Hummer, D.G., & Kunasz, P.B., 1980, ApJ, 236, 609
 * called by hydropesc, return value is escape probability */
double RTesc_lya(
	/* the inward escape probability */
	double *esin, 
	/* the destruction probility */
	double *dest, 
	/* abundance of the species */
	double abund, 
	/* element number,  0 for H */
	long int nelem, 
	realnum DopplerWidth)
{
	double beta, 
	  conopc, 
	  escla_v;
	 realnum dstin, 
	  dstout;

	DEBUG_ENTRY( "RTesc_lya()" );

	/* 
	 * this is one of the three fundamental escape probability functions
	 * the three are esc_CRDwing_1side, esc_PRD_1side, and RTesc_lya
	 * evaluate esc prob for LA
	 * optical depth in outer direction always defined
	 */

	/* check charge */
	ASSERT( nelem >= 0 && nelem < LIMELM );

	if( Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauTot - 
	  Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn < 0. )
	{
		/* this is the case if we overrun the optical depth scale
		 * just leave things as they are */
		escla_v = Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->Pesc;
		rt.fracin = Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->FracInwd;
		*esin = rt.fracin;
		*dest = Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->Pdest;
		return escla_v;
	}

	/* incomplete redistribution */
	conopc = opac.opacity_abs[Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].ipCont-1];
	if( abund > 0. )
	{
		/* the continuous opacity is positive, we have a valid soln */
		beta = conopc/(abund/SQRTPI*Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->opacity/
		  DopplerWidth + conopc);
	}
	else
	{
		/* abundance is zero, set miniumum dest prob */
		beta = 1e-10;
	}

	/* find rt.wayin, the escape prob in inward direction */
	RTesc_lya_1side(
	  Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn,
	  beta,
	  &rt.wayin,
	  &dstin , 
	  /* position of line in energy array on C scale */
	  Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].ipCont-1);

	ASSERT( (rt.wayin <= 1.) && (rt.wayin >= 0.) && (dstin <= 1.) && (dstin >= 0.) );

	/* find rt.wayout, the escape prob in outward direction */
	RTesc_lya_1side(MAX2(Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauTot/100.,
	  Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauTot-Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn),
	  beta,
	  &rt.wayout,
	  &dstout, 
	  Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].ipCont-1);

	ASSERT( (rt.wayout <= 1.) && (rt.wayout >= 0.) && (dstout <= 1.) && (dstout >= 0.) );

	/* esc prob is mean of in and out */
	escla_v = (rt.wayin + rt.wayout)/2.;
	/* the inward escaping part of the line */
	*esin = rt.wayin;

	/* dest prob is mean of in and out */
	*dest = (dstin + dstout)/2.f;
	/* >>chng 02 oct 02, sum of escape and dest prob must be less then unity,
	 * for very thin models this forces dest prob to go to zero, 
	 * rather than the value of DEST0, caught by Jon Slavin */
	*dest = (realnum)MIN2( *dest , 1.-escla_v );
	/* but dest prob can't be negative */
	*dest = (realnum)MAX2(0., *dest );

	/* fraction of line emitted in inward direction */
	rt.fracin = rt.wayin/(rt.wayin + rt.wayout);
	ASSERT( escla_v >=0. && *dest>=0. && *esin>=0. );
	return escla_v;
}

/*esc_PRD escape probability radiative transfer for incomplete redistribution */
double esc_PRD(double tau, 
  double tau_out, 
  double damp )
{
	double escgrd_v, 
	  tt;

	DEBUG_ENTRY( "esc_PRD()" );

	ASSERT( damp > 0. );

	/* find escape prob for incomp redis, average of two 1-sided probs*/

	if( iteration > 1 )
	{
		/*  outward optical depth if defined */
		tt = tau_out - tau;
		/*  help convergence by not letting tau go to zero at back edge of
		 *  when there was a bad guess for the total optical depth
		 *  note that this test is seldom hit since RTMakeStat does check
		 *  for overrun */
		if( tt < 0. )
		{
			tt = tau/SCALE;
		}

		rt.wayin = (realnum)esc_PRD_1side(tau,damp);
		rt.wayout = (realnum)esc_PRD_1side(tt,damp);
		rt.fracin = rt.wayin/(rt.wayin + rt.wayout);
		escgrd_v = 0.5*(rt.wayin + rt.wayout);
	}
	else
	{
		/*  outward optical depth not defined, dont estimate fraction out */
		rt.fracin = 0.;
		rt.wayout = 1.;
		escgrd_v = esc_PRD_1side(tau,damp);
		rt.wayin = (realnum)escgrd_v;
	}

	ASSERT( escgrd_v > 0. );
	return escgrd_v;
}

/*esc_CRDwing escape probability radiative transfer for CRDS in core only */
double esc_CRDwing(double tau_in, 
  double tau_out, 
  double damp)
{
	double escgrd_v, 
	  tt;

	DEBUG_ENTRY( "esc_CRDwing()" );

	/* find escape prob for CRD with damping wings, average of two 1-sided probs*/

	/* crd with wings */
	if( iteration > 1 )
	{
		/*  outward optical depth if defined */
		/* >>chng 03 jun 07, add test for masers here */
		if( tau_out <0 || tau_in < 0. )
		{
			/* we have a maser, use smallest optical depth to damp it out */
			tt = MIN2( tau_out , tau_in );
			tau_in = tt;
		}
		else
		{
			tt = tau_out - tau_in;
			/*  help convergence by not letting tau_in go to zero at back edge of
			 *  when there was a bad guess for the total optical depth
			 *  note that this test is seldom hit since RTMakeStat does check
			 *  for overrun */
			if( tt < 0. )
			{
				tt = tau_in/SCALE;
			}
		}

		rt.wayin = (realnum)esc_CRDwing_1side(tau_in,damp);
		rt.wayout = (realnum)esc_CRDwing_1side(tt,damp);
		rt.fracin = rt.wayin/(rt.wayin + rt.wayout);
		escgrd_v = 0.5*(rt.wayin + rt.wayout);
	}
	else
	{
		/*  outward optical depth not defined, dont estimate fraction out */
		rt.fracin = 0.;
		rt.wayout = 1.;
		escgrd_v = esc_CRDwing_1side(tau_in,damp);
		rt.wayin = (realnum)escgrd_v;
	}

	ASSERT( escgrd_v > 0. );
	return escgrd_v;
}

/*esc_CRDwing escape probability radiative transfer for incomplete redistribution */
double esc_CRDcore(double tau_in, 
  double tau_out)
{
	double escgrd_v, 
	  tt;

	DEBUG_ENTRY( "esc_CRDcore()" );

	/* find escape prob for CRD with damping wings, average of two 1-sided probs*/

	/* crd with wings */
	if( iteration > 1 )
	{
		/*  outward optical depth if defined */
		/* >>chng 03 jun 07, add test for masers here */
		if( tau_out <0 || tau_in < 0. )
		{
			/* we have a maser, use smallest optical depth to damp it out */
			tt = MIN2( tau_out , tau_in );
			tau_in = tt;
		}
		else
		{
			tt = tau_out - tau_in;
			/*  help convergence by not letting tau_in go to zero at back edge of
			 *  when there was a bad guess for the total optical depth
			 *  note that this test is seldom hit since RTMakeStat does check
			 *  for overrun */
			if( tt < 0. )
			{
				tt = tau_in/SCALE;
			}
		}

		rt.wayin = (realnum)esca0k2(tau_in);
		rt.wayout = (realnum)esca0k2(tt);
		rt.fracin = rt.wayin/(rt.wayin + rt.wayout);
		escgrd_v = 0.5*(rt.wayin + rt.wayout);
	}
	else
	{
		/*  outward optical depth not defined, dont estimate fraction out */
		rt.fracin = 0.;
		rt.wayout = 1.;
		escgrd_v = esca0k2(tau_in);
		rt.wayin = (realnum)escgrd_v;
	}

	ASSERT( escgrd_v > 0. );
	return escgrd_v;
}

/*esca0k2 derive Hummer's K2 escape probability for Doppler core only */
double esca0k2(double taume)
{
	double arg, 
	  esca0k2_v, 
	  suma, 
	  sumb, 
	  sumc, 
	  sumd, 
	  tau;
	static const double a[5]={1.00,-0.1117897,-0.1249099917,-9.136358767e-3,
	  -3.370280896e-4};
	static const double b[6]={1.00,0.1566124168,9.013261660e-3,1.908481163e-4,
	  -1.547417750e-7,-6.657439727e-9};
	static const double c[5]={1.000,19.15049608,100.7986843,129.5307533,-31.43372468};
	static const double d[6]={1.00,19.68910391,110.2576321,169.4911399,-16.69969409,
	  -36.664480000};

	DEBUG_ENTRY( "esca0k2()" );

	/* compute Hummer's K2 escape probability function for a=0
	 * using approx from 
	 * >>refer	line	escp	Hummer, D.G., xxxx, JQRST, 26, 187.
	 *
	 * convert to David's opacity */
	tau = taume*SQRTPI;

	if( tau < 0. )
	{
		/* the line mased */
		esca0k2_v = escmase(taume);

	}
	else if( tau < 0.01 )
	{
		esca0k2_v = 1. - 2.*tau;

	}
	else if( tau <= 11. )
	{
		suma = a[0] + tau*(a[1] + tau*(a[2] + tau*(a[3] + a[4]*tau)));
		sumb = b[0] + tau*(b[1] + tau*(b[2] + tau*(b[3] + tau*(b[4] + 
		  b[5]*tau))));
		esca0k2_v = tau/2.5066283*log(tau/SQRTPI) + suma/sumb;

	}
	else
	{
		/* large optical depth limit */
		arg = 1./log(tau/SQRTPI);
		sumc = c[0] + arg*(c[1] + arg*(c[2] + arg*(c[3] + c[4]*arg)));
		sumd = d[0] + arg*(d[1] + arg*(d[2] + arg*(d[3] + arg*(d[4] + 
		  d[5]*arg))));
		esca0k2_v = (sumc/sumd)/(2.*tau*sqrt(log(tau/SQRTPI)));
	}
	return esca0k2_v;
}

/*escmase escape probability for negative (masing) optical depths */
STATIC void FindNeg( void )
{
	long int i;

	DEBUG_ENTRY( "FindNeg()" );

	/* do the level 1 lines */
	for( i=1; i <= nLevel1; i++ )
	{
		/* check if a line was a strong maser */
		if( TauLines[i].Emis->TauIn < -1. )
			DumpLine(&TauLines[i]);
	}

	/* Generic atoms & molecules from databases 
	 * added by Humeshkar Nemala*/
	for( i = 0; i <linesAdded2; i++)
	{
		/* check if a line was a strong maser */
		if(dBaseLines[i].TauIn < -1. )
			DumpLine(dBaseLines[i].tran);
	}

	/* now do the level 2 lines */
	for( i=0; i < nWindLine; i++ )
	{
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			/* check if a line was a strong maser */
			if( TauLine2[i].Emis->TauIn < -1. )
				DumpLine(&TauLine2[i]);
		}
	}

	/* now do the hyperfine structure lines */
	for( i=0; i < nHFLines; i++ )
	{
		/* check if a line was a strong maser */
		if( HFLines[i].Emis->TauIn < -1. )
			DumpLine(&HFLines[i]);
	}

	return;
}

STATIC double escmase(double tau)
{
	double escmase_v;

	DEBUG_ENTRY( "escmase()" );

	/* this is the only routine that computes maser escape probabilities */
	ASSERT( tau <= 0. );

	if( tau > -0.1 )
	{
		escmase_v = 1. - tau*(0.5 + tau/6.);
	}
	else if( tau > -30. )
	{
		escmase_v = (1. - exp(-tau))/tau;
	}
	else
	{
		fprintf( ioQQQ, " DISASTER escmase called with 2big tau%10.2e\n", 
		  tau  );
		fprintf( ioQQQ, " This is zone number%4ld\n", nzone );
		FindNeg();
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( escmase_v >= 1. );
	return escmase_v;
}

/*esccon continuum escape probability */
double esccon(double tau, 
	      double hnukt)
{
	double dinc, 
	  escpcn_v, 
	  sumesc, 
	  sumrec;

	DEBUG_ENTRY( "esccon()" );

	/* computes continuum escape probabilities */
	if( tau < 0.01 )
	{
		escpcn_v = 1.;
		return escpcn_v;
	}

	else if( hnukt > 1. && tau > 100. )
	{
		escpcn_v = 1e-20;
		return escpcn_v;
	}

	my_Integrand_conrec func_conrec;
	func_conrec.chnukt_ContTkt = hnukt;
	Integrator<my_Integrand_conrec,Gaussian32> conrec;

	my_Integrand_escConE2 func_escConE2;
	func_escConE2.chnukt_ContTkt = hnukt;
	func_escConE2.chnukt_ctau = tau;
	Integrator<my_Integrand_escConE2,Gaussian32> escConE2;

	dinc = 10./hnukt;
	sumrec = conrec.sum(1.,1.+dinc,func_conrec);
	sumesc = escConE2.sum(1.,1.+dinc,func_escConE2);

	if( sumrec > 0. )
	{
		escpcn_v = sumesc/sumrec;
	}
	else
	{
		escpcn_v = 0.;
	}
	return escpcn_v;
}

/*RTesc_lya_1side fit Hummer and Kunasz escape probability for hydrogen atom Lya */
STATIC void RTesc_lya_1side(double taume, 
  double beta, 
  realnum *esc, 
  realnum *dest,
  /* position of line in frequency array on c scale */
  long ipLine )
{
	double esc0, 
	  fac, 
	  fac1, 
	  fac2, 
	  tau, 
	  taucon, 
	  taulog;

	/* DEST0 is the smallest destruction probability to return
	 * in high metallicity models, in rt.h
	const double DEST0=1e-8;*/

	DEBUG_ENTRY( "RTesc_lya_1side()" );

	/* fits to numerical results of Hummer and Kunasz Ap.J. 80 */
	tau = taume*SQRTPI;

	/* this is the real escape probability */
	esc0 = 1./((0.6451 + tau)*(0.47 + 1.08/(1. + 7.3e-6*tau)));

	esc0 = MAX2(0.,esc0);
	esc0 = MIN2(1.,esc0);

	if( tau > 0. )
	{
		taulog = log10(MIN2(1e8,tau));
	}
	else
	{
		/* the line mased 
		 *>>chng 06 sep 08, kill xLaMase 
		hydro.xLaMase = MIN2(hydro.xLaMase,(realnum)tau);*/
		taulog = 0.;
		*dest = 0.;
		*esc = (realnum)esc0;
	}

	if( beta > 0. )
	{
		taucon = MIN2(2.,beta*tau);

		if( taucon > 1e-3 )
		{
			fac1 = -1.25 + 0.475*taulog;
			fac2 = -0.485 + 0.1615*taulog;
			fac = -fac1*taucon + fac2*POW2(taucon);
			fac = pow(10.,fac);
			fac = MIN2(1.,fac);
		}
		else
		{
			fac = 1.;
		}

		*esc = (realnum)(esc0*fac);
		/* MIN puts cat at 50 */
		*dest = (realnum)(beta/(0.30972 - MIN2(.28972,0.03541667*taulog)));
	}

	else
	{
		*dest = 0.;
		*esc = (realnum)esc0;
	}

	*dest = MIN2(*dest,1.f-*esc);
	*dest = MAX2(0.f,*dest);

	/* >>chng 99 apr 12, limit destruction prob in case where gas dominated by scattering.
	 * in this case scattering is much more likely than absorption on this event */
	*dest = (realnum)( (1. - opac.albedo[ipLine]) * *dest + opac.albedo[ipLine]*DEST0);
	/* this is for debugging H Lya */
	{
		/*@-redef@*/
		enum {BUG=false};
		/*@+redef@*/
		if( BUG )
		{
			fprintf(ioQQQ,"scatdest tau %.2e beta%.2e 1-al%.2e al%.2e dest%.2e \n",
			taume,
			beta, 
			(1. - opac.albedo[ipLine]), 
			opac.albedo[ipLine] ,
			*dest 
			);
		}
	}
	return;
}

/*RT_DestProb returns line destruction probability due to continuum opacity */
double RT_DestProb(
	  /* abundance of species */
	  double abund, 
	  /* its line absorption cross section */
	  double crsec, 
	  /* pointer to energy within continuum array, to get background opacity,
	   * this is on the f not c scale */
	  long int ipanu, 
	  /* line width */
	  double widl, 
	  /* escape probability */
	  double escp, 
	  /* type of redistribution function */
	  int nCore)
{
	/* this will be the value we shall return */
	double eovrlp_v;

	double conopc, 
	  beta;

	/* DEST0 is the smallest destruction probability to return
	 * in high metallicity models 
	 * this was set to 1e-8 until 99nov18,
	 * in cooling flow model the actual Lya ots dest prob was 1e-16,
	 * and this lower limit of 1e-8 caused energy balance problems,
	 * since dest prob was 8 orders of magnitude too great.  
	 * >>chng 99 nov 18, to 1e-20, but beware that comments indicate that
	 * this will cause problems with high metallicity clouds(?) */
	/* >>chng 00 jun 04, to 0 since very feeble ionization clouds, with almost zero opacity,
	 * this was a LARGE number */
	/*const double DEST0=1e-20;
	const double DEST0=0.;*/

	DEBUG_ENTRY( "RT_DestProb()" );

	/* computes "escape probability" due to continuum destruction of
	 *
	 * if esc prob gt 1 then line is masing - return small number for dest prob */
	/* >>>chng 99 apr 10, return min dest when scattering greater than abs */
	/* no idea of opacity whatsoever, on very first soln for this model */
	/* >>chng 05 mar 20, add test on line being above upper bound of frequency 
	 * do not want to evaluate opacity in this case since off scale */
	if( escp >= 1.0 || !conv.nTotalIoniz || ipanu >= rfield.nflux )
	{
		eovrlp_v = 0.;
		return eovrlp_v;
	}

	/* find continuum opacity */
	conopc = opac.opacity_abs[ipanu-1];

	ASSERT( crsec > 0. );

	/* may be no population, cannot use above test since return 0 not DEST0 */
	if( abund <= 0. || conopc <= 0. )
	{
		/* do not set this to DEST0 since energy not then conserved */
		eovrlp_v = 0.;
		return eovrlp_v;
	}

	/* fac of 1.7 convert to Hummer convention for line opacity */
	beta = conopc/(abund*SQRTPI*crsec/widl + conopc);
	/* >>chng 04 may 10, rm * 1-pesc)
	beta = MIN2(beta,(1.-escp)); */

	if( nCore == ipDEST_INCOM )
	{
		/*  fits to 
		 *  >>>refer	la	esc	Hummer and Kunasz 1980 Ap.J. 236,609.
		 *  the max value of 1e-3 is so that we do not go too far
		 *  beyond what Hummer and Kunasz did, discussed in
		 * >>refer	rt	esc proc	Ferland, G.J., 1999, ApJ, 512, 247 */
		/** \todo	2	this min is because there are no calculations that show what to do
		 * for beta beyound this value */
		eovrlp_v = MIN2(1e-3,8.5*beta);/**/
	}
	else if( nCore == ipDEST_K2 )
	{
		/*  Doppler core only; a=0., Hummer 68 
		eovrlp_v = RT_DestHummer(beta);*/
		eovrlp_v = MIN2(1e-3,8.5*beta);/**/
	}
	else if( nCore == ipDEST_SIMPL )
	{
		/*  this for debugging only 
		eovrlp_v = 8.5*beta;*/
		/* >>chng 04 may 13, use same min function */
		eovrlp_v = MIN2(1e-3,8.5*beta);/**/
	}
	else
	{
		fprintf( ioQQQ, " chCore of %i not understood by RT_DestProb.\n", 
		  nCore );
		cdEXIT(EXIT_FAILURE);
	}

	/* renorm to unity */
	eovrlp_v /= 1. + eovrlp_v;

	/* multiply by 1-escape prob, since no destruction when optically thin */
	eovrlp_v *= 1. - escp;

	/*check results in bounds */
	ASSERT( eovrlp_v >= 0.  );
	ASSERT( eovrlp_v <= 1.  );

	{
		/* debugging code for Lya problems */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			if( rfield.anu[ipanu-1]>0.73 && rfield.anu[ipanu-1]<0.76 &&
			    fp_equal( abund, Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->PopOpc ) )
			{
				fprintf(ioQQQ,"%li RT_DestProb\t%g\n",
					nzone, eovrlp_v  );
			}
		}
	}

	/* >>chng 04 may 10, rm min */
	/* this hack removed since no fundamental reason for it to be here,
	 * this should be added to scattering escape, if included at all */
#	if 0
	/* >>chng 99 apr 12, limit destruction prob in case where gas dominated by scattering.
	 * in this case scattering is much more likely than absorption on this event
	eovrlp_v = (1. - opac.albedo[ipanu-1]) * eovrlp_v + 
		opac.albedo[ipanu-1]*DEST0; */
	/* >>chng 01 aug 11, add factor of 3 for increase in mean free path, and min on 0 */
	/*eovrlp_v = MAX2(DEST0,1. - 3.*opac.albedo[ipanu-1]) * eovrlp_v + 
		opac.albedo[ipanu-1]*DEST0;*/
	eovrlp_v = POW2(1. - opac.albedo[ipanu-1]) * eovrlp_v + 
		opac.albedo[ipanu-1]*0.;
#	endif

	return eovrlp_v;
}

/*RT_LineWidth compute line width (cm/sec), using optical depth array information 
 * this is where effects of wind are done */
double RT_LineWidth(const transition * t, realnum DopplerWidth)
{
	double RT_LineWidth_v, 
	  aa, 
	  atau, 
	  b, 
	  r, 
	  vth;
	realnum tau; 

	DEBUG_ENTRY( "RT_LineWidth()" );

	/* uses line width from 
	 * >>refer	esc	prob	Bonilha et al. Ap.J. (1979) 233 649
	 * return value is half velocity width*(1-ESC PROB) [cm s-1]
	 * this assumes incomplete redistribution, damp.tau^1/3 width */

	/* thermal width */
	vth = DopplerWidth;

	/* optical depth in outer direction is defined
	 * on second and later iterations. 
	 * smaller of inner and outer optical depths is chosen for esc prob */
	if( iteration > 1 )
	{
		/* optical depth to shielded face */
		realnum tauout = t->Emis->TauTot - t->Emis->TauIn;

		/* >>chng 99 apr 22 use smaller optical depth */
		tau = MIN2( t->Emis->TauIn , tauout );
	}
	else
	{
		tau = t->Emis->TauIn;
	}
	/* do not evaluate line width if quite optically thin - will be dominated
	 * by noise in this case */
	if( tau <1e-3 )
		return 0;

	t->Emis->damp = t->Emis->dampXvel / DopplerWidth;
	ASSERT( t->Emis->damp > 0. );

	double Pesc =  esc_PRD_1side( tau , t->Emis->damp);

	/* max optical depth is thermalization length */
	ASSERT( 5.3e16/dense.eden < BIGFLOAT );
	realnum therm = (realnum)(5.3e16/dense.eden);
	if( tau > therm )
	{
		/* \todo 2 this seems to create an inconsistency as it changes tau
		 * for the purposes of this routine (to return the line-width),
		 * but this leaves the actual optical depth unchanged. */
		pressure.lgPradDen = true;
		tau = therm;
	}

	/* >>chng 01 jun 23, use wind vel instead of rt since rt deleted */
	/* >>chng 04 may 13, use thermal for subsonic cases */
	/** \todo	1	dynamics; this test assumes that neg vel are subsonic, so that sobolev length
	 * would overestimate the optical depth, since ion is at most present over computed
	 * slab, and possibly more.  */
	/** \todo	1	rewrite so that this checks on size not sign of windv */
	if( ! wind.lgBallistic() )
	{
		/* static geometry */
		/* esc prob has noise if smaller than FLT_EPSILON, or is masing */
		if( (tau-opac.taumin)/100. < FLT_EPSILON )
		{
			RT_LineWidth_v = 0.;
		}
		else if( tau <= 20. )
		{
			atau = -6.907755;
			if( tau > 1e-3 )
 				atau = log(tau);
			aa = 4.8 + 5.2*tau + (4.*tau - 1.)*atau;
			b = 6.5*tau - atau;
			double escProb = Pesc + t->Emis->Pelec_esc + t->Emis->Pdest;
			RT_LineWidth_v = vth*0.8862*aa/b*(1. - MIN2( 1., escProb ) ); 
			/* small number roundoff can dominate this process */
			if( escProb >= 1. - 100. * FLT_EPSILON )
				RT_LineWidth_v = 0.;
		}
		else
		{
			ASSERT( t->Emis->damp*tau >= 0.);
			atau = log(MAX2(0.0001,tau));
			aa = 1. + 2.*atau/pow(1. + 0.3*t->Emis->damp*tau,0.6667) + pow(6.5*
			  t->Emis->damp*tau,0.333);
			b = 1.6 + 1.5/(1. + 0.20*t->Emis->damp*tau);
			RT_LineWidth_v = vth*0.8862*aa/b*(1. - MIN2( 1. , 
				(Pesc+ t->Emis->Pelec_esc + t->Emis->Pdest)) );
		}

		/* we want full width, not half width */
		RT_LineWidth_v *= 2.;

	}
	else
	{
		/* ballistic wind */
		r = t->Emis->damp*tau/PI;
		if( r <= 1. )
		{
			RT_LineWidth_v = vth*sqrt(log(MAX2(1.,tau))*.2821);
		}
		else
		{
			RT_LineWidth_v = 2.*fabs(wind.windv0);
			if( r*vth <= RT_LineWidth_v )
			{
				RT_LineWidth_v = vth*r*log(RT_LineWidth_v/(r*vth));
			}
		}
	}

	ASSERT( RT_LineWidth_v >= 0. );
	return RT_LineWidth_v;
}

/*RT_DestHummer evaluate Hummer's betaF(beta) function */
double RT_DestHummer(double beta) /* beta is ratio of continuum to mean line opacity,
									 * returns dest prob = beta F(beta) */
{
	double fhummr_v, 
	  x;

	DEBUG_ENTRY( "RT_DestHummer()" );

	/* evaluates Hummer's F(beta) function for case where damping
	 * constant is zero, are returns beta*F(beta)
	 * fit to Table 1, page 80, of Hummer MNRAS 138, 73-108.
	 * beta is ratio of continuum to line opacity; FUMMER is
	 * product of his F() times beta; the total destruction prob
	 * this beta is Hummer's normalization of the Voigt function */

	ASSERT( beta >= 0.);/* non-positive is unphysical */
	if( beta <= 0. )
	{
		fhummr_v = 0.;
	}
	else
	{
		x = log10(beta);
		if( x < -5.5 )
		{
			fhummr_v = 3.8363 - 0.56329*x;
		}
		else if( x < -3.5 )
		{
			fhummr_v = 2.79153 - 0.75325*x;
		}
		else if( x < -2. )
		{
			fhummr_v = 1.8446 - 1.0238*x;
		}
		else
		{
			fhummr_v = 0.72500 - 1.5836*x;
		}
		fhummr_v *= beta;
	}
	return fhummr_v;
}
