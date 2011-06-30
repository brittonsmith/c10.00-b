/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*rt_continuum_shield_fcn computing continuum shielding due to single line */
/*conpmp local continuum pumping rate radiative transfer for all lines */
/*conpmp local continuum pumping rate radiative transfer for all lines */
#include "cddefines.h"
#include "rt.h"
#include "transition.h"

/*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double conpmp(transition * t);
STATIC inline double FITTED( double t );

class my_Integrand
{
public:
	double PumpDamp, PumpTau;

	double operator() (double x)
	{
		double v = vfun(PumpDamp,x);
		double opfun_v = sexp(PumpTau*v)*v;
		return opfun_v;
	}
};


/*rt_continuum_shield_fcn computing continuum shielding due to single line */
double RT_continuum_shield_fcn( transition *t )
{
	double value;

	DEBUG_ENTRY( "rt_continuum_shield_fcn()" );

	value = -1.;

	ASSERT( t->Emis->damp > 0. );

	if( rt.nLineContShield == LINE_CONT_SHIELD_PESC )
	{
		/* set continuum shielding pesc - shieding based on escape probability */
		if( t->Emis->iRedisFun == ipPRD )
		{
			value =  esc_PRD_1side(t->Emis->TauCon,t->Emis->damp);
		}
		else if( t->Emis->iRedisFun == ipCRD )
		{
			value =  esca0k2(t->Emis->TauCon);
		}
		else if( t->Emis->iRedisFun == ipCRDW )
		{
			value =  esc_CRDwing_1side(t->Emis->TauCon,t->Emis->damp);
		}
		else if( t->Emis->iRedisFun == ipLY_A )
		{
			value = esc_PRD_1side(t->Emis->TauCon,t->Emis->damp);
		}
		else
			TotalInsanity();
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_FEDERMAN )
	{
		/* set continuum shielding Federman - this is the default */
		double core, wings;

		/* these expressions implement the appendix of
		 * >>refer	line	shielding	Federman, S.R., Glassgold, A.E., & 
		 * >>refercon	Kwan, J. 1979, ApJ, 227, 466 */
		/* doppler core - equation A8 */
		if( t->Emis->TauCon < 2. )
		{
			core = sexp( t->Emis->TauCon * 0.66666 );
		}
		else if( t->Emis->TauCon < 10. )
		{
			core = 0.638 * pow(t->Emis->TauCon,(realnum)-1.25f );
		}
		else if( t->Emis->TauCon < 100. )
		{
			core = 0.505 * pow(t->Emis->TauCon,(realnum)-1.15f );
		}
		else
		{
			core = 0.344 * pow(t->Emis->TauCon,(realnum)-1.0667f );
		}

		/* do we add damping wings? */
		wings = 0.;
		if( t->Emis->TauCon > 1.f && t->Emis->damp>0. )
		{
			/* equation A6 */
			double t1 = 3.02*pow(t->Emis->damp*1e3,-0.064 );
			double u1 = sqrt(t->Emis->TauCon*t->Emis->damp )/SDIV(t1);
			wings = t->Emis->damp/SDIV(t1)/sqrt( 0.78540 + POW2(u1) );
			/* add very large optical depth tail to converge this with respect
			 * to escape probabilities - if this function falls off more slowly
			 * than escape probability then upper level will become overpopulated.
			 * original paper was not intended for this regime */
			if( t->Emis->TauCon>1e7 )
				wings *= pow( t->Emis->TauCon/1e7,-1.1 );
		}
		value = core + wings;
		/* some x-ray lines have vastly large damping constants, greater than 1.
		 * in these cases the previous wings value does not work - approximation
		 * is for small damping constant - do not let pump efficiency exceed unity
		 * in this case */
		if( t->Emis->TauCon>0. )
			value = MIN2(1., value );
	}
	else if( rt.nLineContShield == LINE_CONT_SHIELD_FERLAND )
	{
		/* set continuum shielding ferland */
		value = conpmp( t );
	}
	else if( rt.nLineContShield == 0 )
	{
		/* set continuum shielding none */
		value = 1.;
	}
	else
	{
		TotalInsanity();
	}

	/* the returned pump shield function must be greater than zero,
	 * and less than 1 if a maser did not occur */
	ASSERT( value>=0 && (value<=1.||t->Emis->TauCon<0.) );

	return value;
}


static const double BREAK = 3.;
/* fit to results for tau less than 10 */
STATIC inline double FITTED( double t )
{
	return (0.98925439 + 0.084594094*t)/(1. + t*(0.64794212 + t*0.44743976));
}

/*conpmp local continuum pumping rate radiative transfer for all lines */
STATIC double conpmp(transition * t)
{
	double a0, 
	  conpmp_v, 
	  tau, 
	  yinc1, 
	  yinc2;

	DEBUG_ENTRY( "conpmp()" );

	/* tau used will be optical depth in center of next zone
	 * >>chng 96 july 6, had been ipLnTauIn, did not work when sphere static set */
	tau = t->Emis->TauCon;
	/* compute pumping probability */
	if( tau <= 10. )
	{
		/* for tau<10 a does not matter, and one fit for all */
		conpmp_v = FITTED(tau);
	}
	else if( tau > 1e6 )
	{
		/* this far in winds line opacity well below electron scattering
			* so ignore for this problem */
		conpmp_v = 0.;
	}
	else
	{
		my_Integrand func;
		func.PumpDamp = t->Emis->damp;
		func.PumpTau = tau;
		Integrator<my_Integrand,Gaussian32> opfun;

		yinc1 = opfun.sum( 0., BREAK, func );
		yinc2 = opfun.sum( BREAK, 100., func );

		a0 = 0.886227*(1. + func.PumpDamp);
		conpmp_v = (yinc1 + yinc2)/a0;
	}

	/* EscProb is escape probability, will not allow conpmp to be greater than it
	 * on second iteration with thick lines, pump prob=1 and esc=0.5
	 * conpmp = MIN( conpmp , t->t(ipLnEscP) )
	 * */
	return( conpmp_v );
}
