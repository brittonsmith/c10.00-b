/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*DrvContPump local continuum pumping rate radiative transfer for all lines */
/*con_pump_op  routine used to get continuum pumping of lines 
 * used in DrvContPump in call to qg32 */
#include "cddefines.h"
#include "rfield.h"
#include "doppvel.h"
#include "radius.h"
#include "continuum.h"
#include "transition.h"

/*con_pump_op  routine used to get continuum pumping of lines 
 * used in DrvContPump in call to qg32 */
class my_Integrand_con_pump_op
{
public:
	/* damping constant used for pumping */
	realnum damp;
	/* variable used for inward optical depth for pumping */
	realnum PumpTau;

	double operator() (double x)
	{
		double v = vfun(damp,x);
		double opfun_v = sexp(PumpTau*v)*v;
		return( opfun_v );
	}
};

/** \todo	2	if used, add damp as arg since calling routine probably evaluated it */
double DrvContPump(transition * t, realnum DopplerWidth)
{
	double a0, 
	  ContPump_v, 
	  tau, 
	  yinc1, 
	  yinc2;

	DEBUG_ENTRY( "DrvContPump()" );

	/* fit to results for tau less than 10 */
#	define FITTED(t)	((0.98925439 + 0.084594094*(t))/(1. + (t)*(0.64794212 + (t)*0.44743976)))

	if( !rfield.lgInducProcess )
	{
		/* option to turn off continuum pumping with no fluorescence */
		ContPump_v = 0.;
	}
	else
	{
		/* tau used will be optical depth in center of next zone */
		tau = t->Emis->TauIn + t->Emis->PopOpc * t->Emis->opacity / DopplerWidth * radius.dRNeff;
		/* compute pumping probability */
		if( tau <= 10. )
		{
			/* for tau<10 a does not matter, and one fit for all */
			ContPump_v = FITTED(tau);
		}
		else if( tau > 1e6 )
		{
			/* this far in winds line opacity well below electron scattering
			 * so ignore for this problem */
			ContPump_v = 0.;
		}
		else
		{
			my_Integrand_con_pump_op func;
			if( t->Emis->iRedisFun > 0 )
			{
				func.damp = t->Emis->damp;
			}
			else
			{
				func.damp = 0.;
			}
			func.PumpTau = (realnum)tau;

			Integrator<my_Integrand_con_pump_op,Gaussian32> con_pump_op;
#			define	BREAK_	3.
			yinc1 = con_pump_op.sum(0.,BREAK_,func);
			yinc2 = con_pump_op.sum(BREAK_,100.,func);

			a0 = 0.886227*(1. + func.damp);
			ContPump_v = (yinc1 + yinc2)/a0;
		}
	}

	/* EscProb is escape probability, will not allow ContPump to be greater than it
	 * on second iteration with thick lines, pump prob=1 and esc=0.5
	 * ContPump = MIN( ContPump , t->t(ipLnEscP) )
	 * */
	return( ContPump_v );
#	undef	FITTED
}
