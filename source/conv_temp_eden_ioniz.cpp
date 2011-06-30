/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvTempEdenIoniz determine  temperature, called by ConPresTempEdenIoniz,
 * calls ConvEdenIoniz to get electron density and ionization */
/*lgConvTemp returns true if heating-cooling is converged */
/*CoolHeatError evaluate ionization, and difference in heating and cooling, for temperature temp */
/*DumpCoolStack helper routine to dump major coolants */
/*DumpHeatStack helper routine to dump major heating agents */
#include "cddefines.h"
#include "hmi.h"
#include "thermal.h"
#include "iso.h"
#include "hydrogenic.h"
#include "colden.h"
#include "h2.h"
#include "pressure.h"
#include "dense.h"
#include "struc.h"
#include "thirdparty.h"
#include "trace.h"
#include "phycon.h"
#include "conv.h"

/*lgConvTemp returns true if heating-cooling is converged */
STATIC bool lgConvTemp(const iter_track& TeTrack);
/*CoolHeatError evaluate ionization, and difference in heating and cooling, for temperature temp */
STATIC double CoolHeatError( double temp );

// debugging routines to print main sources of cooling and heating
STATIC void DumpCoolStack(double thres);
STATIC void DumpHeatStack(double thres);

/*ConvTempEdenIoniz determine  temperature, called by ConPresTempEdenIoniz,
 * calls ConvEdenIoniz to get electron density and ionization 
 * returns 0 if ok, 1 if disaster */
int ConvTempEdenIoniz( void )
{
	DEBUG_ENTRY( "ConvTempEdenIoniz()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "\n  ConvTempEdenIoniz called\n" );
	}
	if( trace.nTrConvg >= 2 )
	{
		fprintf( ioQQQ, "  ConvTempEdenIoniz called, entering temp loop using solver %s.\n",
			 conv.chSolverTemp );
	}
	// this branch uses the van Wijngaarden-Dekker-Brent method
	if( strcmp( conv.chSolverTemp , "vWDB" ) == 0 )
	{
		conv.lgConvTemp = false;

		// deal with special temperature laws first
		if( thermal.lgTLaw )
		{
			double TeNew = phycon.te;
			if( thermal.lgTeBD96 )
			{
				/* Bertoldi & Drain 96 temp law specified by TLAW BD96 command */
				TeNew = thermal.T0BD96 / (1. + thermal.SigmaBD96 * colden.colden[ipCOL_HTOT]);
			}
			else if( thermal.lgTeSN99 )
			{
				/* Sternberg & Neufeld 99 temp law specified by TLAW SN99 command,
				 * this is equation 16 of 
				 * >>refer	H2	temperature law	Sternberg, A., & Neufeld, D.A. 1999, ApJ, 516, 371-380 */
				TeNew = thermal.T0SN99 / 
					(1. + 9.*POW4(2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN]) );
			}
			else
				TotalInsanity();

			TempChange( TeNew, false );
		}

		if( thermal.lgTemperatureConstant || thermal.lgTLaw )
		{
			if( ConvEdenIoniz() ) 
				return 1;
			PresTotCurrent();

			// convergence is automatic...
			conv.lgConvTemp = true;

			if( trace.lgTrace || trace.nTrConvg >= 2 )
			{
				fprintf( ioQQQ, "  ConvTempEdenIoniz: Te %e C %.4e H %.4e\n",
					 phycon.te, thermal.ctot, thermal.htot );
				fprintf( ioQQQ, "  ConvTempEdenIoniz returns ok.\n" );
			}

			return 0;
		}

		// here starts the standard solver for variable temperature
		iter_track TeTrack;
		double t1=0, error1=0, t2, error2;
		bool lgHighTeSearch = false;

		t2 = phycon.te;
		error2 = CoolHeatError( t2 );

		for( int n=0; n < 5 && !lgAbort; ++n )
		{
			const int DEF_ITER = 10;
			const double DEF_FACTOR = 0.02;
			double step, factor = DEF_FACTOR;

			TeTrack.clear();

			// set up an initial guess for the bracket
			// t2 was already initialized outside the main loop, or is copied from the
			// previous iteration. don't record this evaluation, it may be poorly converged
			for( int i=0; i < 2 && !lgAbort; ++i )
			{
				t1 = t2;
				error1 = error2;

				// the factor 1.2 creates 20% safety margin
				step = safe_div( -1.2*error1, conv.dCmHdT, 0. );
				step = sign( min( abs(step), factor*t1 ), step );
				t2 = t1 + step;
				error2 = CoolHeatError( t2 );
				TeTrack.add( t2, error2 );
			}

			int j = 0;

			// now hunt until we have bracketed the solution
			while( error1*error2 > 0. && j++ < 2*DEF_ITER && !lgAbort )
			{
				t1 = t2;
				error1 = error2;
				double deriv = TeTrack.deriv(7);
				// the factor 1.2 creates 20% safety margin
				step = safe_div( -1.2*error1, deriv, 0. );
				step = sign( min( abs(step), factor*t1 ), step );
				t2 = t1 + step;
				error2 = CoolHeatError( t2 );
				TeTrack.add( t2, error2 );
			}

			if( trace.nTrConvg >= 2 && error1*error2 > 0. && !lgAbort )
			{
				fprintf( ioQQQ, "  ConvTempEdenIoniz: bracket1 fails t1: %e %e t2: %e %e\n",
					 t1, error1, t2, error2 );
				TeTrack.print_history();
			}

			// If we reach this point without having bracketed the solution this means
			// that searching using the derivate doesn't work... This typically happens
			// when the heating and cooling curve no longer intersect, which implies
			// that C-H vs Te now has a local maximum, but doesn't change sign any more.
			// Hence we have reached a thermal front and we need to take a discontinuous
			// step up or down. In the first while loop we search upwards, in the second
			// downwards.

			int nUP = 2;

			if( error2 < 0. && !lgHighTeSearch && !lgAbort )
			{
				lgHighTeSearch = true;
				nUP = 25;
			}

			// hunting using the derivative failed, so first start hunting upwards
			//
			// guard against overrunning the maximum temperature limit of the code
			// so that we are guaranteed to get a chance to hunt downwards
			while( error1*error2 > 0. && j++ < (2+nUP)*DEF_ITER && !lgAbort &&
			       t2*(1.+factor) <= 0.99999*phycon.TEMP_LIMIT_HIGH )
			{
				t1 = t2;
				error1 = error2;
				if( lgHighTeSearch && t2 > 4000. )
					factor = 0.05;
				t2 = t1*(1.+factor);
				error2 = CoolHeatError( t2 );
				TeTrack.add( t2, error2 );
			}

			// hunting upwards failed, so start hunting downwards
			// we may need to take a big plunge downwards, so # of iter should be big
			//
			// don't guard against dropping below the minimum temperature limit of the
			// code so that we get an abort if the solver has completely lost the solution.
			// this is a serious condition that needs to be investigated...
			while( error1*error2 > 0. && j++ < (2+2*nUP+9)*DEF_ITER && !lgAbort )
			{
				t1 = t2;
				error1 = error2;
				if( lgHighTeSearch && t2 <= 4000. )
					factor = DEF_FACTOR;
				t2 = t1*(1.-factor);
				error2 = CoolHeatError( t2 );
				TeTrack.add( t2, error2 );
			}

			factor = DEF_FACTOR;

			if( trace.nTrConvg >= 2 && error1*error2 > 0. && !lgAbort )
			{
				fprintf( ioQQQ, "  ConvTempEdenIoniz: bracket2 fails t1: %e %e t2: %e %e\n",
					 t1, error1, t2, error2 );
				TeTrack.print_history();
			}

			// keeping the history up until now has a bad effect on convergence
			// so we wipe the slate clean....
			TeTrack.clear();

			// the bracket should have been found, now set up the Brent solver
			if( TeTrack.init_bracket( t1, error1, t2, error2 ) == 0 )
			{
				// The convergence criterion is based on the relative accuracy of Cool-Heat,
				// combined with a relative accuracy on the temperature itself. We need to
				// keep iterating until both accuracies are reached. Here we set tolerance on
				// Te to 2 ulp. If bracket gets narrower than 3 ulp we declare a convergence
				// failure to avoid changes getting lost in machine precision.
				TeTrack.set_tol(2.*DBL_EPSILON*t2);

				t2 = 0.5*(t1+t2);
				for( int i = 0; i < (1<<(n/2))*DEF_ITER && !lgAbort; i++ )
				{
					// check for convergence, as well as a pathologically narrow bracket
					if( lgConvTemp(TeTrack) || TeTrack.bracket_width() < 3.*DBL_EPSILON*t2 )
						break;

					error2 = CoolHeatError( t2 );
					TeTrack.add( t2, error2 );
					t2 = TeTrack.next_val(factor);
				}

				if( conv.lgConvTemp )
					break;

				if( trace.nTrConvg >= 2 && !lgAbort )
				{
					fprintf( ioQQQ, "  ConvTempEdenIoniz: brent fails\n" );
					TeTrack.print_history();
				}
			}
		}

		if( lgAbort )
			return 1;

		// only declare solution unstable if it is at least at the 2-sigma confidence level
		thermal.lgUnstable = ( conv.dCmHdT + 2.*conv.sigma_dCmHdT < 0. );

		if( trace.lgTrace || trace.nTrConvg >= 2 )
		{
			fprintf( ioQQQ, "  ConvTempEdenIoniz: Te %e C %.4e H %.4e (C-H)/H %.2f%%"
				 " d(C-H)/dT %.2e +/- %.2e\n",
				 phycon.te, thermal.ctot, thermal.htot,
				 (thermal.ctot/thermal.htot-1.)*100.,
				 conv.dCmHdT, conv.sigma_dCmHdT );
			fprintf( ioQQQ, "  ConvTempEdenIoniz returns converged=%c\n", TorF(conv.lgConvTemp) );
		}
	}
	else
	{
		fprintf( ioQQQ, "ConvTempEdenIoniz finds insane solver %s\n", conv.chSolverTemp );
		ShowMe();
	}

	return 0;
}


/* returns true if heating-cooling is converged */
STATIC bool lgConvTemp(const iter_track& TeTrack)
{
	DEBUG_ENTRY( "lgConvTemp()" );

	if( lgAbort )
	{
		/* converging the temperature was aborted */
		conv.lgConvTemp = false;
	}
	else if( ( abs(thermal.htot - thermal.ctot)/thermal.htot <= conv.HeatCoolRelErrorAllowed &&
		   TeTrack.bracket_width()/phycon.te <= conv.HeatCoolRelErrorAllowed/3. ) ||
		 thermal.lgTemperatureConstant || phycon.te <= phycon.TEMP_LIMIT_LOW )
	{
		/* announce that temp is converged if relative heating - cooling mismatch
		 * is less than the relative heating cooling error allowed, or
		 * if this is a constant temperature model */
		conv.lgConvTemp = true;
		// remember numerical derivative to estimate initial stepsize on next call
		conv.dCmHdT = TeTrack.deriv(conv.sigma_dCmHdT);
	}
	else
	{
		/* big mismatch, this has not converged */
		conv.lgConvTemp = false;
	}

	if( trace.nTrConvg >= 2 )
		fprintf( ioQQQ, "  lgConvTemp: C-H abs err %.4e Te err %.4e converged=%c\n",
			 abs(thermal.htot - thermal.ctot)/thermal.htot,
			 TeTrack.bracket_width()/phycon.te,
			 TorF(conv.lgConvTemp) );

	return conv.lgConvTemp;
}

/*CoolHeatError evaluate ionization, and difference in heating and cooling, for temperature temp */
STATIC double CoolHeatError( double temp )
{
	DEBUG_ENTRY( "CoolHeatError()" );

	TempChange( temp, false );

	/* converge the ionization and electron density; 
	 * this calls ionize until lgIonDone is true */
	/* NB should NOT set insanity - but rather return error condition */
	if( ConvEdenIoniz() )
		lgAbort = true;

	/* >>chng 01 mar 16, evaluate pressure here since changing and other values needed */
	/* reevaluate pressure */
	/* this sets values of pressure.PresTotlCurr */
	PresTotCurrent();

	/* keep track of temperature solver in this zone
	 * conv.hist_temp_nzone is reset in ConvInitSolution */
	if( nzone != conv.hist_temp_nzone )
	{
		/* first time in this zone - reset history */
		conv.hist_temp_nzone = nzone;
		conv.hist_temp_temp.clear();
		conv.hist_temp_heat.clear();
		conv.hist_temp_cool.clear();
	}

	conv.hist_temp_temp.push_back( phycon.te );
	conv.hist_temp_heat.push_back( thermal.htot );
	conv.hist_temp_cool.push_back( thermal.ctot );

	// dump major contributors to heating and cooling - for debugging purposes
	if( false )
	{
		DumpCoolStack( conv.HeatCoolRelErrorAllowed/5.*thermal.ctot );
		DumpHeatStack( conv.HeatCoolRelErrorAllowed/5.*thermal.htot );
	}

	if( trace.nTrConvg >= 2 )
		fprintf( ioQQQ, "  CoolHeatError: Te: %.4e C: %.4e H: %.4e (C-H)/H: %.4e\n",
			 temp, thermal.ctot, thermal.htot, thermal.ctot/thermal.htot-1. );

	double error = thermal.ctot - thermal.htot;

	// this can get set if temperature drops below floor temperature -> fake convergence
	if( thermal.lgTemperatureConstant )
		error = 0.;

	return error;
}

STATIC void DumpCoolStack(double thres)
{
	multimap<double,string> output;
	char line[200];

	for( int i=0; i < thermal.ncltot; ++i )
	{
		double fraction;
		if( abs(thermal.heatnt[i]) > thres )
		{
			fraction = thermal.heatnt[i]/thermal.ctot;
			sprintf( line, "heat %s %e: %e %e\n",
				 thermal.chClntLab[i], thermal.collam[i], thermal.heatnt[i], fraction );
			output.insert( pair<const double,string>( fraction, string(line) ) );
		}
		if( abs(thermal.cooling[i]) > thres )
		{
			fraction = thermal.cooling[i]/thermal.ctot;
			sprintf( line, "cool %s %e: %e %e\n",
				 thermal.chClntLab[i], thermal.collam[i], thermal.cooling[i], fraction );
			output.insert( pair<const double,string>( fraction, string(line) ) );
		}
	}

	dprintf( ioQQQ, " >>>>>>> STARTING COOLING DUMP <<<<<<\n" );
	dprintf( ioQQQ, "total cooling %e\n", thermal.ctot );
	// this will produce sorted output in reverse order (largest contributor first)
	for( multimap<double,string>::reverse_iterator i=output.rbegin(); i != output.rend(); ++i )
		dprintf( ioQQQ, "%s", i->second.c_str() );
	dprintf( ioQQQ, " >>>>>>> FINISHED COOLING DUMP <<<<<<\n" );
}

STATIC void DumpHeatStack(double thres)
{
	multimap<double,string> output;
	char line[200];

	for( int nelem=0; nelem < LIMELM; ++nelem )
	{
		for( int i=0; i < LIMELM; ++i )
		{
			double fraction = thermal.heating[nelem][i]/thermal.htot;
			if( abs(thermal.heating[nelem][i]) > thres )
			{
				sprintf( line, "heating[%i][%i]: %e %e\n",
					 nelem, i, thermal.heating[nelem][i], fraction );
				output.insert( pair<const double,string>( fraction, string(line) ) );
			}
		}
	}

	dprintf( ioQQQ, " >>>>>>> STARTING HEATING DUMP <<<<<<\n" );
	dprintf( ioQQQ, "total heating %e\n", thermal.htot );
	// this will produce sorted output in reverse order (largest contributor first)
	for( multimap<double,string>::reverse_iterator i=output.rbegin(); i != output.rend(); ++i )
		dprintf( ioQQQ, "%s", i->second.c_str() );
	dprintf( ioQQQ, " >>>>>>> FINISHED HEATING DUMP <<<<<<\n" );
}
