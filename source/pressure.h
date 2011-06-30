/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PRESSURE_H_
#define PRESSURE_H_

#include "rt.h"
#include "rfield.h"
#include "doppvel.h"
#include "physconst.h"
#include "transition.h"

/**PressureTotal determine the gas and line radiation pressures for current conditions,
 * this sets the variable pressure.PresTotlCurr */
void PresTotCurrent(void);

/**PressureRadiationLine calculate radiation pressure due to a particular line */
inline double PressureRadiationLine( const transition *t, realnum DopplerWidth )
{
	DEBUG_ENTRY( "PressureRadiationLine()" );

	/* return zero if below plasma frequency */	
	if( t->EnergyErg / EN1RYD <= rfield.plsfrq )
		return 0.;

	/* RT_LineWidth gets line width in terms of RT effects */
	double width = RT_LineWidth(t, DopplerWidth);

	double PopOpc = t->Emis->PopOpc/t->Lo->g;
	/* return zero line radiation PressureReturned if line mases or has
	 * zero opacity */
	/* \todo 2 1e-22 is arbtrary but roughly 1/kpc.  Replace with a cloud width if available? */
	if( PopOpc*t->Emis->opacity/ DopplerWidth <= 1.e-22 || width<=0. )
		return 0.;

	double PressureReturned = PI8 * HPLANCK / 3. * POW4(t->EnergyWN) *
		(t->Hi->Pop/t->Hi->g)/PopOpc * width;

	/* this prevents line radiation PressureReturned from being very large when line 
	 * is not optically thick but total opacity at that energy is large 
	 * due to overlapping transitions */
	long int ipLineCenter = t->Emis->ipFine + rfield.ipFineConVelShift;
	if( ipLineCenter > 0 && ipLineCenter < rfield.nfine && rfield.lgOpacityFine &&
		rfield.fine_opac_zone[ipLineCenter] > SMALLFLOAT )
	{
		double FractionThisLine = t->Emis->PopOpc * t->Emis->opacity / DopplerWidth/
			rfield.fine_opac_zone[ipLineCenter];
		if( FractionThisLine<1e-5 )
			FractionThisLine = 0.;
		/* fine opacities are only reevaluated one time per zone due
		 * to the expense - PopOpc is for the current solution - but the two
		 * may be out of step by a few percent, due to the variation in
		 * abundance from zone to zone.  This prevents the change
		 * in solution from increasing the radiation pressure.
		 * This correction is mainly an order of magnitude scaler to prevent
		 * optically thin lines from appearing to be optically thick due to
		 * overlapping lines */
		FractionThisLine = MIN2(1., FractionThisLine);
		ASSERT( FractionThisLine >= 0. && FractionThisLine <= 1.0 );
		PressureReturned *= FractionThisLine;
	}

	return PressureReturned;
}

/** variables dealing with pressure across model */
EXTERN struct t_pressure {

	/** lowest and highest total current pressure to desired
	 * pressures in model, a measure of how well converged it is.  This
	 * includes total pressure, so in a constant gas pressure model will
	 * show a large excursion */
	realnum PresLow, 
	  PresHigh;

	realnum PresPowerlaw;

	/** the ram pressure */
	double PresRamCurr;

	/** current turbulent pressure */
	double PresTurbCurr;

	/** the current pressure, and the correct pressure, the ratio
	 * of the two is the error, PresTotlCurr is set in PressureTotal */
	double PresTotlCurr, 
	  PresTotlCorrect,
	  /** PresGasCurr is gas pressure, nkT, set in PressureTotal */
	  PresGasCurr;

	/** PresTotlInit - total pressure at the illuminated face */
	double PresTotlInit;

	/** option to set pressure at illuminated face as number on 
	 * constant pressure command */
	bool lgPressureInitialSpecified;
	/** linear pressure in nkT units, zero if not set */
	double PressureInitialSpecified;

	/** pres_radiation_lines_curr is line radiation pressure for current zone */
	double pres_radiation_lines_curr;

	/** flag saying whether or not incident continuum should be included
	 * in total pressure, turned off with constant gas pressure */
	bool lgContRadPresOn;

	/** total pressure related variables
	 * PresInteg is integrated radiation pressure */
	realnum PresInteg, 
	  pinzon;

	/** integrated radiative acceleration, if no opacity present,
	 * only dilution due to 1/r^2
	 * this is the acceleration used in Eddington luminosity comparison */
	realnum PresIntegElecThin,
		pinzon_PresIntegElecThin;

	/** (Self-)gravity forces: Yago Ascasibar (UAM, Spring 2009) */
	double RhoGravity_dark;
	double RhoGravity_self;
	double RhoGravity_external;
	double RhoGravity;
	double IntegRhoGravity;
	int gravity_symmetry;
	double self_mass_factor;
	vector<double> external_mass[3];

	/** lgPres_radiation_ON says whether radiation pressure enabled, turned off with
	* constant density, constant gas pressure commands, on with constant pressure */
	bool lgPres_radiation_ON;
	bool lgPres_magnetic_ON;
	bool lgPres_ram_ON;

	realnum 
	  /** RadBetaMax is largest ratio of radiation to gas pressure */
	  RadBetaMax, 
	  /** pbeta is ratio of radition to gas pressure, evaluated in PressureTotal */
	  pbeta, 
	  /** PresMas is largest pressure that occurred in the calculation */
	  PresMax;

	/** pointer to line with greatest radiation pressure */
	long int ipPradMax_line;

	/** zone where greatest radiation pressure occurred */
	long int ipPradMax_nzone;

	/** string with label for line with greatest contributrion to pressure*/
	char chLineRadPres[101];

	/** lgPradCap true if radiation pressure capped on first iteration
	 * lgPradDen capped by thermalization length */
	bool lgPradCap, 
	  lgPradDen;

	/** flag true if radiation pressure is turned on, part of total pressure */
	bool lgLineRadPresOn;

	/** option to not abort on high radiation pressure - set true in initialization,
	 * set false with NO ABORT on constant radiation pressure command */
	bool lgRadPresAbortOK;

	/** option to abort when we reach the sonic point - set according to
	 * the dynamics pressure mode */ 
	bool lgSonicPointAbortOK; 

	/** we hit the sonic point */
	bool lgSonicPoint;

	/** True when we are in limbo while trying to find a strong-D
	 * solution. This is when there is no possible solution for the
	 * current target pressure. We just grit our teeth and plough
	 * onwards, hoping to come out the other side and to sort it all
	 * out on a later iteration of the dynamics */
	bool lgStrongDLimbo;

	}	pressure;


#endif /* PRESSURE_H_ */
