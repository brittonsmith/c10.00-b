/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_one_tauinc increment optical depths for all heavy element lines, zone by zone,
 * mainly called by RT_tau_inc, but also by FeII */
#include "cddefines.h"
#include "doppvel.h"
#include "geometry.h"
#include "rfield.h"
#include "radius.h"
#include "wind.h"
#include "rt.h"
#include "physconst.h"
#include "cosmology.h"
#include "transition.h"

void RT_line_one_tauinc(transition * t ,
	/* following four are flags to generate info if some species has extreme maser */
	long int maser_flag_species,
	long int maser_flag_ion,
	long int maser_flag_hi,
	long int maser_flag_lo,
	realnum DopplerWidth )

{
	DEBUG_ENTRY( "RT_line_one_tauinc()" );

	/* routine increments optical depths for static or expanding atmosphere */

	/* this is line center frequency, including bulk motion of gas */
	long int ipLineCenter = t->Emis->ipFine + rfield.ipFineConVelShift;
	double OpacityEffective, EffectiveThickness;
	realnum dTau_total;

	/* find line center opacity - use fine opacity if array indices are OK */
	if( t->Emis->ipFine>=0 && ipLineCenter>0 && ipLineCenter<rfield.nfine && rfield.lgOpacityFine )
	{
		/* use fine opacities fine grid fine mesh to get optical depth 
		 * to continuum source */
		/* total line center optical depth, all overlapping lines included */
		OpacityEffective = rfield.fine_opac_zone[ipLineCenter];
	}
	else
	{
		OpacityEffective = t->Emis->PopOpc * t->Emis->opacity / DopplerWidth;
	}

#if	0
	if( rfield.anu[ t->ipCont-1 ] < rfield.plsfrq )
	{
		/* transition is below plasma frequency - make optical depth huge */
		dTau_total = 1.e10;

		t->Emis->TauIn = dTau_total;
		t->Emis->TauCon = dTau_total;
		t->Emis->TauTot = dTau_total;
	}
	else 
#endif
	if( cosmology.lgDo )
	{
		/* dv/dr (s-1), equal to dv/dt / v */
		/* in this case, dv/dr is just the Hubble factor */
		wind.dvdr = GetHubbleFactor(cosmology.redshift_current);

		fixit(); //This doppler width is sqrt(2kt/m), but Seager et al use sqrt(3kt/m).  Resolve
		EffectiveThickness = DopplerWidth / wind.dvdr;
		dTau_total = (realnum)(OpacityEffective * EffectiveThickness);
		
		t->Emis->TauIn = dTau_total;
		t->Emis->TauCon = dTau_total;
		t->Emis->TauTot = dTau_total;
	}

	/* use cumulated fine optical depth for both d-critical and static, 
	 * for d-critical speeds are only roughly sonic
	 * optical depth is computed including velocity shift */
	else if( ! wind.lgBallistic() )
	{
		/* static and negative velocity solutions */
		EffectiveThickness = radius.drad_x_fillfac;
		dTau_total = (realnum)(OpacityEffective * EffectiveThickness);
		
		t->Emis->TauIn += dTau_total;
		t->Emis->TauCon += dTau_total;

	}

	else
	{
		/* ballistic outflowing wind
		 * effective length scale for Sobolev or LVG approximation, eqn 3 of
		 * >>refer	RT	wind	Castor, J.I., Abbott, D.C., & Klein, R.I., 1975, ApJ, 195, 157
		 */

		/* dv/dr (s-1), equal to dv/dt / v */
		wind.dvdr = fabs(wind.AccelTotalOutward - wind.AccelGravity) / wind.windv;
		/* depth (cm) over which wind accelerates by one velocity width
		 * include filling factor */
		EffectiveThickness = DopplerWidth / SDIV(wind.dvdr) * geometry.FillFac;

		/* min2 is to not let the physical scale exceed the current depth */
		EffectiveThickness = MIN2( radius.depth, EffectiveThickness );
		dTau_total = (realnum)(OpacityEffective * EffectiveThickness);
		
		t->Emis->TauIn = dTau_total;
		t->Emis->TauCon = dTau_total;
		t->Emis->TauTot = dTau_total;
	}

	/* keep track of any masers */
	if( dTau_total < rt.dTauMase )
	{
		rt.dTauMase = dTau_total;
		rt.mas_species = maser_flag_species;
		rt.mas_ion = maser_flag_ion;
		rt.mas_hi = maser_flag_hi;
		rt.mas_lo = maser_flag_lo;
		if( rt.dTauMase < -1. )
			rt.lgMaserCapHit = true;
	}

	return;
}
