/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_escape do line radiative transfer,
 * evaluates escape and destruction probability */
/* NB no wind distinction - that is done where optical depths are incremented
 * with line opacity - rt_line_one_tauinc and in line width for rad pressure */
/*RT_line_fine_opacity do fine opacities for one line */
/*RT_line_electron_scatter evaluate electron scattering escape probability */
/*RT_line_pumping pumping by external and locally emitted radiation fields */
#include "cddefines.h"
#include "rfield.h"
#include "doppvel.h"
#include "dense.h"
#include "opacity.h"
#include "transition.h"
#include "conv.h"
#include "radius.h"
#include "rt.h"
#include "physconst.h"
#include "cosmology.h"
#include "thirdparty.h"

STATIC void getVoigt(long int nCells_core, long int nCells_damp, 
							realnum *xprofile, realnum *profile, realnum damp);

/*RT_line_pumping pumping by external and locally emitted radiation fields */
STATIC void RT_line_pumping(
		 /* the em line we will work on  */
		 transition* t ,
		 /* this is option to not include line self shielding across this zone.
		 * this can cause pump to depend on zone thickness, and leads to unstable
		 * feedback in some models with the large H2 molecule, due to Solomon
		 * process depending on zone thickness and level populations. */
		 bool lgShield_this_zone,
		 realnum DopplerWidth)
{
	DEBUG_ENTRY( "RT_line_pumping()" );
	
	ASSERT( t->ipCont >= 1 );

	/* pumping by incident and diffuse continua */
	/* option to kill induced processes */
	if( !rfield.lgInducProcess )
	{
		t->Emis->pump = 0.;
	}
	else if( conv.lgFirstSweepThisZone || t->Emis->iRedisFun == ipLY_A )
	{
		double dTau;
		double shield_continuum;

		/* continuum shielding into this function */
		shield_continuum = RT_continuum_shield_fcn( t );

		/* continuum upward pumping rate, A gu/gl abs prob occnum
		* the "no induced" command causes continuum pumping to be set to 0 
		* this includes pumping by diffuse continuum */
		t->Emis->pump = t->Emis->Aul * t->Hi->g / t->Lo->g * shield_continuum *(
			rfield.OccNumbIncidCont[t->ipCont-1] + rfield.OccNumbContEmitOut[t->ipCont-1] );
		
		/* use mean of past few zone thicknesses rather than
		 * current zone thickness since this can induce wild oscillations, see
		 * nova_photos.in and blr_f92.in test caseis */
		dTau =  (t->Emis->PopOpc * t->Emis->opacity / DopplerWidth + 
			opac.opacity_abs[t->ipCont-1])* radius.drad_x_fillfac_mean;
	
		/* \todo 1 -- consider removing this cap.  This happens to late to be capping optical depth here */
		/* Cap optical depth at unity */
		dTau = MIN2(1., dTau);
		/* only limit because smaller just evaluates to 1, not clear why upper limit at all */
		if( lgShield_this_zone && dTau > 1e-3 && dTau < 1e10 )
		{
			/* correction for line optical depth across zone */
			t->Emis->pump *= log(1. + dTau ) / dTau;
		}
		
		/* NB NB line pumping by local diffuse emission is not included. */
	}

	return;
}

/*RT_line_electron_scatter evaluate electron scattering escape probability */
STATIC void RT_line_electron_scatter(
		 /* the em line we will work on  */
		 transition* t ,
		 realnum DopplerWidth)
{

	DEBUG_ENTRY( "RT_line_electron_scatter()" );

	/* escape following scattering off an electron */
	/* this is turned off with the no scattering escape command */
	if( !rt.lgElecScatEscape )
	{
		t->Emis->Pelec_esc = 0.;
		return;
	}

	// in the transition structure PopOpc is the pop relative to the ion,
	// but has been converted to a physical population before this routine
	// was called.  No need to convert here */
	double opac_line = t->Emis->PopOpc * t->Emis->opacity/DopplerWidth;

	if( opac_line > SMALLFLOAT )
	{
		/* the opacity in electron scattering */
		double opac_electron = dense.eden*6.65e-25;
		/* this is equation 5 of 
		 *>>refer	line	desp	Netzer, H., Elitzur, M., & Ferland, G. J. 1985, ApJ, 299, 752*/
		double opacity_ratio = opac_electron/(opac_electron+opac_line) * (1.-t->Emis->Pesc);
		/* keep total probability less than 0.1 */
		t->Emis->Pelec_esc = (realnum)opacity_ratio * MAX2(0.f,1.f-t->Emis->Pesc-t->Emis->Pdest);
	}
	else
		// slowly let escape probability go to zero
		t->Emis->Pelec_esc *= 0.95f;

	return;
}

/*RT_line_escape do line radiative transfer escape and destruction probabilities 
 * this routine sets */
STATIC void RT_line_escape(
	   /* the em line we will work on  */
	   transition* t,
	   /* Stark escape probability to be added to Pesc */
	   realnum pestrk,
	   realnum DopplerWidth,
	   bool lgGoodTau)
{
	int nRedis = -1;

	DEBUG_ENTRY( "RT_line_escape()" );

	/* a runaway maser */
	if( t->Emis->TauIn < -30. )
	{
		fprintf( ioQQQ, "PROBLEM RT_line_escape called with large negative "
			"optical depth, zone %.2f, setting lgAbort true.\n",
			fnzone );
		DumpLine(t);
		/* return busted true instead of exit here, so that code will
		 * not hang under MPI */
		lgAbort = true;
		return;
	}

	/* atomic number of this element on the C scale*/
	long int nelem = t->Hi->nelem-1;

	if( cosmology.lgDo )
	{
		/* Sobolev escape */
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			realnum tau_Sobolev =  t->Emis->TauIn;

			if( tau_Sobolev < 1E-5 )
			{
				t->Emis->Pesc = 1.;
			}
			else
			{
				t->Emis->Pesc = ( 1.f - exp( -1.f * tau_Sobolev ) )/ tau_Sobolev;
			}

			/* inward escaping fraction */
			t->Emis->FracInwd = rt.fracin;
		}
		fixit(); // is this correct?
		nRedis = ipDEST_K2;
	}
	/* static solution - which type of line will determine
	 * which redistribution function */
	/* iRedisFun == 1 - alpha resonance line, partial redistribution,
	 * ipPRD == 1 */
	else if( t->Emis->iRedisFun == ipPRD )
	{
		/* incomplete redistribution with wings */
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			t->Emis->Pesc = (realnum)esc_PRD( t->Emis->TauIn, t->Emis->TauTot, t->Emis->damp );

			/* >>chng 03 jun 07, do not clobber esp prob when line is masing -
			* this had effect of preventing total escape prob from getting larger than 1 */
			if( pestrk > 0.f && t->Emis->Pesc < 1.f )
				t->Emis->Pesc = min( 1.f, t->Emis->Pesc + pestrk );

			/* inward escaping fraction */
			t->Emis->FracInwd = rt.fracin;
		}
		nRedis = ipDEST_INCOM;
	}

	/* complete redistribution without wings - t->ipLnRedis is ipCRD == -1 */
	else if( t->Emis->iRedisFun == ipCRD )
	{
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			/* >>chng 01 mar -6, escsub will call any of several esc prob routines,
			* depending of how core is set.  We always want core-only for this option,
			* so call  esca0k2(tau) directly */
			t->Emis->Pesc = (realnum)esc_CRDcore( t->Emis->TauIn, t->Emis->TauTot );

			if( pestrk > 0.f && t->Emis->Pesc < 1.f )
				t->Emis->Pesc = min( 1.f, t->Emis->Pesc + pestrk );

			/* inward escaping fraction */
			t->Emis->FracInwd = rt.fracin;
		}
		nRedis = ipDEST_K2;
	}

	/* CRD with wings, = 2 */
	else if( t->Emis->iRedisFun == ipCRDW )
	{
		/* complete redistribution with damping wings */
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			t->Emis->Pesc = (realnum)esc_CRDwing( t->Emis->TauIn, t->Emis->TauTot, t->Emis->damp );

			if( pestrk > 0.f && t->Emis->Pesc < 1.f )
				t->Emis->Pesc = min( 1.f, t->Emis->Pesc + pestrk );

			/* inward escaping fraction */
			t->Emis->FracInwd = rt.fracin;
		}
		nRedis = ipDEST_K2;
	}

	/* Lya is special case */
	else if( t->Emis->iRedisFun == ipLY_A )
	{
		/* incomplete redistribution with wings, for special case of Lya
		* uses fits to Hummer & Kunasz numerical results 
		* this routine is different because escape and dest probs
		* are evaluated together, so no test of lgDoEsc */
		if( lgGoodTau )
		{
			double dest , esin;

			/* this will always evaluate escape prob, no matter what lgDoEsc is.
			* The destruction prob comes back as dest */
			t->Emis->Pesc = (realnum)RTesc_lya( &esin, &dest, t->Emis->PopOpc, nelem, DopplerWidth );

			if( pestrk > 0.f && t->Emis->Pesc < 1.f )
				t->Emis->Pesc = min( 1.f, t->Emis->Pesc + pestrk );

			/* this is current destruction rate */
			t->Emis->Pdest = (realnum)dest;

			/*  this is fraction of line which is inward escaping */
			t->Emis->FracInwd = rt.fracin;
		}
	}
	else
	{
		fprintf( ioQQQ, " RT_line_escape called with impossible redistribution function %d\n",
				t->Emis->iRedisFun);
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* only do this if not Lya special case, since dest prob already done */
	if( lgGoodTau && t->Emis->iRedisFun != ipLY_A && t->Emis->opacity > 0. )
	{
		t->Emis->Pdest = (realnum)RT_DestProb(
			/* the abundance of the species */
			t->Emis->PopOpc,
			/* line center opacity in funny units (needs a vel) */
			t->Emis->opacity,
			/* array index for line in continuum array,
			 * on f not c scale */
			t->ipCont,
			/* line width in velocity units */
			DopplerWidth,
			/* escape probability */
			t->Emis->Pesc,
			/* redistribution function */
			nRedis);
	}

	return;
}

/*RT_line_fine_opacity do fine opacities for one line */
STATIC void RT_line_fine_opacity(
	/* the em line we will work on  */
	transition* t ,
	realnum DopplerWidth)

{
	DEBUG_ENTRY( "RT_line_fine_opacity()" );

	/* this is line center frequency, including bulk motion of gas */
	long int ipLineCenter = t->Emis->ipFine + rfield.ipFineConVelShift;

	/* define fine opacity fine grid fine mesh */
	/* rfield.lgOpacityFine flag set false with no fine opacities command */
	if( !conv.lgLastSweepThisZone || ipLineCenter < 0 || t->Emis->PopOpc < SMALLFLOAT ||
		ipLineCenter>rfield.nfine || !rfield.lgOpacityFine )
	{
		return;
	}

	/* number of fine opacity cells corresponding to one doppler width for current
	 * temperature, velocity field, and nuclear mass,
	 * rfield.fine_opac_velocity_width is width per cell, cm/s */
	realnum cells_wide_1x = DopplerWidth/rfield.fine_opac_velocity_width;

	/* line center opacity - type realnum since will add to fine opacity array,
	 * which is realnum */
	realnum opac_line =  (realnum)t->Emis->PopOpc * t->Emis->opacity / DopplerWidth;

	// this is effective optical depth to this point. Do not do line if 
	// this product is less than SMALLFLOAT
	double dTauEffec = opac_line*radius.depth_x_fillfac;
	if( dTauEffec < SMALLFLOAT )
		return;

	/* core width of optically thick line, do 4x with exponential Doppler core,
	 * must be at least one cell, but profile is symmetric */
	const bool doDamp = dTauEffec*t->Emis->damp/9. > 0.1;
	long int nCells_core = (long)(cells_wide_1x*4.f + 1.5f);
	/* core is symmetric - make sure both upper and lower bounds
	 * are within continuum energy bin */
	if( ipLineCenter - nCells_core < 1 )
		nCells_core = ipLineCenter - 1;
	if( ipLineCenter + nCells_core > rfield.nfine )
		nCells_core = ipLineCenter -rfield.nfine - 1;

	/* want core to be at least one cell wide */
	nCells_core = MAX2( 1 , nCells_core );

	long int nCells_damp;
	/* include damping wings if optical depth is large */
	if( doDamp )
	{
		// find number of cells to extend the damping wings - cells_wide_1x is one dop width
		// tests with th85orion, stopping at half -h2 point, showed that 0.01 was
		// needed for outer edge, given the definition of dTauEffec
		realnum x = (realnum)sqrt( dTauEffec * t->Emis->damp*100./SQRTPI ) * cells_wide_1x;
		// test on size of x, which can exceed
		// limits of long in extreme optical depths */
		if( x<LONG_MAX )
		{
			nCells_damp = (long)x;

			if( ipLineCenter-nCells_damp < 1 )
				nCells_damp = ipLineCenter-1;

			if( ipLineCenter+nCells_damp > rfield.nfine )
				nCells_damp = rfield.nfine - ipLineCenter-1;
		}
		else
		{
			/* x was too big, just set to extreme range, which is
			* number of cells to nearest boundary */
			nCells_damp = MIN2( rfield.nfine-ipLineCenter , ipLineCenter )-1;
		}
	}
	else
	{
		nCells_damp = nCells_core;
	}

	static vector<realnum> xprofile, profile;
	xprofile.resize(nCells_damp);
	profile.resize(nCells_damp);

	for( long int i=0; i<nCells_damp; ++i )
	{
		/* distance from line center in units of doppler width */
		xprofile[i] = (realnum) i/cells_wide_1x;
	}
	
	getVoigt(nCells_core,  nCells_damp, &xprofile[0], &profile[0],
				t->Emis->damp);

	/* line center itself, must not double count here */
	rfield.fine_opac_zone[ipLineCenter] += profile[0]*opac_line;
	for( long int i=1; i<nCells_damp; ++i )
	{
		rfield.fine_opac_zone[ipLineCenter+i] += profile[i]*opac_line;
		rfield.fine_opac_zone[ipLineCenter-i] += profile[i]*opac_line;
	}

	return;
}

STATIC void getVoigt(long int nCells_core, long int nCells_damp, 
							realnum *xprofile, realnum *profile, realnum damp)
{
	if (1)
	{
		profile[0] = 1.;
		for( long int i=1; i<nCells_core; ++i )
		{
			realnum xsq = POW2(xprofile[i]);
			
			/* the line opacity at that point, 
				includes doppler core and damping wings */
			profile[i] = sexp(xsq) + damp / MAX2(1.f,((realnum)SQRTPI*xsq));
			
		}
		ASSERT( nCells_core>0 );
		for( long int i=nCells_core; i<nCells_damp; ++i )
		{
			/* distance from line center in doppler units */
			realnum xsq = POW2(xprofile[i]);
			profile[i] = damp/((realnum)SQRTPI*xsq);
		}
	}
	else
	{
		humlik(nCells_damp,xprofile,damp,profile);
	}
	
}

/*RT_line_one do rt for emission line structure - calls RT_line_escape or RT_line_wind */
void RT_line_one(
	/* the em line we will work on  */
	transition* t,
	/* this is option to not include line self shielding across this zone.
	 * this can cause pump to depend on zone thickness, and leads to unstable
	 * feedback in some models with the large H2 molecule, due to Solomon
	 * process depending on zone thickness and level populations. */
	bool lgShield_this_zone,
	/* Stark escape probability to be added to Pesc */
	realnum pestrk,
	realnum DopplerWidth )
{
	DEBUG_ENTRY( "RT_line_one()" );

	// do nothing is population and this is not the very first call 
	// skip line transfer if requested with 'no line transfer' command, but never skip Lya
	if( !rfield.lgDoLineTrans && (t->Emis->iRedisFun != ipLY_A) )
	{
		return;
	}

	/* line damping constant at current temperature  */
	t->Emis->damp = t->Emis->dampXvel / DopplerWidth;
	ASSERT( t->Emis->damp > 0. );

	// do not evaluate if no population 
	if( t->Lo->Pop<=SMALLFLOAT )
	{
		/* zero population, return after setting everything with side effects */
		t->Emis->Pesc = 1.f;

		/* inward escaping fraction */
		t->Emis->FracInwd = 0.5;

		/* pumping rate */
		t->Emis->pump = 0.;

		/* destruction probability */
		t->Emis->Pdest = 0.;
		t->Emis->Pelec_esc = 0.;

		return;
	}

	/* option to keep track of population values during calls,
	 * print out data to make histogram */
	enum {DEBUG_LOC=false};
	if( DEBUG_LOC )
	{
		static long int nTau[100];
		long n;

		if( nzone==0 )
		{
			for(n=0; n<100; ++n )
				nTau[n] = 0;
		}
		if( t->Lo->Pop<=SMALLFLOAT )
			n = 0;
		else
			n = (long)log10( t->Lo->Pop )+37;
		n = MIN2( n , 99 );
		n = MAX2( n , 0 );
		++nTau[n];
		if( nzone > 183 )
		{
			for(n=0; n<100; ++n )
				fprintf(ioQQQ,"%li\t%li\n", n , nTau[n] );
			cdEXIT(EXIT_SUCCESS);
		}
	}
	
	// transition is below plasma frequency - photons not emitted
	if( t->EnergyErg / EN1RYD <= rfield.plsfrq )
	{
		t->Emis->Pesc = SMALLFLOAT;
		t->Emis->Pdest = SMALLFLOAT;
		t->Emis->Pelec_esc = SMALLFLOAT;
		t->Emis->pump = SMALLFLOAT;
	}
	else
	{

		/* this checks if we have overrun the optical depth scale,
		* in which case the inward optical depth is greater than the
		* previous iteration's total optical depth.
		* We do not reevaluate escape probabilities if the optical depth
		* scale has been overrun due to huge bogus change in solution
		* that would result */
		bool lgGoodTau = lgTauGood( t );

		// the last sweep through this zone is to do the fine opacities
		// the populations are not updated on the last sweep so the 
		// line transfer details also should not be updated
		if( conv.lgLastSweepThisZone )
			RT_line_fine_opacity( t , DopplerWidth );
		else
		{
			RT_line_escape( t, pestrk, DopplerWidth , lgGoodTau);
			RT_line_electron_scatter( t , DopplerWidth );
			RT_line_pumping( t , lgShield_this_zone , DopplerWidth );	
		}
	}

	return;
}

