/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atoms_fe2ovr compute FeII overlap with Lya */
/*fe2par evaluate FeII partition function */
#include "cddefines.h"
#include "doppvel.h"
#include "dense.h"
#include "rfield.h"
#include "iso.h"
#include "phycon.h"
#include "hydrogenic.h"
#include "ipoint.h"
#include "opacity.h"
#include "radius.h"
#include "atomfeii.h"
#include "thirdparty.h"
#include "taulines.h"

const double WLAL = 1215.6845;

/* There are 373 transitions:
 * Wavelength (A)
 * absorption oscillator strength
 * Energy of lower level (Ryd)
 * Statistical weight of lower level (g) */
/** t_fe2ovr_la: constructor storing energy levels for Fred's FeII ground */
t_fe2ovr_la::t_fe2ovr_la()
{
	DEBUG_ENTRY( "t_fe2ovr_la()" );

	const long VERSION_MAGIC = 20070717L;
	static const char chFile[] = "fe2ovr_la.dat";

	FILE *io = open_data( chFile, "r" );

	bool lgErr = false;
	long i = -1L;

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	double help=0;
	for( i=0; i < NFEII; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le", &help ) != 1 );
		fe2lam[i] = (realnum)help;
	}
	for( i=0; i < NFEII; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le", &help ) != 1 );
		fe2osc[i] = (realnum)help;
	}
	for( i=0; i < NFEII; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le", &help ) != 1 );
		fe2enr[i] = (realnum)help;
	}
	for( i=0; i < NFEII; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le", &help ) != 1 );
		fe2gs[i] = (realnum)help;
	}
	for( i=0; i < NFE2PR; i++ )
		lgErr = lgErr || ( fscanf( io, "%le", &fe2pt[i] ) != 1 );
	for( i=0; i < NFE2PR; i++ )
		lgErr = lgErr || ( fscanf( io, "%le", &fe2pf[i] ) != 1 );

	fclose( io );

	ASSERT( !lgErr );
	return;
}

void t_fe2ovr_la::zero_opacity()
{
	DEBUG_ENTRY( "zero_opacity()" );

	for( int i=0; i < NFEII; i++ )
	{
		Fe2PopLte[i] = 0.f;
		feopc[i] = 0.f;
		Fe2TauLte[i] = opac.taumin;
	}
	return;
}

void t_fe2ovr_la::init_pointers()
{
	DEBUG_ENTRY( "init_pointers()" );

	for( int i=0; i < NFEII; i++ )
		ipfe2[i] = ipoint(fe2enr[i]);
	return;
}

/** tau_inc: update line opacities */
void t_fe2ovr_la::tau_inc()
{
	DEBUG_ENTRY( "tau_inc()" );

	for( int i=0; i < NFEII; i++ )
		/* optical depths for Feii dest of lya when large feii not used */
		Fe2TauLte[i] += feopc[i]*(realnum)radius.drad_x_fillfac;
	return;
}

/** atoms_fe2ovr compute FeII overlap with Lya */
void t_fe2ovr_la::atoms_fe2ovr(void)
{
	DEBUG_ENTRY( "atoms_fe2ovr()" );

	long int i;

	static long int nZoneEval;

	double Fe2Partn, 
	  displa, 
	  hopc, 
	  rate, 
	  weight;

	static double BigFeWidth, 
	  BigHWidth;

	/* wavelength of Lya in Angstroms */

	/* compute efficiency of FeII emission overlapping with Ly-alpha
	 * implemented with Fred Hamann
	 *
	 * make Ly-a width monotonically increasing to avoid oscillation
	 * in deep regions of x-ray ionized clouds.
	 *
	 * do nothing if large FeII atom is used */
	if( FeII.lgFeIILargeOn )
	{ 
		return;
	}

	if( nzone <= 1 )
	{
		BigHWidth = hydro.HLineWidth;
		BigFeWidth = GetDopplerWidth(dense.AtomicWeight[ipIRON]);
		nZoneEval = nzone;
	}

	/* do not do pumping if no population,line is thin, or turned off */
	if( (dense.xIonDense[ipIRON][1] <= 0. || !FeII.lgLyaPumpOn) || 
		hydro.HLineWidth <= 0. )
	{
		Fe2Partn = 0.;
		hydro.dstfe2lya = 0.;

		for( i=0; i < NFEII; i++ )
		{
			Fe2PopLte[i] = 0.;
		}
		return;
	}

	/* only evaluate this one time per zone to avoid oscillations
	 * deep in x-ray ionized clouds */
	if( nZoneEval == nzone && nzone > 1 )
	{ 
		return;
	}

	BigHWidth = MAX2(BigHWidth,(double)hydro.HLineWidth);
	BigFeWidth = MAX2(BigFeWidth,(double)GetDopplerWidth(dense.AtomicWeight[ipIRON]) );
	nZoneEval = nzone;

	/* check that data is linked in */
	ASSERT( fe2lam[0] > 0. );

	rate = 0.;
	Fe2Partn = fe2par(phycon.te);
	for( i=0; i < NFEII; i++ )
	{
		/* this is displacement from line center in units of Lya width */
		displa = fabs(fe2lam[i]-WLAL)/WLAL*3e10/BigHWidth;
		if( displa < 1.333 )
		{
			/* have variable weighting factor depending on distance away
			 * this comes form the Verner's fits to Adam's results */
			weight = ( displa <= 0.66666 ) ? 1. : MAX2(0.,1.-(displa-0.666666)/0.66666);

			Fe2PopLte[i] = (realnum)(fe2gs[i]/Fe2Partn*rfield.ContBoltz[ipfe2[i]-1]*
					       dense.xIonDense[ipIRON][1]);

			feopc[i] = (realnum)(Fe2PopLte[i]*fe2osc[i]*0.0150*(fe2lam[i]*1e-8)/BigFeWidth);

			/* Ly-alpha line-center opacity */
			if( StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop > 0. )
			{
				hopc = 7.60e-8*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop/
					GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]);
			}
			else
			{
				hopc = 7.60e-8*dense.xIonDense[ipHYDROGEN][0]/GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]);
			}
			rate += (feopc[i]/SDIV((feopc[i] + hopc)))*(BigFeWidth/GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]))*
				(1. - 1./(1. + 1.6*Fe2TauLte[i]))*weight;
		}
	}

	/* dstfe2lya is total Lya deexcitation probability due to line overlap */
	hydro.dstfe2lya = (realnum)rate;
	return;
}

/** fe2par evaluate FeII partition function */
double t_fe2ovr_la::fe2par(double te)
{
	DEBUG_ENTRY( "fe2par()" );

	double fe2par_v;

	/* function to evaluate partition function for FeII
	 *
	 * Temperature (K) */

	if( te <= fe2pt[0] )
		fe2par_v = fe2pf[0];
	else if( te >= fe2pt[NFE2PR-1] )
		fe2par_v = fe2pf[NFE2PR-1];
	else
	{
		long i = hunt_bisect(fe2pt,NFE2PR,te);
		double slope = (fe2pf[i+1] - fe2pf[i])/(fe2pt[i+1] - fe2pt[i]);
		double du = slope*(te - fe2pt[i]);
		fe2par_v = fe2pf[i] + du;
	}
	return( fe2par_v );
}
