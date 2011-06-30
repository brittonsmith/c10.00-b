/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitSimPostparse initialization at start of simulation, 
 * called from cloudy after parser, 
 * called one time per sim in grid 
 * not called for every iteration */
#include "cddefines.h" 
#include "dense.h"
#include "struc.h"
#include "rfield.h"
#include "elementnames.h"
#include "atoms.h"
#include "physconst.h"
#include "iterations.h"
#include "pressure.h"
#include "mole.h"
#include "trace.h"
#include "radius.h"
#include "thermal.h"
#include "heavy.h"
#include "wind.h"
#include "init.h" 
#include "iso.h"

/*InitSimPostparse initialization after parser, 
 * called one time for each simulation in a grid, 
 * only before first iteration
 * sets initial or zeros values before start of a simulation */
void InitSimPostparse( void )
{
	DEBUG_ENTRY( "InitSimPostparse()" );

	struc.nzonePreviousIteration = -1;

	thermal.thist = 0.;
	thermal.tlowst = 1e20f;
	thermal.nUnstable = 0;
	thermal.lgUnstable = false;

	CO_Init();
	CO_zero();

	/* >> chng 06 Nov 24 rjrw: Move reaction definitions here so they can be controlled by parsed commands */
	CO_create_react();

	/* zero out diffuse Lya emission */
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		for( long ion=0; ion<=nelem; ++ion )
			Heavy.xLyaHeavy[nelem][ion] = 0;

	/* convert STOP RADIUS command to STOP THICKNESS command now that inner radius is known */
	if( radius.StopRadius[0] > 0. )
	{
		for( long j=0; j < iterations.iter_malloc; j++ )
		{
			if( radius.StopRadius[j] > radius.Radius )
				radius.StopThickness[j] = radius.StopRadius[j] - radius.Radius;
			else
			{
				fprintf( ioQQQ, " PROBLEM stopping radius is <= inner radius. Bailing out.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* this term applies centrifugal acceleration if DISK option on wind command */
	wind.DiskRadius = 0;
	if( wind.lgDisk )
		wind.DiskRadius = radius.Radius;

	iterations.lgLastIt = false;

	rfield_opac_zero( 0 , rfield.nupper );

	rfield.lgUSphON = false;

	/* where cloud is optically thick to brems */
	rfield.ipEnergyBremsThin = 0;
	rfield.EnergyBremsThin = 0.;
	rfield.comtot = 0.;
	rfield.comoff = 1.;
	rfield.cmcool = 0.;
	rfield.cinrat = 0.;
	rfield.extin_mag_B_point = 0.;
	rfield.extin_mag_V_point = 0.;
	rfield.extin_mag_B_extended = 0.;
	rfield.extin_mag_V_extended = 0.;
	rfield.EnerGammaRay = 7676.;

	/* usually computed in pressure_change, but not called before first zone */
	wind.AccelGravity = (realnum)(GRAV_CONST*wind.comass*SOLAR_MASS/POW2(radius.Radius)*
		(1.-wind.DiskRadius/radius.Radius) );
	if( trace.lgTrace && trace.lgWind )
	{
		fprintf(ioQQQ," InitSimPostparse sets AccelGravity %.3e lgDisk?%c\n",
			wind.AccelGravity , 
			TorF(wind.lgDisk) );
	}

	/* pressure related variables */
	pressure.PresInteg = 0.;
	pressure.PresIntegElecThin = 0.;
	pressure.pinzon = 0.;

	/* update iso sequence level information */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* these elements that are turned off */
			if( nelem>ipHELIUM && !dense.lgElmtOn[nelem] )
			{
				iso.numLevels_max[ipISO][nelem] = 0;
				iso.nCollapsed_max[ipISO][nelem] = 0;
				iso.n_HighestResolved_max[ipISO][nelem] = 0;

				iso.numLevels_local[ipISO][nelem] = 0;
				iso.nCollapsed_local[ipISO][nelem] = 0;
				iso.n_HighestResolved_local[ipISO][nelem] = 0;
			}
			else
			{
				iso_update_num_levels( ipISO, nelem );
				/* must always have at least one collapsed level, unless compiling recomb data file. */
				ASSERT( iso.nCollapsed_max[ipISO][nelem] > 0 || iso.lgCompileRecomb[ipISO] );
			}
		}
	}

	if( iso.lgPrintNumberOfLevels )
	{
		fprintf(ioQQQ,"\n\n Number of levels in ions treated by iso sequences.\n");
		fprintf(ioQQQ," ISO   Element  hi-n(l-resolved) #(l-resolved)  n(collapsed)\n");
		/* option to print number of levels for each element */
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( long nelem=ipISO; nelem<LIMELM; ++nelem )
			{
				/* print number of levels */
				fprintf(ioQQQ," %s  %s    %4li            %4li           %4li \n",
					iso.chISO[ipISO] ,
					elementnames.chElementSym[nelem],
					iso.n_HighestResolved_max[ipISO][nelem],
					iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem],
					iso.nCollapsed_max[ipISO][nelem]);
			}
		}
	}
	atoms.p2nit = 0.;
	atoms.d5200r = 0.;

	return;
}
