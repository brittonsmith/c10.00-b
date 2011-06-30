/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_driving derive radiative acceleration due to line absorption of incident continuum,
 * return value is line radiative acceleration */
#include "cddefines.h"
#include "physconst.h"
#include "iso.h"
#include "dense.h"
#include "taulines.h"
#include "h2.h"
#include "atomfeii.h"
#include "rt.h"

/*RT_line_driving derive radiative acceleration due to line absorption of incident continuum,
 * return value is line radiative acceleration */
double RT_line_driving(void)
{
	long int i, 
	  ipHi, 
	  nelem, 
	  ipLo,
	  ipISO;

	double AllHeavy, 
	  AllRest, 
	  OneLine, 
	  fe2drive, 
	  forlin_v, 
	  h2drive,
	  accel_iso[NISO];

	/* following used for debugging */
	/* double 
	  RestMax, 
	  HeavMax, 
	  hydromax;
	  long int 
	  ipRestMax, 
	  ihmax; */

	DEBUG_ENTRY( "RT_line_driving()" );

	/* this function finds the total rate the gas absorbs energy
	 * this result is divided by the calling routine to find the
	 * momentum absorbed by the gas, and eventually the radiative acceleration
	 *
	 * the energy absorbed by the line is
	 * Abundance * energy * A *(g_up/g_lo) * occnum * escape prob
	 * where occnum is the photon occupation number, and the g's are
	 * the ratios of statistical weights */

	/* total energy absorbed in this zone, per cubic cm
	 * do hydrogen first */

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		accel_iso[ipISO] = 0;
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( (dense.IonHigh[nelem] >= nelem + 1-ipISO)  )
			{
				for( ipHi=1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
				{
					/* do not put in highest level since its not real */
					for( ipLo=0; ipLo < ipHi - 1; ipLo++ )
					{
						/* do not include bogus lines */
						if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont > 0 )
						{
							OneLine = Transitions[ipISO][nelem][ipHi][ipLo].Emis->pump*
							Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg*
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc;

							accel_iso[ipISO] += OneLine;
						}
					}
				}
				
				// SatelliteLines are indexed by lower level, summed over satellite levels
				for( ipLo=0; ipLo < iso.numLevels_local[ipISO][nelem]; ipLo++ )
				{
					/* do not include bogus lines */
					if( SatelliteLines[ipISO][nelem][ipLo].ipCont > 0 )
					{
						OneLine = SatelliteLines[ipISO][nelem][ipLo].Emis->pump*
						SatelliteLines[ipISO][nelem][ipLo].EnergyErg*
						SatelliteLines[ipISO][nelem][ipLo].Emis->PopOpc;

						accel_iso[ipISO] += OneLine;
					}

				}

				for( ipHi=StatesElemNEW[nelem][nelem-ipISO][iso.numLevels_max[ipISO][nelem]-1].n; ipHi < iso.nLyman[ipISO]; ipHi++ )
				{
					/* do not include bogus lines */
					if( ExtraLymanLines[ipISO][nelem][ipHi].ipCont > 0 )
					{
						OneLine = ExtraLymanLines[ipISO][nelem][ipHi].Emis->pump*
						ExtraLymanLines[ipISO][nelem][ipHi].EnergyErg*
						ExtraLymanLines[ipISO][nelem][ipHi].Emis->PopOpc;

						accel_iso[ipISO] += OneLine;
					}

				}
			}
		}
	}

	/* all heavy element lines in calculation of cooling 
	 * these are the level 1 lines */
	AllHeavy = 0.;
	for( i=1; i <= nLevel1; i++ )
	{
		OneLine = 
		  TauLines[i].Emis->pump*
		  TauLines[i].EnergyErg*
		  TauLines[i].Emis->PopOpc;
		AllHeavy += OneLine;
	}

	/* all heavy element lines treated with g-bar 
	 * these are the level 2 lines, f should be ok */
	AllRest = 0.;
	for( i=0; i < nWindLine; i++ )
	{
		OneLine = 
			TauLine2[i].Emis->pump*
			TauLine2[i].EnergyErg*
			TauLine2[i].Emis->PopOpc;
		AllRest += OneLine;
	}
	for( i=0; i < nUTA; i++ )
	{
		OneLine = 
			UTALines[i].Emis->pump*
			UTALines[i].EnergyErg*
			UTALines[i].Emis->PopOpc;
		AllRest += OneLine;
	}
	for( i=0; i < nHFLines; i++ )
	{
		OneLine = 
			HFLines[i].Emis->pump*
			HFLines[i].EnergyErg*
			HFLines[i].Emis->PopOpc;
		AllRest += OneLine;
	}
	for( i=0; i < linesAdded2; i++)
	{
		OneLine =
			dBaseLines[i].pump*
			dBaseLines[i].tran->EnergyErg*
			dBaseLines[i].PopOpc;
		AllRest += OneLine;
	}

	/* the H2 molecule */
	h2drive = H2_Accel();

	/* The large model FeII atom */
	fe2drive = 0.;
	FeIIAccel(&fe2drive);

	forlin_v = AllHeavy + accel_iso[ipH_LIKE] + accel_iso[ipHE_LIKE] + 
		fe2drive + h2drive + AllRest;

	/*fprintf(ioQQQ," wind te %e %e %e %e %e\n", 	
		AllHeavy , HydroAccel , fe2drive , he1l , AllRest );*/
	return( forlin_v );
}
