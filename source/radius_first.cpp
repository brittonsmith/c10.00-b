/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*radius_first derive thickness of first zone, called after conditions in first zone 
 * are established, sets 
radius.drad_x_fillfac
radius.drad
 */
#include "cddefines.h"
#define	Z	1.0001
#include "wind.h"
#include "stopcalc.h"
#include "thermal.h"
#include "dynamics.h"
#include "trace.h"
#include "save.h"
#include "pressure.h"
#include "iso.h"
#include "h2.h"
#include "rfield.h"
#include "dense.h"
#include "hmi.h"
#include "geometry.h"
#include "opacity.h"
#include "ipoint.h"
#include "radius.h"

void radius_first(void)
{
	long int i ,
		ip;

	bool lgDoPun;

	int indexOfSmallest = 0;

	const int NUM_DR_TYPES = 13;

	struct t_drValues{
		double dr;
		char whatToSay[40];
	} drValues[NUM_DR_TYPES];

	double accel, 
	  BigOpacity, 
	  change,
	  dr912, 
	  drH2 ,
	  drContPres ,
	  drOpacity ,
	  drStromgren, /* used for Stromgren length */
	  drTabDen, 
	  dradf, 
	  drcol, 
	  dr_time_dep,
	  drthrm, 
	  factor, 
	  winddr;
	static double drad_last_iteration=-1.;

	DEBUG_ENTRY( "radius_first()" );

	/***********************************************************************
	 *
	 * wind model, use acceleration length                          
	 *
	 ***********************************************************************/

	if( wind.lgBallistic() )
	{
		/* evaluate total pressure, although value not used (so stuffed into dr912) */
		/* >>chng 01 nov 02, remove call, confirm vals defined with assert */
		ASSERT( dense.pden > 0. && dense.wmole > 0. );
		accel = 1.3e-10*dense.pden*dense.wmole;
		winddr = POW2(wind.windv)/25./accel;
	}
	else
	{
		winddr = 1e30;
	}

	/* key off of Lyman continuum optical depth */
	if( StopCalc.taunu > 0.99 && StopCalc.taunu < 3. )
	{
		dr912 = StopCalc.tauend/6.3e-18/(dense.xIonDense[ipHYDROGEN][0]*geometry.FillFac)*Z/50.;
	}
	else
	{
		dr912 = 1e30;
	}

	if( dynamics.lgTimeDependentStatic && iteration > 2 )
	{
		/* when time dependent case on do not let dr change since current continuum
		 * is not good indicator of conditions */
		dr_time_dep = drad_last_iteration;
	}
	else
	{
		dr_time_dep = 1e30;
	}

	/***********************************************************************
	 *
	 * key off of column density; total, neutral, or ionized                          
	 *
	 ***********************************************************************/

	if( StopCalc.HColStop < 5e29 )
	{
		/* this is useful for very thin columns, normally larger than 1st DR */
		drcol = log10(StopCalc.HColStop) - log10(dense.gas_phase[ipHYDROGEN]*geometry.FillFac* 20.);
	}
	else if( StopCalc.colpls < 5e29 )
	{
		/* ionized column density */
		drcol = log10(StopCalc.colpls) - log10(dense.xIonDense[ipHYDROGEN][1]*geometry.FillFac* 20.);
	}
	else if( StopCalc.colnut < 5e29 )
	{
		/* neutral column denisty */
		drcol = log10(StopCalc.colnut) - log10(dense.xIonDense[ipHYDROGEN][0]*geometry.FillFac*50.);
	}
	else
	{
		/* not used */
		drcol = 30.;
	}
	/* finally convert the drived column density to linear scale */
	drcol = pow(10.,MIN2(35.,drcol));

	/***********************************************************************
	 *
	 * key off of density or abundance fluctuations, must be small part of wavelength                          
	 *
	 ***********************************************************************/

	if( dense.flong != 0. )
	{
		/* flong set => density fluctuations */
		dradf = 6.283/dense.flong/10.;
		dradf = MIN4(dradf,radius.StopThickness[iteration-1]*Z,drcol,dr912);
	}
	else
	{
		dradf = FLT_MAX;
	}

	/* >>>chng 99 nov 18, add check on stromgren length */
	/* estimate Stromgren length, but only if there are ionizing photons
	 * and not constant temperature model */
	if( (rfield.qhtot>0.) && (rfield.qhtot> rfield.qbal*0.01) && (rfield.uh>1e-10) )
	{
		/* >>chng 99 dec 23, double to allow lte.in to work on alphas */
		/* >>chng 03 mar 15, density to double to avoid overflow, PvH */
		drStromgren = (double)(rfield.qhtot)/iso.RadRec_caseB[ipH_LIKE][ipHYDROGEN]/
			POW2((double)dense.gas_phase[ipHYDROGEN]);

		// Hubble radius is largest value
		drStromgren = MIN2(1e28 , drStromgren );

		/* different logic if this is a sphere */
		if( drStromgren/radius.rinner > 1. )
		{
			/* >>chng 03 mar 15, to double to avoid FP overflow, PvH */
			drStromgren = (double)rfield.qhtot*3./(radius.rinner*
			    iso.RadRec_caseB[ipH_LIKE][ipHYDROGEN]*POW2((double)dense.gas_phase[ipHYDROGEN]) );
			drStromgren += 1.;
			/* this results in r_out / r_in */
			drStromgren = pow( drStromgren , 0.33333);
			/* make it a physics thickness in cm */
			drStromgren *= radius.rinner;
		}

		/* remember the Stromgren thickness */
		radius.thickness_stromgren = (realnum)drStromgren;

		/* take one hundredth of this, do nothing if near underflow */
		if( drStromgren > SMALLFLOAT *100.)
			drStromgren /= 100.;
	}
	else
	{
		drStromgren = FLT_MAX;
		radius.thickness_stromgren = FLT_MAX;
	}

	/***********************************************************************
	 *
	 * find largest opacity, to keep the first zone optical depth 1
	 * this is usually the physics that sets the first zone thickness
	 *
	 ***********************************************************************/

	/* >>>chng 99 jun 25, this is to simulate behavior of code before extension
	 * of continuum array to 1e-8 Ryd */
	ip = ipoint(1e-5);

	/* find largest opacity */
	BigOpacity = 0.;
	for( i=ip; i < rfield.nflux; i++ )
	{
		/* remember largest opacity, and energy where this happened,
		 * make sure flux at energy is gt 0, can be zero for case where
		 * nflux increased to include emission from some ions */
		if( rfield.flux[0][i]>0. && opac.opacity_abs[i] > BigOpacity )
		{
			BigOpacity = opac.opacity_abs[i];
		}
	}
	/* BigOpacity may be zero on very first call */

	/* drChange set with set didz command, is only number set with this command,
	 * default in zerologic is 0.15 
	 * set drad to small part of*/
	if( BigOpacity > SMALLFLOAT )
	{
		drOpacity = (radius.drChange/100.)/BigOpacity/geometry.FillFac;
	}
	else
	{
		drOpacity = 1e30;
	}

	/***********************************************************************
	 *
	 * thermalization length of typical lines
	 *
	 ***********************************************************************/

	drthrm = 1.5e31/MAX2(1.,POW2((double)dense.gas_phase[ipHYDROGEN]));
	/* thermalization length is irrelevant above critical density, which we
	 * regard as 1e16 cm-3.  Put a floor at the corresponding drthrm. */
	drthrm = MAX2( 0.15, drthrm );

	/***********************************************************************
	 *
	 * make sure we resolve initial structure in dense_tabden command
	 * if interpolated table we need to make sure that we resolve the 
	 * initial changes in the structure
	 *
	 ***********************************************************************/

	if( strcmp(dense.chDenseLaw,"DLW2") == 0 )
	{
		drTabDen = 1.;
		i = 1;
		factor = 0.;
		while( i < 100 && factor < 0.05 && radius.Radius+drTabDen*2.<radius.StopThickness[0] )
		{
			/* check densities at ever larger dr's, until factor becomes more than 5% */
			factor = dense.gas_phase[ipHYDROGEN]/
				dense_tabden(radius.Radius+drTabDen, drTabDen );
			/* density change can be positive or negative sign */
			factor = fabs(factor-1.);
			drTabDen *= 2.;
			i += 1;
		}
		drTabDen /= 2.;
	}
	else
	{
		drTabDen = 1e30;
	}

	/* >>chng 03 mar 20, add check on lyman band optical depth - want first zone
	 * to be thin in H2 bands */
	/* some tests are fully molecular with solomon process turned off,
	 * do not sense this when already almost fully molecular */
	if( hmi.H2_total/dense.gas_phase[ipHYDROGEN] < 0.1 )
	{
		change = 0.1;
	}
	else
	{
		/* >>chng 04 mar 14, this branch, H is quite molecular,
		 * still do not want large changes in solomon rate since linearization
		 * would not work in hmole network, bu do not need such fine steps */
		change = 1.;
	}

	/* >>chng 04 mar 14 go back to original logic since molecular
	 * pdr's had big jump in conditions from
	 * first to second zon even when most H in H2
	change = 0.1; */
	/* >>chng 04 apr 18, change from 0.1 to 0.001, inital zones too large in 
	 * leiden test case f1 */
	change = 0.001;
	/* >>chng 04 mar 13, not too large when big H2 is on */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated )
	{
		if( fabs(hmi.HeatH2Dexc_BigH2)/thermal.ctot > 0.05 )
		{
			/* changes in H2 heating caused by changes in solomon rate
			 * would drive temperature failures */
			/* >>chng 04 apr 18, change from 0.001 to 0.0001, inital zones too large in 
			 * leiden test case f1 */
			change = 0.0001;
		}
		else
		{
			/* >>chng 04 apr 18, change from 0.01 to 0.001, inital zones too large in 
			 * leiden test case f1 */
			change = 0.001;
		}
	}
	drH2 = change / SDIV( 
		hmi.H2_total * geometry.FillFac * hmi.H2Opacity );

	/* >>chng 06 feb 01, very high U ulirg models had dramatic increase in
	 * cont pre in first few zones,
	 * in constant total pressure case, don't want acceleration across first zone to
	 * be large compared with current gas pressure */
	if( (strcmp( dense.chDenseLaw, "CPRE" )==0) && pressure.lgContRadPresOn )
	{
		/* radiative acceleration was evaluated in PressureTotal */
		drContPres = 0.05 * pressure.PresTotlCurr / 
			((double)wind.AccelTotalOutward*dense.xMassDensity*geometry.FillFac*geometry.DirectionalCosin);
	}
	else if( !wind.lgStatic() )
	{
		/* acceleration and change in v in wind */
		double g = fabs(wind.AccelTotalOutward-wind.AccelGravity);
		/* wind - do not let velocity change by too much */
		drContPres = 0.05*POW2(wind.windv)/(2.*SDIV(g));
	}
	else
		drContPres = 1e30;

	drValues[0].dr  = drOpacity;
	drValues[1].dr  = radius.Radius/20.;
	drValues[2].dr  = drStromgren;
	drValues[3].dr  = radius.StopThickness[iteration-1]/10.;
	drValues[4].dr  = drcol;
	drValues[5].dr  = dr912;
	drValues[6].dr  = drthrm;
	drValues[7].dr  = winddr;
	drValues[8].dr  = dradf;
	drValues[9].dr  = drTabDen;
	drValues[10].dr = drH2;
	drValues[11].dr = drContPres;
	drValues[12].dr = dr_time_dep;

	strcpy( drValues[0].whatToSay,  "drOpacity" );
	strcpy( drValues[1].whatToSay,  "radius.Radius/20.");
	strcpy( drValues[2].whatToSay,  "drStromgren");
	strcpy( drValues[3].whatToSay,  "radius.StopThickness[iteration-1]/10.");
	strcpy( drValues[4].whatToSay,  "drcol");
	strcpy( drValues[5].whatToSay,  "dr912");
	strcpy( drValues[6].whatToSay,  "drthrm");
	strcpy( drValues[7].whatToSay,  "winddr");
	strcpy( drValues[8].whatToSay,  "dradf");
	strcpy( drValues[9].whatToSay,  "drTabDen");
	strcpy( drValues[10].whatToSay, "drH2");
	strcpy( drValues[11].whatToSay, "drContPres");
	strcpy( drValues[12].whatToSay, "dr_time_dep");

	for( i=0; i<NUM_DR_TYPES; i++ )
	{
		if( drValues[i].dr < drValues[indexOfSmallest].dr )
		{
			indexOfSmallest = i;
		}
	}

	radius.drad = drValues[indexOfSmallest].dr;

	double rfacmin = radius.lgSdrminRel ? radius.Radius : 1.;
	/* reset if radius.drad is less than radius.sdrmin */
	if( rfacmin*radius.sdrmin >= radius.drad )
	{
		radius.drad = rfacmin*radius.sdrmin;
		/* set flag for comment if the previous line forced a larger dr than
		 * would otherwise have been chosen.  will cause comment to be generated
		 * in PrtComment if set true*/
		radius.lgDR2Big = true;
	}
	else
	{
		radius.lgDR2Big = false;
	}

	/* this min had been in the big min set above, but caused a false alarm
	 * on the lgDR2Big test above since the set dr command sets both in and max */
	// lgSdrmaxRel true if sdrmax is relative to current radius, false if limit in cm
	double rfacmax = radius.lgSdrmaxRel ? radius.Radius : 1.;
	radius.drad = MIN2( rfacmax*radius.sdrmax, radius.drad );
	radius.drad_mid_zone = radius.drad/2.;

#if	0
	/***********************************************************************
	 *
	 * we have now generated range of estimates of first thickness,
	 * now choose smallest of the group
	 *
	 ***********************************************************************/

	/* radius div by 20, to prevent big change in iron ionization for high ioniz gas,
	 * this is also the ONLY place that sphericity comes in */
	radius.drad = MIN4( MIN3( drOpacity, radius.Radius/20., drStromgren ),
			    MIN3( radius.StopThickness[iteration-1]/10., drcol, dr912 ),
			    MIN4( drthrm, winddr, dradf, drTabDen ),
			    MIN3( drH2, drContPres, dr_time_dep ) );

	/* option to set lower limit to zone thickness, with set drmin command*/
	radius.drad = MAX2( radius.drad, rfacmin*radius.sdrmin );

	/* set flag for comment if the previous line forced a larger dr than
	 * would otherwise have been chosen.  will cause comment to be generated
	 * in PrtComment if set true*/
	if( fp_equal( radius.drad, rfacmin*radius.sdrmin ) )
	{
		radius.lgDR2Big = true;
	}
	else
	{
		radius.lgDR2Big = false;
	}

	/* this min had been in the big min set above, but caused a false alarm
	 * on the lgDR2Big test above since the set dr command sets both in and max */
	radius.drad = MIN2( rfacmax*radius.sdrmax, radius.drad );
#endif

	radius.drad_x_fillfac = radius.drad * geometry.FillFac;

	/* save dr for this iteration */
	drad_last_iteration = radius.drad;

	/* drMinimum is smallest acceptable DRAD, and is 1/100 OF DRAD(1) */
	/* this can be turned off by GLOB command */
	if( radius.lgDrMnOn )
	{
		/* >>chng 05 mar 05, drMinimum is now drad * hden, to make propro to optical depth
		 * avoid false trigger across thermal fronts 
		 * add * dense.gas_phase */
		/* NB - drMinimum not used in code - delete? */
		radius.drMinimum = (realnum)(radius.drad * dense.gas_phase[ipHYDROGEN]/1e7);
	}
	else
	{
		radius.drMinimum = 0.;
	}

	/* if set drmin is used, make sure drMinimum (which will cause an abort) is
	 * smaller than drmin */
	if( radius.lgSMinON )
	{
		/* >>chng 05 mar 05, drMinimum is now drad * hden, to make propro to optical depth
		 * avoid false trigger across thermal fronts 
		 * add * dense.gas_phase */
		/* NB - drMinimum not used in code - delete? */
		radius.drMinimum = MIN2(radius.drMinimum * dense.gas_phase[ipHYDROGEN],
					(realnum)(rfacmin*radius.sdrmin/10.f) );
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" radius_first called, finds dr=%13.5e drMinimum=%12.3e sdrmin=%10.2e sdrmax=%10.2e\n", 
			 radius.drad, radius.drMinimum/ dense.gas_phase[ipHYDROGEN],
			 rfacmin*radius.sdrmin, rfacmax*radius.sdrmax );
	}

	if( radius.drad < SMALLFLOAT*1.1 )
	{
		fprintf( ioQQQ, 
			" PROBLEM radius_first detected likely insanity, found dr=%13.5e \n", radius.drad);
		fprintf( ioQQQ, 
			" radius_first: calculation continuing but crash is likely. \n");
		/* this sets flag that insanity has occurred */
		TotalInsanity();
	}

	/* all this is to only save on last iteration
	 * the save dr command is not really a save command, making this necessary
	 * lgDRon is set true if "save dr" entered */
	if( save.lgDROn )
	{
		lgDoPun = true;
	}
	else
	{
		lgDoPun = false;
	}

	/* save what we decided up? */
	if( lgDoPun )
	{
		/* create hash marks on second and later iterations */
		if( iteration > 1 && save.lgDRHash )
		{
			static int iter_punch=-1;
			if( iteration !=iter_punch )
				fprintf( save.ipDRout, "%s\n",save.chHashString );
			iter_punch = iteration;
		}
		/* this is common part of each line, the zone count, depth, chosen dr, and depth2go */
		/* >>chng 05 aug 15, had printed drNext, always zero, rather the drad, which is set here */
		fprintf( save.ipDRout , "%ld\t%.5e\t%.3e\t%.3e\t", nzone, radius.depth, radius.drad, radius.Depth2Go );

		if( radius.lgDR2Big )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from radius.sdrmin\n");

		}
		else if( fp_equal( radius.drad, rfacmax*radius.sdrmax ) )
		{

			fprintf( save.ipDRout, 
				"radius_first keys from radius.sdrmax\n");
		}
		else
		{
			ASSERT( indexOfSmallest < NUM_DR_TYPES - 1 );
			fprintf( save.ipDRout, "radius_first keys from %s\n", 
				drValues[indexOfSmallest].whatToSay);
		}

		/* \todo	1 improve this printout and the drValues treatment above. */

#if	0
		if( fp_equal( radius.drad, drOpacity ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from drOpacity, opac was %.2e at %.2e Ryd\n", 
			   BigOpacity , BigOpacityAnu );
		}
		else if( fp_equal( radius.drad, radius.Radius/20. ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from radius.Radius\n" );
		}
		else if( fp_equal( radius.drad, drStromgren ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from drStromgren\n");
		}
		else if( fp_equal( radius.drad, dr_time_dep ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from time dependent\n");
		}
		else if( fp_equal( radius.drad, radius.StopThickness[iteration-1]/10. ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from radius.StopThickness[iteration-1]\n");
		}
		else if( fp_equal( radius.drad, drcol ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from drcol\n");
		}
		else if( fp_equal( radius.drad, rfacmin*radius.sdrmin ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from radius.sdrmin\n");
		}
		else if( fp_equal( radius.drad, dr912 ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from dr912\n");
		}
		else if( fp_equal( radius.drad, rfacmax*radius.sdrmax ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from radius.sdrmax\n");
		}
		else if( fp_equal( radius.drad, drthrm ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from drthrm\n");
		}
		else if( fp_equal( radius.drad, winddr ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from winddr\n");
		}
		else if( fp_equal( radius.drad, drH2 ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from H2 lyman lines\n");
		}
		else if( fp_equal( radius.drad, dradf ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from dradf\n");
		}
		else if( fp_equal( radius.drad, drTabDen ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from drTabDen\n");
		}
		else if( fp_equal( radius.drad, drContPres ) )
		{
			fprintf( save.ipDRout, 
				"radius_first keys from radiative acceleration across zone\n");
		}
		else
		{
			fprintf( save.ipDRout,  "radius_first insanity\n" );
			fprintf( ioQQQ, "radius_first insanity, radius is %e\n" ,
				radius.drad);
			ShowMe();
		}
#endif

	}
	return;
}
