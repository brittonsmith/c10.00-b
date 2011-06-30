/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PressureChange called by ConvPresTempEdenIoniz
 * evaluate the current pressure, change needed to get it to converge,
 * the global static variable pressure_change_factor
 * applies this correction factor to all gas constituents,
 * sets conv.lgConvPres true if good pressure, false if pressure change capped */
/*lgConvPres finds amount by which density should change to move towards pressure equilibrium
 * returns true if pressure is converged */
#include "cddefines.h"
#include "abund.h"
#include "hmi.h"
#include "struc.h"
#include "trace.h"
#include "wind.h"
#include "phycon.h"
#include "thermal.h"
#include "dense.h"
#include "geometry.h"
#include "radius.h"
#include "mole.h"
#include "dynamics.h"
#include "pressure.h"
#include "colden.h"
#include "conv.h"
#include "cosmology.h"
#include "dark_matter.h"

/* this is the pressure change pressure_change_factor, needed to move current pressure to correct pressure */
static double pressure_change_factor;

/*PressureChange evaluate the current pressure, and change needed to
 * get it to PresTotlInit,
 * return value is true is density was changed, false if no changes were necessary */
int PressureChange( 
	/* this is change factor, 1 at first, becomes smaller as oscillations occur */
	double dP_chng_factor )
{
	long int ion, 
	  nelem,
	  lgChange,
	  mol;

	double abun, 
	  edensave, 
	  hold;

	/* biggest multiplicative change in the pressure that we will allow */
	/* allowed error in the pressure is conv.PressureErrorAllowed*/
	double pdelta;

	static double FacAbun, 
	  FacAbunSav, 
	  OldAbun;

	DEBUG_ENTRY( "PressureChange()" );

	edensave = dense.eden;

	/* first evaluate total pressure for this location, and current conditions
	 * CurrentPressure is just sum of gas plus local line radiation pressure */
	/* this sets values of pressure.PresTotlCurr */
	PresTotCurrent();

	/* this will save the  history of density - pressure relationship
	 * for the current zone */
	if( nzone != conv.hist_pres_nzone )
	{
		/* first time in this zone - reset history */
		conv.hist_pres_nzone = nzone;
		conv.hist_pres_density.clear();
		conv.hist_pres_current.clear();
		conv.hist_pres_correct.clear();
	}

	/* >>chng 04 feb 11, add option to remember current density and pressure */
	conv.hist_pres_density.push_back( dense.gas_phase[ipHYDROGEN] );
	conv.hist_pres_current.push_back( pressure.PresTotlCurr );
	conv.hist_pres_correct.push_back( pressure.PresTotlCorrect );

	/* remember old hydrogen density */
	hold = dense.gas_phase[ipHYDROGEN];

	/* this will be set true if density or abundances change in this zone */
	lgChange = false;

	/* this evaluates current pressure, sets pressure_change_factor and 
	 * pressure.PresTotlCorrect, updates velocity,
	 * and returns true if pressure has converged 
	 * sets pressure.PresTotlCorrect */
	conv.lgConvPres = lgConvPres();
	/*fprintf(ioQQQ,"DEBUG pressure nH %.3e init %.3e correct %.3e currnt %.3e \n",
		dense.gas_phase[ipHYDROGEN] ,
		pressure.PresTotlInit,
		pressure.PresTotlCorrect ,
		pressure.PresTotlCurr);*/

	/* if convergence is OK at present state, so no change reqd, simply return 
	 * >>chng 05 feb 04, cannot do this test here since variable abundances, test has
	 * not been done, so variable abundances did not work,
	 * caught by Marcelo Castellanos and David Valls-Gabaud 
	if( conv.lgConvPres )
		return false;*/

	/* >> chng 02 dec 13 rjrw: short-circuit if nothing changes */
	if( pressure_change_factor != 1. )
	{
		lgChange = true;
	}

	/* allow 3 percent changes,
	 * dP_chng_factor is initially 1, becomes smaller if sign of pressure change, changes */
	pdelta = 0.03 * dP_chng_factor;

	/* make sure that change is not too extreme */
	pressure_change_factor = MIN2(pressure_change_factor,1.+pdelta);
	pressure_change_factor = MAX2(pressure_change_factor,1.-pdelta);

	{
		/*@-redef@*/
		enum{DEBUG_LOC=false};
		static long int nsave=-1;
		/*@+redef@*/
		if( DEBUG_LOC /*&& nzone > 150 && iteration > 1*/ )
		{
			if( nsave-nzone ) fprintf(ioQQQ,"\n");
			nsave = nzone;
			fprintf(ioQQQ,"nnzzone\t%li\t%.2f%%\t%.3f\t%.2e\t%.2e\t%.2e\t%.2e\n", 
				nzone,
				pressure_change_factor,
				/* when this is negative we need to raise the density */
				(pressure.PresTotlCorrect-pressure.PresTotlCurr)*100./pressure.PresTotlCorrect, 
				pressure.PresTotlCorrect,
			    pressure.PresTotlCurr, 
				pressure.PresGasCurr,
				pressure.pres_radiation_lines_curr );
		}
	}

	/* >>chng 97 jun 03, added variable abundances for element table command
	 * option to get abundances off a table with element read command */
	if( abund.lgAbTaON )
	{
		lgChange = true;
		for( nelem=1; nelem < LIMELM; nelem++ )
		{
			if( abund.lgAbunTabl[nelem] )
			{
				abun = AbundancesTable(radius.Radius,radius.depth,nelem+1)*
				  dense.gas_phase[ipHYDROGEN];

				hold = abun/dense.gas_phase[nelem];
				dense.gas_phase[nelem] = (realnum)abun;

				for( ion=0; ion < (nelem + 2); ion++ )
				{
					dense.xIonDense[nelem][ion] *= (realnum)hold;
				}
			}
		}
	}

	/* this is set false if fluctuations abundances command entered,
	 * when density variations are in place this is true,
	 * and dense.chDenseLaw is "SINE" */
	if( !dense.lgDenFlucOn )
	{
		/* abundance variations are in place */
		lgChange = true;
		if( nzone <= 1 )
		{
			OldAbun = 1.;
			FacAbun = 1.;
			if( dense.lgDenFlucRadius )
			{
				/* cycle over radius */
				FacAbunSav = dense.cfirst*cos(radius.depth*dense.flong+
				dense.flcPhase) + dense.csecnd;
			}
			else
			{
				/* cycle over column density */
				FacAbunSav = dense.cfirst*cos( colden.colden[ipCOL_HTOT]*dense.flong+
				dense.flcPhase) + dense.csecnd;
			}
		}
		else
		{
			OldAbun = FacAbunSav;
			/* rapid abundances fluctuation */
			if( dense.lgDenFlucRadius )
			{
				/* cycle over radius */
				FacAbunSav = dense.cfirst*cos(radius.depth*dense.flong+
				dense.flcPhase) + dense.csecnd;
			}
			else
			{
				/* cycle over column density */
				FacAbunSav = dense.cfirst*cos( colden.colden[ipCOL_HTOT]*dense.flong+
				dense.flcPhase) + dense.csecnd;
			}
			FacAbun = FacAbunSav/OldAbun;
		}
	}
	else
	{
		/* abundance variations are NOT in place */
		FacAbun = 1.;
	}

	/* chng 02 dec 11 rjrw -- remove test, saves marginal time could generate nasty intermittent bug when 
	 * ( pressure_change_factor*FacAbun == 1 ) && (pressure_change_factor != 1) */
	if( lgChange )
	{
		/* H, He not affected by abundance fluctuations, so only change these
		 * by the pressure change factor */
		for( nelem=ipHYDROGEN; nelem <= ipHELIUM; ++nelem )
		{
			dense.xMolecules[nelem] *= (realnum)pressure_change_factor;
			dense.gas_phase[nelem] *= (realnum)pressure_change_factor;
			for( ion=0; ion < (nelem + 2); ion++ )
			{
				/* >>chng 96 jul 12 had only multiplied total abun, not ions */
				dense.xIonDense[nelem][ion] *= (realnum)pressure_change_factor;
			}
		}

		/* Li on up is affect by both pressure and abundance variations,
		 * so multiply by both factors */
		for( nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
		{
			dense.xMolecules[nelem] *= (realnum)(pressure_change_factor*FacAbun);
			dense.gas_phase[nelem] *= (realnum)(pressure_change_factor*FacAbun);
			for( ion=0; ion < (nelem + 2); ion++ )
			{
				/* >>chng 96 jul 12 had only multiplied total abun, not ions */
				dense.xIonDense[nelem][ion] *= (realnum)(pressure_change_factor*FacAbun);
			}
		}

		dense.eden *= pressure_change_factor;

		/* must call TempChange since ionization has changed, there are some
		 * terms that affect collision rates (H0 term in electron collision) */
		TempChange(phycon.te , false);

		/* molecules done in hmole, only change pressure, not abundances */
		for(mol = 0; mol < N_H_MOLEC; mol++) 
		{
			hmi.Hmolec[mol] *= (realnum)pressure_change_factor;
		}
		hmi.H2_total *= (realnum)pressure_change_factor;

		/* molecules done in comole */
		/* >>chng 03 sep 15, upper limit had not included the C, O atoms/ions */
		/*for( ion=0; ion < NUM_HEAVY_MOLEC; ion++ )*/
		for( ion=0; ion < mole.num_comole_calc; ion++ )
		{
			COmole[ion]->hevmol *= (realnum)(pressure_change_factor*FacAbun);

			/* check for NaN */
			ASSERT( !isnan( COmole[ion]->hevmol ) );
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
		  " PressureChange called, changing HDEN from %10.3e to %10.3e Set fill fac to %10.3e\n", 
		  hold, dense.gas_phase[ipHYDROGEN], geometry.FillFac );

		if( trace.lgNeBug )
		{
			fprintf( ioQQQ, " EDEN change PressureChange from to %10.3e %10.3e %10.3e\n", 
			  edensave, dense.eden, edensave/dense.eden );
		}
	}

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && (nzone>215)/**/ )
		{
			fprintf( ioQQQ, 
				"%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%c\n", 
			  radius.depth, 
			  pressure.PresTotlCurr, 
			  pressure.PresTotlInit + pressure.PresInteg, 
			  pressure.PresTotlInit, 
			  pressure.PresGasCurr, 
			  pressure.PresRamCurr, 
			  pressure.pres_radiation_lines_curr, 
			  /* subtract continuum rad pres which has already been added on */
			  pressure.PresInteg - pressure.pinzon, 
			  wind.windv/1e5,
			  sqrt(5.*pressure.PresGasCurr/3./dense.xMassDensity)/1e5,
			  TorF(conv.lgConvPres) );
		}
	}

	return lgChange;
}

/*lgConvPres finds amount by which density should change to move towards pressure equilibrium
 * sets pressure.PresTotlCorrect
 * returns true if pressure is converged */
bool lgConvPres(void)
{
	double dnew, 
	  term;
	bool lgRet;

	DEBUG_ENTRY( "lgConvPres()" );

	/* make sure this is set by one of the following branches - set to zero here, 
	 * then assert that it is greater than zero at end */
	pressure.PresTotlCorrect = 0.;

	/* evaluate a series of possible pressure options, and set the file static variable
	 * pressure_change_factor */
	/* inside out globule */
	if( strcmp(dense.chDenseLaw,"GLOB") == 0 )
	{
		/* GLBDST is distance from globule, or glbrad-DEPTH */
		if( radius.glbdst < 0. )
		{
			fprintf( ioQQQ, " Globule distance is negative, internal overflow has occured,  sorry.\n" );
			fprintf( ioQQQ, " This is routine lgConvPres, GLBDST is%10.2e\n", 
			  radius.glbdst );
			cdEXIT(EXIT_FAILURE);
		}
		pressure_change_factor = (radius.glbden*pow(radius.glbrad/(radius.glbdst),radius.glbpow))/
		  dense.gas_phase[ipHYDROGEN];
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}
	else if( cosmology.lgDo )
	{
		/* cosmological - density varies because of expansion of universe */
		dnew = GetDensity( cosmology.redshift_current );
		pressure_change_factor = dnew/dense.gas_phase[ipHYDROGEN];
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}
	else if( (strcmp(dense.chDenseLaw,"WIND") == 0) ) {

		/* this is impossible - wind with identically zero velocity */
		if ( wind.lgStatic() )
		{
			/* declare insanity  */
			fprintf( ioQQQ, " PROBLEM WIND called with zero velocity - this is impossible.\n Sorry.\n" );
			/* TotalInsanity announces fatal problem, ShowMe, then cdEXIT with failure */
			TotalInsanity();
		}
		
		/* this is positive wind velocity the outflowing wind beyond sonic point */
		else if( wind.lgBallistic() )
		{
			
			/* this is logic for supersonic outflowing wind solution,
			 * which assumes positive velocity, well above sonic point.
			 * following makes sure wind v only updated once per zone */
			if( /*! fp_equal( radius.depth, rsave ) &&*/ nzone > 1 )
			{
				/* Wind model */
				
				if( trace.lgTrace && trace.lgWind )
				{
					fprintf(ioQQQ," lgConvPres sets AccelGravity %.3e lgDisk?%c\n",
							  wind.AccelGravity , 
							  TorF(wind.lgDisk) );
				}
				
				/* following is form of energy equation
				 * struc.windv[nzone-2] is velocity of previous zone
				 * this increments that velocity to form square of new wind 
				 * velocity for outer edge of this zone */
				term = POW2(struc.windv[nzone-2]) + 2.*(wind.AccelTotalOutward - wind.AccelGravity)* radius.drad;
				
				/* increment velocity if it is substantially positive */
				if( term <= 1e3 )
				{
					/* wind velocity is well below sonic point, give up, 
					 * do not change velocity */
					wind.lgVelPos = false;
				}
				else
				{
					/* struc.windv[nzone-2] is velocity of previous zone
					 * this increments that velocity to form square of new wind 
					 * velocity for outer edge of this zone 
					 double windnw = (double)POW2(struc.windv[nzone-2]) + 
					 (double)(2.*(wind.AccelTotalOutward-wind.AccelGravity))*radius.drad;*/
					
					/* wind.windv is velocity at OUTER edge of this zone */
					term = sqrt(term);
					if (wind.windv > 0)
					{
						wind.windv = (realnum) term;
					}
					else
					{
						wind.windv = -(realnum) term;
					}
					wind.lgVelPos = true;
				}
				
				if( trace.lgTrace && trace.lgWind )
				{
					fprintf(ioQQQ," lgConvPres new wind V zn%li %.3e AccelTotalOutward %.3e AccelGravity %.3e\n",
							  nzone,wind.windv, wind.AccelTotalOutward, wind.AccelGravity );
				}
			}
			else
			{
				pressure_change_factor = 1.;
			}
			
			/* conservation of mass sets density here */
			pressure_change_factor = wind.emdot/(wind.windv*dense.gas_phase[ipHYDROGEN])/radius.r1r0sq;
			pressure.PresTotlCorrect = pressure.PresTotlCurr * pressure_change_factor;
		}
		
		/* this is negative wind velocity the new dynamics */
		else 
		{
			/* sets pressure.PresTotlCorrect  */
			pressure_change_factor = DynaPresChngFactor();
		}
	}

	else if( strcmp(dense.chDenseLaw,"SINE") == 0 )
	{
		/* rapid density fluctuation */
		if( dense.lgDenFlucRadius )
		{
			pressure_change_factor = (dense.cfirst*cos(radius.depth*dense.flong+dense.flcPhase) + 
			dense.csecnd)/dense.gas_phase[ipHYDROGEN];
		}
		else
		{
			pressure_change_factor = (dense.cfirst*cos(colden.colden[ipCOL_HTOT]*dense.flong+dense.flcPhase) + 
			dense.csecnd)/dense.gas_phase[ipHYDROGEN];
		}
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}

	else if( strcmp(dense.chDenseLaw,"POWR") == 0 )
	{
		/* power law function of radius */
		dnew = dense.den0*pow(radius.Radius/radius.rinner,(double)dense.DensityPower);
		pressure_change_factor = dnew/dense.gas_phase[ipHYDROGEN];
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}

	else if( strcmp(dense.chDenseLaw,"POWD") == 0 )
	{
		/* power law function of depth */
		dnew = dense.den0*pow(1. + radius.depth/dense.rscale,(double)dense.DensityPower);
		pressure_change_factor = dnew/dense.gas_phase[ipHYDROGEN];
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}

	else if( strcmp(dense.chDenseLaw,"POWC") == 0 )
	{
		/* power law function of column density */
		dnew = dense.den0*pow(1.f + colden.colden[ipCOL_HTOT]/
		  dense.rscale,dense.DensityPower);
		pressure_change_factor = dnew/dense.gas_phase[ipHYDROGEN];
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}

	else if( strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		/* constant pressure */
		if( pressure.lgContRadPresOn )
		{
			/* >>chng 01 oct 31, replace pneed with CorretPressure */
			/* this has pressure due to incident continuum */
			pressure.PresTotlCorrect = pressure.PresTotlInit + pressure.PresInteg;
		}
		else
		{
			/* this does not have pressure due to incident continuum*/
			pressure.PresTotlCorrect = pressure.PresTotlInit*
				/* following term normally unity, power law set with option par on cmmnd*/
				pow(radius.Radius/radius.rinner,(double)pressure.PresPowerlaw);
		}

		if( dark.lgNFW_Set || pressure.gravity_symmetry >=0 )
		{
			fixit();  // doesn't this exclude current zone gravity pressure.RhoGravity?
				  // and doesn't the above exclude current rad press, pressure.pinzon?
			pressure.PresTotlCorrect += pressure.IntegRhoGravity;
		}

		/* ratio of correct to current pressures */
		pressure_change_factor = pressure.PresTotlCorrect / pressure.PresTotlCurr;
	}

	else if( strncmp( dense.chDenseLaw ,"DLW" , 3) == 0 )
	{
		if( strcmp(dense.chDenseLaw,"DLW1") == 0 )
		{
			/* call ACF sub */
			pressure_change_factor = dense_fabden(radius.Radius,radius.depth)/dense.gas_phase[ipHYDROGEN];
		}
		else if( strcmp(dense.chDenseLaw,"DLW2") == 0 )
		{
			/* call table interpolation subroutine
			 * >>chng 96 nov 29, added dense_tabden */
			pressure_change_factor = dense_tabden(radius.Radius,radius.depth)/dense.gas_phase[ipHYDROGEN];
		}
		else if( strcmp(dense.chDenseLaw,"DLW3") == 0 )
		{
			/* call parametrized wind subroutine */
			pressure_change_factor = dense_parametric_wind(radius.Radius)/dense.gas_phase[ipHYDROGEN];
		}
		else
		{
			fprintf( ioQQQ, " Insanity, lgConvPres gets chCPres=%4.4s\n", 
			  dense.chDenseLaw );
			cdEXIT(EXIT_FAILURE);
		}
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
		/* >>chng 06 feb 21, from above to below.  we want to look at change,
		 * not the value.  Bug found by Valentina Luridiana */
		/*if( pressure.PresTotlCorrect > 3. || pressure.PresTotlCorrect< 1./3 )*/
		if( pressure_change_factor > 3. || pressure_change_factor< 1./3 )
		{
			static bool lgWARN2BIG=false;
			if( !lgWARN2BIG )
			{
				lgWARN2BIG = true;
				fprintf(ioQQQ,"\n\n >========== Warning!  The tabulated or functional change in density as a function of depth was VERY large. This is zone %li.\n",nzone);
				fprintf(ioQQQ," >========== Warning!  This will cause convergence problems. \n");
				fprintf(ioQQQ," >========== Warning!  The current radius is %.3e. \n",radius.Radius);
				fprintf(ioQQQ," >========== Warning!  The current depth is %.3e. \n",radius.depth);
				fprintf(ioQQQ," >========== Warning!  Consider using a more modest change in density vs radius. \n\n\n");
			}
		}
	}

	else if( strcmp(dense.chDenseLaw,"CDEN") == 0 )
	{
		/* this is the default, constant density */
		pressure_change_factor = 1.;
		pressure.PresTotlCorrect = pressure.PresTotlCurr*pressure_change_factor;
	}

	else
	{
		fprintf( ioQQQ, " Unknown pressure law=%s= This is a major internal error.\n", 
		  dense.chDenseLaw );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* one of the branches above must have reset this variable,
	 * and it was init to 0 at start.  Confirm that non-zero */
	ASSERT( pressure.PresTotlCorrect > FLT_MIN );

	/* now see whether current pressure is within error bounds */
	if( pressure_change_factor > 1. + conv.PressureErrorAllowed || pressure_change_factor < 1. - conv.PressureErrorAllowed )
	{
		lgRet = false;
		conv.lgConvPres = false;
	}
	else
	{
		lgRet = true;
		conv.lgConvPres = true;
	}

	return lgRet;
}
