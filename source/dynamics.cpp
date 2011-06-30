/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* DynaIterEnd called at end of iteration when advection is turned on */
/* DynaStartZone called at start of zone calculation when advection is turned on */
/* DynaEndZone called at end of zone calculation when advection is turned on */
/* DynaIonize, called from ionize to evaluate advective terms for current conditions */
/* DynaPresChngFactor, called from PressureChange to evaluate new density needed for
 * current conditions and wind solution, returns ratio of new to old density */
/* timestep_next - find next time step in time dependent case */
/* DynaZero zero some dynamics variables, called from zero.c */
/* DynaCreateArrays allocate some space needed to save the dynamics structure variables, 
 * called from DynaCreateArrays */
/* DynaPrtZone - called to print zone results */
/* DynaSave save info related to advection */
/* DynaSave, save output for dynamics solutions */
/* ParseDynaTime parse the time command, called from ParseCommands */
/* ParseDynaWind parse the wind command, called from ParseCommands */
#include "cddefines.h"
#include "cddrive.h"
#include "struc.h"
#include "input.h"
#include "colden.h"
#include "radius.h"
#include "thirdparty.h"
#include "stopcalc.h"
#include "hextra.h"
#include "rfield.h"
#include "iterations.h"
#include "trace.h"
#include "conv.h"
#include "timesc.h"
#include "dense.h"
#include "mole.h"
#include "thermal.h"
#include "pressure.h"
#include "phycon.h"
#include "wind.h"
#include "hmi.h"
#include "iso.h"
#include "dynamics.h"
#include "cosmology.h"
#include "taulines.h"
#include "parser.h"
static int ipUpstream=-1,iphUpstream=-1,ipyUpstream=-1;

/* 
 * >>chng 01 mar 16, incorporate advection within dynamical solutions
 * this file contains the routines that incorporate effects of dynamics and advection upon the
 * thermal and ionization solutions.  
 *
 * This code was originally developed in March 2001 during
 * a visit to UNAM Morelia, in collaboration with Robin Williams, Jane Arthur, and Will Henney.
 * Development was continued by email, and in a meeting in July/Aug 2001 at U Cardiff
 * Further development June 2002, UNAM Morelia (Cloudy and the Duendes Verdes)
 */

/* the interpolated upstream densities of all ionization stages,
 * [element][ion] */
static double **UpstreamIon;
static double ***UpstreamStatesElemNEW;
/* total abundance of each element per hydrogen */
static double *UpstreamElem;

/* hydrogen molecules */
static double *Upstream_H2_molec;

/* CO molecules */
static double *Upstream_CO_molec;

/* space to save continuum when time command is used 
static realnum *dyna_flux_save;*/

/* array of times and continuum values we will interpolate upon */
static double *time_elapsed_time , 
	*time_flux_ratio ,
	*time_dt,
	*time_dt_scale_factor;
bool lgtime_dt_specified;
int *lgtime_Recom;
#define	NTIME	200

/* number of time steps actually read in */
static long int nTime_flux=0;

/* routine called at end of iteration to determine new step sizes */
STATIC void DynaNewStep(void);

/* routine called at end of iteration to save values in previous iteration */
STATIC void DynaSaveLast(void);

/* routine called to determine mass flux at given distance */
/* static realnum DynaFlux(double depth); */

/* lookahead distance, as derived in DynaIterEnd */
static double Dyn_dr;

/* advected work */
static double AdvecSpecificEnthalpy;

/* HI ionization structure from previous iteration */
static realnum *Old_histr/*[NZLIM]*/ ,
	/* Lyman continuum optical depth from previous iteration */
	*Old_xLyman_depth/*[NZLIM]*/,
	/* depth of position in previous iteration */
	*Old_depth/*[NZLIM]*/,
	/* old n_p density from previous iteration */
	*Old_hiistr/*[NZLIM]*/,
	/* old pressure from previous iteration */
	*Old_pressure/*[NZLIM]*/,
	/* H density - particles per unit vol */
	*Old_hden/*[NZLIM]*/ ,
	/* density - total grams per unit vol */
	*Old_DenMass/*[NZLIM]*/ ,
	/* sum of enthalpy and kinetic energy per gram */
	*EnthalpyDensity/*[NZLIM]*/,
	/* old electron density from previous iteration */
	*Old_ednstr/*[NZLIM]*/,
	/* sum of enthalpy and kinetic energy per gram */
	*Old_EnthalpyDensity/*[NZLIM]*/;

static realnum **Old_H2_molec;
static realnum **Old_CO_molec;

/* the ionization fractions from the previous iteration */
static realnum ***Old_xIonDense;

/* the gas phase abundances from the previous iteration */
static realnum **Old_gas_phase;

/* the iso levels from the previous iteration */
static realnum ****Old_StatesElemNEW;

/* the number of zones that were saved in the previous iteration */
static long int nOld_zone;

/* the number of zones that were saved in the previous iteration */
static realnum DivergePresInteg;

/*timestep_next - find next time step in time dependent case */
STATIC double timestep_next( void )
{
	static double te_old=-1;
	double timestep_Hp_temp , timestep_return;

	DEBUG_ENTRY( "timestep_next()" );

	timestep_return = dynamics.timestep;

	if( dynamics.lgRecom )
	{
		double te_new;
		if( cdTemp(
			/* four char string, null terminzed, giving the element name */
			"hydr", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number */
			2 ,
			/* will be temperature */
			&te_new, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"VOLUME" ) )
			TotalInsanity();

		if( te_old>0 )
		{
			double dTdStep = fabs(te_new-te_old)/te_new;
			/* this is the change factor to get 0.1 frac change in mean temp */
			double dT_factor = 0.04/SDIV(dTdStep);
			dT_factor = MIN2( 2.0 , dT_factor );
			dT_factor = MAX2( 0.01 , dT_factor );
			timestep_Hp_temp = dynamics.timestep * dT_factor;
		}
		else
		{
			timestep_Hp_temp = -1.;
		}
		te_old = te_new;
		if( timestep_Hp_temp > 0. )
			timestep_return = timestep_Hp_temp;
	}
	else
	{
		timestep_return = dynamics.timestep_init;
	}
	fprintf(ioQQQ,"DEBUG timestep_next returns %.3e, old temp %.2e\n" , timestep_return, te_old );

	return( timestep_return );
}

/* ============================================================================== */
/* DynaPresChngFactor, called from PressureChange to evaluate factor needed
 * to find new density needed for 
 * current conditions and wind solution, returns ratio of new to old density,
 * called when wind velocity is negative */

/* object is to get the local ram pressure
 * RamNeeded = pressure.PresTotlInit + pressure.PresGasCurr + pressure.PresInteg;
 * to equal the sum of the inital pressur at the illuminated face, the local gas pressure,
 * and the integrate radiative acceleration from the incident continuum 
 *
 * the local gas pressure is linear in density if the temperature is constant,
 *
 * the local ram pressure is inversely linear in density because of the relationship
 * between velocity and density introduced by the mass flux conservation 
 */
#define SUBSONIC   1
#define SUPERSONIC 2
/*#define FREE       3*/
#define STRONGD    4
#define ORIGINAL   5
#define SHOCK      6
#define ANTISHOCK  7
#define ANTISHOCK2 8

double DynaPresChngFactor(void)
{
	double factor,
		er,
		width;
	static double dp = -1.,
		dpp = -1.,
		erp = -1.,
		erpp = -1.;
		/** \todo	2	this should be setable at run time as per Robin Williams email
		rshock = 4e16; */
	static int lastzone = -1,
		zonepresmode,
		globalpresmode;
	int ipPRE;

	DEBUG_ENTRY( "DynaPresChngFactor()" );

	/* update the current pressure and radiative acceleration */
	PresTotCurrent();

	/* update the current desired pressure */
	pressure.PresTotlCorrect = pressure.PresTotlInit + pressure.PresInteg*pressure.lgContRadPresOn
		+ DivergePresInteg;
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" DynaPresChngFactor set PresTotlCorrect=%.3e PresTotlInit=%.3e PresInteg=%.3e DivergePresInteg=%.3e\n", 
			pressure.PresTotlCorrect , pressure.PresTotlInit , pressure.PresInteg*pressure.lgContRadPresOn,
			DivergePresInteg );
	}

	if( dynamics.lgTracePrint)
		fprintf(ioQQQ,"Pressure: init %g +rad %g +diverge %g = %g cf %g\n",
			pressure.PresTotlInit, pressure.PresInteg*pressure.lgContRadPresOn,
			DivergePresInteg, pressure.PresTotlCorrect, pressure.PresTotlCurr);
	/*fprintf(ioQQQ,"DEBUG\t%.2f\thden\t%.4e\tPtot\t%.4e\tPgas\t%.4e\n",
		fnzone, dense.gas_phase[ipHYDROGEN],pressure.PresTotlCurr,pressure.PresGasCurr );*/

	/* dynamics.lgSetPresMode is flag to indicate sane value of dynamics.chPresMode.  
	 * If set true with SET DYNAMICS PRESSURE MODE
	 * then do not override that choice */
	if( !dynamics.lgSetPresMode )
	{
		/* above set true if pressure mode was set with
		 * SET DYNAMICS PRESSURE MODE - if we got here
		 * it was not set, and must make a guess */
		if(pressure.PresGasCurr < pressure.PresRamCurr)
			strcpy( dynamics.chPresMode , "supersonic" );
		else
			strcpy( dynamics.chPresMode , "subsonic" );
		/* clear the flag - pressure mode has been set */
		dynamics.lgSetPresMode = true;
	}

	/* if globally looking for transonic solution, then locally sometimes
	 * supersonic, sometimes subsonic - this branch sets global flag,
	 * which can also be set with SET DYNAMICS PRESSURE MODE.
	 * Under default conditions, ChPresMode was set in previous branch
	 * to sub or super sonic depending on current velocity on first time*/
	if( strcmp( dynamics.chPresMode , "original" ) == 0 )
	{
		globalpresmode = ORIGINAL;
		pressure.lgSonicPointAbortOK = true;
	}
	else if( strcmp( dynamics.chPresMode , "subsonic" ) == 0 )
	{
		globalpresmode = SUBSONIC;
		pressure.lgSonicPointAbortOK = true;
	}
	/* supersonic */
	else if( strcmp( dynamics.chPresMode , "supersonic" ) == 0 )
	{
		globalpresmode = SUPERSONIC;
		pressure.lgSonicPointAbortOK = true;
	}
	/* strong d */
	else if( strcmp( dynamics.chPresMode , "strongd" ) == 0 )
	{
		globalpresmode = STRONGD;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "shock" ) == 0 )
	{
		globalpresmode = SHOCK;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "antishock-by-mach" ) == 0 )
	{
		globalpresmode = ANTISHOCK2;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "antishock" ) == 0 )
	{
		globalpresmode = ANTISHOCK;
		pressure.lgSonicPointAbortOK = false;
	}

	/* this sets pressure mode for the current zone based on global mode
	 * and local conditions */
	if(globalpresmode == ORIGINAL)
	{
		/* in this mode use comparison between ram and gas pressure to determine
		 * whether sub or super sonic */
		if(pressure.PresGasCurr > pressure.PresRamCurr)
		{
			zonepresmode = SUBSONIC;
		}
		else
		{
			zonepresmode = SUPERSONIC;
		}
	}
	else if(globalpresmode == STRONGD)
	{
		if(nzone <= 1)
			zonepresmode = SUPERSONIC;
	}
	else if(globalpresmode == SUBSONIC)
	{
		zonepresmode = SUBSONIC;
	}
	else if(globalpresmode == SUPERSONIC)
	{
		zonepresmode = SUPERSONIC;
	}
	else if(globalpresmode == SHOCK)
	{
		if(radius.depth < dynamics.ShockDepth)
		{
			zonepresmode = SUBSONIC;
		}
		else
		{
			zonepresmode = SUPERSONIC;
		}
	}
	else if(globalpresmode == ANTISHOCK)
	{
		if(radius.depth < dynamics.ShockDepth)
		{
			zonepresmode = SUPERSONIC;
		}
		else
		{
			zonepresmode = SUBSONIC;
		}
	}
	else if(globalpresmode == ANTISHOCK2)
	{
		/* WJH 19 May 2004: This version of the antishock mode will
		  insert the antishock at the point where the isothermal Mach
		  number falls below a certain value, dynamics.ShockMach */
		if( fabs(wind.windv) > dynamics.ShockMach*sqrt(pressure.PresGasCurr/dense.xMassDensity))
		{
			zonepresmode = SUPERSONIC;
		}
		else
		{
			zonepresmode = SUBSONIC;
		}
	}
	else
	{
		printf("Need to set global pressure mode\n");
		cdEXIT( EXIT_FAILURE );
	}

	er = pressure.PresTotlCurr-pressure.PresTotlCorrect;
	/* fprintf(ioQQQ,"Ds %ld: %ld; %g %g %g %g %g %g %g %g\n",iteration,nzone,dense.gas_phase[ipHYDROGEN],er,pressure.PresTotlCorrect,pressure.PresTotlCurr,phycon.te,thermal.ctot,thermal.htot,phycon.EnthalpyDensity); */
	/* fprintf(ioQQQ,"Ds %ld: %ld; %g %g %g %g\n",iteration,nzone,dense.gas_phase[ipHYDROGEN],er,pressure.PresTotlCurr,phycon.te); */
	if(globalpresmode == ORIGINAL || lastzone != nzone || fabs(er-erp) < SMALLFLOAT) 
	{
		/* First time in this zone, or when last step made no change, 
		 * take step hopefully in the right direction...
		 * ...or at least somewhere */
		if( zonepresmode == SUBSONIC ) 
		{
			/* when gas pressure dominated need to increase density to increase pressure */
			factor = pressure.PresTotlCorrect / pressure.PresTotlCurr;
			ipPRE = 0;
		}
		else
		{
			/* when ram pressure dominated need to decrease density to increase pressure */
			factor = pressure.PresTotlCurr / pressure.PresTotlCorrect;
			ipPRE = 1;
		}
		if(fabs(factor-1.) > 0.01)
		{
			factor = 1.+sign(0.01,factor-1.);
		}
		erp = er;
		dp = dense.gas_phase[ipHYDROGEN];
		erpp = -1.;
		dpp = -1.;
	}
	else
	{
#if 0
		printf("Ds: %d; %g %g %g; %g %g %g tot %g\n",nzone,dense.gas_phase[ipHYDROGEN],dp,dpp,er,erp,erpp,
					 pressure.PresTotlCorrect);
#endif
		if(1 || dpp < 0. || fabs((dense.gas_phase[ipHYDROGEN]-dp)*(dp-dpp)*(dpp-dense.gas_phase[ipHYDROGEN])) < SMALLFLOAT) 
		{
			/* second iteration on this zone, do linerar fit to previous Pres - rho curve
			 * Linear approximation to guess root with two independent samples */
			factor = (dense.gas_phase[ipHYDROGEN]*erp-dp*er)/(erp-er)/dense.gas_phase[ipHYDROGEN];
			/* Limit step length to `reasonable' extrapolation */
			width = fabs(1.-dp/dense.gas_phase[ipHYDROGEN]);
			if(width > 1e-2)
				width = 1e-2;

			/* Subsonic case: pressure ^ with density ^ => increase density further */
			/* Super "  case: pressure ^ with density v => decrease density further */

			/* printf("Presmode %d flag %g factor %g\n",zonepresmode,(er-erp)*(dense.gas_phase[ipHYDROGEN]-dp),factor); */
			/* WJH 21 May 2004: I think that this part is to force the solution
			 * towards the required branch when it appears to have the "wrong" value
			 * of dP/drho */
			if(zonepresmode == SUBSONIC && (er-erp)*(dense.gas_phase[ipHYDROGEN]-dp) < 0)
			{				
				factor = 1+3*width;
			}
			else if(zonepresmode == SUPERSONIC && (er-erp)*(dense.gas_phase[ipHYDROGEN]-dp) > 0)
			{
				factor = 1-3*width;
			}

			if(fabs(factor-1.) > 3*width)
			{
				factor = 1.+sign(3*width,factor-1.);
			}
			ipPRE = 2;
			if(fabs(dp-dense.gas_phase[ipHYDROGEN]) > SMALLFLOAT) 
			{
				dpp = dp;
				erpp = erp;
				dp = dense.gas_phase[ipHYDROGEN];
				erp = er;
			}
		}
		else
		{
			/* 3rd or more solutions in this zone so have extensive data
			 * on pressure - density relation 
			 * Use quadratic fit to last three errors to estimate optimum */
			double a, b, c, q, dmin, emin, dsol, f1 , f2;
			int iroot;
			a = er/(dense.gas_phase[ipHYDROGEN]-dp)/(dense.gas_phase[ipHYDROGEN]-dpp) +
				erp/(dp-dense.gas_phase[ipHYDROGEN])/(dp-dpp)+
				erpp/(dpp-dense.gas_phase[ipHYDROGEN])/(dpp-dp);
			b = (erp-erpp)/(dp-dpp) - a * (dp+dpp);
			c = erp - dp*(a*dp+b);
			dmin = (-0.5*b/a);
			emin = (a*dmin+b)*dmin+c;
			if(a < 0) 
			{
				a = -a;
				b = -b;
				c = -c;
			}
#if 0
			fprintf(ioQQQ,"Check 1: %g %g\n",er,(a*dense.gas_phase[ipHYDROGEN]+b)*dense.gas_phase[ipHYDROGEN]+c);
			fprintf(ioQQQ,"Check 2: %g %g\n",erp,(a*dp+b)*dp+c);
			fprintf(ioQQQ,"Check 3: %g %g\n",erpp,(a*dpp+b)*dpp+c);
#endif
			q = b*b-4*a*c;
			if(q < 0) 
			{
				/* Imaginary root, search for local minimum */
				/* printf("No root at %d (%g cf %g) => %g\n",nzone,q,b*b,dmin); */
				factor = dmin/dense.gas_phase[ipHYDROGEN];

				/** \todo	2	Switch to supersonic when bad enough failure for STRONGD --
				 * need to improve this logic when the p(rho) data is cleaner
				 * */
				if(globalpresmode == STRONGD && -q > 1e-3*b*b)
				{
					zonepresmode = SUBSONIC;
				}
			} 
			else
			{
				/* Look for nearest root */
				/* if(zonepresmode == SUPERSONIC || (zonepresmode != SUBSONIC && (dense.gas_phase[ipHYDROGEN]-dmin) < 0)) */

				if(zonepresmode == SUPERSONIC)
					iroot = 1;
				else
					iroot = 0;

				if(emin > 0)
					iroot = ! iroot;

				if(iroot)
				{
					/* Look for smaller root */
					if(b > 0) 
					{
						dsol = -(b+sqrt(q))/(2*a);
					} 
					else 
					{
						dsol = 2*c/(-b+sqrt(q));
					}
				}
				else
				{
					if(b < 0)
					{
						dsol = (-b+sqrt(q))/(2*a);
					}
					else
					{
						dsol = -2*c/(b+sqrt(q));
					}
				}
				factor = dsol/dense.gas_phase[ipHYDROGEN];
#if 0
				fprintf(ioQQQ,"Check 4: %g %g %d %g %g\n",dsol,(a*dsol+b)*dsol+c,iroot,dmin,emin);
#endif
			}
			/* Limit step length */
			f1 = fabs(1.-dpp/dense.gas_phase[ipHYDROGEN]);
			f2 = fabs(1.- dp/dense.gas_phase[ipHYDROGEN]);
			/*width = MAX2(fabs(1.-dpp/dense.gas_phase[ipHYDROGEN]),fabs(1.-dp/dense.gas_phase[ipHYDROGEN]));*/
			width = MAX2(f1,f2);
			/* width = MAX2(width,1e-2); */
			if(fabs(factor-1.) > 3*width)
			{
				factor = 1.+sign(3*width,factor-1.);
			}
			ipPRE = 3;
			if(fabs(dp-dense.gas_phase[ipHYDROGEN]) > SMALLFLOAT) 
			{
				dpp = dp;
				erpp = erp;
				dp = dense.gas_phase[ipHYDROGEN];
				erp = er;
			}
		}
	}		

#if 0
	printf("Out: %d; %g; %g %g; %g %g\n",nzone,factor*dense.gas_phase[ipHYDROGEN],dp,dpp,erp,erpp);
#endif
	lastzone = nzone;

	if( dynamics.lgTracePrint )
		fprintf(ioQQQ,"windv %li r %g v %g f %g\n",
			nzone,radius.depth,wind.windv,DynaFlux(radius.depth));

	/* convergence trace at this level */
	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			" DynaPresChngFactor mode %s scale factor %.3f vel %.3e\n",
			dynamics.chPresMode , factor , wind.windv );
	}

	{
		/*@-redef@*/
		enum{DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			char chPRE[][4] = {"gas" , "ram", "sec", "par" };
			fprintf(ioQQQ,
				"pre %s\tfac\t%.5f\n",
				chPRE[ipPRE],
				factor
				);
		}
	}

	/* identically zero velocities cannot occur */
	ASSERT( wind.windv != 0. || (wind.windv == 0. && dynamics.lgTimeDependentStatic) );

	return factor;
}

/* ============================================================================== */
/* DynaIonize, called from ConvBase to evaluate advective terms for current conditions,
 * calculates terms to add to ionization balance equation */
void DynaIonize(void)
{
	long int nelem, 
		ion, 
		mol ,
		i;

	DEBUG_ENTRY( "DynaIonize()" );

	/* the time (s) needed for gas to move dynamics.Dyn_dr  */
	/* >>chng 02 dec 11 rjrw increase dynamics.dynamics.timestep when beyond end of previous zone, to allow -> eqm */
	if( !dynamics.lgTimeDependentStatic )
	{
		/* in time dependent model dynamics.timestep only changes once at end of iteration
		 * and is constant across a model */
		/* usual case - full dynamics */
		dynamics.timestep = -Dyn_dr/wind.windv;
	}
	/* printf("%d %g %g %g %g\n",nzone,radius.depth,Dyn_dr,radius.depth-Old_depth[nOld_zone-1],dynamics.timestep); */

	ASSERT(nzone<struc.nzlim );
	if( nzone > 0 )
		EnthalpyDensity[nzone-1] = (realnum)phycon.EnthalpyDensity;

	/* do nothing on first iteration or when looking beyond previous iteration */
	/* >>chng 05 jan 27, from hardwired "iteration < 2" to more general case,
	 * this is set with SET DYNAMICS RELAX command with the default of 2 */
	//For advection cases, switch to local equilibrium when depth is outside
	//the region of previous iteration.
	//Possibly should limit range of dynamical sources further by adding Dyn_dr???
	double depth = radius.depth; 
	if( iteration < dynamics.n_initial_relax+1 || 
			( ! dynamics.lgTimeDependentStatic && 
				( depth < 0 || depth > dynamics.oldFullDepth ) ) )
	{
		/* first iteration, return zero */
		dynamics.Cool_r = 0.;
		dynamics.Heat_v = 0.;
		dynamics.dHeatdT = 0.;

		dynamics.Rate = 0.;

		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			for( ion=0; ion<nelem+2; ++ion )
			{
				dynamics.Source[nelem][ion] = 0.;
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_local[ipISO][nelem]; ++level )
					{
						dynamics.StatesElemNEW[nelem][nelem-ipISO][level] = 0.;
					}
				}
			}
		}
		for(mol=0;mol<N_H_MOLEC;mol++)
		{
			dynamics.H2_molec[mol] = 0.;
		}
		for(mol=0;mol<mole.num_comole_calc;mol++)
		{
			dynamics.CO_molec[mol] = 0.;
		}
		return;
	}

	if( dynamics.lgTracePrint )
	{
		fprintf(ioQQQ,"workwork\t%li\t%.3e\t%.3e\t%.3e\n",
			nzone,
			phycon.EnthalpyDensity,
			0.5*POW2(wind.windv)*dense.xMassDensity ,
			5./2.*pressure.PresGasCurr
			); /**/
	}

	/* net cooling due to advection */
	/* >>chng 01 aug 01, removed hden from dilution, put back into here */
	/* >>chng 01 sep 15, do not update this variable the very first time this routine
	 * is called at the new zone. */
	dynamics.Cool_r = 1./dynamics.timestep*dynamics.lgCoolHeat;
	dynamics.Heat_v = AdvecSpecificEnthalpy/dynamics.timestep*dynamics.lgCoolHeat;
	dynamics.dHeatdT = 0.*dynamics.lgCoolHeat;

#if 0
		/* printf("DynaCool %g %g %g\n",	dynamics.Cool,	phycon.EnthalpyDensity/dynamics.timestep,AdvecSpecificEnthalpy*dense.gas_phase[ipHYDROGEN]/dynamics.timestep); */
		if(dynamics.Cool > 0) {
			dynamics.Heat = 0.;
			/* temp deriv of cooling minus heating */
			dynamics.dCooldT = 5./2.*pressure.PresGasCurr/phycon.te/dynamics.timestep*dynamics.lgCoolHeat;
			dynamics.dHeatdT = 0.*dynamics.lgCoolHeat;
		} else {
			dynamics.Heat = -dynamics.Cool;
			dynamics.Cool = 0.;
			/* temp deriv of cooling minus heating */
			dynamics.dCooldT = 0.*dynamics.lgCoolHeat;
			dynamics.dHeatdT = -5./2.*pressure.PresGasCurr/phycon.te/dynamics.timestep*dynamics.lgCoolHeat;
		}
#endif

#	if 0
	if( dynamics.lgTracePrint || nzone>17 && iteration == 10)
	{
		fprintf(ioQQQ,
			"dynamics cool-heat\t%li\t%.3e\t%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			nzone,
			phycon.te, 
			dynamics.lgCoolHeat,
			thermal.htot,
			phycon.EnthalpyDensity/dynamics.timestep,
			AdvecSpecificEnthalpy*dense.gas_phase[ipHYDROGEN]/dynamics.timestep,
			phycon.EnthalpyDensity,
			AdvecSpecificEnthalpy*dense.gas_phase[ipHYDROGEN],
			dense.gas_phase[ipHYDROGEN],
			dynamics.timestep);
	}
#	endif

	/* second or greater iteration, have advective terms */
	/* this will evaluate advective terms for current physical conditions */

	/* the rate (s^-1) atoms drift in from upstream, a source term for the ground */

	/* dynamics.Hatom/dynamics.timestep is the source (cm^-3 s^-1) of neutrals,
		 normalized to (s^-1) by the next higher ionization state as for all 
		 other recombination terms */

	/* dynamics.xIonDense[ipHYDROGEN][1]/dynamics.timestep is the sink (cm^-3 s^-1) of
		 ions, normalized to (s^-1) by the ionization state as for all other
		 ionization terms */

	dynamics.Rate = 1./dynamics.timestep;

	for(mol=0;mol<N_H_MOLEC;mol++)
	{
		dynamics.H2_molec[mol] = Upstream_H2_molec[mol]* dense.gas_phase[ipHYDROGEN]*dynamics.Rate;
	}

	for(mol=0;mol<mole.num_comole_calc;mol++)
	{
		dynamics.CO_molec[mol] = Upstream_CO_molec[mol]* dense.gas_phase[ipHYDROGEN]*dynamics.Rate;
	}

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* check that the total number of each element is conserved along the flow */
			if(fabs(UpstreamElem[nelem]*dense.gas_phase[ipHYDROGEN] -dense.gas_phase[nelem])/dense.gas_phase[nelem]>=1e-3) 
			{
				/* sum of all H in standard H chem */
				realnum sumh = 0.;
				for(i=0;i<N_H_MOLEC;i++) 
				{
					sumh += hmi.Hmolec[i]*hmi.nProton[i];
				}
				fprintf(ioQQQ,
					"PROBLEM conservation error: zn %li elem %li upstream %.8e abund %.8e (up-ab)/up %.2e f1(H n CO) %.2e f2(H n CO) %.2e\n",
					nzone ,
					nelem,
					UpstreamElem[nelem]*dense.gas_phase[ipHYDROGEN],
					dense.gas_phase[nelem] ,
					(UpstreamElem[nelem]*dense.gas_phase[ipHYDROGEN]-dense.gas_phase[nelem]) /
						(UpstreamElem[nelem]*dense.gas_phase[ipHYDROGEN]) ,
					(dense.gas_phase[ipHYDROGEN]-sumh) / dense.gas_phase[ipHYDROGEN] ,
					dense.H_sum_in_CO / dense.gas_phase[ipHYDROGEN] );
			}
			/* ASSERT( fabs(UpstreamElem[nelem]*dense.gas_phase[ipHYDROGEN] -dense.gas_phase[nelem])/dense.gas_phase[nelem]<1e-3); */
			for( ion=0; ion<dense.IonLow[nelem]; ++ion )
			{
				dynamics.Source[nelem][ion] = 0.;
			}
			for( ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
			{
				/* Normalize to next higher state in current zone, except at first iteration
				where upstream version may be a better estimate (required for
				convergence in the small dynamics.timestep limit) */

				dynamics.Source[nelem][ion] = 
					/* UpstreamIon is ion number per unit hydrogen because dilution is 1/hden 
					 * this forms the ratio of upstream atom over current ion, per dynamics.timestep,
					 * so Source has units cm-3 s-1 */
					UpstreamIon[nelem][ion] * dense.gas_phase[ipHYDROGEN] / dynamics.timestep;

			}
			for( ion=dense.IonHigh[nelem]+1;ion<nelem+2; ++ion )
			{
				dynamics.Source[nelem][ion] = 0.;
			}
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long level=0; level < iso.numLevels_local[ipISO][nelem]; ++level )
				{
					dynamics.StatesElemNEW[nelem][nelem-ipISO][level] = 
						UpstreamStatesElemNEW[nelem][nelem-ipISO][level] * dense.gas_phase[ipHYDROGEN] / dynamics.timestep;
				}
			}
		}
	}

#	if 0
	fprintf(ioQQQ,"dynamiccc\t%li\t%.2e\t%.2e\t%.2e\t%.2e\n",
		nzone,
		dynamics.Rate,
		dynamics.Source[ipHYDROGEN][0],
		dynamics.Rate,
		dynamics.Source[ipCARBON][3]);
#	endif
#if 0
	nelem = ipCARBON;
	ion = 3;
	/*if( nzone > 160 && iteration > 1 )*/
	fprintf(ioQQQ,"dynaionizeeee\t%li\t%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
		nzone, 
		ipUpstream,
		radius.depth ,
		Old_depth[ipUpstream],
		dense.xIonDense[nelem][ion], 
		UpstreamIon[nelem][ion]* dense.gas_phase[ipHYDROGEN],
		Old_xIonDense[ipUpstream][nelem][ion] ,
		dense.xIonDense[nelem][ion+1], 
		UpstreamIon[nelem][ion+1]* dense.gas_phase[ipHYDROGEN],
		Old_xIonDense[ipUpstream][nelem][ion+1] ,
		dynamics.timestep,
		dynamics.Source[nelem][ion]
		);
#endif
	if( dynamics.lgTracePrint )
	{
		fprintf(ioQQQ,"    DynaIonize, %4li photo=%.2e , H recom= %.2e \n",
			nzone,dynamics.Rate , dynamics.Source[0][0]  );
	}
	return;
}

/* ============================================================================== */
/* DynaStartZone called at start of zone calculation when advection is turned on */
void DynaStartZone(void)
{
	/* this routine is called at the start of a zone calculation, by ZoneStart:
	 *
	 * it sets deduced variables to zero if first iteration,
	 *
	 * if upstream depth is is outside the computed structure on previous iteration,
	 * return value at shielded face 
	 *
	 * Also calculates discretization_error, an estimate of the accuracy of the source terms.
	 *
	 * */

	/* this is index of upstream cell in structure stored from previous iteration */
	double upstream, dilution, dilutionleft, dilutionright, frac_next;

	/* Properties for cell half as far upstream, used to converge dynamics.timestep */
	double hupstream, hnextfrac=-BIGFLOAT, hion, hmol, hiso;

	/* Properties for cell at present depth, used to converge dynamics.timestep */
	double ynextfrac=-BIGFLOAT, yion, ymol, yiso;

	long int nelem , ion, mol;

	DEBUG_ENTRY( "DynaStartZone()" );

	/* do nothing on first iteration */
	if( iteration < 2 )
	{
		dynamics.Upstream_hden = 0.;
		AdvecSpecificEnthalpy = 0.;
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			for( ion=0; ion<nelem+2; ++ion )
			{
				UpstreamIon[nelem][ion] = 0.;
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_max[ipISO][nelem]; ++level )
					{
						UpstreamStatesElemNEW[nelem][nelem-ipISO][level] = 0.;
					}
				}
			}
		}
		/* hydrogen molecules */
		for(mol=0; mol<N_H_MOLEC;mol++) 
		{
			Upstream_H2_molec[mol] = 0;
		}
		for(mol=0; mol<mole.num_comole_calc;mol++) 
		{
			Upstream_CO_molec[mol] = 0;
		}

		ipUpstream = -1;
		iphUpstream = -1;
		ipyUpstream = -1;
		return;
	}

	/* radius.depth is distance from illuminated face of cloud to outer edge of
	 * current zone, which has thickness radius.drad */

	/* find where the down stream point is, in previous iteration,
	 * NB neg sign since velocity is negative, we are looking for point
	 * where current gas cell was, in previous iteration */

	/* true, both advection and wind solution */
	/* don't interpolate to the illuminated side of the first cell */
	upstream = MAX2(Old_depth[0] , radius.depth + Dyn_dr);
	hupstream = 0.5*(radius.depth + upstream);

	/* now locate upstream point in previous stored structure,
	 * will be at the same point or higher than we found previously */
	while( (Old_depth[ipUpstream+1] < upstream ) && 
			 ( ipUpstream < nOld_zone-1 ) )
	{
		++ipUpstream;
	}
	ASSERT( ipUpstream <= nOld_zone-1 );

	/* above loop will give ipUpstream == nOld_zone-1 if computed structure has been overrun */

	if( (ipUpstream != -1) && (ipUpstream != nOld_zone-1)&& ((Old_depth[ipUpstream+1] - Old_depth[ipUpstream])> SMALLFLOAT) )
	{
		/* we have not overrun radius scale of previous iteration */
		frac_next = ( upstream - Old_depth[ipUpstream])/
			(Old_depth[ipUpstream+1] - Old_depth[ipUpstream]);
		dynamics.Upstream_hden = (realnum)(Old_hden[ipUpstream] +
			(Old_hden[ipUpstream+1] - Old_hden[ipUpstream])*
			frac_next);
		dilutionleft = 1./Old_hden[ipUpstream];
		dilutionright = 1./Old_hden[ipUpstream+1];

		/* fractional changes in density from passive advection */
		/* >>chng 01 May 02 rjrw: use hden for dilution */
		/* >>chng 01 aug 01, remove hden here, put back into vars when used in DynaIonize */
		dilution = 1./dynamics.Upstream_hden;

		/* the advected enthalphy */
		AdvecSpecificEnthalpy = (Old_EnthalpyDensity[ipUpstream]*dilutionleft +
			(Old_EnthalpyDensity[ipUpstream+1]*dilutionright - Old_EnthalpyDensity[ipUpstream]*dilutionleft)*
			frac_next);

		ASSERT( Old_depth[ipUpstream] <= upstream && upstream <= Old_depth[ipUpstream+1] );

		// we have a mix of realnum and double here, so make sure realnum is used for the test...
		realnum lo = (realnum)(Old_EnthalpyDensity[ipUpstream]*dilutionleft);
		realnum x = (realnum)AdvecSpecificEnthalpy;
		realnum hi = (realnum)(Old_EnthalpyDensity[ipUpstream+1]*dilutionright);
		if( ! fp_bound(lo,x,hi) )
		{
			fprintf(ioQQQ,"PROBLEM interpolated enthalpy density is not within range  %.16e\t%.16e\t%.16e\t%e\t%e\n",
				lo, x, hi, (hi-x)/(hi-lo), (x-lo)/(hi-lo));
		}

		ASSERT( fp_bound(lo,x,hi) );

		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			UpstreamElem[nelem] = 0.;
			for( ion=0; ion<nelem+2; ++ion )
			{
				/* Robin - I made several changes like the following - it seems easier to
				 * bring out the division by the old hydrogen density rather than putting in
				 * dilution and then looking for how diluation is defined.  I think the code
				 * is equivalent */
				/* UpstreamIon is ion number per unit hydrogen, both at the upstream position */
				UpstreamIon[nelem][ion] = 
					Old_xIonDense[ipUpstream][nelem][ion]/Old_hden[ipUpstream] +
					(Old_xIonDense[ipUpstream+1][nelem][ion]/Old_hden[ipUpstream+1] - 
					Old_xIonDense[ipUpstream][nelem][ion]/Old_hden[ipUpstream])*
					frac_next;

				UpstreamElem[nelem] += UpstreamIon[nelem][ion];
			}
		}

		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_max[ipISO][nelem]; ++level )
					{
						UpstreamStatesElemNEW[nelem][nelem-ipISO][level] = 
							Old_StatesElemNEW[ipUpstream][nelem][nelem-ipISO][level]/Old_hden[ipUpstream] +
							(Old_StatesElemNEW[ipUpstream+1][nelem][nelem-ipISO][level]/Old_hden[ipUpstream+1] - 
							 Old_StatesElemNEW[ipUpstream][nelem][nelem-ipISO][level]/Old_hden[ipUpstream])*
							frac_next;
						/* check for NaN */
						ASSERT( !isnan( UpstreamStatesElemNEW[nelem][nelem-ipISO][level] ) );
					}
				}
			}
		}

		for(mol=0;mol<N_H_MOLEC;mol++)
		{
			Upstream_H2_molec[mol] = Old_H2_molec[ipUpstream][mol]/Old_hden[ipUpstream] +
				(Old_H2_molec[ipUpstream+1][mol]/Old_hden[ipUpstream+1] - 
				 Old_H2_molec[ipUpstream][mol]/Old_hden[ipUpstream])*
				frac_next;
			/* test on mol > 1, first two elements are H0 and H+, which are alread
			 * counted in upstreamElem */
			if(mol > 1)
				UpstreamElem[ipHYDROGEN] += Upstream_H2_molec[mol]*hmi.nProton[mol];
		}
		/* loop is only up to mole.num_heavy_molec, not mole.num_comole_calc, since do not want
		 * to consider atoms and ions that have also been done */
		for(mol=0;mol<mole.num_comole_calc;mol++)
		{
			Upstream_CO_molec[mol] = Old_CO_molec[ipUpstream][mol]/Old_hden[ipUpstream] +
				(Old_CO_molec[ipUpstream+1][mol]/Old_hden[ipUpstream+1] - 
				 Old_CO_molec[ipUpstream][mol]/Old_hden[ipUpstream])*
				frac_next;

			if(COmole[mol]->n_nuclei > 1)
			{
				for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
				{
					if( mole.lgElem_in_chemistry[nelem] )
					{
						UpstreamElem[nelem] +=
							COmole[mol]->nElem[nelem]*Upstream_CO_molec[mol];
					}
				}
			}
		}
	}
	else
	{
		/* SPECIAL CASE - we have overrun the previous iteration's radius */
		long ipBound = ipUpstream;
		if (ipBound == -1)
			ipBound = 0;
		dynamics.Upstream_hden = Old_hden[ipBound];
		/* fractional changes in density from passive advection */
		dilution = 1./dynamics.Upstream_hden;
		/* AdvecSpecificEnthalpy enters as heat term */
		AdvecSpecificEnthalpy = Old_EnthalpyDensity[ipBound]/Old_hden[ipBound];
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			UpstreamElem[nelem] = 0.;
			for( ion=0; ion<nelem+2; ++ion )
			{
				/* UpstreamIon is ion number per unit hydrogen */
				UpstreamIon[nelem][ion] = 
					Old_xIonDense[ipBound][nelem][ion]/Old_hden[ipBound];
				UpstreamElem[nelem] += UpstreamIon[nelem][ion];
			}
		}

		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_max[ipISO][nelem]; ++level )
					{
						UpstreamStatesElemNEW[nelem][nelem-ipISO][level] = 
							Old_StatesElemNEW[ipBound][nelem][nelem-ipISO][level]/Old_hden[ipBound];
						/* check for NaN */
						ASSERT( !isnan( UpstreamStatesElemNEW[nelem][nelem-ipISO][level] ) );
					}
				}
			}
		}

		for(mol=0;mol<N_H_MOLEC;mol++) 
		{
			Upstream_H2_molec[mol] = Old_H2_molec[ipBound][mol]/Old_hden[ipBound];
			if(mol > 1)
				UpstreamElem[ipHYDROGEN] += Upstream_H2_molec[mol]*hmi.nProton[mol];
		}
		for(mol=0;mol<mole.num_comole_calc;mol++) 
		{
			Upstream_CO_molec[mol] = Old_CO_molec[ipBound][mol]/Old_hden[ipBound];
			if(COmole[mol]->n_nuclei > 1)
			{
				for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
				{
					if( mole.lgElem_in_chemistry[nelem] )
					{
						UpstreamElem[nelem] +=
							Upstream_CO_molec[mol]*COmole[mol]->nElem[nelem];
					}
				}
			}
		}
	}

	/* Repeat enough of the above for half-step and no-step to judge convergence:
	 * the result of this code is the increment of discretization_error */

	while( (Old_depth[iphUpstream+1] < hupstream ) && 
		( iphUpstream < nOld_zone-1 ) )
	{
		++iphUpstream;
	}
	ASSERT( iphUpstream <= nOld_zone-1 );

	while( (Old_depth[ipyUpstream+1] < radius.depth ) && 
		( ipyUpstream < nOld_zone-1 ) )
	{
		++ipyUpstream;
	}
	ASSERT( ipyUpstream <= nOld_zone-1 );

	dynamics.dRad = BIGFLOAT;

	if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
		hnextfrac = ( hupstream - Old_depth[iphUpstream])/
			(Old_depth[iphUpstream+1] - Old_depth[iphUpstream]);
	else
		hnextfrac = 0.;

	if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
		ynextfrac = ( radius.depth - Old_depth[ipyUpstream])/
			(Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]);
	else
		ynextfrac = 0.;

	// Shouldn't be jumping over large numbers of upstream cells
	if(ipUpstream != -1 && ipUpstream < nOld_zone-1)
		dynamics.dRad = MIN2(dynamics.dRad,
												 5*(Old_depth[ipUpstream+1] - Old_depth[ipUpstream]));
	

	// Value for scaling zonal changes to set zone width
	const double STEP_FACTOR=0.05;

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		for( ion=0; ion<nelem+2; ++ion )
		{
			double f1;
			if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
			{
				yion = 
					Old_xIonDense[ipyUpstream][nelem][ion]/Old_hden[ipyUpstream] +
					(Old_xIonDense[ipyUpstream+1][nelem][ion]/Old_hden[ipyUpstream+1] - 
					 Old_xIonDense[ipyUpstream][nelem][ion]/Old_hden[ipyUpstream])*
					ynextfrac;
			}
			else
			{		
				long ipBound = ipyUpstream;
				if (-1 == ipBound)
					ipBound = 0;
				yion = Old_xIonDense[ipBound][nelem][ion]/Old_hden[ipBound];
			}
			if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
			{
				hion = 
					Old_xIonDense[iphUpstream][nelem][ion]/Old_hden[iphUpstream] +
					(Old_xIonDense[iphUpstream+1][nelem][ion]/Old_hden[iphUpstream+1] - 
					 Old_xIonDense[iphUpstream][nelem][ion]/Old_hden[iphUpstream])*
					hnextfrac;
			}
			else
			{		
				long ipBound = iphUpstream;
				if (-1 == ipBound)
					ipBound = 0;
				hion = Old_xIonDense[ipBound][nelem][ion]/Old_hden[ipBound];
			}

			/* the proposed thickness of the next zone */
			f1 = fabs(yion - UpstreamIon[nelem][ion] );
			if( f1 > SMALLFLOAT )
			{
				dynamics.dRad = MIN2(dynamics.dRad,STEP_FACTOR*fabs(Dyn_dr) * 
				/* don't pay attention to species with abundance relative to H less than 1e-6 */
					MAX2(yion + UpstreamIon[nelem][ion],1e-6 ) / f1);
			}

			/* Must be consistent with convergence_error below */
			/* this error is error due to the advection length not being zero - a finite
			 * advection length.  no need to bring convergence error to below
			 * discretization error.  when convergece error is lower than a fraction of this
			 * errror we reduce the advection length. */
			dynamics.discretization_error += POW2(yion-2.*hion+UpstreamIon[nelem][ion]); 
			dynamics.error_scale2 += POW2(UpstreamIon[nelem][ion]-yion);
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long level=0; level < iso.numLevels_max[ipISO][nelem]; ++level )
				{
					double f1;
					if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && 
						(Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
					{
						yiso = 
							Old_StatesElemNEW[ipyUpstream][nelem][nelem-ipISO][level]/Old_hden[ipyUpstream] +
							(Old_StatesElemNEW[ipyUpstream+1][nelem][nelem-ipISO][level]/Old_hden[ipyUpstream+1] - 
							 Old_StatesElemNEW[ipyUpstream][nelem][nelem-ipISO][level]/Old_hden[ipyUpstream])*
							ynextfrac;
					}
					else
					{		
						long ipBound = ipyUpstream;
						if (-1 == ipBound)
							ipBound = 0;
						yiso = Old_StatesElemNEW[ipBound][nelem][nelem-ipISO][level]/Old_hden[ipBound];
					}
					if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && 
						(Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
					{
						hiso = 
							Old_StatesElemNEW[iphUpstream][nelem][nelem-ipISO][level]/Old_hden[iphUpstream] +
							(Old_StatesElemNEW[iphUpstream+1][nelem][nelem-ipISO][level]/Old_hden[iphUpstream+1] - 
							 Old_StatesElemNEW[iphUpstream][nelem][nelem-ipISO][level]/Old_hden[iphUpstream])*
							hnextfrac;
					}
					else
					{		
						long ipBound = iphUpstream;
						if (-1 == ipBound)
							ipBound = 0;
						hiso = Old_StatesElemNEW[ipBound][nelem][nelem-ipISO][level]/Old_hden[ipBound];
					}

					/* the proposed thickness of the next zone */
					f1 = fabs(yiso - UpstreamStatesElemNEW[nelem][nelem-ipISO][level] );
					if( f1 > SMALLFLOAT )
					{
						dynamics.dRad = MIN2(dynamics.dRad,fabs(Dyn_dr*STEP_FACTOR) * 
						/* don't pay attention to species with abundance relative to H less tghan 1e-6 */
						  MAX2(yiso + UpstreamStatesElemNEW[nelem][nelem-ipISO][level],1e-6 ) / f1);
					}

					/* Must be consistent with convergence_error below */
					/* this error is error due to the advection length not being zero - a finite
					 * advection length.  no need to bring convergence error to below
					 * discretization error.  when convergece error is lower than a fraction of this
					 * errror we reduce the advection length. */
					dynamics.discretization_error += POW2(yiso-2.*hiso+UpstreamStatesElemNEW[nelem][nelem-ipISO][level]); 
					dynamics.error_scale2 += POW2(UpstreamStatesElemNEW[nelem][nelem-ipISO][level]);
				}
			}
		}
	}


	for(mol=0; mol < N_H_MOLEC; mol++) 
	{
		double f1;
		if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
		{
			ymol = 
				Old_H2_molec[ipyUpstream][mol]/Old_hden[ipyUpstream] +
				(Old_H2_molec[ipyUpstream+1][mol]/Old_hden[ipyUpstream+1] - 
				 Old_H2_molec[ipyUpstream][mol]/Old_hden[ipyUpstream])*
				ynextfrac;
		}
		else
		{		
			long ipBound = ipyUpstream;
			if (-1 == ipBound)
				ipBound = 0;
			ymol = Old_H2_molec[ipBound][mol]/Old_hden[ipBound];
		}
		if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
		{
			hmol = 
				Old_H2_molec[iphUpstream][mol]/Old_hden[iphUpstream] +
				(Old_H2_molec[iphUpstream+1][mol]/Old_hden[iphUpstream+1] - 
				 Old_H2_molec[iphUpstream][mol]/Old_hden[iphUpstream])*
				hnextfrac;
		}
		else
		{		
			long ipBound = iphUpstream;
			if (-1 == ipBound)
				ipBound = 0;
			hmol = Old_H2_molec[ipBound][mol]/Old_hden[ipBound];
		}

		/* the proposed thickness of the next zone */
		f1 = fabs(ymol - Upstream_H2_molec[mol] );
		if( f1 > SMALLFLOAT )
		{
			dynamics.dRad = MIN2(dynamics.dRad,fabs(Dyn_dr*STEP_FACTOR) * 
				/* don't pay attention to species with abundance relative to H less tghan 1e-6 */
				MAX2(ymol + Upstream_H2_molec[mol],1e-6 ) / f1 );
		}

		/* Must be consistent with convergence_error below */
		dynamics.discretization_error += POW2(ymol-2.*hmol+Upstream_H2_molec[mol]); 
		dynamics.error_scale2 += POW2(Upstream_H2_molec[mol]-ymol);
	}

	for(mol=0; mol < mole.num_comole_calc; mol++) 
	{
		double f1;
		if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
		{
			ymol = 
				Old_CO_molec[ipyUpstream][mol]/Old_hden[ipyUpstream] +
				(Old_CO_molec[ipyUpstream+1][mol]/Old_hden[ipyUpstream+1] - 
				 Old_CO_molec[ipyUpstream][mol]/Old_hden[ipyUpstream])*
				ynextfrac;
		}
		else
		{		
			long ipBound = ipyUpstream;
			if (-1 == ipBound)
				ipBound = 0;
			ymol = Old_CO_molec[ipBound][mol]/Old_hden[ipBound];
		}
		if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
		{
			hmol = 
				Old_CO_molec[iphUpstream][mol]/Old_hden[iphUpstream] +
				(Old_CO_molec[iphUpstream+1][mol]/Old_hden[iphUpstream+1] - 
				 Old_CO_molec[iphUpstream][mol]/Old_hden[iphUpstream])*
				hnextfrac;
		}
		else
		{		
			long ipBound = iphUpstream;
			if (-1 == ipBound)
				ipBound = 0;
			hmol = Old_CO_molec[ipBound][mol]/Old_hden[ipBound];
		}

		/* the proposed thickness of the next zone */
		f1 = fabs(ymol - Upstream_CO_molec[mol] );
		if( f1 > SMALLFLOAT )
		{
			dynamics.dRad = MIN2(dynamics.dRad,fabs(Dyn_dr*STEP_FACTOR) * 
				/* don't pay attention to species with abundance relative to H less than 1e-6 */
				MAX2(ymol + Upstream_CO_molec[mol],1e-6 ) / f1 );
		}

		/* Must be consistent with convergence_error below */
		dynamics.discretization_error += POW2(ymol-2.*hmol+Upstream_CO_molec[mol]); 
		dynamics.error_scale2 += POW2(Upstream_CO_molec[mol]-ymol);
	}

	if( dynamics.lgTracePrint )
	{
		fprintf(ioQQQ," DynaStartZone, %4li photo=%.2e , H recom= %.2e dil %.2e \n",
			nzone,dynamics.Rate , dynamics.Source[0][0] , dilution*dense.gas_phase[ipHYDROGEN] );
	}
	return;
}

/* DynaEndZone called at end of zone calculation when advection is turned on */
void DynaEndZone(void)
{
	DEBUG_ENTRY( "DynaEndZone()" );

	/* this routine is called at the end of a zone calculation, by ZoneEnd */

	DivergePresInteg += wind.windv*(DynaFlux(radius.depth)-DynaFlux(radius.depth-radius.drad)); 

	if(dynamics.lgTracePrint)
		fprintf(ioQQQ,"Check dp: %g %g mom %g %g mas %g\n",
			wind.windv*(DynaFlux(radius.depth)-DynaFlux(radius.depth-radius.drad)),
			2*wind.windv*DynaFlux(radius.depth)*radius.drad/(1e16-radius.depth),
			wind.windv*DynaFlux(radius.depth),
			wind.windv*DynaFlux(radius.depth)*(1e16-radius.depth)*(1e16-radius.depth),
			DynaFlux(radius.depth)*(1e16-radius.depth)*(1e16-radius.depth));
	return;
}


/* ============================================================================== */
/* DynaIterEnd called at end of iteration when advection is turned on */
void DynaIterEnd(void)
{
	/* this routine is called by IterRestart at the end of an iteration 
	 * when advection is turned on.  currently it only derives a 
	 * dynamics.timestep by looking at the spatial derivative of
	 * some stored quantities */
	long int i;
	static long int nTime_dt_array_element = 0;

	DEBUG_ENTRY( "DynaIterEnd()" );

	/* set stopping outer radius to current outer radius once we have let
	 * solution relax dynamics.n_initial_relax times
	 * note the off by one confusion - relax is 2 by default,
	 * want to do two static iterations then start dynamics 
	 * iteration was incremented before this call so iteration == 2 at
	 * end of first iteration */
	if( iteration == dynamics.n_initial_relax+1)
	{
		fprintf(ioQQQ,"DYNAMICS DynaIterEnd sets stop radius to %.2e after "
			"dynamics.n_initial_relax=%li iterations.\n",
			radius.depth,
			dynamics.n_initial_relax);
		for( i=iteration-1; i<iterations.iter_malloc; ++i )
			/* set stopping radius to current radius, this stops
			 * dynamical solutions from overrunning the upstream scale
			 * and extending to very large radius due to unphysical heat
			 * appearing to advect into region */
			radius.StopThickness[i] = radius.depth;
	}

	DivergePresInteg = 0.;

	/* This routine is only called if advection is turned on at the end of
	 * each iteration.  The test 
	 * is to check whether wind velocity is also set by dynamics code */

	/* !dynamics.lgStatic true - a true dynamical model */
	if( !dynamics.lgTimeDependentStatic )
	{
		if(iteration == dynamics.n_initial_relax+1 )
		{
			/* let model settle down for n_initial_relax iterations before we begin
			 * dynamics.n_initial_relax set with "set dynamics relax" 
			 * command - it gives the first iteration where we do dynamics 
			 * note the off by one confusion - relax is 2 by default,
			 * want to do two static iterations then start dynamics 
			 * iteration was incremented before this call so iteration == 2 
			 * at end of first iteration */
			if( dynamics.AdvecLengthInit> 0. )
			{
				Dyn_dr =  dynamics.AdvecLengthInit;
			}
			else
			{
				/* -ve dynamics.adveclengthlimit sets length as fraction of first iter total depth */
				Dyn_dr = -dynamics.AdvecLengthInit*radius.depth;
			}

			if (wind.windv0 > 0)
				Dyn_dr = -Dyn_dr;

			if( dynamics.lgTracePrint )
			{
				fprintf(ioQQQ," DynaIterEnd, dr=%.2e \n",
					Dyn_dr );
			}
		}
		else if(iteration > dynamics.n_initial_relax+1 )
		{
			/* evaluate errors and update Dyn_dr */
			DynaNewStep();
		}
	}
	else
	{
		/* this is time-dependent static model */
		static double HeatInitial=-1. , HeatRadiated=-1. ,
			DensityInitial=-1;
		/* n_initial_relax is number of time-steady models before we start 
		 * to evolve, set with "set dynamics relax" command */
		Dyn_dr =  0.;
		fprintf(ioQQQ,
			"DEBUG times enter dynamics.timestep %.2e elapsed_time %.2e iteration %li relax %li \n",
			dynamics.timestep ,
			dynamics.time_elapsed,
			iteration , dynamics.n_initial_relax);
		if( iteration > dynamics.n_initial_relax )
		{
			/* evaluate errors */
			DynaNewStep();

			/* this is set true on CORONAL XXX TIME INIT command, says to use constant
			 * temperature for first n_initial_relax iterations, then let run free */
			if( dynamics.lg_coronal_time_init  )
			{
				thermal.lgTemperatureConstant = false;
				thermal.ConstTemp = 0.;
			}

			/* time variable branch, now without dynamics */
			/* elapsed time - don't increase dynamics.time_elapsed during first two
			 * two iterations since this sets static model */
			dynamics.time_elapsed += dynamics.timestep;
			/* dynamics.timestep_stop is -1 if not set with explicit stop time */
			if( dynamics.timestep_stop > 0 && dynamics.time_elapsed > dynamics.timestep_stop )
			{
				dynamics.lgStatic_completed = true;
			}

			/* stop lowest temperature time command */
			if( (phycon.te < StopCalc.TempLoStopIteration) ||
				(phycon.te > StopCalc.TempHiStopIteration ) )
				dynamics.lgStatic_completed = true;

			/* this is heat radiated, after correction for change of H density in constant
			 * pressure cloud */
			HeatRadiated += (thermal.ctot-dynamics.Cool()) * dynamics.timestep *
				(DensityInitial / dense.gas_phase[ipHYDROGEN]);
		}
		else
		{
			/* this branch, during initial relaxation of solution */
			HeatInitial = 1.5 * pressure.PresGasCurr;
			HeatRadiated = 0.;
			DensityInitial = dense.gas_phase[ipHYDROGEN];
			fprintf(ioQQQ,"DEBUG relaxing times requested %li this is step %li\n", 
				dynamics.n_initial_relax , iteration);
		}
		fprintf(ioQQQ,"DEBUG heat conser HeatInitial=%.2e HeatRadiated=%.2e\n",
			HeatInitial , HeatRadiated );

		/* at this point dynamics.time_elapsed is the time at the end of the 
		 * previous iteration.  We need dt for the next iteration */
		if( dynamics.time_elapsed > time_elapsed_time[nTime_dt_array_element] )
		{
			/* time has advanced to next table point, 
				* set dynamics.timestep to specified value */
			++nTime_dt_array_element;
			/* this is an assert since we already qualified the array
				* element above */
			ASSERT( nTime_dt_array_element < nTime_flux );

			/* option to set flag for recombination logic */
			if( lgtime_Recom[nTime_dt_array_element] )
			{
				fprintf(ioQQQ,"DEBUG dynamics turn on recombination logic\n");
				dynamics.lgRecom = true;
				/* set largest possible zone thickness to value on previous
					* iteration when light was on - during recombination conditions
					* become much more homogeneous and dr can get too large,
					* crashing into H i-front */
				radius.sdrmax = radius.dr_max_last_iter/3.;
				radius.lgSdrmaxRel = false;
			}

			if( lgtime_dt_specified )
			{
				/* this is the new dynamics.timestep */
				fprintf(ioQQQ,"DEBUG lgtimes increment Time to %li %.2e\n" ,nTime_dt_array_element,
					dynamics.timestep);
				dynamics.timestep = time_dt[nTime_dt_array_element];
				/* option to change time step factor - default is 1.2 set in DynaZero */
				if( time_dt_scale_factor[nTime_dt_array_element] > 0. )
					dynamics.timestep_factor = time_dt_scale_factor[nTime_dt_array_element];
			}
		}
		else if( lgtime_dt_specified )
		{
			/* we are between two points in the table, increase dynamics.timestep */
			dynamics.timestep *= dynamics.timestep_factor;
			fprintf(ioQQQ,"DEBUG lgtimes increment Timeby dynamics.timestep_factor  to %li %.2e\n" ,
				nTime_dt_array_element,
				dynamics.timestep );
			if(dynamics.time_elapsed+dynamics.timestep > time_elapsed_time[nTime_dt_array_element] )
			{
				fprintf(ioQQQ,"DEBUG lgtimes but reset to %.2e\n" ,dynamics.timestep );
				dynamics.timestep = 1.0001*(time_elapsed_time[nTime_dt_array_element]-dynamics.time_elapsed);
			}
		}
		else
		{
			/* time not specified in third column, so use initial */
			dynamics.timestep = timestep_next();
		}

		if( cosmology.lgDo )
		{
			cosmology.redshift_step = -(1.f + cosmology.redshift_current ) *
				GetHubbleFactor( cosmology.redshift_current ) * (realnum)dynamics.timestep;
			cosmology.redshift_current += cosmology.redshift_step;
		}

		fprintf(ioQQQ,"DEBUG times exit dynamics.timestep %.2e elapsed_time %.2e scale %.2e ",
			dynamics.timestep ,
			dynamics.time_elapsed,
			rfield.time_continuum_scale);

		if( cosmology.lgDo )
		{
			fprintf(ioQQQ,"redshift %.3e ", cosmology.redshift_current );
		}

		fprintf(ioQQQ,"\n" );
	}

	/* Dyn_dr == 0 is for static time dependent cloud */
	ASSERT( (iteration < dynamics.n_initial_relax+1) ||
		Dyn_dr != 0. || (Dyn_dr == 0. && wind.lgStatic()) );

	/* reset the upstream counters */
	ipUpstream = iphUpstream = ipyUpstream = -1;
	dynamics.discretization_error = 0.;
	dynamics.error_scale2 = 0.;

	/* save results from previous iteration */
	DynaSaveLast();
	return;
}

/*DynaNewStep work out convergence errors */
STATIC void DynaNewStep(void)
{
	long int ilast = 0,
		i,
		nelem,
		ion,
		mol;

	double frac_next=-BIGFLOAT,
		Oldi_hden,
		Oldi_ion,
		Oldi_iso,
		Oldi_mol;

	DEBUG_ENTRY( "DynaNewStep()" );

	/*n = MIN2(nzone, NZLIM-1);*/
	dynamics.convergence_error = 0;
	dynamics.error_scale1 = 0.;

	ASSERT( nzone < struc.nzlim);
	for(i=0;i<nzone;++i) 
	{
		/* Interpolate for present position in previous solution */
		while( (Old_depth[ilast] < struc.depth[i] ) && 
			( ilast < nOld_zone-1 ) )
		{
			++ilast;
		}
		ASSERT( ilast <= nOld_zone-1 );

		if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
		{
			frac_next = ( struc.depth[i] - Old_depth[ilast])/
				(Old_depth[ilast+1] - Old_depth[ilast]);
			Oldi_hden = Old_hden[ilast] +
				(Old_hden[ilast+1] - Old_hden[ilast])*
				frac_next;
		}
		else
		{
			Oldi_hden = Old_hden[ilast];
		}
		/* Must be consistent with discretization_error above */
		/* >>chngf 02 aug 01, multiply by cell width */
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			for( ion=0; ion<nelem+2; ++ion )
			{
				if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
				{
					Oldi_ion = (Old_xIonDense[ilast][nelem][ion] +
						(Old_xIonDense[ilast+1][nelem][ion]-Old_xIonDense[ilast][nelem][ion])*
						frac_next);
				}
				else
				{
					Oldi_ion = Old_xIonDense[ilast][nelem][ion];
				}
				dynamics.convergence_error += POW2(Oldi_ion/Oldi_hden-struc.xIonDense[nelem][ion][i]/struc.hden[i]) /* *struc.dr[i] */;  

				/* >>chng 02 nov 11, add first error scale estimate from Robin */
				dynamics.error_scale1 += POW2(struc.xIonDense[nelem][ion][i]/struc.hden[i]);			
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_local[ipISO][nelem]; ++level )
					{
						if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
						{
							Oldi_iso = (Old_StatesElemNEW[ilast][nelem][nelem-ipISO][level] +
							  (Old_StatesElemNEW[ilast+1][nelem][nelem-ipISO][level]-Old_StatesElemNEW[ilast][nelem][nelem-ipISO][level])*
													frac_next);
						}
						else
						{
							Oldi_iso = Old_StatesElemNEW[ilast][nelem][nelem-ipISO][level];
						}
						dynamics.convergence_error += POW2(Oldi_iso/Oldi_hden-struc.StatesElemNEW[nelem][nelem-ipISO][level][i]/struc.hden[i]) /* *struc.dr[i] */;  

						/* >>chng 02 nov 11, add first error scale estimate from Robin */
						dynamics.error_scale1 += POW2(struc.StatesElemNEW[nelem][nelem-ipISO][level][i]/struc.hden[i]);			
					}
				}
			}
		}
		for( mol=0; mol < N_H_MOLEC; mol++)
		{
			if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
			{
				Oldi_mol = (Old_H2_molec[ilast][mol] +
										(Old_H2_molec[ilast+1][mol]-Old_H2_molec[ilast][mol])*
										frac_next);
			}
			else
			{
				Oldi_mol = Old_H2_molec[ilast][mol];
			}
			dynamics.convergence_error += POW2(Oldi_mol/Oldi_hden-struc.H2_molec[mol][i]/struc.hden[i]) /* *struc.dr[i] */;  

			/* >>chng 02 nov 11, add first error scale estimate from Robin
			 * used to normalize the above convergence_error */
			dynamics.error_scale1 += POW2(struc.H2_molec[mol][i]/struc.hden[i]);			
		}
		for( mol=0; mol < mole.num_comole_calc; mol++)
		{
			if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
			{
				Oldi_mol = (Old_CO_molec[ilast][mol] +
					(Old_CO_molec[ilast+1][mol]-Old_CO_molec[ilast][mol])*
					frac_next);
			}
			else
			{
				Oldi_mol = Old_CO_molec[ilast][mol];
			}
			dynamics.convergence_error += POW2(Oldi_mol/Oldi_hden-struc.CO_molec[mol][i]/struc.hden[i]);  

			/* >>chng 02 nov 11, add first error scale estimate from Robin
			 * used to normalize the above convergence_error */
			dynamics.error_scale1 += POW2(struc.CO_molec[mol][i]/struc.hden[i]);
		}
	}

	/* convergence_error is an estimate of the convergence of the solution from its change during the last iteration,
		 discretization_error is an estimate of the accuracy of the advective terms, calculated in DynaStartZone above:
		 if dominant error is from the advective terms, need to make them more accurate.
	*/

	/* report properties of previous iteration */
	fprintf(ioQQQ,"DYNAMICS DynaNewStep: Dyn_dr %.2e convergence_error %.2e discretization_error %.2e error_scale1 %.2e error_scale2 %.2e\n",
		Dyn_dr, dynamics.convergence_error , dynamics.discretization_error ,
		dynamics.error_scale1 , dynamics.error_scale2 
		);

	/* >>chng 02 nov 29, dynamics.convergence_tolerance is now set to 0.1 in init routine */
	if( dynamics.convergence_error < dynamics.convergence_tolerance*dynamics.discretization_error ) 
		Dyn_dr /= 1.5;
	return;
}

/*DynaSaveLast save results from previous iteration */
STATIC void DynaSaveLast(void)
{
	long int i,
		ion,
		nelem,
		mol;

	DEBUG_ENTRY( "DynaSaveLast()" );

	/* Save results from previous iteration */
	nOld_zone = nzone;
	dynamics.oldFullDepth = struc.depth[nzone-1];
	ASSERT( nzone < struc.nzlim );
	for( i=0; i<nzone; ++i )
	{
		Old_histr[i] = struc.histr[i];
		Old_depth[i] = struc.depth[i];
		Old_xLyman_depth[i] = struc.xLyman_depth[i];
		/* old n_p density from previous iteration */
		Old_hiistr[i] = struc.hiistr[i];
		/* old pressure from previous iteration */
		Old_pressure[i] = struc.pressure[i];
		/* old electron density from previous iteration */
		Old_ednstr[i] = struc.ednstr[i];
		/* energy term */
		Old_EnthalpyDensity[i] = EnthalpyDensity[i];
		/* >>chng 02 May 2001 rjrw: add hden for dilution */
		Old_hden[i] = struc.hden[i];
		Old_DenMass[i] = struc.DenMass[i];

		for(mol=0;mol<N_H_MOLEC;mol++)
		{
			Old_H2_molec[i][mol] = struc.H2_molec[mol][i];
		}
		for(mol=0;mol<mole.num_comole_calc;mol++)
		{
			Old_CO_molec[i][mol] = struc.CO_molec[mol][i];
		}

		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			Old_gas_phase[i][nelem] = dense.gas_phase[nelem];
			for( ion=0; ion<nelem+2; ++ion )
			{
				Old_xIonDense[i][nelem][ion] = struc.xIonDense[nelem][ion][i];
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_max[ipISO][nelem]; ++level )
					{ 
						Old_StatesElemNEW[i][nelem][nelem-ipISO][level] = struc.StatesElemNEW[nelem][nelem-ipISO][level][i];
						ASSERT( !isnan( Old_StatesElemNEW[i][nelem][nelem-ipISO][level] ) );
					}
				}
			}
		}		
	}
	return;
}

realnum DynaFlux(double depth)

{
	realnum flux;

	DEBUG_ENTRY( "DynaFlux()" );

	if(dynamics.FluxIndex == 0) 
	{
		flux = (realnum)dynamics.FluxScale;
	}	
	else 
	{
		flux = (realnum)(dynamics.FluxScale*pow(fabs(depth-dynamics.FluxCenter),dynamics.FluxIndex));
		if(depth < dynamics.FluxCenter)
			flux = -flux;
	}
	if(dynamics.lgFluxDScale) 
	{
		/*flux *= struc.DenMass[0]; */
		/* WJH 21 may 04, changed to use dense.xMassDensity0, which should be strictly constant */
		flux *= dense.xMassDensity0; 
	}
	return flux;
}

/* ============================================================================== */
/*DynaZero zero some dynamics variables, called from zero.c, 
 * before parsing commands */
void DynaZero( void )
{
	int ipISO;

	DEBUG_ENTRY( "DynaZero()" );

	/* the number of zones in the previous iteration */
	nOld_zone = 0;

	/* by default advection is turned off */
	dynamics.lgAdvection = false;
	/*dynamics.Velocity = 0.;*/
	AdvecSpecificEnthalpy = 0.;
	dynamics.Cool_r = 0.;
	dynamics.Heat_v = 0.;
	dynamics.dHeatdT = 0.;
	dynamics.HeatMax = 0.;
	dynamics.CoolMax = 0.;
	dynamics.Rate = 0.;

	/* sets recombination logic, keyword RECOMBINATION on a time step line */
	dynamics.lgRecom = false;

	/* don't force populations to equilibrium levels */
	dynamics.lgEquilibrium = false;

	/* set true if time dependent calculation is finished */
	dynamics.lgStatic_completed = false;

	/* vars that determine whether time dependent soln only - set with time command */
	dynamics.lgTimeDependentStatic = false;
	dynamics.timestep_init = -1.;
	/* this factor multiplies the time step */
	dynamics.timestep_factor = 1.2;
	dynamics.time_elapsed = 0.;

	/* set the first iteration to include dynamics rather than constant pressure */
	/* iteration number, initial iteration is 1, default is 2 - changed with SET DYNAMICS FIRST command */
	dynamics.n_initial_relax = 2;

	/* set initial value of the advection length,
	 * neg => fraction of depth of init model, + length cm */
	dynamics.AdvecLengthInit = -0.1;

	/* this is a tolerance for determining whether dynamics has converged */
	dynamics.convergence_tolerance = 0.1;

	/* this says that set dynamics pressure mode was set */
	dynamics.lgSetPresMode = false;

	/* set default values for uniform mass flux */
	dynamics.FluxScale = 0.;
	dynamics.lgFluxDScale = false;
	dynamics.FluxCenter = 0.;
	dynamics.FluxIndex = 0.;
	dynamics.dRad = BIGFLOAT;

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* factor to allow turning off advection for one of the iso seq, 
		 * this is done with command "no advection h-like" or "he-like" 
		 * only for testing */
		dynamics.lgISO[ipISO] = true;
	}
	/* turn off advection for rest of ions, command "no advection metals" */
	dynamics.lgMETALS = true;
	/* turn off thermal effects of advection, command "no advection cooling" */
	dynamics.lgCoolHeat = true;
	DivergePresInteg = 0.;

	dynamics.discretization_error = 0.;
	dynamics.error_scale2 = 0.;
	return;
}


/* ============================================================================== */
/* DynaCreateArrays allocate some space needed to save the dynamics structure variables, 
 * called from DynaCreateArrays */
void DynaCreateArrays( void )
{
	long int nelem,
		ns,
		i,
		ion,
		mol;

	DEBUG_ENTRY( "DynaCreateArrays()" );

	Upstream_H2_molec = (double*)MALLOC((size_t)N_H_MOLEC*sizeof(double) );
	Upstream_CO_molec = (double*)MALLOC((size_t)mole.num_comole_calc*sizeof(double) );

	dynamics.H2_molec = (double*)MALLOC((size_t)N_H_MOLEC*sizeof(double) );
	dynamics.CO_molec = (double*)MALLOC((size_t)mole.num_comole_calc*sizeof(double) );

	UpstreamElem = (double*)MALLOC((size_t)LIMELM*sizeof(double) );

	dynamics.Source = ((double**)MALLOC( (size_t)LIMELM*sizeof(double *) ));
	UpstreamIon = ((double**)MALLOC( (size_t)LIMELM*sizeof(double *) ));
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		dynamics.Source[nelem] = ((double*)MALLOC( (size_t)(nelem+2)*sizeof(double) ));
		UpstreamIon[nelem] = ((double*)MALLOC( (size_t)(nelem+2)*sizeof(double) ));
		for( ion=0; ion<nelem+2; ++ion )
		{
			dynamics.Source[nelem][ion] = 0.;
		}
	}

	UpstreamStatesElemNEW = ((double***)MALLOC( (size_t)LIMELM*sizeof(double **) ));
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		if( dense.lgElmtOn[nelem] )
		{
			UpstreamStatesElemNEW[nelem] = (double**)MALLOC(sizeof(double*)*(unsigned)(nelem+1) );
			for( long ion=0; ion<nelem+1; ion++ )
			{
				long ipISO = nelem-ion;
				if( ipISO < NISO )
				{
					UpstreamStatesElemNEW[nelem][nelem-ipISO] = (double*)MALLOC(sizeof(double)*(unsigned)iso.numLevels_max[ipISO][nelem]);
				}
				else
				{
					fixit();  // for now, point non-iso ions to NULL
					UpstreamStatesElemNEW[nelem][nelem-ipISO] = NULL;
				}
			}
		}
	}


	dynamics.StatesElemNEW = ((double***)MALLOC( (size_t)LIMELM*sizeof(double **) ));
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		if( dense.lgElmtOn[nelem] )
		{
			dynamics.StatesElemNEW[nelem] = (double**)MALLOC(sizeof(double*)*(unsigned)(nelem+1) );
			for( long ion=0; ion<nelem+1; ion++ )
			{
				long ipISO = nelem-ion;
				if( ipISO < NISO )
				{
					dynamics.StatesElemNEW[nelem][nelem-ipISO] = (double*)MALLOC(sizeof(double)*(unsigned)iso.numLevels_max[ipISO][nelem]);
				}
				else
				{
					fixit();  // for now, point non-iso ions to NULL
					dynamics.StatesElemNEW[nelem][nelem-ipISO] = NULL;
				}
			}
		}
	}

	dynamics.Rate = 0.;

	Old_EnthalpyDensity = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_ednstr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	EnthalpyDensity = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_DenMass = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_hden = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_pressure = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_histr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_hiistr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_depth = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_xLyman_depth = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	Old_xIonDense = (realnum ***)MALLOC(sizeof(realnum **)*(unsigned)(struc.nzlim) );

	Old_StatesElemNEW = (realnum ****)MALLOC(sizeof(realnum ***)*(unsigned)(struc.nzlim) );

	Old_gas_phase = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)(struc.nzlim) );

	Old_H2_molec = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)(struc.nzlim) );

	Old_CO_molec = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)(struc.nzlim) );

	/* now create diagonal of space for ionization arrays */
	for( ns=0; ns < struc.nzlim; ++ns )
	{
		Old_xIonDense[ns] = 
			(realnum**)MALLOC(sizeof(realnum*)*(unsigned)(LIMELM+3) );

		Old_StatesElemNEW[ns] = 
			(realnum***)MALLOC(sizeof(realnum**)*(unsigned)(LIMELM) );

		Old_gas_phase[ns] = 
			(realnum*)MALLOC(sizeof(realnum)*(unsigned)(LIMELM+3) );

		Old_H2_molec[ns] = 
 		 (realnum*)MALLOC(sizeof(realnum)*(unsigned)(N_H_MOLEC) );

		Old_CO_molec[ns] = 
 		 (realnum*)MALLOC(sizeof(realnum)*(unsigned)(mole.num_comole_calc) );

		for( nelem=0; nelem< (LIMELM+3);++nelem )
		{
			Old_xIonDense[ns][nelem] = 
				(realnum*)MALLOC(sizeof(realnum)*(unsigned)(LIMELM+1) );
		}

		for( nelem=0; nelem< LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				Old_StatesElemNEW[ns][nelem] = 
					(realnum**)MALLOC(sizeof(realnum*)*(unsigned)(nelem+1) );
				for( ion=0; ion<nelem+1; ion++ )
				{
					long ipISO = nelem-ion;
					if( ipISO < NISO )
					{
						Old_StatesElemNEW[ns][nelem][ion] = 
							(realnum*)MALLOC(sizeof(realnum)*(unsigned)iso.numLevels_max[ipISO][nelem]);
					}
					else
					{
						fixit();  // for now, point non-iso ions to NULL
						Old_StatesElemNEW[ns][nelem][ion] = NULL;
					}
				}
			}
		}
	}

	for( i=0; i < struc.nzlim; i++ )
	{
		/* these are values if H0 and tau_912 from previous iteration */
		Old_histr[i] = 0.;
		Old_xLyman_depth[i] = 0.;
		Old_depth[i] = 0.;
		dynamics.oldFullDepth = 0.;
		/* old n_p density from previous iteration */
		Old_hiistr[i] = 0.;
		/* old pressure from previous iteration */
		Old_pressure[i] = 0.;
		/* old electron density from previous iteration */
		Old_ednstr[i] = 0.;
		Old_hden[i] = 0.;
		Old_DenMass[i] = 0.;
		Old_EnthalpyDensity[i] = 0.;
		for( nelem=0; nelem< (LIMELM+3);++nelem )
		{
			for( ion=0; ion<LIMELM+1; ++ion )
			{
				Old_xIonDense[i][nelem][ion] = 0.;
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso.numLevels_max[ipISO][nelem]; ++level )
					{
						Old_StatesElemNEW[i][nelem][nelem-ipISO][level] = 0.;
					}
				}
			}
		}
		for(mol=0;mol<N_H_MOLEC;mol++)
		{
			Old_H2_molec[i][mol] = 0.;
		}
		for(mol=0;mol<mole.num_comole_calc;mol++)
		{
			Old_CO_molec[i][mol] = 0.;
		}
	}
	return;
}

/*advection_set_default - called to set default conditions
 * when time and wind commands are parsed,
 * lgWind is true if dynamics, false if time dependent */
STATIC void advection_set_default( bool lgWind )
{

	DEBUG_ENTRY( "advection_set_default()" );

	/* turn on advection */
	dynamics.lgAdvection = true;

	/* turn off prediction of next zone's temperature, as guessed in ZoneStart,
	 * also set with no tepredictor */
	thermal.lgPredNextTe = false;

	/* turn off both CO and H2 networks since advection not included, also say physical
	 * conditions are not ok*/
	/* the heavy element molecules are turned off for the moment,
	 * until I am confident that the H-only models work robustly */
	/* co.lgNoH2Mole = true; */
	/** \todo	1	get following to work */
	/* >>chng 06 jun 29, add advective terms to CO solver */
	co.lgNoCOMole = true; 
#	if 0
	co.lgNoCOMole = true; /* >>chng 04 apr 23, tried to rm this line - problems */
	phycon.lgPhysOK = false;/* >>chng 04 apr 23, rm this line */
#	endif
	/* >>chng 06 nov 29, there is a conservation problem in the ionization -
	 * molecular solvers that is demonstrated by the */
	/**\todo	0	document this problem with molecules - time dependent */

	/* use the new temperature solver
	strcpy( conv.chSolverEden , "new" ); */

	/*  constant total pressure, gas+rad+incident continuum
	 *  turn on radiation pressure */
	pressure.lgPres_radiation_ON = true;
	pressure.lgPres_magnetic_ON = true;
	pressure.lgPres_ram_ON = true;

	/* we need to make the solvers much more exact when advection is in place */
	if( lgWind )
	{
		/* increase precision of solution */
		conv.EdenErrorAllowed = 1e-3;
		/* the actual relative error is relative to the total heating and cooling,
		 * which include the dynamics.heat and .cool, which are the advected heating/cooling.
		 * the two terms can be large and nearly cancel, what is written to the .heat and cool files
		 * by save files has had the smaller of the two subtracted, leaving only the net advected 
		 * heating and cooling */
		conv.HeatCoolRelErrorAllowed = 3e-4f;
		conv.PressureErrorAllowed = 1e-3f;

		if( cosmology.lgDo )
		{
			conv.EdenErrorAllowed = 1e-5;
			conv.PressureErrorAllowed = 1e-5f;
		}
	}
	return;
}

/* ============================================================================== */
/* ParseDynaTime parse the time command, called from ParseCommands */
void ParseDynaTime( Parser &p )
{
	DEBUG_ENTRY( "ParseDynaTime()" );

	/*flag set true when time dependent only */
	dynamics.lgTimeDependentStatic = true;

	dynamics.timestep_init = p.getNumberCheckAlwaysLogLim("dynamics.timestep",30.);

	dynamics.timestep = dynamics.timestep_init;
	if( p.nMatch( "TRAC" ) )
		dynamics.lgTracePrint = true;

	/* this is the stop time and is optional */
	dynamics.timestep_stop = p.getNumberDefaultAlwaysLog("stop time", -1.);

	/* set default flags - false says that time dependent, not dynamical solution */
	advection_set_default(false);

	wind.windv0 = 0.;
	wind.setStatic();
	wind.windv = wind.windv0;

	/* this is used in convpres to say wind solution - both cases use this*/
	//strcpy( dense.chDenseLaw, "WIND" );

	/* create time step and flux arrays */
	time_elapsed_time = (double*)MALLOC((size_t)NTIME*sizeof(double));
	time_flux_ratio = (double*)MALLOC((size_t)NTIME*sizeof(double));
	time_dt = (double*)MALLOC((size_t)NTIME*sizeof(double));
	time_dt_scale_factor = (double*)MALLOC((size_t)NTIME*sizeof(double));
	lgtime_Recom = (int*)MALLOC((size_t)NTIME*sizeof(int));

	/* number of lines we will save */
	nTime_flux = 0;

	/* get the next line, and check for eof */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
			" Hit EOF while reading time-continuum list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* third column might set dt - if any third column is missing then
	 * this is set false and only time on command line is used */
	lgtime_dt_specified = true;

	while( p.strcmp("END") != 0 )
	{
		if( nTime_flux >= NTIME )
		{
			fprintf( ioQQQ, 
				" Too many time points have been entered; the limit is %d.  Increase variable NTIME in dynamics.c.\n", 
			  NTIME );
			cdEXIT(EXIT_FAILURE);
		}

		if( p.nMatch("CYCLE") )
		{
			double period = p.getNumberCheckAlwaysLog("log time");
			ASSERT( period > time_elapsed_time[nTime_flux-1] );
			long pointsPerPeriod = nTime_flux;
			while( nTime_flux < NTIME - 1 )
			{
				time_elapsed_time[nTime_flux] = period + time_elapsed_time[nTime_flux-pointsPerPeriod];
				time_flux_ratio[nTime_flux] = time_flux_ratio[nTime_flux-pointsPerPeriod];
				time_dt[nTime_flux] = time_dt[nTime_flux-pointsPerPeriod];
				time_dt_scale_factor[nTime_flux] = time_dt_scale_factor[nTime_flux-pointsPerPeriod];
				nTime_flux++;
			}
			//Tell the code to continue cyclically by equating two named time points
			fprintf( ioQQQ, " Adding cycles with period = %e s.\n", period );

			/* get next line and check for eof */
			p.getline();
			if( p.m_lgEOF )
			{
				fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			continue;
		}

		time_elapsed_time[nTime_flux] = p.getNumberCheckAlwaysLog("log time");
		time_flux_ratio[nTime_flux] = p.getNumberCheckAlwaysLog("log flux ratio");

		/* this is optional dt to set time step - if not given then initial
		 * time step is always used */
		time_dt[nTime_flux] =  p.getNumberDefaultAlwaysLog("log time step",-1.);

		/* if any of these are not specified then do not use times array */
		if( time_dt[nTime_flux] < 0.0 )
			lgtime_dt_specified = false;

		/* this is optional scale factor to increase time */
		time_dt_scale_factor[nTime_flux] = p.getNumberDefaultAlwaysLog(
			"scale factor to increase time",-1.);

		/* turn on recombination front logic */
		if( p.nMatch("RECOMBIN") )
		{
			/* this sets flag dynamics.lgRecom true so that all of code knows recombination
			 * is taking place */
			lgtime_Recom[nTime_flux] = true;
		}
		else
		{
			lgtime_Recom[nTime_flux] = false;
		}

		/* this is total number stored so far */
		++nTime_flux;

		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}

	}

	if( nTime_flux < 2 )
	{
		fprintf( ioQQQ, " At least two instances of time must be specified.  There is an implicit instance at t=0.\n" \
			" The user must specify at least one additional time.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	for( long i=0; i < nTime_flux; i++ )
	{
		fprintf( ioQQQ, "DEBUG time dep %.2e %.2e %.2e %.2e\n",
			time_elapsed_time[i],
			time_flux_ratio[i] ,
			time_dt[i],
			time_dt_scale_factor[i]);
	}
	fprintf( ioQQQ, "\n" );
	return;
}
/* ============================================================================== */
/* ParseDynaWind parse the wind command, called from ParseCommands */
void ParseDynaWind( Parser &p)
{
	int iVelocity_Type;
	bool lgModeSet=false;
	/* compiler flagged possible paths where dfdr used but not set -
	 * this is for safety/keep it happy */
	double dfdr=-BIGDOUBLE;

	DEBUG_ENTRY( "ParseDynaWind()" );

	if( p.nMatch( "TRAC" ) )
		dynamics.lgTracePrint = true;

	/* Flag for type of velocity law: 
	 * 1 is original, give initial velocity at illuminated face
	 * 2 is face flux gradient (useful if face velocity is zero),
	 * set to zero, but will be reset if velocity specified */
	iVelocity_Type = 0;
	/* wind structure, parameters are initial velocity and optional mass
	 * v read in in km s-1 and convert to cm s-1, mass in solar masses */
	if( p.nMatch( "VELO" ) )
	{
		wind.windv0 = (realnum)(p.getNumberPlain("velocity")*1e5);
		wind.windv = wind.windv0;
		wind.setDefault();
		iVelocity_Type = 1;
	}

	if( p.nMatch( "BALL" ) )
	{
		wind.setBallistic();
		lgModeSet = true;
	}

	if( p.nMatch( "STAT" ) )
	{
		wind.windv0 = 0.;
		wind.setStatic();
		lgModeSet = true;
		iVelocity_Type = 1;
	}
	
 	if ( 1 == iVelocity_Type && !lgModeSet)
	{
		if (wind.windv0 > 0.)
		{
			fprintf(ioQQQ,"Warning, BALListic option needed to switch off pressure gradient terms\n");
		}
		else if (wind.windv0 == 0.)
		{
			fprintf(ioQQQ,"Warning, STATic option needed for zero speed solutions\n");
		}	
	}

	if( p.nMatch("DFDR") )
	{
		/* velocity not specified, rather mass flux gradient */
		dfdr = p.getNumberPlain("flux gradient");
		iVelocity_Type = 2;
	}

	/* center option, gives xxx */
	if( p.nMatch("CENT") )
	{
		/* physical length in cm, can be either sign */
		dynamics.FluxCenter = p.getNumberPlain(
			"centre of mass flux distribution");
	}

	/* flux index */
	if( p.nMatch("INDE") )
	{
		/* power law index */
		dynamics.FluxIndex = p.getNumberPlain(
			"power law index of mass flux distribution");
	}

	/* the case where velocity was set */
	if(iVelocity_Type == 1) 
	{
		/* was flux index also set? */
		if(dynamics.FluxIndex == 0)
		{
			dynamics.FluxScale = wind.windv0;
			dynamics.lgFluxDScale = true;
			/* Center doesn't mean much in this case -- make sure it's
			 * in front of grid so DynaFlux doesn't swap signs where
			 * it shouldn't */
			dynamics.FluxCenter = -1.;
		}
		else
		{
			/** \todo	2	Need to include divergence terms in pressure balance 
			 * if flux index is != 0 */
			/* velocity was set but flux index was not set - estimate it */
			dynamics.FluxScale = wind.windv0*
				pow(fabs(dynamics.FluxCenter),-dynamics.FluxIndex);

			dynamics.lgFluxDScale = true;
			if(dynamics.FluxCenter > 0)
			{
				dynamics.FluxScale = -dynamics.FluxScale;
			}
		}
	} 
	/* the case where flux gradient is set */
	else if(iVelocity_Type == 2)
	{
		if(dynamics.FluxIndex == 0)
		{
			fprintf(ioQQQ,"Can't specify gradient when flux is constant!\n");
			/* use this exit handler, which closes out MPI when multiprocessing */
			cdEXIT(EXIT_FAILURE);
		}
		/** \todo	2	Need to include divergence terms in pressure balance 
		 * if flux index is != 0 */
		/* Can't specify FluxScale from dvdr rather than dfdr, as
		 * d(rho)/dr != 0 */ 
		dynamics.FluxScale = dfdr/dynamics.FluxIndex*
			pow(fabs(dynamics.FluxCenter),1.-dynamics.FluxIndex);
		if(dynamics.FluxCenter > 0)
		{
			dynamics.FluxScale = -dynamics.FluxScale;
		}
		dynamics.lgFluxDScale = false;

		/* put in bogus value simply as flag -- assume that surface velocity
		 * is small or we wouldn't be using this to specify. */
		wind.windv0 = -0.01f;
		wind.setDefault();
	}
	else
	{
		/* assume old-style velocity-only specification */ 
		/* wind structure, parameters are initial velocity and optional mass
		 *  v in km/sec, mass in solar masses */
		wind.windv0 = (realnum)(p.getNumberCheck("wind velocity")*1e5);
		if (wind.windv0 < 0.) 
		{
			wind.setDefault();
		}
		else if (wind.windv0 > 0.)
		{
			wind.setBallistic();
		}
		else
		{
			wind.setStatic();
		}

		dynamics.FluxScale = wind.windv0;
		dynamics.FluxIndex = 0.;
		dynamics.lgFluxDScale = true;
		/* Center doesn't mean much in this case -- make sure it's
		 * in front of grid so DynaFlux doesn't swap signs where
		 * it shouldn't */
		dynamics.FluxCenter = -1.;
	}

	wind.windv = wind.windv0;

#	ifdef FOO
	fprintf(ioQQQ,"Scale %g (*%c) Index %g Center %g\n",
		dynamics.FluxScale,(dynamics.lgFluxDScale)?'D':'1',
		dynamics.FluxIndex,dynamics.FluxCenter);
#	endif

	/* option to include advection */
	if( p.nMatch( "ADVE" ) )
	{
		/* set default flags - true says dynamical solution */
		advection_set_default(true);
	}

	else
	{
		/* this is usual hypersonic outflow */
		if( wind.windv0 <= 1.e6 )
		{
			/* speed of sound roughly 10 km/s */
			fprintf( ioQQQ, " >>>>Initial wind velocity should be greater than speed of sound; calculation only valid above sonic point.\n" );
			wind.lgWindOK = false;
		}

		/* set the central object mass, in solar masses */
		wind.comass = (realnum)p.getNumberDefault("central object mass",1.);
		/* default is one solar mass */

		/* option for rotating disk, keyword is disk */
		wind.lgDisk = false;
		if( p.nMatch( "DISK") )
			wind.lgDisk = true;

	}

	/* this is used in convpres to say wind solution - both cases use this*/
	strcpy( dense.chDenseLaw, "WIND" );

	/*  option to turn off continuum radiative acceleration */
	if( p.nMatch("NO CO") )
	{
		pressure.lgContRadPresOn = false;
	}
	else
	{
		pressure.lgContRadPresOn = true;
	}
	return;
}

/*DynaPrtZone called to print zone results */
void DynaPrtZone( void )
{

	DEBUG_ENTRY( "DynaPrtZone()" );

	ASSERT( nzone>0 && nzone<struc.nzlim );

	if( nzone > 0 )
	{
		fprintf(ioQQQ," DYNAMICS Advection: Uad %.2f Uwd%.2e FRCcool: %4.2f Heat %4.2f\n",
			timesc.sound_speed_adiabatic/1e5 ,
			wind.windv/1e5 ,
			dynamics.Cool()/thermal.ctot,
			dynamics.Heat()/thermal.ctot);
	}

	ASSERT( EnthalpyDensity[nzone-1] > 0. );

	fprintf(ioQQQ," DYNAMICS Eexcit:%.4e Eion:%.4e Ebin:%.4e Ekin:%.4e ET+pdv %.4e EnthalpyDensity/rho%.4e AdvSpWork%.4e\n",
		phycon.EnergyExcitation,
		phycon.EnergyIonization,
		phycon.EnergyBinding,
  		0.5*POW2(wind.windv)*dense.xMassDensity,
		5./2.*pressure.PresGasCurr ,
		EnthalpyDensity[nzone-1]/dense.gas_phase[ipHYDROGEN] , AdvecSpecificEnthalpy
	);
	return;
}

/*DynaPunchTimeDep - save info about time dependent solution */
void DynaPunchTimeDep( FILE* ipPnunit , const char *chJob )
{

	DEBUG_ENTRY( "DynaPunchTimeDep()" );

	if( strcmp( chJob ,  "END" ) == 0 )
	{
		double te_mean,
			H2_mean,
			H0_mean,
			Hp_mean,
			Hep_mean;
		/* save info at end */
		if( cdTemp(
			/* four char string, null terminated, giving the element name */
			"HYDR", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 means that chLabel is a special case */
			2, 
			/* will be temperature */
			&te_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"HYDR", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			2, 
			/* will be fractional ionization */
			&Hp_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"HYDR", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			1, 
			/* will be fractional ionization */
			&H0_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"H2  ", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			0, 
			/* will be fractional ionization */
			&H2_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"HELI", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			2, 
			/* will be fractional ionization */
			&Hep_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		fprintf( ipPnunit , 
			"%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n" , 
			dynamics.time_elapsed , 
			dynamics.timestep ,
			rfield.time_continuum_scale , 
			dense.gas_phase[ipHYDROGEN],
			te_mean , 
			Hp_mean , 
			H0_mean , 
			H2_mean , 
			Hep_mean ,
			/* ratio of CO to total H column densities */
			findspecies("CO")->hevcol / SDIV( colden.colden[ipCOL_HTOT] ),
			cosmology.redshift_current,
			dense.eden/dense.gas_phase[ipHYDROGEN]
			);
	}
	else
		TotalInsanity();
	return;
}

/*DynaSave save dynamics - info related to advection */
void DynaSave(FILE* ipPnunit , char chJob )
{
	DEBUG_ENTRY( "DynaSave()" );

	if( chJob=='a' )
	{
		/* this is save dynamics advection, the only save dynamics */
		fprintf( ipPnunit , "%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			radius.depth_mid_zone,
			thermal.htot , 
			dynamics.Cool() , 
			dynamics.Heat() , 
			dynamics.dCooldT() ,
			dynamics.Source[ipHYDROGEN][ipHYDROGEN],
			dynamics.Rate,
			phycon.EnthalpyDensity/dense.gas_phase[ipHYDROGEN] ,
			AdvecSpecificEnthalpy
			);
	}
	else
		TotalInsanity();
	return;
}

double t_dynamics::Heat()
{
	return Heat_v*dense.gas_phase[ipHYDROGEN];
}

double t_dynamics::Cool()
{
	return Cool_r*phycon.EnthalpyDensity;
}

double t_dynamics::dCooldT()
{
	return Cool_r*5./2.*pressure.PresGasCurr/phycon.te;
}

void DynaIterStart(void)
{
	DEBUG_ENTRY( "DynaIterStart()" );
	
	if( 0 == nTime_flux )
	{
		rfield.time_continuum_scale = 1.;
		return;
	}
	else if( dynamics.time_elapsed <= time_elapsed_time[0] )
	{
		/* if very early times not specified assume no flux variation yet */
		rfield.time_continuum_scale = (realnum)time_flux_ratio[0];
	}
	else if( dynamics.time_elapsed > time_elapsed_time[nTime_flux-1] )
	{
		fprintf( ioQQQ, 
					" PROBLEM - DynaIterStart - I need the continuum at time %.2e but the table ends at %.2e.\n" ,
					dynamics.time_elapsed ,
					time_elapsed_time[nTime_flux-1]);
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		rfield.time_continuum_scale = (realnum)linint(
			/* the times in seconds */
			time_elapsed_time, 
			/* the rfield.time_continuum_scale factors */
			time_flux_ratio, 
			/* the number of rfield.time_continuum_scale factors */
			nTime_flux,
			/* the desired time */
			dynamics.time_elapsed);
	}
	
	fprintf(ioQQQ,"DEBUG time dep reset continuum iter %ld dynamics.timestep %.2e elapsed time %.2e scale %.2e",
			  iteration,
			  dynamics.timestep ,
			  dynamics.time_elapsed, 
			  rfield.time_continuum_scale);
	if( dynamics.lgRecom )
	{
		fprintf(ioQQQ," recom");
	}
	fprintf(ioQQQ,"\n");
	
	/* make sure that at least one continuum source is variable */
	long int nTimeVary = 0;
	for( long int i=0; i < rfield.nShape; i++ )
	{
		/* this is set true if particular continuum source can vary with time 
		 * set true if TIME appears on intensity / luminosity command line */
		if( rfield.lgTimeVary[i] )
			++nTimeVary;
	}
	
	if( hextra.lgTurbHeatVaryTime )
	{
		/* vary extra heating */
		hextra.TurbHeat = hextra.TurbHeatSave * rfield.time_continuum_scale;
		fprintf(ioQQQ,"DEBUG TurbHeat vary new heat %.2e\n",
				  hextra.TurbHeat);
	}
	else if( !nTimeVary )
	{
		fprintf(ioQQQ," DISASTER - there were no variable continua "
				  "or heat sources - put TIME option on at least one "
				  "luminosity or hextra command.\n");
		cdEXIT(EXIT_FAILURE);
	}
}

