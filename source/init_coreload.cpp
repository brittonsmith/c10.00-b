/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitCoreload one time initialization of core load, called from cdDrive, this sets
* minimum set of values needed for the code to start - called after
* input lines have been read in and checked for VARY or GRID - so
* known whether single or multiple sims will be run */
#include "cddefines.h"
#include "optimize.h"
#include "parse.h"
#include "struc.h"
#include "input.h"
#include "dense.h"
#include "hcmap.h"
#include "h2.h"
#include "version.h"
#include "hextra.h"
#include "mole.h"
#include "heavy.h"
#include "grid.h"
#include "ionbal.h"
#include "iso.h"
#include "taulines.h"
#include "cosmology.h"
#include "physconst.h"
#include "broke.h"

// //////////////////////////////////////////////////////////////////////////
//
//
// NB DO NOT ADD VARIABLES TO THIS FILE!  THE GOAL IS TO REMOVE THIS FILE
// initialization of variables done one time per coreload should be done in
// a constructor for the data
//
//
// //////////////////////////////////////////////////////////////////////////

/*InitCoreload one time initialization of core load, called from cdDrive, this sets
* minimum set of values needed for the code to start - called after
* input lines have been read in and checked for VARY or GRID - so
* known whether single or multiple sims will be run */
void InitCoreload( void )
{
	static int nCalled=0;
	long int nelem;

	DEBUG_ENTRY( "InitCoreload()" );

	/* return if already called */
	if( nCalled )
		return;

	++nCalled;

	hcmap.lgMapOK = true;
	hcmap.lgMapDone = false;

	// subdirectory delimiter character
#	ifdef _WIN32
	strcpy( input.chDelimiter, "\\" );
#	else
	strcpy( input.chDelimiter, "/" );
#	endif

	/* will be reset to positive number when map actually done */
	hcmap.nMapAlloc = 0;
	hcmap.nmap = 0;
	hcmap.lgMapBeingDone = false;

	/* name of output file from optimization run */
	strncpy( chOptimFileName , "optimal.in" , sizeof( chOptimFileName ) );

	/* number of electrons in valence shell that can Compton recoil ionize */
	long int nCom[LIMELM] = 
	{
		1 , 2 ,								/* K 1s shell */
		1 , 2 ,								/* L 2s shell */
		1 , 2 , 3 , 4 , 5 , 6 ,				/* L 2p shell */
		1 , 2 ,								/* M 3s shell */
		1 , 2 , 3 , 4 , 5 , 6 ,				/* M 3p shell */
		1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 ,		/* M 3d shell */
		1 , 2 ,								/* N 4s shell */
		1 , 2								/* N 4p shell */
	};

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		ionbal.nCompRecoilElec[nelem] = nCom[nelem];
	}

	/* list of shells, 1 to 7 */
	strcpy( Heavy.chShell[0], "1s" );
	strcpy( Heavy.chShell[1], "2s" );
	strcpy( Heavy.chShell[2], "2p" );
	strcpy( Heavy.chShell[3], "3s" );
	strcpy( Heavy.chShell[4], "3p" );
	strcpy( Heavy.chShell[5], "3d" );
	strcpy( Heavy.chShell[6], "4s" );

	/* variables for H-like sequence */
	/* default number of levels for hydrogen iso sequence */
	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		iso.n_HighestResolved_max[ipH_LIKE][nelem] = 5;
		iso.nCollapsed_max[ipH_LIKE][nelem] = 2;
	}

	iso.n_HighestResolved_max[ipH_LIKE][ipHYDROGEN] = 10;
	iso.n_HighestResolved_max[ipH_LIKE][ipHELIUM] = 10;

	iso.nCollapsed_max[ipH_LIKE][ipHYDROGEN] = 15;
	iso.nCollapsed_max[ipH_LIKE][ipHELIUM] = 15;

	/* variables for He-like sequence */
	/* "he-like" hydrogen (H-) is treated elsewhere */
	iso.n_HighestResolved_max[ipHE_LIKE][ipHYDROGEN] = -LONG_MAX;
	iso.numLevels_max[ipHE_LIKE][ipHYDROGEN] = -LONG_MAX;
	iso.nCollapsed_max[ipHE_LIKE][ipHYDROGEN] = -LONG_MAX;

	for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
	{
		/* put at least three resolved and 1 collapsed in every element for he-like */
		iso.n_HighestResolved_max[ipHE_LIKE][nelem] = 3;
		iso.nCollapsed_max[ipHE_LIKE][nelem] = 1;
	}

	iso.n_HighestResolved_max[ipHE_LIKE][ipHELIUM] = 6;
	iso.nCollapsed_max[ipHE_LIKE][ipHELIUM] = 20;

	/* And n=5 for these because they are most abundant */
	iso.n_HighestResolved_max[ipHE_LIKE][ipCARBON] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipNITROGEN] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipOXYGEN] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipNEON] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipSILICON] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipMAGNESIUM] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipSULPHUR] = 5;
	iso.n_HighestResolved_max[ipHE_LIKE][ipIRON] = 5;
	/* also set this, for exercising any possible issues with biggest charge models */
	iso.n_HighestResolved_max[ipHE_LIKE][LIMELM-1] = 5;

	iso.chISO[ipH_LIKE] = "H-like ";
	iso.chISO[ipHE_LIKE] = "He-like";

	max_num_levels = 0;
	for( long ipISO = ipH_LIKE; ipISO < NISO; ipISO++ )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* set this to LONG_MAX, reduce to actual number later,
			 * then verify number of levels is not increased after initial coreload */
			iso.numLevels_malloc[ipISO][nelem] = LONG_MAX;
			iso_update_num_levels( ipISO, nelem );
		}
	}

	statesAdded = 0;
	lgStatesAdded = false;
	linesAdded = 0;
	lgLinesAdded = false;

	/** \todo	3	these should be const since cannot change, are flags */
	/* these are used to set trace levels of output by options on atom h2 trace command 
	 * first is minimum level of trace, keyword is FINAL */
	mole.nH2_trace_final = 1;
	/* intermediate level of trace, info per iteration, key ITERATION */
	mole.nH2_trace_iterations = 2;
	/* full trace, keyword is FULL */
	mole.nH2_trace_full = 3;
	/* print matrices used in solving X */
	mole.nH2_trace_matrix = 4;

	/* turn element on if this is first call to this routine,
	* but do not reset otherwise since would clobber how elements were set up */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		/* never turn on element if turned off on first iteration */
		dense.lgElmtOn[nelem] = true;

		/* option to set ionization with element ioniz cmnd,
		* default (false) is to solve for ionization */
		dense.lgSetIoniz[nelem] = false;
		
		// are we using Chianti for this ion?
		for( int ion=0; ion<=nelem+1; ++ion )
			dense.lgIonChiantiOn[nelem][ion] = false;
	}

	/* smallest density to permit in any ion - if smaller then set to zero */
	dense.density_low_limit = SMALLFLOAT * 1e3;
	dense.density_low_limit = MAX2( dense.density_low_limit, 1e-50 );

	/* default cosmic ray background */
	/* >>chng 99 jun 24, slight change to value
	* quoted by 
	* >>refer	cosmic ray	ionization rate	McKee, C.M., 1999, astro-ph 9901370
	* this will produce a total
	* secondary ionization rate of 2.5e-17 s^-1, as tested in 
	* test suite cosmicray.in.  If each ionization produces 2.4 eV of heat,
	* the background heating rate should be 9.6e-29 * n*/
	/* >>chng 04 jan 26, update cosmic ray ionization rate for H0 to
	* >>refer	cosmic ray	ionization	Williams, J.P., Bergin, E.A., Caselli, P., 
	* >>refercon	Myers, P.C., & Plume, R. 1998, ApJ, 503, 689,
	* H0 ionization rate of 2.5e-17 s-1 and a H2 ionization rate twice this
	* >>chng 04 mar 15, comment said 2.5e-17 which is correct, but code produce 8e-17,
	* fix back to correct value 
	*/
	/* NB - the rate is derived from the density.  these two are related by the secondary
	* ionization efficiency problem.  background_rate is only here to provide the relationship
	* for predominantly neutral gas.  the background_density is the real rate. 
	hextra.background_density = 1.99e-9f;*/
	/* >>chng 05 apr 16, to get proper ionization rate in ism_set_cr_rate, where
	* H is forced to be fully atomic, no molecules, density from 1.99 to 2.15 */
	hextra.background_density = 2.15e-9f;
	hextra.background_rate = 2.5e-17f;

	/* initialization for save files - must call after input commands have
	 * been scanned for grid and vary options.  So known if grid or single run 
	 * default save is different for these */
	grid.lgGridDone = false;
	grid.lgStrictRepeat = false;

	/* these are energy range... if not changed with command, 0. says just use energy limits of mesh */
	grid.LoEnergy_keV = 0.;
	grid.HiEnergy_keV = 0.;
	grid.ipLoEnergy = 0;
	grid.ipHiEnergy = 0;

	for( long i=0; i < LIMPAR; ++i )
		optimize.lgOptimizeAsLinear[i] = false;

	/* limit on ionization we check for zoning and prtcomment */
	struc.dr_ionfrac_limit = 1e-3f;

	SaveFilesInit();

	H2_init_coreload();

	/* initialize cosmological information */
	cosmology.redshift_current = 0.f;
	cosmology.redshift_start = 0.f;
	cosmology.redshift_step = 0.f;
	cosmology.omega_baryon = 0.04592f;
	cosmology.omega_rad = 8.23e-5f;
	cosmology.omega_lambda = 0.7299177f;
	cosmology.omega_matter = 0.27f;
	cosmology.omega_k = 0.f;
	/* the Hubble parameter in 100 km/s/Mpc */
	cosmology.h = 0.71f;
	/* the Hubble parameter in km/s/Mpc */
	cosmology.H_0 = 100.f*cosmology.h;
	cosmology.lgDo = false;

	// the code is ok at startup; only init here so that code remains broken
	// in a grid if any single model announces broken
	broke.lgBroke = false;
	broke.lgFixit = false;
	broke.lgCheckit = false;

	return;
}
