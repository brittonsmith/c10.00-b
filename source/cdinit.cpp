/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdInit routine to initialize variables, called at start of calculation */
/*cdPrepareExit prepare termination of the code, but do not terminate yet */
/* unset EXTERN so that everything is defined here */
#define EXTERN 
#include "cddefines.h"

/* used for saving map*/
FILE *ioMAP = NULL;

/* external ZeroNum used to div by zero 
 * ok here since never changed*/
const double ZeroNum = 0.;

/* this must go here since it defines NTA needed for other lines*/
#include "taulines.h"

/* following is true extern in taulines.h */
long nWindLine = NWINDDIM;

#include "abund.h"
#include "atmdat.h"
#include "atoms.h"
#include "atomfeii.h"
#include "monitor_results.h"
#include "broke.h"
#include "ca.h"
#include "called.h"
#include "carb.h"
#include "cddrive.h"
/* this will be set true when cdInit is called.  The definition is in cdInit.
* Other routines will check that this is true when they are called, 
* to verify that cdInit was called first */
bool lgcdInitCalled=false;
#include "colden.h"
#include "conv.h"
#include "continuum.h"
#include "coolheavy.h"
#include "cosmology.h"
#include "dark_matter.h"
#include "dense.h"
#include "doppvel.h"
#include "dynamics.h"
#include "elementnames.h"
#include "embesq.h"
#include "fe.h"
#include "fudgec.h"
#include "geometry.h"
#include "grainvar.h"
#include "grid.h"
#include "h2.h"
#include "h2_priv.h"
bool lgH2_READ_DATA = false;
#include "he.h"
#include "heavy.h"
#include "hextra.h"
#include "hmi.h"
#include "hydrogenic.h"
/* this is set true once space malloced, then never change
* number of levels again with hydrogenic command, 
* also to make sure MALLOC only happens one time  */
bool lgHydroMalloc = false;
/*  */
#include "hyperfine.h"
#include "input.h"
#include "ionbal.h"
#include "iso.h"
#include "iterations.h"
#include "lines.h"
/* these are the definitions of the line save arrays in lines.h */
LinSv *LineSv=NULL;
LinSv *LineSvSortWL=NULL;
#include "magnetic.h"
#include "hcmap.h"
#include "mean.h"
#include "mewecoef.h"
#include "mpi_utilities.h" 
#include "mole.h" 
#include "nitro.h"
#include "noexec.h"
#include "numderiv.h"
#include "oxy.h"
#include "parse.h"
#include "peimbt.h"
#include "phycon.h"
#include "plot.h"
#include "sil.h"
#include "version.h"
/* this is set true when space is allocated for the FeII arrays,
* once this happens the total number of levels cannot be changed with the atom feii levels command */
bool lgFeIIMalloc=false;
/* */
#include "pressure.h"
#include "prt.h"
#include "save.h"
#include "radius.h"
#include "rfield.h"
/* set true when malloced, init to false */
bool lgRfieldMalloced=false;
#include "opacity.h"
bool lgOpacMalloced=false;
#include "rt.h"
#include "secondaries.h"
#include "state.h"
#include "stopcalc.h"
#include "struc.h"
#include "thermal.h"
#include "timesc.h"
#include "trace.h"
#include "warnings.h"
#include "wind.h"
#include "init.h"


/* =================================================================== */
void cdInit(void)
{
	long i;
	double vtest;

	DEBUG_ENTRY( "cdInit()" );

	/* set flag saying that cdInit has been called */
	lgcdInitCalled = true;

	/*test if the following integer types have the correct width*/
	if( sizeof(int16) != 2 || sizeof(uint16) != 2 || sizeof(int32) != 4 || sizeof(uint32) != 4 )
		TotalInsanity();

	/*********************************************************
	 *  on a VAX compile with /G_FLOATING option on FORTRAN; *
	 *  following makes sure this happened.                  *
	 *********************************************************/
	vtest = 1e-35;
	vtest /= 1e35;
	if( vtest == 0. )
	{
		fprintf( ioQQQ, " Something is wrong with the double precision.  Use /g_floating on a VAX\n" );
	}

	/* initialize some variables dealing with cloudy's interaction with machine environment */
	/* if TALK is true then do standard printout
	 * if false then never say anything */
	/* only the master rank produces output */
	called.lgTalk = cpu.lgMPI_talk();
	/* this flag is needed to turn print on to have effect */
	called.lgTalkIsOK = cpu.lgMPI_talk();
	/* means talk not forced off by call to cdTalk*/
	called.lgTalkForcedOff = false;

	optimize.lgNoVary = false;
	optimize.lgVaryOn = false;
	optimize.lgOptimr = false;
	grid.lgGrid = false;
	grid.nGridCommands = 0;

	for( i=0; i<NUM_OUTPUT_TYPES; i++ )
	{
		grid.lgOutputTypeOn[i] = false;
	}

	/* this is a global variable in monitor_results.h, and can be checked by
	 * other routines to see if asserts are ok - (most calculations will not use asserts,
	 * and this will be the only place values are set, although they will be checked in maincl) */
	lgMonitorsOK = true;
	lgBigBotch = false;
	lgPrtSciNot = false;

	/* number of lines entered with cdLine
	 * both check that number less than NKRD, the limit
	 * the line save array is defined from 0 through input.nSave */
	input.nSave = -1;

	/* nRead is the number of the command in the input stream - many optimize options
	 * point to it to refer to the original command.  it is incremented before
	 * it is used, so will become 0.  it is the array element within the stack
	 * of emission lines */
	input.nRead = -1;

	/* this is number of init lines read in */
	input.nSaveIni = 0;
	input.lgUnderscoreFound = false;
	input.lgBracketFound = false;

	/* this is sanity check that lines are read in ok */
	for( i=0; i < NKRD; i++ )
	{
		strcpy( input.chCardSav[i], "error! - no line image input" );
	}

	/* start the timer to log execution time */
	cdSetExecTime();

	/* zero out lots of variables */
	zero();
	return;
}


/* =================================================================== */
/* cdPrepareExit prepare termination of the code, but do not terminate yet
 * this routine should only be called by exception handlers, never from the main code */
void cdPrepareExit()
{
	enum {DEBUG_LOC=false};
	if( DEBUG_LOC )
		fprintf(ioQQQ," cdExit called\n");

	// make sure file descriptors are closed in case they were redirected
	cdInput( "", "" );
	cdOutput( "", "" );

	/* close any open units */
	CloseSaveFiles( true );
}

