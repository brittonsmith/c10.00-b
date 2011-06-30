/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveLineData punches selected line data for all lines transferred in code */
/*Save1LineData save data for one line */
#include "cddefines.h"
#include "taulines.h"
#include "hmi.h"
#include "iso.h"
#include "phycon.h"
#include "physconst.h"
#include "elementnames.h"
#include "hydrogenic.h"
#include "lines_service.h"
#include "dense.h"
#include "atomfeii.h"
#include "lines.h"
#include "atmdat.h"
#include "prt.h"
#include "h2.h"
#include "thermal.h"
#include "cooling.h"
#include "save.h"

NORETURN void SaveLineData(FILE * ioPUN)
{

	long int i, 
	  j, 
	  limit ,
	  nelem ,
	  ipHi , 
	  ipLo;

	const long nskip=2; /* number of emission lines per line of output */
	double tot;
	bool lgElemOff=false;

	DEBUG_ENTRY( "SaveLineData()" );

	/* routine punches out (on unit ioPUN) line data
	 * for all recombination lines, and all transitions that are transferred */

	/* say what is happening so we know why we stopped */
	fprintf( ioQQQ, " punching line data, then stopping\n" );

	/* first check that all lines are turned on */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( !dense.lgElmtOn[nelem] )
		{
			fprintf(ioQQQ," WARNING - I am punching line data but element %s is turned off.\n",
			elementnames.chElementName[nelem]);
			lgElemOff = true;
		}
	}
	if( lgElemOff )
	{
		fprintf(ioQQQ,"Some elements are turned off and save line data requested.\n");
		fprintf(ioQQQ,"Code is now designed to do save line data only with all elements on.\n");
		fprintf(ioQQQ,"Please try again with all elements on.\n");
		fprintf(ioQQQ,"Please try again with all elements on.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* evaluate rec coefficient at constant temperature if this is set, else
	 * use 10,000K */
	double TeNew;
	if( thermal.lgTemperatureConstant )
	{
		TeNew = thermal.ConstTemp;
	}
	else
	{
		TeNew = 1e4;
	}
	TempChange(TeNew , false);

	/* this is set of Dima's recombination lines */
	t_ADfA::Inst().rec_lines(phycon.te,LineSave.RecCoefCNO);
	fprintf( ioPUN, "\n       Recombination lines of C, N, O\n" );
	fprintf( ioPUN, "    Ion  WL(A)   Coef     Ion   WL(A)  Coef\n" );
	for( i=0; i<471; i+=nskip)
	{
		/* nskip is set to 2 above */
		limit = MIN2(471,i+nskip);
		fprintf( ioPUN, "    " );
		for( j=i; j < limit; j++ )
		{
			fprintf( ioPUN, "%2.2s%2ld%6ld%8.3f    ", 
				elementnames.chElementSym[(long)(LineSave.RecCoefCNO[0][j])-1], 
			  (long)(LineSave.RecCoefCNO[0][j]-LineSave.RecCoefCNO[1][j]+1.01), 
			  (long)(LineSave.RecCoefCNO[2][j]+0.5), 
			  log10(SDIV(LineSave.RecCoefCNO[3][j]) ) );
		}
		fprintf( ioPUN, "    \n" );
	}
	fprintf( ioPUN, "\n\n" );

	dense.eden = 1.;
	dense.gas_phase[ipHYDROGEN] = 1.;
	dense.EdenHCorr = 1.;

	/* want very small neutral fractions so get mostly e- cs */
	dense.xIonDense[ipHYDROGEN][1] = 1.e-5f;
	hmi.Hmolec[ipMH2g] = 0.;
	dense.xIonDense[ipHYDROGEN][1] = 1.;
	for( i=1; i <= nLevel1; i++ )
	{
		TauLines[i].Lo->Pop = 1.;
	}

	for( i=0; i < nWindLine; i++ )
	{
		TauLine2[i].Lo->Pop = 1.;
	}

	for( i=0; i < nUTA; i++ )
	{
		UTALines[i].Lo->Pop = 1.;
	}

	for( i=0; i < LIMELM; i++ )
	{
		for( j=0; j < LIMELM+1; j++ )
		{
			dense.xIonDense[i][j] = 1.;
		}
	}

	/* evaluate cooling, this forces evaluation of collision strengths */
	CoolEvaluate(&tot);

	fprintf( ioPUN, "       Level 1 transferred lines\n" );

	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );

	for( i=1; i <= nLevel1; i++ )
	{
		/* chLineLbl generates a 1 char string from the line transfer array info -
		 * "Ne 2  128" the string is null terminated,
		 * in following printout the first 4 char is used first, followed by
		 * an integer, followed by the rest of the array*/
		Save1LineData( &TauLines[i] , ioPUN , true);
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end level 1, start level 2\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );
	for( i=0; i < nWindLine; i++ )
	{
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			Save1LineData( &TauLine2[i] , ioPUN , true);
		}
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end level 2, start inner shell\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );

	for( i=0; i < nUTA; i++ )
	{
		Save1LineData( &UTALines[i] , ioPUN , true);
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end inner shell, start H-like iso seq\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );

	/* h-like iso sequence */
	/* the hydrogen like iso-sequence */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		iso_collide( ipH_LIKE, nelem );
		if( nelem < 2 || dense.lgElmtOn[nelem] )
		{
			/* arrays are dim'd iso.numLevels_max[ipH_LIKE][nelem]+1 */
			/* keep this limit to iso.numLevels_max, instead of _local.  */
			for( ipLo=ipH1s; ipLo < iso.numLevels_max[ipH_LIKE][nelem]-1; ipLo++ )
			{
				for( ipHi=ipLo+1; ipHi < iso.numLevels_max[ipH_LIKE][nelem]; ipHi++ )
				{
					Save1LineData( &Transitions[ipH_LIKE][nelem][ipHi][ipLo] , ioPUN , false );
				}
			}
		}
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end H-like iso seq, start He-like iso seq\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );
	for( nelem=1; nelem < LIMELM; nelem++ )
	{
		if( nelem < 2 || dense.lgElmtOn[nelem] )
		{
			/* arrays are dim'd iso.numLevels_max[ipH_LIKE][nelem]+1  */
			for( ipLo=ipHe1s1S; ipLo < iso.numLevels_max[ipHE_LIKE][nelem]-1; ipLo++ )
			{
				for( ipHi=ipLo+1; ipHi < iso.numLevels_max[ipHE_LIKE][nelem]; ipHi++ )
				{
					Save1LineData( &Transitions[ipHE_LIKE][nelem][ipHi][ipLo] , ioPUN , false );
				}
			}
		}
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end he-like iso seq, start hyperfine structure lines\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );
	/* fine structure lines */
	for( i=0; i < nHFLines; i++ )
	{
		Save1LineData( &HFLines[i] , ioPUN , true);
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end hyperfine, start database lines\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );
	/* Databases: Atoms & Molecules*/
	for( i=0; i < linesAdded2; i++ )
	{
		Save1LineData( dBaseLines[i].tran , ioPUN , true);	
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end database, start satellite lines\n" );
	fprintf( ioPUN, "Ion         WL  gl gu    gf       A       CS   n(crt)\n" );
	for( long ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
			{
				for( i=0; i<iso.numLevels_max[ipISO][nelem]; i++ )
				{
					Save1LineData( &SatelliteLines[ipISO][nelem][i], ioPUN , true );
				}
			}
		}
	}

	/* want very small ionized fractions so get mostly H2 cs */
	dense.eden = 1e-6;
	dense.gas_phase[ipHYDROGEN] = 1e-6f;
	dense.EdenHCorr = 1e-6f;
	dense.xIonDense[ipHYDROGEN][1] = 1.;
	hmi.Hmolec[ipMH2g] = 1.;
	hmi.Hmolec[ipMH2s] = 1.;
	dense.xIonDense[ipHYDROGEN][1] = 1e-6f;

	/* H2 molecule */
	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end satellite, start H2 lines\n" );
	fprintf( ioPUN, "Eu Vu Ju El Vl Jl        WL    gl gu    gf       A       CS   n(crt)\n" );

	/* ioPUN unit, and option to print all possible lines - false indicates
	 * save only significant lines */
	H2_LevelPops();
	H2_Punch_line_data( ioPUN , false );

	/* save FeII data if atom is currently enabled */
	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end H2, start FeII lines\n" );
	fprintf( ioPUN, " Lo  Hi  Ion  label   WL  gl gu    gf       A       CS   n(crt)\n" );

	/* ioPUN unit, and option to print all possible lines - false indicates
	 * save only significant lines */
	FeIIPunData( ioPUN , false );

	/* ChkMonitorstring is searched for by one of the scripts in the nightly run
	 * this run will terminate with no asserts but that is the correct behavior */
	fprintf( ioQQQ , "\n The code is left in a disturbed state after creating the SAVE LINE DATA file.\n"
			" No calculation is actually performed, only the SAVE LINE DATA file is produced.\n"
			" Remove the SAVE LINE DATA command to do the calculation.\n\n ChkMonitorend is ok.\n" );

	/* stop when done, we have done some serious damage to the code */
	cdEXIT(EXIT_SUCCESS);
}

/*Save1LineData save data for one line */
void Save1LineData( transition * t , FILE * ioPUN  , 
	/* flag saying whether to give collision strength too - in multi level atoms
	 * it will be not valid without a great deal more work */
	 bool lgCS_2 )
{

	double CritDen;
	/* these are used for parts of the line label */
	char chLbl[11];

	DEBUG_ENTRY( "Save1LineData()" );

	if( t->ipCont <= 0 )
	{
		return;
	}

	/** \todo	1	define lifetime and collision rate for multi-level species so that the
	 * critical density is derived correctly in this routine.  For now the flag lgCS_2
	 * being true means to save critical den and is only true for two-level systems 
	 * all places where this routine is called with lgCS_2 false need to be fixed */

	/*iWL = iWavLen( t , &chUnits , &chShift );*/
	/* ion label, like C  1 */
	chIonLbl(chLbl , t );
	fprintf(ioPUN,"%s\t", chLbl );

	/* this is the second piece of the line label, pick up string after start */

	/* the wavelength */
	prt_wl(ioPUN, t->WLAng );

	fprintf( ioPUN, " %3ld%3ld",
	  /* lower and upper stat weights */
	  (long)(t->Lo->g), 
	  (long)(t->Hi->g) );

	/* oscillator strength */
	fprintf( ioPUN,PrintEfmt("%9.2e",  t->Emis->gf));

	/* Einstein A for transition */
	fprintf( ioPUN,PrintEfmt("%9.2e",  t->Emis->Aul));

	/* next collision strengths, use different formats depending on size 
	 * of collision strength */
	if( t->Coll.col_str > 100. )
	{
		fprintf( ioPUN, "%7.1f", t->Coll.col_str );
	}
	else if( t->Coll.col_str > 10. )
	{
		fprintf( ioPUN, "%7.2f", t->Coll.col_str );
	}
	else if( t->Coll.col_str > 1. )
	{
		fprintf( ioPUN, "%7.3f", t->Coll.col_str );
	}
	else if( t->Coll.col_str > .01 )
	{
		fprintf( ioPUN, "%7.4f", t->Coll.col_str );
	}
	else if( t->Coll.col_str > 0.0 )
	{
		fprintf( ioPUN, " %.3e", t->Coll.col_str );
	}
	else
	{
		fprintf( ioPUN, "%7.4f", 0. );
	}

	/* now print critical density but only if cs is positive 
	 * >>chng 06 mar 24, add flag lgCS_2 - in multi-level systems do not want
	 * to save cs since not computed properly */
	if( lgCS_2 && t->Coll.col_str> 0. )
	{
		CritDen = t->Emis->Aul * t->Hi->g*phycon.sqrte / (t->Coll.col_str*COLL_CONST);
		CritDen = log10(CritDen);
	}
	else
	{
		CritDen = 0.;
	}
	fprintf( ioPUN, "%7.3f\n",CritDen );
	return;
}
