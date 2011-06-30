/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
* others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTest parse the test command */
#include "cddefines.h"
#include "atomfeii.h"
#include "input.h"
#include "monitor_results.h"
#include "parser.h"

void ParseTest(Parser &p)
{
	char chStuff[INPUT_LINE_LENGTH];	

	DEBUG_ENTRY( "ParseTest()" );

	/* do a smoke test - 
	* generate test input stream for rapid testing of code */

	/* optional keyword PRINT will cause all the commands to be printed */
	int nPrintTest = p.nMatch("PRIN"  );

	bool lgH2 = p.nMatch(" H2 ");
	bool lgFeII = p.nMatch("FEII") || p.nMatch("FE II");
	bool lgLARG = p.nMatch("LARG");
	bool lgMOLE = p.nMatch("MOLE");

	/* option to turn on the H2 molecule */
	if( lgH2 )
	{
		/* this will both enable the molecule and also set the threshold for actually
		* computing it to a very low value, so that it is actually used in ionized gas */
		sprintf( chStuff , "ATOM H2 LIMIT -20 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseAtomH2(p);
	}

	/* this is option to also turn on large FeII atom */
	if( lgFeII )
	{
		sprintf( chStuff , "ATOM FEII " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseAtomFeII(p);
	}

	/* use largest possible hydrogen atom */
	if( lgLARG )
	{
		sprintf( chStuff , "ATOM H-LIKE ELEMENT HYDROGEN LEVELS LARGER  " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		p.set_point(4);
		ParseAtomISO(ipH_LIKE,p);
	}

	if( lgMOLE )
	{
#		if 0
		{
			/* hydrogen density */
			sprintf( chStuff , "TRACE TEMPERATURE CONVERGENCE " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseTrace(chStuff);
		}
#		endif

		{
			/* hydrogen density */
			sprintf( chStuff , "HDEN 5 " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseHDEN(p);
		}

		{
			/* make a constant temperature model */
			sprintf( chStuff , "CONSTANT TEMPER 50K  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseConstant(p);
		}

		{
			/* continuum to include full energy range */
			sprintf( chStuff , "TABLE ISM  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseTable( p);
		}

		{
			/* extinguish this continuum */
			sprintf( chStuff , "EXTINGUISH 23  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseExtinguish( p );
		}
		
		{
			/* stop in second zone, so we do use the zone increment logic */
			sprintf( chStuff , "STOP ZONE 2  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseStop(p);
		}

		{
			/* set thickness */
			sprintf( chStuff , "SET DR 0  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseSet(p);
		}
		
		{
			/* do Case B so lyman line pumping of H is not important */
			sprintf( chStuff , "CASE B  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseCaseB( p );
		}

		{
			/* set cosmic rays */
			sprintf( chStuff , "COSMIC RAY BACKGROUND  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseCosmicRays( p );
		}

		{
			/* add grains and ism abundances */
			sprintf( chStuff , "ABUNDANCES ISM  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseAbundances( p );
		}

		{
			/* add grains and ism abundances */
			sprintf( chStuff , "CONSTANT GRAIN TEMPERATURE 20K  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseConstant( p );
		}

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR EDEN 0.528 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR MOLECULAR FRACTION H2 -1.193 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR COLUMN CO 1.214 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);			

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR EDEN 0.528 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR HYDROGEN 1 TEMPERATURE 50K " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);
	}

	else
	{
		{
			/* hydrogen density */
			sprintf( chStuff , "HDEN 4 " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseHDEN( p );
		}

		{
			/* make a constant temperature model */
			sprintf( chStuff , "CONSTANT TEMPER 4  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseConstant( p );
		}

		{
			/* continuum to include full energy range */
			sprintf( chStuff , "TABLE AGN  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseTable( p);
		}

		{
			/* set ionization parameter */
			sprintf( chStuff , "IONIZATION PARAMETER -2  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseIonParI( p);
		}

		{
			/* use old abundances */
			sprintf( chStuff , "ABUNDANCES OLD SOLAR 84  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseAbundances(p);
		}

		{
			/* >>chng 02 apr 19, add this */
			/* set total Lyman continuum depth - this is to prevent caution
			 * that Lyman continuum was thin but expected to be thick */
			sprintf( chStuff , "STOP LYMAN OPTICAL -4  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseStop(p);
		}
		
		{
			/* stop in second zone, so we do use the zone increment logic */
			sprintf( chStuff , "STOP ZONE 2  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseStop(p);
		}

		{
			/* set thickness */
			sprintf( chStuff , "SET DR 0  " );
			if( nPrintTest )
				fprintf(ioQQQ , "%s\n" , chStuff );
			p.setline(chStuff);
			p.set_point(4);
			ParseSet(p);
		}

		/* create series of monitor commands */
		sprintf( chStuff , "MONITOR HYDROGEN 1 IONIZATION -3.040 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		sprintf( chStuff , "MONITOR HELIUM 2 IONIZATION -1.067 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		sprintf( chStuff , "MONITOR CARBON 2 IONIZATION -2.301 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.653 to -0.560. Badnell DR by default */
		sprintf( chStuff , "MONITOR CARBON 3 IONIZATION -0.560 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.348 to -0.373. Badnell DR by default */
		sprintf( chStuff , "MONITOR CARBON 4 IONIZATION -0.373 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.490 to -0.530. Badnell DR by default */
		sprintf( chStuff , "MONITOR CARBON 5 IONIZATION -0.530 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.800 to -0.861. Badnell DR by default */
		sprintf( chStuff , "MONITOR OXYGEN 3 IONIZATION  -0.861 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.180 to -0.157. Badnell DR by default */
		sprintf( chStuff , "MONITOR OXYGEN 4 IONIZATION -0.157 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.770 to -0.808. Badnell DR by default */
		sprintf( chStuff , "MONITOR OXYGEN 5 IONIZATION -0.808 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* >>chng 02 apr 19, from 0.7258 to 0.946, due to adding Lyman cont depth */
		/* >>chng 07 oct 22, from 0.946  to 1.108, resolve l-levels of h-like sequence */
		sprintf( chStuff , "MONITOR LINE \"CA B\" 4861.36 1.108 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		/* must have copy of this in chOrgCard, which is used by the routine to get lab */
		p.setline(chStuff);
		ParseMonitorResults(p);

		/* >>chng 02 apr 19, from 2.4603 to 3.25, due to adding Lyman cont depth 
		 * >>chng 06 nov 17, asserted value from 3.25 to 3.11, drift over last year */
		/* >>chng 06 dec 01, from 3.11 to 2.72. Badnell DR by default */
		/* >>chng 07 oct 22, from 2.72 to 3.18, resolve l-levels of h-like sequence */
		sprintf( chStuff , "MONITOR LINE \"O  3\" 5007 3.18 " );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		/* must have copy of this in chOrgCard, which is used by the routine to get lab */
		p.setline(chStuff);
		ParseMonitorResults(p);

		sprintf( chStuff , "MONITOR HTOT -15.011" );
		if( nPrintTest )
			fprintf(ioQQQ , "%s\n" , chStuff );
		p.setline(chStuff);
		ParseMonitorResults(p);
	}

	return;
}
