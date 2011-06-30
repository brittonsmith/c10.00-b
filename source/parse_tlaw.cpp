/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTLaw parse parameters on the tlaw command to set temperature as function of depth,
 * currently only does Bertoldi & Draine simple T law */
#include "cddefines.h"
#include "thermal.h"
#include "parser.h"

void ParseTLaw( Parser &p )
{
	DEBUG_ENTRY( "ParseTLaw()" );

	/* this says that some type of temperature law has been specified */
	thermal.lgTLaw = true;
	thermal.lgTemperatureConstant = true;
	thermal.lgTemperatureConstantCommandParsed = true;

	if( p.nMatch("DB96") )
	{
		/* this is to simulate the temperature law given by equation 41 in 
		 * >>refer	H2	temperature law	Draine, B.T., & Bertoldi, Frank, 1996, ApJ, 468, 269-289 */
		thermal.lgTeBD96 = true;

		/* this is the initial temperature for the BD96 temperature law */
		thermal.T0BD96 = 500.f;
		TempChange(thermal.T0BD96 , false);

		/* the coefficient on column density for temp dropoff */
		thermal.SigmaBD96 = 6e-22f;
	}
	else if( p.nMatch("SN99") )
	{
		/* this is to simulate the temperature law given by equation 16 in 
		 * >>refer	H2	temperature law	Sternberg, A., & Neufeld, D.A. 1999, ApJ, 516, 371-380 */
		thermal.lgTeSN99 = true;

		/* this is the inital temperature for the BD96 temperature law */
		thermal.T0SN99 = 500.f;
		TempChange(thermal.T0SN99 , false);
	}
	else
	{
		fprintf(ioQQQ," There must be a keyword on this command.  The one I know about is BD96\n");
		cdEXIT(EXIT_FAILURE);
	}

#if 0
#include "dense.h"
#include "optimize.h"
#include "input.h"
	bool lgEnd;
	long int j;
	/* all remainder is currently dead code, a copy of DLAW command,
	 * which could be activated if needs arose */
	/* call fcn dense_fabden(RADIUS) which uses the ten parameters
	 * N.B.; existing version of dense_fabden must be deleted
	 * >>chng 96 nov 29, added table option */
	if( p.nMatch("TABL") )
	{
		/* when called, read in densities from input stream */
		strcpy( dense.chDenseLaw, "DLW2" );
		if( p.nMatch("DEPT") )
		{
			dense.lgDLWDepth = true;
		}
		else
		{
			dense.lgDLWDepth = false;
		}

		p.getline();
		dense.frad[0] = (realnum)p.FFmtRead();
		dense.fhden[0] = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " No pairs entered - can\'t interpolate.\n Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		dense.nvals = 2;
		lgEnd = false;

		/* read pairs of numbers until we find line starting with END */
		while( !lgEnd && dense.nvals < LIMTABDLAW )
		{
			p.getline();
			lgEnd = p.m_lgEOF;
			if( !lgEnd )
			{
				if( p.strcmp("END") == 0 )
					lgEnd = true;
			}

			if( !lgEnd )
			{
				dense.frad[dense.nvals-1] = (realnum)p.FFmtRead();
				dense.fhden[dense.nvals-1] = (realnum)p.FFmtRead();
				dense.nvals += 1;
			}
		}
		--dense.nvals;

		for( i=1; i < dense.nvals; i++ )
		{
			/* the radius values are assumed to be strictly increasing */
			if( dense.frad[i] <= dense.frad[i-1] )
			{
				fprintf( ioQQQ, " density.in radii must be in increasing order\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	else
	{
		/* this is usual case, call dense_fabden to get density */
		for( j=0; j < 10; j++ )
		{
			dense.DensityLaw[j] = p.FFmtRead();
		}

		/* set flag so we know which law to use later */
		strcpy( dense.chDenseLaw, "DLW1" );

		/* vary option */
		if( optimize.lgVarOn )
		{
			/* NB - there are 5 = LIMEXT numbers on this line - if LIMEXT ever changes,
			 * chnage this too */
			strcpy( optimize.chVarFmt[optimize.nparm], "DLAW %f %f %f %f %f " );

			/* index for where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			for( j=0; j<LIMEXT; ++j )
			{
				optimize.vparm[j][optimize.nparm] = (realnum)dense.DensityLaw[j];
			}
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = LIMEXT;
			++optimize.nparm;
		}
	}
#	endif
	return;
}
