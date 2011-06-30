/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*molcol generate and print molecular column densities */
#include "cddefines.h"
#include "radius.h"
#include "colden.h"
#include "h2.h"
#include "mole.h"
#include "atomfeii.h"
#include "molcol.h"

void molcol(
	const char *chLabel,
	/* file for printout */
	FILE *ioMEAN )
{
	long int i;

	DEBUG_ENTRY( "molcol()" );

	if( strcmp(chLabel,"PRIN") == 0 )
	{
		/* total hydrogen column density, all forms */
		fprintf( ioMEAN, "\n                                                     Log10 Column density (cm^-2)\n");
		fprintf( ioMEAN, "   Htot  :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_HTOT])));
		fprintf( ioMEAN, "   HII   :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_Hp])));
		fprintf( ioMEAN, "   HI    :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_H0])));
		fprintf( ioMEAN, "   H-    :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_HMIN])));
		fprintf( ioMEAN, "   H2g   :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_H2g])));
		fprintf( ioMEAN, "   H2*   :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_H2s])));
		fprintf( ioMEAN, "   H2+   :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_H2p])));
		fprintf( ioMEAN, "   HeH+  :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_HeHp] )));
		fprintf( ioMEAN, "\n");
		fprintf( ioMEAN, "   H3+   :");
		fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,colden.colden[ipCOL_H3p] )));
		fprintf( ioMEAN, "\n");
	}

	/* call large H2 and CO column density routines which will do their jobs */
	FeII_Colden( chLabel);
	H2_Colden( chLabel);

	if( strcmp(chLabel,"ZERO") == 0 )
	{
		/*  zero out the column densities */
		for( i=0; i < mole.num_comole_calc; i++ )
		{
			COmole[i]->hevcol = 0.;
		}
	}

	else if( strcmp(chLabel,"ADD ") == 0 )
	{
		/*  add together column densities */
		for( i=0; i < mole.num_comole_calc; i++ )
		{
			COmole[i]->hevcol += COmole[i]->hevmol*(realnum)radius.drad_x_fillfac;
		}
	}

	else if( strcmp(chLabel,"PRIN") == 0 )
	{
		/*  print the molecular column densities
		 * want to print all the molecules, not including the atoms/ions
		 * that are part of the co solver.  use first to print them all */
		/*for( i=0; i < mole.num_comole_calc; i++ )*/
		int j=0;
		for( i=0; i < mole.num_comole_calc; i++ )
		{
			if(COmole[i]->n_nuclei <= 1)
				continue;			
			/* print 7 column densities per line */
			if( j!=0 && j%8==0 )
				fprintf( ioMEAN, "\n" );
			fprintf( ioMEAN, "   %-6.6s:", COmole[i]->label );
			fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,COmole[i]->hevcol )));
			j++;
		}
		fprintf( ioMEAN, "\n" );
	}

	else
	{
		fprintf( ioMEAN, " molcol does not understand the label %4.4s\n", 
		  chLabel );
		cdEXIT(EXIT_FAILURE);
	}
	return;

}
