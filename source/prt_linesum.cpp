/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtLineSum parse print line sum command to enter set of lines into sum  */
#include "cddefines.h"
#include "cddrive.h"
#include "radius.h"
#include "lines.h"
#include "input.h"
#include "parser.h"
/* this is the limit to the number of lines we can save */
#define	NRDSUM	300L
#include "prt.h"

static char **chSMLab;
static long int *ipLine;
static bool lgFirst=true;
static long nlsum;
static realnum *wavelength;

void ParsePrtLineSum(Parser &p)
{
	/*static char chSMLab[NRDSUM][5];*/
	bool lgEND;
	long int i;

	/* remember whether we have been called before */

	DEBUG_ENTRY( "ParsePrtLineSum()" );

#		if 0
	if( !lgFirst )
	{
		/* error - more than one read in input stream */
		fprintf(ioQQQ," more than one print line sum has appeared - only first one is used.\n");
		fprintf(ioQQQ," Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}
	else
#		endif
		/* >>chng 03 jan 23, if not first call, do not allocate space, 
		 * had aborted, which was bad in optized runs, or in a grid. 
		 * Bug caught by Melekh Bohdan */
		if( lgFirst )
		{
			/* do not malloc space again */
			lgFirst = false;
			wavelength = ((realnum *)MALLOC( sizeof(realnum )*NRDSUM ));
			
			/* create space for the array of array indices for lines*/
			ipLine = ((long int *)MALLOC(NRDSUM*sizeof(long)));
			
			/* create space for the array of labels*/
			chSMLab = ((char **)MALLOC(NRDSUM*sizeof(char *)));
			
			for( i=0; i<NRDSUM; ++i )
			{
				chSMLab[i] = ((char *)MALLOC(5*sizeof(char )));
			}
		}
	
	/* now read in lines */
	nlsum = 0;
	lgEND = false;
	while( !lgEND )
	{
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
		if( p.strcmp("END" ) != 0 )
		{
			if( nlsum >= NRDSUM )
			{
				fprintf( ioQQQ, 
							" Too many lines have been entered; the limit is %li.  Increase NRDSUM in PrtLineSum.\n", 
							NRDSUM );
				cdEXIT(EXIT_FAILURE);
			}
			
			/*  order on line is label (col 1-4), wavelength */
			strncpy( chSMLab[nlsum], p.getCommand(4).c_str() , 4 );
			chSMLab[nlsum][4] = 0;
			wavelength[nlsum] = (realnum)p.getWaveOpt();
			++nlsum;
		}
		else
		{
			lgEND = true;
		}
	}
}
double PrtLineSum(void)
{
	long int i;

	/* remember whether we have been called before */

	double absint, 
	  relint ,
	  sum=-1.;

	DEBUG_ENTRY( "PrtLineSum()" );

	sum = 0.;
	/* this can be called during setup mode, in which case we do nothing */
	if( LineSave.ipass <= 0 )
	{ 
		return sum;
	}
	
	if( nzone == 1 )
	{
		for( i=0; i < nlsum; i++ )
		{
			/* save the array index for each line */
			if( (ipLine[i] = cdLine((char*)chSMLab[i],wavelength[i],&relint,&absint) ) <=0 )
			{
				fprintf( ioQQQ, " PrtLineSum could not fine line %4.4s %5f\n", 
							chSMLab[i], wavelength[i] );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	
	/* now sum the line */
	for( i=0; i < nlsum; i++ )
	{
		/* this version of chLine uses index, does not search*/
		cdLine_ip(ipLine[i],&relint,&absint);
		absint = pow(10.,absint - radius.Conv2PrtInten);
		sum += absint;
	}
	return sum;
}

