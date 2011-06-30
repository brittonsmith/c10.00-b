/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseNorm parse parameters on the normalize command */
#include "cddefines.h"
#include "lines.h"
#include "input.h"
#include "parser.h"
#include "lines_service.h"

void ParseNorm(Parser &p)
{
	char chLabel[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "ParseNorm()" );

	/* these are flags saying that normalization line has been set */
	LineSave.lgNormSet = true;

	/* >>chng 01 aug 23, insist on a line label */
	/* 
	 * get possible label - must do first since it can contain a number.*/
	/* is there a double quote on the line?  if so then this is a line label */
	if( p.nMatch(  "\"" ) )
	{

		/* GetQuote does the following -
		 * first copy original version of name into chLabel, 
		 * string does include null termination.
		 * set label in OrgCard and second parameter to spaces so 
		 * that not picked up below as keyword */
		p.GetQuote( chLabel , true );
		if( chLabel[4] != '\0' || strlen(chLabel) != 4 )
		{
			fprintf( ioQQQ, " The label identifying the line on the normalize command must be exactly 4 char long.\n" );
			fprintf( ioQQQ, " The command line was as follows:\n %s\n", input.chCardSav[input.nRead] );
			fprintf( ioQQQ, " The label I found was: \"%s\", where were not 4 characters between the quotes.\n", chLabel );
			fprintf( ioQQQ, "Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* copy first four char of label into caps, and null terminate*/
		cap4( LineSave.chNormLab, chLabel);
	}
	else
	{
		fprintf( ioQQQ, "The normalize command does not have a valid line.\n" );
		fprintf( ioQQQ, "A 4 char long line label must also be specified, between double quotes, like \"H  1\" 4861.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* normalise lines to this rather than h-b, sec number is scale factor */
	LineSave.WavLNorm = (realnum)p.getWave();
	
	if( LineSave.WavLNorm < 0 )
	{
		fprintf( ioQQQ, "A negative wavelength does not make sense to me.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* get the error assocated with the 4 significant figures that are visible,
	 * wavelength of 0 (a continuum) has error of zero */
	LineSave.errorwave = WavlenErrorGet( LineSave.WavLNorm );

	LineSave.ScaleNormLine = p.FFmtRead();

	if( p.lgEOL() )
		LineSave.ScaleNormLine = 1.;

	/* confirm that scale factor is positive */
	if( LineSave.ScaleNormLine <= 0. )
	{
		fprintf( ioQQQ, " The scale factor for relative intensities must be greater than zero.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
