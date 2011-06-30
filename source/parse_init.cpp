/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseInit bring an initialization file into input stream before parse  */
#include "cddefines.h"
#include "input.h"
#include "trace.h"
#include "parser.h"

void ParseInit(Parser &p)
{
	char *ipEndL;
	char chName[FILENAME_PATH_LENGTH_2];
	long int ip, 
	  k;
	FILE *ioInitFile; /* will use this as pointer to ini file */

	DEBUG_ENTRY( "ParseInit()" );

	/*bring an initialization file into input stream before parsing  */

	/* check whether single quote on line, this was used in c90 */
	if( p.nMatch( "\'" ) )
	{
		fprintf( ioQQQ, 
			" ParseInit found a single quote on this line.  This was used for file names in C90, but double quotes are used now.\n");
		fprintf( ioQQQ, " The single quote has been ignored.\n");
	}

	if( p.nMatch( "\"" ) )
	{
		/* 
		 * if a quote occurs on the line then get the ini file name 
		 * this will also set the name in chCard and OrgCard to spaces
		 * so later keywords do not key off it
		 */
		p.GetQuote( chName, true );
	}
	else
	{
		/* no quote appeared, so this is the default name, cloudy.ini */
		strcpy( chName, "cloudy.ini" );
	}

	/* at this point we have init file name, now make full name 
	 * this can be a local file, or on the path if the key path appears */

	/* option to get cloudy.ini from a path */
	if( p.nMatch("PATH") )
	{
		ioInitFile = open_data( chName, "r" );
	}
	else
	{
		/* just use file name, and try to open file in current directory first */
		ioInitFile = open_data( chName, "r", AS_LOCAL_DATA );
	}

	/* at this point the init file is open, now bring it into the command stack */
	input.nSaveIni = 1;
	ip = NKRD + 1 - input.nSaveIni;
	while( (read_whole_line( input.chCardSav[ip-1],(int)sizeof(input.chCardSav[ip-1]),ioInitFile)!=NULL ) )
	{
		/* add extra space to be trailing space, needed for commands that end with space */
		ipEndL = strrchr( input.chCardSav[ip-1] , '\n' );
		/* make sure that we found the newline */
		if(ipEndL == NULL )
		{
			fprintf(ioQQQ," ParseInit read in a init file line that did not end with a newline\n");
			fprintf(ioQQQ," line was the following=>%s<=\n",input.chCardSav[ip-1]);
			cdEXIT(EXIT_FAILURE);
		}
		/* >>chng 01 oct 22, add cast */
		/* find offset to end of line, the cr */
		k = (long)(ipEndL - input.chCardSav[ip-1]);
		/* replace cr with space */
		input.chCardSav[ip-1][k] = ' ';
		/* add extra space */
		input.chCardSav[ip-1][k+1] = ' ';
		/* finally null terminate the line */
		input.chCardSav[ip-1][k+2] = '\0';
		/* line starting with space is one way to end input stream */
		if( input.chCardSav[ip-1][0]==' ' ) break; 
		/* totally ignore these lines 
		 * >>chng 06 sep 04 use routine to check for comments */
		if( lgInputComment(input.chCardSav[ip-1]) /*input.chCardSav[ip-1][0]=='#' || input.chCardSav[ip-1][0]=='*' ||
			input.chCardSav[ip-1][0]=='%' || input.chCardSav[ip-1][0]=='/'*/ ) 
			continue;

		/* print input lines if trace specified */
		if( trace.lgTrace )
		{
			fprintf( ioQQQ,"initt=%s=\n",input.chCardSav[ip-1] );
		}

		input.nSaveIni += 1;
		ip = NKRD + 1 - input.nSaveIni;
		if( ip <= input.nSave )
		{
			fprintf( ioQQQ, 
				" Too many ini lines.  Total of all input and ini lines cannot exceed NKRD, presently%4i\n", 
			  NKRD );
			cdEXIT(EXIT_FAILURE);
		}
	}
	fclose(ioInitFile);
	/* last one with real data is NKRD+1-nSaveIni */
	input.nSaveIni -= 1;
	return;
}
