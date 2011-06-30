/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_colden parse save column density command, or actually do the save lines output */
#include "cddefines.h"
#include "cddrive.h"
#include "input.h"
#include "save.h"
#include "parser.h"
/* this is the limit to the number of lines we can store */
#define	NPUNLM	100L

static char chElement[NPUNLM][5];
static long int nColdenEntered;
static long int ionstage[NPUNLM];

/*save_colden parse save column density command, or actually do the save lines output */
void parse_save_colden(
	Parser &p,
  /* the header for the file, a list of identifications */
  char chHeader[] )
{
	char chTemp[INPUT_LINE_LENGTH];

	long int i;

	DEBUG_ENTRY( "parse_save_colden()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine */
	
	/* number of lines we will save */
	nColdenEntered = 0;
	
	/* get the next line, and check for eof */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	while( p.strcmp( "END" ) != 0 )
	{
		if( nColdenEntered >= NPUNLM )
		{
			fprintf( ioQQQ, 
						" Too many lines have been entered; the %ld limit is.  Increase variable NPUNLM in routine save_colden.\n", 
						NPUNLM );
			cdEXIT(EXIT_FAILURE);
		}
		
		/* order on line is label (col 1-4), ionstage */
		strncpy( chElement[nColdenEntered], p.getCommand(4).c_str() , 4 );
		
		/* null terminate the string*/
		chElement[nColdenEntered][4] = 0;
		
		/* now get ionstage - 1 for atom, 2 for first ion, etc */
		ionstage[nColdenEntered] = (long)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("ion stage");
		
		/* this is total number stored so far */
		++nColdenEntered;
		
		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
	}
	
	/*fprintf( ioPUN, "%li lines were entered, they were;\n", 
	  nColdenEntered );*/
	/* give header line */
	
	sprintf( chHeader , "#colden %s %3li", chElement[0] , ionstage[0] );
	for( i=1; i < nColdenEntered; i++ )
	{
		sprintf( chTemp, "\t%s %3li", chElement[i] , ionstage[i] );
		strcat( chHeader, chTemp );
	}
	strcat( chHeader, "\n" );
}

/*save_colden parse save column density command, or actually do the save lines output */
void save_colden(
  /* the file we will write to */
  FILE * ioPUN )
{
	long int i;

	DEBUG_ENTRY( "save_colden()" );

	/* save some column densities  */
	double colden;
	/* save some column column densities command */
	for( i=0; i < nColdenEntered; i++ )
	{
		if( i )
			fprintf(ioPUN,"\t");
		/* get column density, returns 0 if all ok */
		if( cdColm(
				 /* four char string, null terminated, giving the element name */
				 chElement[i], 
				 /* IonStage is ionization stage */
				 ionstage[i], 
				 /* will be column density */
				 &colden) )
		{
			fprintf( ioQQQ, 
						"\n PROBLEM save_colden could not find a column density for "
						"the species with label %s %li \n\n",
						chElement[i] , ionstage[i] );
			colden = 1.;
		}
		fprintf( ioPUN, "%.4f", log10( MAX2(SMALLFLOAT , colden ) ) );
	}
	fprintf( ioPUN, "\n" );
}

