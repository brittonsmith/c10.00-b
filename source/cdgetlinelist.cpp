/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdGetLineList routine to read in master list of emission line wavelengths and ids, for
 * generating loc grids */
#include "cddefines.h"
#include "cddrive.h"

/* return value is number of lines, -1 if file could not be opened */
long int cdGetLineList( 
	/* chFile is optional filename, if void then use BLRLineList,
	 * if not void then use file specified */
	const char chFile[] ,
	/* 2d array of null term strings giving line labels char chLabels[nLines][10] */
	char ***chLabels ,
	/* a 1-d array of line wavelengths */
	realnum **wl )
{
	long int i ,
		nLines;
	bool lgDONE;
	FILE *ioData;

	char chLine[FILENAME_PATH_LENGTH_2];
	const char* chFilename;

	DEBUG_ENTRY( "cdGetLineList()" );

	/* first check that cdInit has been called, since we may have to write
	 * error output */
	if( !lgcdInitCalled )
	{
		fprintf(stderr," cdInit must be called before cdGetLineList.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* use default filename LineList_BLR.dat if void string, else use file specified */
	chFilename = ( strlen(chFile) == 0 ) ? "LineList_BLR.dat" : chFile;

	/* we will check local space first, then on path if not present */
	ioData = open_data( chFilename, "r", AS_LOCAL_DATA_TRY );

	if( ioData == NULL )
	{
		/* did not find file, return -1 */
		return -1;
	}

	/* count how many lines are in the file, ignoring all lines
	 * starting with '#' */
	nLines = 0;
	lgDONE = false;
	while( (read_whole_line( chLine , (int)sizeof(chLine) , ioData ) != NULL) && !lgDONE )
	{
		if( chLine[0] == '\n')
		{
			lgDONE = true;
			continue;
		}

		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( (chLine[0] != '#') )
			++nLines;
	}

	*wl = (realnum *)MALLOC( (size_t)(nLines+1)*sizeof(realnum ) );

	/* create 1-d array of string labels */
	*chLabels = (char**)MALLOC((size_t)(nLines+1)*sizeof(char *) );

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioData , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " cdGetLineList could not rewind line list.\n");
		return( -1 );
	}

	/* actually read and save the lines */
	i = 0;
	lgDONE = false;
	while( (read_whole_line( chLine , (int)sizeof(chLine) , ioData ) != NULL) && !lgDONE)
	{
		long j;
		bool lgEOL;

		if( chLine[0] == '\n')
		{
			lgDONE = true;
			continue;
		}
		/* skip lines that begin with # */
		if( chLine[0] == '#')
			continue;

		/* create second dim of space for labels */
		(*chLabels)[i] = (char*)MALLOC(5*sizeof(char) );

		strncpy( (*chLabels)[i] , chLine , 4);
		(*chLabels)[i][4] = 0;

		/* get and save the wavelength */
		j = 5;
		(*wl)[i] = (realnum)FFmtRead(chLine,&j,INPUT_LINE_LENGTH,&lgEOL);

		/* check for optional micron or cm units, else interpret as Angstroms */
		if( chLine[j-1] == 'M' || chLine[j-1] == 'm')
		{
			/* microns */
			(*wl)[i] *= 1e4;
		}
		else if( chLine[j-1] == 'C'  || chLine[j-1] == 'c')
		{
			/* centimeters */
			(*wl)[i] *= 1e8;
		}

		++i;
	}

	fclose( ioData );

	/* return number of lines we found */
	return nLines;
}
