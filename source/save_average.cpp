/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_average routine to read in list of species to output as averages */
#include "cddefines.h"
#include "cddrive.h"
#include "input.h"
#include "elementnames.h"
#include "save.h"
#include "parser.h"

/* return value is number of lines, -1 if file could not be opened */
void parse_save_average( 
	Parser &p,
	/* the file we will write to */
	long int ipPun, 
	char *chHeader)
{
	long int i;
	long int nLine;

	DEBUG_ENTRY( "parse_save_average()" );

	char chCap[INPUT_LINE_LENGTH];
	char chTemp[INPUT_LINE_LENGTH];

	/* this is index for first line we will read.  use this to count 
	 * total number of species */
	nLine = input.nRead;
	/* very first time this routine is called, chJob is "READ" and we read
	 * in lines from the input stream.  The species labels and other information
	 * are store in the save structure.  These are output in a later call
	 * to this routine with argument PUNCH  */
	
	/* number of lines we will save */
	save.nAverageList[ipPun] = 0;
	
	/* get the next line, and check for eof */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Punch average hit EOF while reading list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
		
	/* keep reading until we hit END */
	while( p.strcmp( "END" ) != 0 )
	{
		/* count number of species we will save */
		++save.nAverageList[ipPun];
		
		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " save averages hit EOF while reading species list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}		
	}
/*#		define PADEBUG*/
#		ifdef PADEBUG
	fprintf(ioQQQ , "DEBUG save_average %li species read in.\n", 
			  save.nAverageList[ipPun] );
#		endif
	
	/* now create space that will be needed to hold these arrays */
	
	if( save.ipPnunit[ipPun] == NULL )
	{
		/* nAverageIonList is set of ions for averages */
		save.nAverageIonList[ipPun] = (int *)MALLOC(
			(size_t)(save.nAverageList[ipPun]*sizeof(int)) );
		
		/* nAverage2ndPar is set of second parameters for averages */
		save.nAverage2ndPar[ipPun] = (int *)MALLOC((size_t)
																  (save.nAverageList[ipPun]*sizeof(int)) );

		/* chAverageType is label for type of average */
		save.chAverageType[ipPun] = (char **)MALLOC((size_t)
																	(save.nAverageList[ipPun]*sizeof(char *)) );
		
		/* chAverageSpeciesLabel is label for species */
		save.chAverageSpeciesLabel[ipPun] = (char **)MALLOC((size_t)
																			  (save.nAverageList[ipPun]*sizeof(char *)) );
		for( i=0; i<save.nAverageList[ipPun]; ++i )
		{
			/* create space for labels themselves */
			save.chAverageType[ipPun][i] = (char *)MALLOC((size_t)
																		  (5*sizeof(char)) );
			
			/* chAverageSpeciesLabel is label for species */
			save.chAverageSpeciesLabel[ipPun][i] = (char *)MALLOC((size_t)
																					 (5*sizeof(char)) );
		}
	}
	
	/* reset array input read to first of the species lines */
	// Evil input rewind CodeSmell -- grrw
	input.nRead = nLine;
	
#		ifdef PADEBUG
	fprintf(ioQQQ , "DEBUG save_average %li species read in.\n", 
			  save.nAverageList[ipPun] );
#		endif
	
	/* reread the input lines and save the data */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Punch average hit EOF while reading list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	/* use this to count number of species, and will assert equal to
	 * total malloced above */
	nLine = 0;
	/* keep reading until we hit END */
	while( p.strcmp("END" ) != 0 )
	{
		/* count number of species we will save */
		++nLine;
		if( p.nMatch("TEMP" ))
		{
			/* temperature */
				strcpy( save.chAverageType[ipPun][nLine-1] , "TEMP" );
		}
		else if( p.nMatch("COLU" ))
		{
			/* column density */
			strcpy( save.chAverageType[ipPun][nLine-1] , "COLU" );
		}
		else if( p.nMatch("IONI" ))
		{
			/* ionization fraction */
			strcpy( save.chAverageType[ipPun][nLine-1] , "IONI" );
		}
		else
		{
			fprintf(ioQQQ,"PROBLEM one of the jobs TEMPerature, COLUmn density, or IONIzation, must appear.\n");
			cdEXIT(EXIT_FAILURE);
		}
		
		/* get element name, a string we shall pass on to the routine
		 * that computes final quantities */
		if( (i = p.GetElem( )) < 0 )
		{
			/* no name found */
			fprintf(ioQQQ, "save average did not an element on this line, sorry %s\n",
					  chCap );
			cdEXIT(EXIT_FAILURE);
		}
		strcpy( save.chAverageSpeciesLabel[ipPun][nLine-1] , elementnames.chElementNameShort[i]);
		
		/* now get ionization stage */
		save.nAverageIonList[ipPun][nLine-1] = (int) p.FFmtRead();
		if( p.lgEOL() )
		{
			/* error - needed that ionization stage */
			p.NoNumb("ionization stage" );
		}
		
		/* look for volume keyword, otherwise will be radius 
		 * only used for some options */
		if( p.nMatch( "VOLU" ) )
		{
			/* volume */
			save.nAverage2ndPar[ipPun][nLine-1] = 1;
		}
		else
		{
			/* radius */
			save.nAverage2ndPar[ipPun][nLine-1] = 0;
		}
		
		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " save averages hit EOF while reading species list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	
	/* these must equal or we have a major logical error */
	ASSERT( nLine == save.nAverageList[ipPun]);
	
#		ifdef PADEBUG
	for( i=0; i<nLine ; ++i )
	{
		fprintf(ioQQQ, "PDDEBUG %s %s %i %i\n",
				  save.chAverageType[ipPun][i],
				  save.chAverageSpeciesLabel[ipPun][i] ,
				  save.nAverageIonList[ipPun][i] ,
				  save.nAverage2ndPar[ipPun][i] );
	}
#		endif
	
	/* save headers */
	sprintf(chHeader, "#averages");
	for( i=0; i<nLine ; ++i )
	{
		sprintf(chTemp, "\t %s %s %i %i",
				  save.chAverageType[ipPun][i],
				  save.chAverageSpeciesLabel[ipPun][i] ,
				save.nAverageIonList[ipPun][i] ,
				  save.nAverage2ndPar[ipPun][i] );
		strcat( chHeader, chTemp );
	}
	strcat( chHeader, "\n");
}

void save_average( 
	/* the file we will write to */
	long int ipPun)
{
	long int i;

	DEBUG_ENTRY( "save_average()" );

		/* do the output */
		for( i=0; i<save.nAverageList[ipPun] ; ++i )
		{
			double result;
			char chWeight[7];
			if( save.nAverage2ndPar[ipPun][i] == 0 )
				strcpy( chWeight , "RADIUS");
			else
				strcpy( chWeight , "VOLUME");

			if( strncmp( save.chAverageType[ipPun][i] , "TEMP" , 4 ) == 0)
			{
				/* temperature */
				if( cdTemp( 
					save.chAverageSpeciesLabel[ipPun][i] ,
					save.nAverageIonList[ipPun][i] ,
					&result , 
					chWeight ) )
				{
					fprintf( ioQQQ, " save average temperature could not identify the species.\n" );
					cdEXIT(EXIT_FAILURE);
				}
				/* will report log of temperature */
				result = log10( result );
			}
			else if( strncmp( save.chAverageType[ipPun][i] , "IONI" , 4 ) == 0)
			{
				/* ionization fraction 
				 * H2 is a special case, HYDRO 0 requests
				 * the H2 fraction, n(H2)/n(H) */
				if( strncmp( "HYDR" , 
					save.chAverageSpeciesLabel[ipPun][i] , 4)==0 &&
					save.nAverageIonList[ipPun][i]== 0 )
					strncpy( save.chAverageSpeciesLabel[ipPun][i],
					"H2  " , 4 );
				if( cdIonFrac( 
					save.chAverageSpeciesLabel[ipPun][i] ,
					save.nAverageIonList[ipPun][i] ,
					&result , 
					chWeight ,
					false 
					) )
				{
					fprintf( ioQQQ, " save average ionization fraction could not identify the species.\n" );
					cdEXIT(EXIT_FAILURE);
				}
				/* will report log of ionization fraction */
				result = log10( result );
			}
			else if( strncmp( save.chAverageType[ipPun][i] , "COLU" , 4 ) == 0)
			{
				/* column density */
				if( cdColm( 
					save.chAverageSpeciesLabel[ipPun][i] ,
					save.nAverageIonList[ipPun][i] ,
					&result ) )
				{
					fprintf( ioQQQ, " save average column density fraction could not identify the species.\n" );
					cdEXIT(EXIT_FAILURE);
				}
				/* will report log of column density */
				result = log10( result );
			}
			else
				TotalInsanity();

			fprintf(save.ipPnunit[ipPun], "\t %.3f", result );
		}
		fprintf(save.ipPnunit[ipPun], "\n");
}
