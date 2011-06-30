/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_line parse save lines command, or actually do the save lines output */
/*Save_Line_RT parse the save line rt command - read in a set of lines */
#include "cddefines.h"
#include "cddrive.h"
#include "radius.h"
#include "taulines.h"
#include "opacity.h"
#include "phycon.h"
#include "dense.h"
#include "lines.h"
#include "h2.h"
#include "lines_service.h"
#include "input.h"
#include "prt.h"
#include "save.h"
#include "iso.h"
#include "parser.h"
/* this is the limit to the number of emission lines we can store */
#define	NPUNLM	200L

/* implement the save line xxx command.  cumulative, structure, and
 * emissivity all use same code base and variables, so only one can be used
 * at present */

static char chPLab[NPUNLM][5];
static long int nLinesEntered;
static realnum wavelength[NPUNLM];
static long int ipLine[NPUNLM];
static bool lgRelativeIntensity;

void parse_save_line(Parser &p, 
  /* true, return rel intensity, false, log of luminosity or intensity I */
  bool lgLog3,
  char *chHeader)
{
	char chTemp[INPUT_LINE_LENGTH];

	// save return value of cdLine, 0 for success, -number of lines for fail
	long int i;

	DEBUG_ENTRY( "parse_save_line()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine
	 * following is flag saying whether to do relative intensity or
	 * absolute emissivity */
	lgRelativeIntensity = lgLog3;
	
	/* number of lines we will save */
	nLinesEntered = 0;
	
	/* get the next line, and check for eof */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	/* convert line to caps */
	
	while( p.strcmp("END") != 0 )
	{
		if( nLinesEntered >= NPUNLM )
		{
			fprintf( ioQQQ, 
						" Too many lines have been entered; the limit is %ld.  Increase variable NPUNLM in routine save_line.\n", 
						NPUNLM );
			cdEXIT(EXIT_FAILURE);
		}
		
		/* order on line is label (col 1-4), wavelength */
		strncpy( chPLab[nLinesEntered], p.getCommandRaw(4).c_str() , 4 );
		
		/* null terminate the string*/
		chPLab[nLinesEntered][4] = 0;
		
		/* now get wavelength */
		wavelength[nLinesEntered] = (realnum)p.getWaveOpt();
		
		/* this is total number stored so far */
		++nLinesEntered;
		
		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
	}
	
	sprintf( chHeader, "depth");
	for( i=0; i < nLinesEntered; i++ )
	{
		sprintf( chTemp, "\t%s ", chPLab[i] );
		strcat( chHeader, chTemp );
		sprt_wl(  chTemp, wavelength[i] );
		strcat( chHeader, chTemp );
	}
	strcat( chHeader, "\n" );
}

void save_line(FILE * ioPUN, /* the file we will write to */
  const char *chDo, 
  // intrinsic or emergent line emission?
  bool lgEmergent)
{
	// save return value of cdLine, 0 for success, -number of lines for fail
	long int nCdLineReturn;
	long int i;
	double a[32], 
	  absint, 
	  emiss, 
	  relint;
	double dlum[NPUNLM];

	DEBUG_ENTRY( "save_line()" );

	if( strcmp(chDo,"PUNS") == 0 )
	{
		/* save lines emissivity command */
		static bool lgMustGetLines=true,
			lgBadLine=false;
		lgBadLine = false;
		static bool lgBadH2Line;
		/* it is possible that we will get here after an initial temperature
		 * too high abort, and the line arrays will not have been defined.
		 * do no lines in this case.  must still do save so that there
		 * is not a missing line in the grid save output */
		if( LineSave.nsum >0 )
		{
			lgBadH2Line = false;
			lgBadLine = false;
			/* save lines structure command */
			for( i=0; i < nLinesEntered; i++ )
			{
				if( nzone <= 1 && lgMustGetLines )
				{
					if( (ipLine[i] = cdEmis((char*)chPLab[i],wavelength[i],
						&emiss,lgEmergent)) <= 0 )
					{
						// missed line - ignore if H2 line since large model may be off
						if( !h2.lgH2ON && strncmp( chPLab[i] , "H2  " , 4 )==0 )
						{
							static bool lgMustPrintFirstTime = true;
							if( lgMustPrintFirstTime )
							{
								/* it's an H2 line and H2 is not being done - ignore it */
								fprintf( ioQQQ,"\nPROBLEM Did not find an H2 line, the large model is not "
									"included, so I will ignore it.  Log intensity set to -30.\n" );
								fprintf( ioQQQ,"I will totally ignore any future missed H2 lines\n\n");
								lgMustPrintFirstTime = false;
							}
							/* flag saying to ignore this line */
							ipLine[i] = -2;
							lgBadH2Line = true;
						}
						else
						{
							fprintf( ioQQQ, " save_line could not find line: %s %f\n", 
							  chPLab[i], wavelength[i] );
							ipLine[i] = -1;
							lgBadLine = true;
						}
					}
				}
				/* 0th line is dummy, can't be used, so this is safe */
				ASSERT( ipLine[i] > 0 || lgBadLine || lgBadH2Line || 
					/* this is case where we did not find line on previous iteration,
					 * perhaps because H2 is not turned on, and -1 or -2 was
					 * stored */
					(ipLine[i]<0&&!lgMustGetLines) );
				/* this version of cdEmis uses index, does not search, do not call if line could not be found */
				/* test on case where we abort before first zone is done
				 * this happens in grid when temperature bounds of code
				 * are exceeded.  In this case return small float */
				if( lgAbort && nzone <=1 )
					dlum[i] = SMALLFLOAT;
				else if( ipLine[i]>0 )
					cdEmis_ip(ipLine[i],&dlum[i],lgEmergent);
				else
					dlum[i] = 1e-30f;
			}
			if( lgBadLine )
			{
				cdEXIT(EXIT_FAILURE);
			}
		}
		lgMustGetLines = false;

		/* print depth */
		fprintf( ioPUN, "%.5e", radius.depth_mid_zone );

		/* then print all line emissivity */
		for( i=0; i < nLinesEntered; i++ )
		{
			/*lint -e644 dlum not initialized */
			fprintf( ioPUN, "\t%.4f", log10( MAX2( SMALLFLOAT , dlum[i] ) ) );
			/*lint +e644 dlum not initialized */
		}
		fprintf( ioPUN, "\n" );
	}

	else if( strcmp(chDo,"PUNC") == 0 )
	{
		/* save lines cumulative command */
		fprintf( ioPUN, "%.5e", radius.depth_mid_zone );

		/* it is possible that we will get here after an initial temperature
		 * too high abort, and the line arrays will not have been defined.
		 * do no lines in this case.  must still do save so that there
		 * is not a missing line in the grid save output */
		if( LineSave.nsum >0 )
		{
			for( i=0; i < nLinesEntered; i++ )
			{
				nCdLineReturn = cdLine((char*)chPLab[i],wavelength[i],
					&relint,&absint,lgEmergent);
				if( lgRelativeIntensity )
					/* relative intensity case */
					a[i] = relint;
				else
					/* emissivity or luminosity case */
					a[i] = absint;

				if( nCdLineReturn<=0 )
				{
					/* missed line - ignore if H2 line */
					if( !h2.lgH2ON && strncmp( chPLab[i] , "H2  " , 4 )==0 )
					{
						static bool lgMustPrintFirstTime = true;
						if( lgMustPrintFirstTime )
						{
							/* it's an H2 line and H2 is not being done - ignore it */
							fprintf( ioQQQ,"Did not find an H2 line, the large model is not "
								"included, so I will ignore it.  Log intensity set to -30.\n" );
							fprintf( ioQQQ,"I will totally ignore any future missed H2 lines\n");
							lgMustPrintFirstTime = false;
						}
						/* flag saying to ignore this line */
						a[i] = -30.;
						absint = -30.;
						relint = -30.;
					}
					else
					{
						fprintf( ioQQQ, " save_line could not fine line: %s %f\n", 
							chPLab[i], wavelength[i] );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}

			for( i=0; i < nLinesEntered; i++ )
			{
				fprintf( ioPUN, "\t%.4e", a[i] );
			}
		}
		fprintf( ioPUN, "\n" );
	}

	else
	{
		fprintf( ioQQQ, 
			" unrecognized key for save_line=%4.4s\n", 
		  chDo );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

#define LIMLINE 10
static long int line_RT_type[LIMLINE] = 
  {LONG_MIN , LONG_MIN ,LONG_MIN , LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
	line_RT_ipISO[LIMLINE] =  
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
		line_RT_nelem[LIMLINE] =  
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
			line_RT_ipHi[LIMLINE] =  
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN },
				line_RT_ipLo[LIMLINE] = 
  {LONG_MIN , LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,
	LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN ,LONG_MIN };
static bool lgMustPrintHeader=true;
static long int nLine=-1;

/*Save_Line_RT parse the save line rt command - read in a set of lines */
void Parse_Save_Line_RT(Parser &p)
{
	/* save line rt parameters */
	DEBUG_ENTRY( "Parse_Save_Line_RT()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine */
	
	/* say that we must print the header */
	lgMustPrintHeader = true;
	
	/* get the next line, and check for eof */
	nLine = 0;
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	do
	{
		if(nLine>=LIMLINE )
		{
			fprintf(ioQQQ," PUNCH RT has too many lines - increase LIMLINE in save_line.cpp\n");
			cdEXIT(EXIT_FAILURE);
		}
		
		/* right now it just does lines in the iso sequences */
		line_RT_type[nLine] = (long int)p.FFmtRead();
		line_RT_ipISO[nLine] = (long int)p.FFmtRead();
		line_RT_nelem[nLine] = (long int)p.FFmtRead();
		line_RT_ipHi[nLine] = (long int)p.FFmtRead();
		line_RT_ipLo[nLine] = (long int)p.FFmtRead();
		
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, 
						" there must be five numbers on this line\n");
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}
		
		/* now increment number of lines */
		++nLine;
		
		/* now get next line until we hit eof or the word END */
		p.getline();
	} while( !p.m_lgEOF && !p.nMatch( "END") );
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Save_Line_RT hit end of file looking for END of RT lines\n");
		p.PrintLine(ioQQQ);
		cdEXIT(EXIT_FAILURE);
	}
}

void Save_Line_RT( 
	FILE * ioPUN )
{
	/* save line rt parameters */

	DEBUG_ENTRY( "Save_Line_RT()" );


	static char chLabel[LIMLINE][30];
	long int n;
	if( lgMustPrintHeader )
	{
		fprintf( ioPUN , "Line\tP(con,inc)\tAul\tgl\tgu\n");
		for( n=0; n<nLine; ++n )
		{
			/* print info for header of file, line id and pump rate */
			sprintf( chLabel[n] , "%s ", 
						chLineLbl(&Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]]) );
			fprintf( ioPUN , "%s ", chLabel[n] );
			fprintf( ioPUN , "%.4e ",
						Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Emis->pump);
			fprintf( ioPUN , "%.4e ",
						Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Emis->Aul);
			fprintf( ioPUN , "%.0f ",
						Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Lo->g);
			fprintf( ioPUN , "%.0f ",
						Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Hi->g);
			fprintf( ioPUN , "\n");
			
			if( line_RT_type[n]!=0. )
			{
				/* for now code only exists for H He like iso - assert this */
				fprintf( ioQQQ, 
							" PunchLine_RT only H, He like allowed for now\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		fprintf( ioPUN , "Line\tTauIn\tPopLo\tPopHi\tCul\tk(line)\tk(con,abs)\tk(con,scat)\n");
		lgMustPrintHeader = false;
	}
	
	fprintf(ioPUN, "RADIUS\t%e\tDEPTH\t%e\tTe\t%e\tNe\t%e\n",
			  radius.Radius_mid_zone ,
			  radius.depth_mid_zone ,
			  phycon.te ,
			  dense.eden );
	for( n=0; n<nLine; ++n )
	{
		/* index for line within continuum array */
		long int ipCont = Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].ipCont;
		fprintf( ioPUN , "%s ", chLabel[n] );
		fprintf( ioPUN , "\t%e\t%e\t%e",
					Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Emis->TauIn ,
					StatesElemNEW[line_RT_nelem[n]][line_RT_nelem[n]-line_RT_ipISO[n]][line_RT_ipLo[n]].Pop ,
					StatesElemNEW[line_RT_nelem[n]][line_RT_nelem[n]-line_RT_ipISO[n]][line_RT_ipHi[n]].Pop
			);
		fprintf( ioPUN , "\t%e",
					Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Coll.ColUL 
			);
		
		fprintf( ioPUN , "\t%e\t%e\t%e\n",
					Transitions[line_RT_ipISO[n]][line_RT_nelem[n]][line_RT_ipHi[n]][line_RT_ipLo[n]].Emis->PopOpc ,
					opac.opacity_abs[ipCont-1] ,
					opac.opacity_sct[ipCont-1]
			);
	}
}
 
#	undef LIMELM

