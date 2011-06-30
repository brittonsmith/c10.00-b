/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*prt_LineLabels save all labels and wavelengths for emission line array */
/*sprt_wl write wavelength to string - must be kept parallel with prt_wl */
/*prt_wl - print floating wavelength in Angstroms, in output format */
#include "cddefines.h"
#include "lines.h"
#include "prt.h"

/*prt_wl print floating wavelength in Angstroms, in output format */
void prt_wl( FILE *ioOUT , realnum wl )
{
	char chString[100];
	DEBUG_ENTRY( "prt_wl()" );

	sprt_wl( chString , wl );

	fprintf(ioOUT, "%s", chString );
	return;
}

/* write wavelength to string */
void sprt_wl( char *chString , realnum wl )
{
	char chUnits[10];

	DEBUG_ENTRY( "sprt_wl()" );

	/* print in A unless > 1e4, then use microns */
	if( wl > 1e8 )
	{
		/* centimeters */
		strcpy( chUnits , "c" );
		wl /= 1e8;
	}
	else if( wl > 1e4 )
	{
		/* microns */
		strcpy( chUnits , "m" );
		wl /= 1e4;
	}
	else if( wl == 0. )
	{
		strcpy( chUnits , " " );
	}
	else
	{
		/* Angstroms units */
		strcpy( chUnits , "A" );
	}

	/* want total of four sig figs */
	if( LineSave.sig_figs == 4 )
	{
		if( wl==0. )
		{
			sprintf(chString, "%5i", 0 );
		}
		else if( wl<10. )
		{
			sprintf(chString, "%5.3f", wl );
		}
		else if( wl<100. )
		{
			sprintf(chString, "%5.2f", wl );
		}
		else if( wl < 1e3 )
		{
			sprintf(chString, "%5.1f", wl );
		}
		else if( wl < 1e4 )
		{
			sprintf(chString, "%5.0f", wl );
		}
		else if( wl < 1e5 )
		{
			sprintf(chString, "%5i", (int)wl );
		}
		else
		{
			TotalInsanity();
		}
	}
	else if( LineSave.sig_figs == 5 )
	{
		/* this branch five sig figs */
		if( wl==0. )
		{
			sprintf(chString, "%5i", 0 );
		}
		else if( wl<10. )
		{
			sprintf(chString, "%5.4f", wl );
		}
		else if( wl<100. )
		{
			sprintf(chString, "%5.3f", wl );
		}
		else if( wl < 1e3 )
		{
			sprintf(chString, "%5.2f", wl );
		}
		else if( wl < 1e4 )
		{
			sprintf(chString, "%5.1f", wl );
		}
		else if( wl < 1e5 )
		{
			sprintf(chString, "%5.0f", wl );
		}
		else if( wl < 1e6 )
		{
			sprintf(chString, "%5i", (int)wl );
		}
		else
		{
			TotalInsanity();
		}
	}
	else
	{
		ASSERT( LineSave.sig_figs == 6 );
		/* this branch five sig figs */
		if( wl==0. )
		{
			sprintf(chString, "%6i", 0 );
		}
		else if( wl<10. )
		{
			sprintf(chString, "%6.5f", wl );
		}
		else if( wl<100. )
		{
			sprintf(chString, "%6.4f", wl );
		}
		else if( wl < 1e3 )
		{
			sprintf(chString, "%6.3f", wl );
		}
		else if( wl < 1e4 )
		{
			sprintf(chString, "%6.2f", wl );
		}
		else if( wl < 1e5 )
		{
			sprintf(chString, "%6.1f", wl );
		}
		else if( wl < 1e6 )
		{
			sprintf(chString, "%6.0f", wl );
		}
		else if( wl < 1e7 )
		{
			sprintf(chString, "%6i", (int)wl );
		}
		else
		{
			TotalInsanity();
		}
	}
	strcat( chString , chUnits );
	return;
}

/*prt_LineLabels save all labels and wavelengths for emission line array */
void prt_LineLabels(
	/* io file handle */
	FILE * ioOUT ,
	/* print all if true, if false then do not print parts of 
	 * transferred lines */
	bool lgPrintAll )
{
	long int i;

	DEBUG_ENTRY( "prt_LineLabels()" );

	for( i=0; i < LineSave.nsum; i++ )
	{
		if( strcmp( LineSv[i].chALab , "####" )==0 )
		{
			/*fprintf( ioOUT, "%s ", LineSv[i].chALab );*/
			fprintf( ioOUT, "####\t%s",LineSave.chHoldComments[(int)LineSv[i].wavelength] ); 
		}
		else
		{
			if( !lgPrintAll &&
				(strcmp( LineSv[i].chALab , "Inwd" )==0 ||
				 strcmp( LineSv[i].chALab , "Coll" )==0 ||
				 strcmp( LineSv[i].chALab , "Pump" )==0 ||
				 strcmp( LineSv[i].chALab , "Heat" )==0)
				)
				/* option to do not print lots of redundant labels 
				 * lgPrintAll is false by default set true with LONG option
				 * on save line labels command */
				continue;
			/* this format chosen to be identical to that used by final */
			fprintf( ioOUT, "%li\t%s\t", 
				i,
				LineSv[i].chALab );
			/* wavelength as given in printout */
			prt_wl( ioOUT, LineSv[i].wavelength );
			/* skip over leading spaces - a formatting problem */
			long int j = 0;
			while( LineSv[i].chComment[j]!='\0' && LineSv[i].chComment[j]==' ')
				++j;
			/* comment entered when line intensity generated */
			fprintf( ioOUT , "\t%s" , &LineSv[i].chComment[j] );
		}
		fprintf( ioOUT, "\n" );
	}
	return;
}
