/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*wcnint initialize stack or warnings, cautions, notes */
/*warnin enter warnings at the end of the calculations into large stack */
/*notein enter a note about calculation into comment array */
/*bangin called by routine comment to enter surprise into comment stack */
/*caunin called by comment to enter caution into comment stack */
#include "cddefines.h"
#include "warnings.h"

void wcnint(void)
{

	DEBUG_ENTRY( "wcnint()" );

	/* this sub is called first, to initialize the variables */
	warnings.nwarn = 0;
	warnings.ncaun = 0;
	warnings.nnote = 0;
	warnings.nbang = 0;
	return;
}

/*warnin enter warnings at the end of the calculations into large stack */
void warnin(char *chLine)
{

	DEBUG_ENTRY( "warnin()" );

	if( warnings.nwarn >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many warnings have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( warnings.chWarnln[warnings.nwarn], chLine  );
	}

	++warnings.nwarn;
	return;
}

/*notein enter a note about calculation into comment array */
void notein(char *chLine)
{

	DEBUG_ENTRY( "notein()" );

	if( warnings.nnote >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many notes have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( warnings.chNoteln[warnings.nnote], chLine );
	}

	++warnings.nnote;
	return;
}

/*bangin called by routine comment to enter surprise into comment stack */
void bangin(char *chLine)
{

	DEBUG_ENTRY( "bangin()" );

	if( warnings.nbang >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many surprises have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( warnings.chBangln[warnings.nbang], chLine );
	}

	++warnings.nbang;
	return;
}

/*caunin called by comment to enter caution into comment stack */
void caunin(char *chLine)
{

	DEBUG_ENTRY( "caunin()" );

	if( warnings.ncaun >= LIMWCN )
	{
		static bool lgFirst=true;
		if( lgFirst )
			fprintf( ioQQQ, 
				" Too many cautions have been entered; increase the value of LIMWCN everywhere in the code.\n" );
		lgFirst = false;
	}
	else
	{
		strcpy( warnings.chCaunln[warnings.ncaun], chLine );
	}

	++warnings.ncaun;
	return;
}
