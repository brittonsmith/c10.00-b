/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseMap parse map command to produce map of heating and cooling,
 * map is produced by calling punt(" map") */
#include "cddefines.h"
#include "hcmap.h"
#include "parser.h"

void ParseMap(Parser &p )
{
	bool
	  lgLogOn;

	DEBUG_ENTRY( "ParseMap()" );

	/* say output goes to stdout */
	ioMAP = ( ioQQQ == NULL ) ? stdout : ioQQQ;

	/* do cooling space map for specified zones
	 * if no number, or <0, do map and punch out without doing first zone */
	hcmap.MapZone = (long)p.FFmtRead();
	if( p.lgEOL() )
	{
		hcmap.MapZone = 0;
		return;
	}

	if( p.nMatch("RANG") )
	{
		hcmap.RangeMap[0] = (realnum)p.FFmtRead();
		if( hcmap.RangeMap[0] <= 10. )
		{
			hcmap.RangeMap[0] = (realnum)pow((realnum)10.f,hcmap.RangeMap[0]);
			lgLogOn = true;
		}
		else
		{
			lgLogOn = false;
		}
		hcmap.RangeMap[1] = (realnum)p.FFmtRead();
		if( lgLogOn )
			hcmap.RangeMap[1] = (realnum)pow((realnum)10.f,hcmap.RangeMap[1]);

		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be a zone number, followed by two temperatures, on this line.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		return;
	}
	return;
}
