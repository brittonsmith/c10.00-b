/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolCoba compute cobalt cooling */
#include "cddefines.h"
#include "taulines.h"
#include "lines_service.h"
#include "cooling.h"
#include "atoms.h"

void CoolCoba(void)
{

	DEBUG_ENTRY( "CoolCoba()" );

	/* [Co XI] 5168.A
	 * Y(ik) from 
	 *  >>refer	co11	as	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(1.36, &TauLines[ipCo11527]);
	atom_level2( &TauLines[ipCo11527]);
	return;
}
