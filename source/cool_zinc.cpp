/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolZinc compute zinc cooling */
#include "cddefines.h"
#include "taulines.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolZinc(void)
{

	DEBUG_ENTRY( "CoolZinc()" );

	/* zinc iv 3.625 microns, 
	 *  >>refer	Zn 4	wl	Dinerstein, H.L.,  Gaballe, T.R> 2001, ApJ, 562, 515-520 */
	PutCS(1., &TauLines[ipZn04363]);
	atom_level2( &TauLines[ipZn04363]);
	return;
}
