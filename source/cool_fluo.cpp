/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolFluo evaluate total cooling due to fluorine */
#include "cddefines.h"
#include "taulines.h"
#include "lines_service.h"
#include "phycon.h"
#include "atoms.h"
#include "cooling.h"

void CoolFluo(void)
{
	double cs;

	DEBUG_ENTRY( "CoolFluo()" );

	/* [F II] 29.33 micron, 67.2 micron
	 * collision strength transition prob
	 * >>refer	f2	cs	Galavis, M.E., et al. 1997, A&AS 123, 159
	 * >>refer	f2	as	Buttler, K., & Zeippen, C.J., 1994, A&AS 108, 1 */
	PutCS(0.60,&TauLines[ipF0229]);
	PutCS(0.206,&TauLines[ipF0267]);
	PutCS(0.160,&TauDummy);

	/* subroutine atom_level3( t10,t21,t20) */
	atom_level3(&TauLines[ipF0229],&TauLines[ipF0267],&TauDummy);

	/* collision strength 
	 * >>refer	f4	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273
	 * [F IV] 44.07 microns */
	cs = MIN2(0.711,0.1245*phycon.te10*phycon.te05*phycon.te01*
	  phycon.te001*phycon.te001);
	PutCS(cs,&TauLines[ipF444]);

	/* [F IV] 25.83 microns */
	cs = MIN2(1.89,0.2023*phycon.te20*phycon.te003*phycon.te003);
	PutCS(cs,&TauLines[ipF425]);
	cs = MIN2(0.451,0.02922*phycon.te20*phycon.te05);
	PutCS(cs,&TauDummy);

	/* subroutine atom_level3( t10,t21,t20) */
	atom_level3(&TauLines[ipF444],&TauLines[ipF425],&TauDummy);
	return;
}
