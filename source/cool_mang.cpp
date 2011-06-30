/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolMang compute manganese cooling */
#include "cddefines.h"
#include "taulines.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"
#include "phycon.h"

void CoolMang(void)
{
	realnum cs;

	DEBUG_ENTRY( "CoolMang()" );

	/* [Mn IX] 7968.5A
	 * Y(ik) from 
	 *  >>refer	mn9	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(2.48,&TauLines[ipxMn0979]);
	atom_level2(&TauLines[ipxMn0979]);

	/* [Mn 17] 1169.59, cs from 
	 *  >>referold	mn17	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/*  >>refer	mn17	cs	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 129, 161 */
	/*>>chng 06 jul 19 Changes made-Humeshkar Nemala*/
	if(phycon.te < 1.151E6)
	{
		cs = (realnum)(0.107);
	}
	else if(phycon.te < 4.58E6)
	{
		cs = (realnum)((6.252E-03)*(phycon.te20*phycon.te003*phycon.te0005));
	}
	else
	{
		cs = (realnum)((81.365)/(phycon.te40*phycon.te01*phycon.te004*phycon.te0002));
	}
	/*PutCS(0.12,&TauLines[ipxMn1712]);*/
	PutCS(cs,&TauLines[ipxMn1712]);
	atom_level2(&TauLines[ipxMn1712]);
	return;
}
