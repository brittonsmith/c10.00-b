/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolVana compute vanadium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "dense.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"
#include "phycon.h"

void CoolVana(void)
{
	realnum a21, 
	  a31, 
	  a32, 
	  g1, 
	  g2, 
	  g3, 
	  p2, 
	  p3,
	  cs;

	DEBUG_ENTRY( "CoolVana()" );

	/* V Vanadium cooling - element 23
	 *
	 * V III 8823 */
	a21 = 0.05f;
	a31 = 0.10f;
	a32 = 0.00;
	g1 = 28.;
	g2 = 12.;
	g3 = 18.;

	/* [V III] 8823, multiplet average */
	p3 = (realnum)atom_pop3(g1,g2,g3,g1,g2,g3,a21,a31,a32,16303.,606.,&p2,
		dense.xIonDense[ipVANADIUM][2],  0.,0.,0.);

	CoolHeavy.V38830 = p2*a21*2.25e-12;
	CoolHeavy.V38507 = p3*a31*2.34e-12;
	CoolAdd("V  3",8823,CoolHeavy.V38830);
	CoolAdd("V  3",8507,CoolHeavy.V38507);

	/* V IV */
	a21 = 0.054f;
	a31 = 0.039f;
	a32 = 0.007f;
	g1 = 21.;
	g2 = 5.;
	g3 = 9.;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * energies are in kelvin */
	p3 = (realnum)atom_pop3(g1,g2,g3,g1,g2,g3,a21,a31,a32,15159.,3437.,&p2,
		dense.xIonDense[ipVANADIUM][3], 0.,0.,0.);
	/* 7735 ang - 3=>1 */
	CoolHeavy.V47741 = p3*a31*2.57e-12;
	/* 9489 - 2=>1 */
	CoolHeavy.V49496 = p2*a21*2.09e-12;
	/* 4.19 microns 3=>2 */
	CoolHeavy.V44p2m = p3*a32*4.74e-13;
	CoolAdd("V  4",7735,CoolHeavy.V47741);
	CoolAdd("V  4",9489,CoolHeavy.V49496);
	CoolAdd("V  4",42,CoolHeavy.V44p2m);

	/* [V VII] 1.3038 mic
	 * Y(ik) from 
	 * >>refer	v7	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(2.39,&TauLines[ipVa07130]);
	atom_level2(&TauLines[ipVa07130]);

	/* [V 15] 1721.38, cs from 
	 * >>referold	v15	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/* >>refer	v15	cs	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 129, 161 */
	/*>>chng 06 jul 19 Changes made-Humeshkar Nemala*/
	if(phycon.te < 3.566E6)
	{
		cs = (realnum)(0.0149*phycon.te10*phycon.te05*phycon.te004*phycon.te0003);
	}
	else
	{
		cs = (realnum)(47.350653/((phycon.te40/phycon.te02)*phycon.te0002));
	}

	/*PutCS(0.10,&TauLines[ipVa15172]);*/
	PutCS(cs,&TauLines[ipVa15172]);
	atom_level2(&TauLines[ipVa15172]);
	return;
}
