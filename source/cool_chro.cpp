/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolChro compute chromium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "lines_service.h"
#include "dense.h"
#include "atoms.h"
#include "cooling.h"
#include "phycon.h"


void CoolChro(void)
{
	double a21, 
	  a31, 
	  a32;
	realnum p2, 
	  p3,
	  cs;

	DEBUG_ENTRY( "CoolChro()" );

	/* Cr Chromium cooling
	 *
	 * POPEXC( O12,g1,g2,A21,excit,abund); result already*a21
	 * [Cr III] 5828, multiplet average */
	CoolHeavy.Cr3l21 = atom_pop2(25.,25.,9.,0.05,2.47e4,dense.xIonDense[ipCHROMIUM][2])*
	  3.41e-12;
	CoolAdd("Cr 3",5828,CoolHeavy.Cr3l21);

	/* Cr IV
	 * these are 2 lines estimated by Jim Kingdon
	 * a's are bad, collision strengths just one */
	a21 = 0.053;
	a31 = 0.102;
	a32 = 0.00;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * energies are in kelvin */
	p3 = (realnum)atom_pop3(28.,12.,18.,28.,12.,18.,a21,a31,a32,19795.,1356.,&p2,
	  dense.xIonDense[ipCHROMIUM][3],0.,0.,0.);
	/* multiplet at roughly 6801 A */
	CoolHeavy.Cr4l31 = p3*a31*2.92e-12;
	/* multiplet at roughly 7267 A */
	CoolHeavy.Cr4l21 = p2*a21*2.74e-12;
	CoolAdd("Cr 4",6801,CoolHeavy.Cr4l31);
	CoolAdd("Cr 4",7267,CoolHeavy.Cr4l21);

	/* Cr V
	 * these are 3 lines estimated by Jim Kingdon
	 * a's are bad, collision strengths just one */
	a21 = 0.157;
	a31 = 0.048;
	a32 = 0.016;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * energies are in kelvin */
	p3 = (realnum)atom_pop3(21.,5.,9.,21.,5.,9.,a21,a31,a32,18028.,3842.,&p2,dense.xIonDense[ipCHROMIUM][4],
	  0.,0.,0.);
	/* multiplet at roughly 6577 A */
	CoolHeavy.Cr5l31 = p3*a31*3.02e-12;
	/* multiplet at roughly 7979 A */
	CoolHeavy.Cr5l21 = p2*a21*2.49e-12;
	/* multiplet at roughly 3.74 microns */
	CoolHeavy.Cr5l32 = p2*a32*5.31e-13;
	CoolAdd("Cr 5",6577,CoolHeavy.Cr5l31);
	CoolAdd("Cr 5",37,CoolHeavy.Cr5l32);
	CoolAdd("Cr 5",7979,CoolHeavy.Cr5l21);

	/* [Cr VIII] 1.0106 mic
	 * Y(ik) from 
	 * >>refer	Cr8	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(1.76, &TauLines[ipCr08101]);
	atom_level2( &TauLines[ipCr08101]);

	/* [Cr 16] 1410.60, cs from 
	 * >>referold	Cr16	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/* >>refer	Cr16	cs	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 129, 161 */
	/*>>chng 06 jul 19 Changes made-Humeshkar Nemala*/
	if(phycon.te < 4.057E6)
	{
		cs = (realnum)(0.0195*phycon.te10*phycon.te03*phycon.te0001);
	}
	else
	{
		cs = (realnum)(59.141956/((phycon.te40/phycon.te01)*(phycon.te007/phycon.te0002)));
	}
	/*PutCS(0.11, &TauLines[ipCr16141]);*/
	PutCS(cs, &TauLines[ipCr16141]);
	atom_level2( &TauLines[ipCr16141]);

	return;
}
