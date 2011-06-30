/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolScan compute scandium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "dense.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolScan(void)
{
	double a21, 
	  a31, 
	  a32, 
	  g1, 
	  g2, 
	  g3, 
	  p3;
	realnum p2;

	DEBUG_ENTRY( "CoolScan()" );

	/* Sc scandium cooling
	 *
	 * Sc II
	 * these are 3 lines estimated by Jim Kingdon
	 * a's are bad, collision strengths just one */
	a21 = 0.004;
	a31 = 6.95e-3;
	a32 = 5.0e-4;
	g1 = 15.;
	g2 = 5.;
	g3 = 21.;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * energies are in kelvin */
	p3 = atom_pop3(g1,g2,g3,g1,g2,g3,a21,a31,a32,3504.,3407.,&p2,
		dense.xIonDense[ipSCANDIUM][1], 0.,0.,0.);

	/* 2.08 microns */
	CoolHeavy.Sc22p08m = p3*a31*9.56e-13;

	/* 4.1 micron */
	CoolHeavy.Sc24p1m = p2*a21*4.85e-13;

	/* 4.22 micron */
	CoolHeavy.Sc24p2m = p3*a32*4.71e-13;

	CoolAdd("Sc 2",41,CoolHeavy.Sc24p1m);
	CoolAdd("Sc 2",42,CoolHeavy.Sc24p2m);
	CoolAdd("Sc 2",21,CoolHeavy.Sc22p08m);

	/* Sc III 3933
	 * POPEXC( O12,g1,g2,A21,excit,abund); result already*a21
	 * [Sc III] 3933, multiplet average */
	g1 = 10.;
	g2 = 2.;
	CoolHeavy.Sc33936 = atom_pop2(g1,g1,g2,0.03,3.66e4,
		dense.xIonDense[ipSCANDIUM][2])* 5.06e-12;
	CoolAdd("Sc 3",3933,CoolHeavy.Sc33936);

	/* [Sc V] 2.31 mic
	 * Y(ik) from 
	 * >>refer	sc5	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(6.00,&TauLines[ipSc05231]);
	atom_level2(&TauLines[ipSc05231]);

	/* Sc VI */
	a21 = 4.94;
	a31 = 49.24;
	a32 = 4.3;
	g1 = 9.;
	g2 = 5.;
	g3 = 1.;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * energies are in kelvin */
	p3 = atom_pop3(g1,g2,g3,g1,g2,g3,a21,a31,a32,28464.,40045.,&p2,
		dense.xIonDense[ipSCANDIUM][5], 0.,0.,0.);
	/* 2100 ang - 3=>1 */
	CoolHeavy.Sc42100 = p3*a31*9.47e-12;
	/* 5054 - 2=>1 */
	CoolHeavy.Sc45058 = p2*a21*3.93e-12;
	/* 3592 3=>2 */
	CoolHeavy.Sc43595 = p3*a32*5.53e-12;
	CoolAdd("Sc 6",2100,CoolHeavy.Sc42100);
	CoolAdd("Sc 6",5054,CoolHeavy.Sc45058);
	CoolAdd("Sc 6",3592,CoolHeavy.Sc43595);

	/* [Sc 13] 2637.97A, cs from 
	 * >>refer	sc13	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	PutCS(0.182,&TauLines[ipSc13264]);
	atom_level2(&TauLines[ipSc13264]);

	return;
}
