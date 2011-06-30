/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolTita compute titanium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "dense.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolTita(void)
{
	realnum a21, 
	  a31, 
	  a32, 
	  p2, 
	  p3;

	DEBUG_ENTRY( "CoolTita()" );

	/* Ti Titanium cooling
	 *
	 * these are 3 lines estimated by Jim Kingdon
	 * a's are bad, collision strengths just one */
	a21 = 0.015f;
	a31 = 0.032f;
	a32 = 0.002f;

	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * energies are in kelvin */
	p3 = (realnum)atom_pop3(21.,5.,9.,21.,5.0,9.0,a21,a31,a32,11844.,3148.4,&p2,
	  dense.xIonDense[ipTITANIUM][2],0.,0.,0.);

	/* multiplet at roughly 9594 */
	CoolHeavy.Ti3l31 = p3*a31*2.07e-12;
	/* multiplet at roughly 4.57 microns */

	CoolHeavy.Ti3l32 = p3*a32*4.35e-13;
	/* multiplet at roughly 1.21 microns */

	CoolHeavy.Ti3l21 = p2*a21*1.64e-12;
	CoolAdd("Ti 3",9594,CoolHeavy.Ti3l31);
	CoolAdd("Ti 3",4,CoolHeavy.Ti3l32);
	CoolAdd("Ti 3",1,CoolHeavy.Ti3l21);

	/* [Ti VI] 1.7150 mic
	 * Y(ik) from 
	 * >>refer	ti6	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	PutCS(3.48,&TauLines[ipTi06172]);
	atom_level2(&TauLines[ipTi06172]);

	/* [Ti XIV] 2117.79, cs from 
	 * >>refer	ti14	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	PutCS(0.23,&TauLines[ipTi14212]);
	atom_level2(&TauLines[ipTi14212]);
	return;
}
