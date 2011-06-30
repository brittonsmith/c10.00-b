/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSodi compute sodium cooling */
#include "cddefines.h"
#include "coolheavy.h"
#include "taulines.h"
#include "dense.h"
#include "phycon.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolSodi(void)
{
	double a12, 
	  a13, 
	  a23, 
	  cs, 
	  cs13, 
	  cs23, 
	  p3;
	realnum p2;

	DEBUG_ENTRY( "CoolSodi()" );

	/* NaD Na D Na1 Na 1 lines
	 * transition is 3s-3p, ^2P^o - ^2P^o
	 * change 26 Feb 96 cs to Verner's value */
	cs = 2.12e-2*(phycon.te70/phycon.te02);
	PutCS(cs,&TauLines[ipT5895]);
	atom_level2(&TauLines[ipT5895]);

	/* [NaIII] 7.319 microns
	 * cs 
	 * >>refer	na3	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	cs = MIN2(0.40,0.198*phycon.te05*phycon.te01);
	cs = MAX2(0.35,cs);
	PutCS(cs,&TauLines[ipfsNa373]);
	atom_level2(&TauLines[ipfsNa373]);

	/* collision data from 
	 * >>refer	na4	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1
	 * [NaIV] 9.048 microns */
	PutCS(0.802,&TauLines[ipfsNa490]);
	/* [NaIV] 21.29 mic */
	PutCS(0.273,&TauLines[ipfsNa421]);
	/* following from pradhan review */
	PutCS(.111,&TauDummy);
	atom_level3(&TauLines[ipfsNa490],&TauLines[ipfsNa421],&TauDummy);

	/* >>chng 97 mar 19, added NaV lines
	 * NaV lines, collision data from 
	 * >>refer	na5	cs	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143
	 * A's from 
	 * >>refer	na5	as	Kaufman, V., & Sugar, J. 1986, J Phys Chem Ref Dat, 15, 321
	 * >>chng 97 jul 25, cs had been 1, changed to Me */
	cs = 0.919;
	cs13 = 0.359;
	cs23 = 1.915;
	a12 = 8.16e-3;
	a13 = 3.81;
	a23 = 2.68;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	p3 = atom_pop3(4.,10.,6.,cs,cs13,cs23,a12,a13,a23,6.96e4,3.58e4,&p2,
	  dense.xIonDense[ipSODIUM][4],0.,0.,0.);
	CoolHeavy.c1365 = p3*a13*1.46e-11;
	CoolHeavy.c4017 = p3*a23*4.95e-12;
	CoolHeavy.c2067 = p2*a12*9.63e-12;
	CoolAdd("Na 5",1365,CoolHeavy.c1365);
	CoolAdd("Na 5",4017,CoolHeavy.c4017);
	CoolAdd("Na 5",2067,CoolHeavy.c2067);

	/* [Na VI] 14.32 mic, 8.62 mic */
	cs = MIN2(0.77,2.346/(phycon.te10*phycon.te02*phycon.te001));
	PutCS(cs,&TauLines[ipxNa6143]);
	cs = MIN2(2.15,6.934/(phycon.te10*phycon.te03/phycon.te001/
	  phycon.te001));
	PutCS(cs,&TauLines[ipxNa6862]);
	cs = MIN2(0.53,1.518/(phycon.te10*phycon.te01*phycon.te003*
	  phycon.te003));
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipxNa6143],&TauLines[ipxNa6862],&TauDummy);

	/* [Na VI] UV lines, 2971.9, 2872.7 doublet, 2578.9, 1356.6
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	cs = MIN2(0.2876,2.603e-3/(phycon.sqrte/phycon.te10));
	cs = MAX2(0.1,cs);

	CoolHeavy.c2569 = atom_pop3(9.,5.,1.,1.38,0.173,cs,1.68,16.9,5.27,5.01e4,
	  5.60e4,&p2,dense.xIonDense[ipSODIUM][5],0.,0.,0.)*5.27*6.70e-12;

	CoolHeavy.c1357 = CoolHeavy.c2569*(16.9/5.27)*(2972./1357.);
	CoolHeavy.c2972 = p2*1.68*6.70e-12;

	CoolAdd("Na 6",2569,CoolHeavy.c2569);
	CoolAdd("Na 6",1357,CoolHeavy.c1357);
	CoolAdd("Na 6",2972,CoolHeavy.c2972);

	/* [Na VII] 4.675 microns, no CS (interpolated), A NIST
	 * collision strength interpolated from 
	 * >>refer	na7	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	PutCS(1.5,&TauLines[ipxNa0746]);
	atom_level2(&TauLines[ipxNa0746]);
	return;
}
