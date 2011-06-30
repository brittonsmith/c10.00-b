/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolAlum compute aluminum cooling */
#include "cddefines.h"
#include "embesq.h"
#include "taulines.h"
#include "phycon.h"
#include "coolheavy.h"
#include "dense.h"
#include "ligbar.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolAlum(void)
{
	double cs ,
	  cs2s2p, 
	 cs2s3p;
	realnum p2;

	DEBUG_ENTRY( "CoolAlum()" );

	/* Al I 3957 */
	MakeCS( &TauLines[ipAlI3957] );
	atom_level2( &TauLines[ipAlI3957]);

	/* Al I 3090 */
	MakeCS(&TauLines[ipAlI3090]);
	atom_level2(&TauLines[ipAlI3090]);

	/* AlII  1670.787
	 * >>chng 96 may 09 put back into level 2 for array processor speed up
	 * cs from 
	 * >>refer	al2	cs	Tayal, S.S., Burke, P.G., Kingston, A.E. 1985, J.Phys. B, 18, 4321
	 * >>refer	al2	cs	Tayal, S.S., Burke, P.G., Kingston, A.E. 1984, J.Phys. B, 17, 3847
	 * cs = MIN( 5.0 , 0.0125 * sqrte*te10*te003 )
	 * call PutCS( cs , al1671 )
	 * call atom_level2( al1671 )
	 *
	 * Al 2, 3P in 
	 * >>refer	al2	cs	Keenan, F.P., Harra, L.K., Aggarwal, K.M., Feibelman, W.A. 1992,
	 * >>refercon	ApJ, 385, 375
	 * doublet at 2660, 2669 */
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	PutCS(3.56,&TauLines[ipT2670]);

	/* C2670 = AtomSeqBeryllium( 1.67,2.00,6.54, T2670 , 3.67E-3 )  * 7.45E-12 */
	AtomSeqBeryllium(1.67,2.00,6.54,&TauLines[ipT2670],3.67e-3);
	embesq.em2669 = (realnum)(atoms.PopLevels[3]*3.67e-3*7.45e-12);

	/* Aluminum al 3, 1854, 1862 doublet, 3s ^2 S gnd, ^2P^o 1/2 3/5 exc
	 * f=0.854 from 
	 * >>refer	al3	as	Dufton, P.L., Brown, P.J.F., Lennon, D.J., Lynas-Gray, A.E. 1986,
	 * >>refercon	MNRAS, 222, 713
	 * cs from 
	 * >>refer	al3	cs	Dufton, P.L., & Kingston, A.E. 1987, J.Phys. B, 20, 3899 */
	cs = 4.407*phycon.te10*phycon.te03*phycon.te01;
	cs = MIN2(25.0,cs);
	PutCS(cs*0.667,&TauLines[ipT1855]);
	PutCS(cs*0.333,&TauLines[ipT1863]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipT1863],&TauDummy,&TauLines[ipT1855]);

	/* Al V 2.905mm, cs 
	 * >>referold	al5	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/*cs = MIN2(0.524,1.113/(phycon.te10/phycon.te02/phycon.te003));*/
	/* >>refer	al5	cs	Berrington, K.A., Saraph, H.E. & Tully, J.A. 1998, A&AS, 129,161 */
	/* >>chng 06 jul 11 - Humeshkar Nemala*/
	if(phycon.te < 1.58E5)
	{
		cs = (realnum)(0.893/(phycon.te05*phycon.te005*phycon.te001*phycon.te0002));
	}
	else
	{
		cs = (realnum)(3.1991/((phycon.te20/phycon.te04)*(phycon.te003/phycon.te0002)));
	}

	PutCS(cs,&TauLines[ipAl529]);

	atom_level2(&TauLines[ipAl529]);

	/* Al VI 3.66, 9.12 microns */
	cs = 639.1/(phycon.sqrte*pow(phycon.te03,phycon.te003)*phycon.te001);
	cs = MIN2(5.5 , cs);
	PutCS(cs,&TauLines[ipAl6366]);

	cs = MIN2(1.10,49.37/(phycon.sqrte/phycon.te10*phycon.te02/
	  phycon.te001));
	PutCS(cs,&TauLines[ipAl6912]);

	cs = MIN2(2.0,319.11/(phycon.sqrte*phycon.te10/phycon.te02/
	  phycon.te001));
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipAl6366],&TauLines[ipAl6912],&TauDummy);

	/* [Al VI] 2428.4, 2601, 1169.8, 2124.9
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * cs from 
	 * >>refer	al6	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1
	 */
	CoolHeavy.c1170 = atom_pop3(9.,5.,1.,1.044,0.145,0.463,6.63,72.9,7.79,
	  5.92e4,6.767e4,&p2,dense.xIonDense[ipALUMINIUM][5],0.,0.,0.)*72.9*1.70e-11;

	CoolHeavy.c2428 = p2*6.63*8.19e-12;
	CoolHeavy.c2125 = CoolHeavy.c1170*(7.79/72.9)*(1169.5/2124.9);
	CoolAdd("Al 6",1170,CoolHeavy.c1170);
	CoolAdd("Al 6",2428,CoolHeavy.c2428);
	CoolAdd("Al 6",2125,CoolHeavy.c2125);

	/* Al VIII 5.85, 3.72 microns
	 * collision strength 
	 * >>refer	al8	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	cs = MIN2(0.39,0.0459*phycon.te20/phycon.te003/phycon.te003);
	PutCS(cs,&TauLines[ipAl8575]);
	cs = MIN2(1.062,0.0407*phycon.te30/phycon.te003/phycon.te003);
	PutCS(cs,&TauLines[ipAl8370]);
	cs = MIN2(0.27,2.694e-3*phycon.te20*phycon.te20*phycon.te01*
	  phycon.te003);
	PutCS(cs,&TauDummy);
	atom_level3(&TauLines[ipAl8575],&TauLines[ipAl8370],&TauDummy);

	/* [Al IX] 2.04 micron, no collision strength, A NIST */
	PutCS(1.,&TauLines[ipAl09204]);
	atom_level2(&TauLines[ipAl09204]);

	/* Al 10 639, CS 
	 * >>refer	al10	cs	Keenan, F.P. Berrington, K.A., Burke, P.G., Dufton, P.L.,
	 * >>refercon	Kingston, A.E. 1986, PhyS 34, 216
	 * A is extrapolation along iso seq */
	cs = 0.73492 - 0.16964*phycon.alogte + 0.0096631*POW2(phycon.alogte);
	cs = MAX2(0.01,cs);
	PutCS(cs,&TauLines[ipT639]);
	atom_level2(&TauLines[ipT639]);

	/* Al 11 Li seq 2s2p 556
	 * >>refer	al11	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(13,&TauLines[ipTAl550],&TauLines[ipTAl48],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipTAl550]);
	PutCS(cs2s2p*0.5,&TauLines[ipTAl568]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipTAl568],&TauDummy,&TauLines[ipTAl550]);

	PutCS(cs2s3p,&TauLines[ipTAl48]);
	atom_level2(&TauLines[ipTAl48]);
	return;
}
