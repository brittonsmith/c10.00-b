/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolPhos compute phosphorus cooling */
#include "cddefines.h"
#include "taulines.h"
#include "dense.h"
#include "lines_service.h"
#include "phycon.h"
#include "coolheavy.h"
#include "atoms.h"
#include "cooling.h"

void CoolPhos(void)
{
	double cs, cs21 , cs31 , cs32;
	double a21 , a31 , a32;
	realnum p2, p3;

	DEBUG_ENTRY( "CoolPhos()" );

	/* [P II] 60.64, 32.87 microns
	 * cs, As from
	 * >>refer	p2	as	Mendoza, C., & Zeippen, C.J., 1982, MNRAS 199, 1025
	 * >>refer	p2	cs	Krueger, T.K., and Czyzak, S.J., 1970, Proc Roy Soc London A 318, 531 */
	/** \todo	2	update to Tayal data, email of April 22 2003, must be published */
	PutCS(1.587,&TauLines[ipP0260]);
	PutCS(3.566,&TauLines[ipP0233]);
	PutCS(1.0,&TauDummy);

	/* atom_level3(  t10,t21,t20) */
	atom_level3(&TauLines[ipP0260],&TauLines[ipP0233],&TauDummy);

	/* >>chng 01 may 15, add these three lines, discussed in 
	 * Oliva, E., Marconi, A., et al. A&A 2001, 369, L5 */
	/* the 1D-3P and 1S-1d forbidden lines of [PII] */
	/* >>refer	p2	as	Mendoza, C., & Zeippen, C.J., 1982, MNRAS 199, 1025*/
	a21 = 1.952e-2;
	a31 = 0.2025;
	a32 = 1.64;
	/* these are just a guess */
	cs21 = 1.;
	cs31 = 1.;
	cs32 = 1.;
	p3 = (realnum)(atom_pop3(9.,5.,1.,cs21,cs31,cs32,
	  a21,a31,a32,12534.,7877.9,&p2,dense.xIonDense[ipPHOSPHORUS][1], 0.,0.,0.));
	CoolHeavy.p2_32 = p3*a32*1.21e-12;
	CoolHeavy.p2_31 = p3*a31*4.23e-12;
	CoolHeavy.p2_21 = p2*a21*1.72e-12;
	/* 3-2 1.64 mic */
	CoolAdd("p  2",16400,CoolHeavy.p2_32);
	/* 3-1 4670, 4738 */
	CoolAdd("p  2",4700,CoolHeavy.p2_31);
	/* 2-1 1.147, 1.189 mic */
	CoolAdd("p  2",11600,CoolHeavy.p2_21);

	/* [P III] 17.885 microns
	 * cs, A from
	 * >>refer	p3	as	Kaufman, V., & Sugar, J., 1986, J Phys Chem Ref Data 15, 321
	 * >>refer	p3	cs	Krueger, T.K., and Czyzak, S.J., 1970, Proc Roy Soc London A 318, 531 */
	PutCS(1.859,&TauLines[ipP0318]);
	atom_level2(&TauLines[ipP0318]);

	/* [P VII] 1.374 microns
	 * cs from 
	 * >>referold	p7	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */

	/* >>refer	p7	cs	Berrington,K.A., Saraph, H.E. & Tully, J.A. 1998, A&AS, 129, 161 */
	/*>>chng 06 jul 18 Changes made-Humeshkar Nemala*/
	/*There are two fits to the cs:Above and below 7.77e5*/
	if(phycon.te < 7.77E5)
	{
		cs = (realnum)(0.0986*(phycon.te10/(phycon.te01*phycon.te002)));
	}
	else
	{
		cs = (realnum)(12.2273/((phycon.te30/phycon.te04)*phycon.te007*phycon.te0004));
	}
	/*PutCS(0.27,&TauLines[ipP713]);*/
	PutCS(cs,&TauLines[ipP713]);
	atom_level2(&TauLines[ipP713]);

	/* [P VIII] 4.85, 1.735 microns */
	PutCS(0.30,&TauLines[ipP848]);
	PutCS(0.97,&TauLines[ipP817]);
	PutCS(0.26,&TauDummy);

	/* j levels inverted */
	atom_level3(&TauLines[ipP817],&TauLines[ipP848],&TauDummy);

	/* [P X] 2.708, 1.868 microns */
	cs = MIN2(0.33,0.892/(phycon.te10/phycon.te001));
	PutCS(cs,&TauLines[ipP1027]);

	cs = MIN2(1.082,5.949/(phycon.te20/phycon.te03*phycon.te003));
	PutCS(cs,&TauLines[ipP1018]);

	cs = MIN2(0.33,3.054/(phycon.te20*phycon.te02*phycon.te005));
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipP1027],&TauLines[ipP1018],&TauDummy);
	return;
}
