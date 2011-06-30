/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolMagn compute magnesium cooling */
#include "cddefines.h"
#include "coolheavy.h"
#include "taulines.h"
#include "phycon.h"
#include "dense.h"
#include "ligbar.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolMagn(void)
{
	realnum cs, 
	  csdum, 
	  csoi;
	double cs2s2p, 
	  cs2s3p;

	DEBUG_ENTRY( "CoolMagn()" );

	/* Mg I 2853 
	 * fit to Dima's  integration of
	 * >>refer	mg1	cs	Leep, D., & Gallagher, A. 1976, Phys Rev A, 13, 148 */
	cs = (realnum)(5.21e-4*phycon.te/phycon.te10);
	PutCS(cs,&TauLines[ipMgI2853]);
	atom_level2(&TauLines[ipMgI2853]);

	/* Mg I 2026,  */
	MakeCS(&TauLines[ipMgI2026]);
	atom_level2(&TauLines[ipMgI2026]);

	/* Mg 1 4571, data from Mendoza
	 * >>refer	mg1	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143
	 * cs set to O I */
	csoi = (realnum)(2.68e-5*phycon.te*(1. + 1.67e-6*phycon.te - 2.95e-10*phycon.te*
	  phycon.te));
	csoi = (realnum)MAX2(0.1,csoi);
	PutCS(csoi/9.,&TauLines[ipT4561]);
	atom_level2(&TauLines[ipT4561]);

	/* Mg II 2798
	 * cs from 
	 * >>refer	mg2	cs	Sigut, A., & Pradhan, A.K., 1994, J Phys B sub
	 * refer not publ as of '97, cs agrees fairly well with Harrington et al.
	 * previous reference (~'82) */
	cs2s2p = 4.50*phycon.te10;
	PutCS(cs2s2p,&TauLines[ipT2796]);
	PutCS(cs2s2p*0.5,&TauLines[ipT2804]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipT2804],&TauDummy,&TauLines[ipT2796]);
	/* call PutCS( cs , t2800 )
	 * call atom_level2( t2800 )
	 *
	 * following used in MAGNES for photo destruction rate */
	if( atoms.PopLevels[0] > 0. )
	{
		atoms.popmg2 = (realnum)((atoms.PopLevels[2] + atoms.PopLevels[1])/
		  atoms.PopLevels[0]);
	}
	else
	{
		atoms.popmg2 = 0.;
	}

	/* MG IV 4.487 MIC 
	 * cs 
	 * >>referold	mg4	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/*cs = (realnum)MIN2(0.425,0.180*phycon.te05*phycon.te02);
	cs = (realnum)MAX2(0.356,cs);*/
	/* >>chng 06 jul 06-Humeshkar Nemala*/
	/*>> refer Mg IV cs Berrington,K.A., Saraph,H. E., & Tully, J.A. 1998,A&AS,129,161 */
	/*This is the cs of the transition between the levels of 2P^o term(J=1/2 - J=3/2)*/
	if(phycon.te < 4E5)
	{
		cs = (realnum)(0.155*phycon.te07*phycon.te01*
			phycon.te003*phycon.te0005*phycon.te0001);
	}
	else
	{
		cs = (realnum)(5.124/((phycon.te20/phycon.te02)*
			phycon.te007*phycon.te0005*phycon.te0001));
	}
	PutCS(cs,&TauLines[ipTMg4]);

	atom_level2(&TauLines[ipTMg4]);

	/* MG V 5.61, 13.54 micron, 
	 * >>refer	mg5	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1
	 * >>chng 96 jul 16 had been constant 0.3 */
	cs = (realnum)MIN2(0.311,0.11*phycon.te10);
	PutCS(cs,&TauLines[ipTMg14]);

	cs = (realnum)MIN2(1.06,0.339*phycon.te10);
	PutCS(cs,&TauLines[ipTMg6]);
	/* >>chng 96 jul 16 had been constant 0.3 */
	cs = (realnum)MIN2(0.297,0.0745*phycon.te10*phycon.te02);
	PutCS(cs,&TauDummy);
	atom_level3(&TauLines[ipTMg6],&TauLines[ipTMg14],&TauDummy);

	/* [Mg V] 2751+2893- Ne III-like, cs 
	 * >>refer	mg5	as	Mendoza, C., & Zeippen, C.J. 1987, MNRAS, 224, 7p
	 * >>chng 96 aug 5 to three level atom
	 * c2751 = atom_pop2(1.33,9.,5.,2.4,5.14e4,xmg(5))*7.11e-12
	 * dCooldT = dCooldT + c2751*5.14e4*tsq1
	 * call CoolAdd( 'Mg 5' , 2751 , C2751 )
	 *
	 * following is 2-1 transition, both 2928 and 2783, 
	 * >>refer	mg5	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1 */
	PutCS(1.187,&TauLines[ipxMg52855]);

	cs = (realnum)MIN2(0.278,0.0171*phycon.te20*phycon.te05/
	  phycon.te005/phycon.te003);
	cs = (realnum)MAX2(0.182,cs);

	/* 3-2 transition, 2417.5 */
	PutCS(cs,&TauLines[ipxMg52417]);

	/* 3-1 transition, 1324.58 */
	PutCS(0.153,&TauLines[ipxMg51325]);

	atom_level3(&TauLines[ipxMg52855],&TauLines[ipxMg52417],&TauLines[ipxMg51325]);

	/* Mg VI, 1806- OII like, data 
	 * >>refer	mg6	all	Kafatos, M., & Lynch, J.P. 1980, ApJS, 42, 611 */
	/* >>refer	mg6	as	Becker, Butler, Zeippen, 1989, A&A 221, 375 
	 * >>refer	mg6	cs	Ramsbottom & Bell 1997, A&AS 125, 543 */
	CoolHeavy.c1806 = atom_pop2(0.6,4.,10.,0.1,7.974e4,dense.xIonDense[ipMAGNESIUM][5])*
	  1.11e-11;
	CoolAdd("Mg 6",1806,CoolHeavy.c1806);

	/* [Mg VII] IR lines at 5.517 and 9.03 microns, 
	 * carbon-like, 
	 * >>refer	mg7	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	if( phycon.alogte < 4.4 )
	{
		cs = (realnum)(0.027*phycon.te30/phycon.te03*phycon.te003*phycon.te001);
	}
	else
	{
		cs = 0.44f;
	}
	PutCS(cs,&TauLines[ipfsMg790]);
	if( phycon.alogte < 4.6 )
	{
		cs = (realnum)(MIN2(1.456,0.0577*phycon.te30*phycon.te02/phycon.te001/
		  phycon.te001));
		csdum = (realnum)(8.275e-3*phycon.sqrte/phycon.te10/phycon.te001);
	}
	else
	{
		cs = (realnum)(3.257/(phycon.te05*phycon.te02*phycon.te003*phycon.te003));
		csdum = (realnum)(1.456/(phycon.te10*phycon.te01*phycon.te005));
	}
	PutCS(cs,&TauLines[ipfsMg755]);
	PutCS(csdum,&TauDummy);
	/* atom_level3(  t10,t21,t20) */
	atom_level3(&TauLines[ipfsMg790],&TauLines[ipfsMg755],&TauDummy);

	/* [mg vii] 2510, 2629
	 * c2596 = atom_pop2( 1.7,9.,5.,10.,5.76e4,dense.xIonDense(12,7))*7.96e-12
	 * dCooldT = dCooldT + c2596 * 5.76e4*tsq1
	 * call CoolAdd( 'Mg 7' , 2596 , C2596 )
	 *
	 * >>chng 96 aug 5, converted to 3 level atom */
	cs = (realnum)MIN2(0.22,0.3622/(phycon.te05*phycon.te02*phycon.te003));
	PutCS(cs,&TauLines[ipxMg71190]);

	/* 2-1 transitions, 2509.2A+2629.1A together */
	cs = (realnum)MIN2(1.067,0.247*phycon.te10*phycon.te03*phycon.te005);
	PutCS(cs,&TauLines[ipxMg72569]);

	/* 3-2 transition, 2261.5A */
	cs = (realnum)MIN2(0.542,3.863/(phycon.te20*phycon.te03*
	  phycon.te01/phycon.te003));
	cs = (realnum)MAX2(0.3735,cs);
	PutCS(cs,&TauLines[ipxMg72261]);

	atom_level3(&TauLines[ipxMg72569],&TauLines[ipxMg72261],&TauLines[ipxMg71190]);
	/* atom_level3(  t10,t21,t20)
	 *
	 * Mg VIII 3.03 micron, data from 
	 * >>refer	mg8	as	Chandra, S. 1982, SoPh, 75, 133
	 * cs from 
	 * >>refer	mg8	cs	Zhang, H.L., Graziani, M., Pradhan, A.K. 1994, A&A, 283, 319 */
	PutCS(1.0,&TauLines[ipxMg08303]);
	atom_level2(&TauLines[ipxMg08303]);
	/* fs303 = atom_pop2(0.26,2.,4.,0.324,4752.,xmg(8))*6.58E-13
	 * call CoolAdd( 'Mg 8' , 303 , FS303 )
	 *
	 * Mg IX 704.5, 1909-like, A from
	 * >>refer	mg9	as	Muhlethaler, H.P., & Nussbaumer, H. 1976, A&A 48, 109
	 * AtomSeqBeryllium line, cs data from 
	 * >>refer	mg9	cs	Keenan, F.P. Berrington, K.A., Burke, P.G., Dufton, P.L.,
	 * >>refercon Kingston, A.E. 1986, PhyS 34, 216
	 * A's 
	 * >>refer	mg9	as	Fleming, J., Bell, K.L, Hibbert, A., Vaeck, N., Godefroid, M.R.
	 * >>refercon 1996, MNRAS, 279, 1289 */
	/** \todo	2	use AtomSeqBeryllium here */
	cs = (realnum)(0.98288 - 0.23766*phycon.alogte + 0.014334*POW2(phycon.alogte));
	cs = (realnum)MAX2(0.01,cs);
	PutCS(cs,&TauLines[ipT705]);
	atom_level2(&TauLines[ipT705]);

	/* Mg X 610
	 * >>refer	mg10	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(12,&TauLines[ipTMg610],&TauLines[ipT58],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipTMg610]);
	PutCS(cs2s2p*0.5,&TauLines[ipTMg625]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipTMg625],&TauDummy,&TauLines[ipTMg610]);

	/* Mg X 58A */
	PutCS(cs2s3p,&TauLines[ipT58]);
	atom_level2(&TauLines[ipT58]);
	return;
}
