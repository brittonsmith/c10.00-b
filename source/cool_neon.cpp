/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolNeon evaluate total cooling due to neon */
#include "cddefines.h"
#include "coolheavy.h"
#include "taulines.h"
#include "dense.h"
#include "phycon.h"
#include "embesq.h"
#include "ligbar.h"
#include "thermal.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolNeon(void)
{
	double a21, 
	  a31, 
	  a32,
	  cs2s2p, 
	  cs2s3p,
	  /*Introduce the following three variables for use in Ne V-Humeshkar Nemala*/
	  cs12,
	  cs13,
	  cs23; 
	realnum 
	  cs, 
	  pop2;

	DEBUG_ENTRY( "CoolNeon()" );

	/* Neon II 12.8 micron
	 * >>referold	ne2	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/*These statements were commented out-Humeshkar Nemala*/
	/*cs = (realnum)MIN2(0.4,0.1294*phycon.te10/phycon.te02*phycon.te005);*/
	/*cs = (realnum)MAX2(0.272,cs);*/
	/*>>chng 06 jun 30 Changes made by Humeshkar Nemala-Jun 2006*/
	/*>>refer	ne2	cs Griffin,D.C., Mitnik,D.M. & Badnell, N. R.,2001,JPhB,34,4401*/
	/*Temperature range from 1000 TO 4E5 */

	/*>>chng 06 jul 02, adjust precoef to get exact result at 1e4 K */
	cs = (realnum)(0.132*phycon.te07*phycon.te02 *phycon.te004 *phycon.te0001 );

	PutCS(cs,&TauLines[ipTNe13]);
	atom_level2(&TauLines[ipTNe13]);

	/***********************************************************************************
	**************************************Ne III*****************************************
	************************************************************************************/
	double neiii_cs3P13P0,
		neiii_cs3P23P1,
		neiii_cs3P23P0,
		neiii_cs3P1D2,
		neiii_cs3P1S0;

	//neiii_cs calculates the collision strenghts for Ne III
	neiii_cs(neiii_cs3P13P0,
		neiii_cs3P23P1,
		neiii_cs3P23P0,
		neiii_cs3P1D2,
		neiii_cs3P1S0);

	PutCS(neiii_cs3P13P0,&TauLines[ipTNe36]);
	
	PutCS(neiii_cs3P23P1,&TauLines[ipTNe16]);
	
	PutCS(neiii_cs3P23P0,&TauDummy);

	/* now do the level populations */
	atom_level3(&TauLines[ipTNe16],&TauLines[ipTNe36],&TauDummy);
	
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	>>refer	Ne3	as	Froese Fischer, C., & Tachiev, G. 2004, At. Data Nucl. Data Tables, 87, 1
	 */
	a21=0.231;
	a31=2.069;
	a32=2.545;
	CoolHeavy.c3343 = atom_pop3(9.,5.,1.,neiii_cs3P1D2,neiii_cs3P1S0,0.32,a21,a31,a32,
	  3.583e4,4.301e4,&pop2,dense.xIonDense[ipNEON][2],0.,0.,0.)*a32*5.954e-12;
	CoolHeavy.c3869 = pop2*a21*5.14e-12;
	thermal.dCooldT += CoolHeavy.c3869*(3.72e4*thermal.tsq1 - thermal.halfte);
	CoolAdd("Ne 3",3342,CoolHeavy.c3343);
	CoolAdd("Ne 3",3869,CoolHeavy.c3869);
	CoolAdd("Ne 3",1793,CoolHeavy.c3343*1.38);

	/***********************************************************************************
	**************************************Ne IV*****************************************
	************************************************************************************/
	/* Ne IV 2425.4+2422.8  A'S from 
	 * >>refer	ne4	as	Zeippen, C.J. 1982, MNRAS 198 111
	 * Ne IV CS from 
	 * >>referold	ne4	cs	Giles, K. 1981, MNRAS, 195, 63
	 * above gave table, actually used, but following is most recent
	 * calculation, in great agreement, but only gives figures
	 * >>refer	ne4	cs	Ramsbottom, C.A., Bell, K.L., & Keenan, F.P. 1998, MNRAS, 293, 233 
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	a21 = 2.46e-3;
	a31 = 1.02;
	a32 = 0.693;
	CoolHeavy.c4720 = atom_pop3(4.,10.,6.,1.37,0.464,2.14,a21,a31,a32,5.94e4,
	  3.05e4,&pop2,dense.xIonDense[ipNEON][3],0.,0.,0.)*a32*4.22e-12;
	CoolHeavy.c2424 = pop2*8.21e-12*a21;
	CoolAdd("Ne 4",2424,CoolHeavy.c2424);
	CoolAdd("Ne 4",4720,CoolHeavy.c4720);
	CoolAdd("Ne 4",1602,CoolHeavy.c4720*4.1);

	/* Ne V 3426, CS data from 
	 * >>referold	ne5	cs	Lennon, D.J. & Burke, V.M. 1991, MNRAS 251, 628
	 * revised from 
	 * >>refer	ne5	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273
	 * A's from 
	 * >>refer	ne5	as	Baluja, K.L. 1985, J.Phys. B, 18, L413
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	/*>>chng 06 jun 30 The collision strengths are changed-Humeshkar Nemala*/
	/*Temperature dependence for the cs is included*/
	/*>>refer ne5 cs Griffin,D.C., & Badnell,N.R.2000,JPhB,33,4389*/
	/*The cs for the 1-2 transition was obtained by summing the cs
	for transitions 3P0-1D2,3P1-1D2,3P2-1D2*/
	cs12=1.2172*(phycon.te04*phycon.te003*phycon.te0004);
	/*The cs for the 1-3 i.e., the 3P-1S0 transition was obtained by summing the cs 
	for the 3P0-1S0,3P1-1S0,3P2-1S0transitions */
	cs13=1.2598/(phycon.te10*phycon.te04*phycon.te007*phycon.te002*phycon.te0004);
	/*Temperature dependent data avilable for the cs between states 2 and 3,i.e.,1D2-1S0
	on the internet at the Oak Ridge National Laboratory CFADC.
	Fit was done for data obtained through private communication*/
	cs23 = 0.2524*(phycon.te07*phycon.te01*phycon.te007*phycon.te0003);
	/*CoolHeavy.c2975 = atom_pop3(9.,5.,1.,2.18,0.254,0.688,0.521,4.35,2.76,
	  4.297e4,4.835e4,&pop2,dense.xIonDense[ipNEON][4],0.,0.,0.)*2.76*6.69e-12;*/
	CoolHeavy.c2975 = atom_pop3(9.,5.,1.,cs12,cs13,cs23,0.521,4.35,2.76,
	  4.297e4,4.835e4,&pop2,dense.xIonDense[ipNEON][4],0.,0.,0.)*2.76*6.69e-12;
	/* following are old values
	 * C2975 = POP3( 9.,5.,1., 1.8,0.25,0.52, 0.521,4.34,2.76,
	 *  1 4.297E4,4.835E4, POP2 , ANEON(5),0.) * 2.76*6.69E-12 */
	CoolHeavy.c1565 = CoolHeavy.c2975*1.901*1.572;
	CoolHeavy.c3426 = pop2*0.521*5.81e-12;
	thermal.dCooldT += CoolHeavy.c3426*(4.20e4*thermal.tsq1 - thermal.halfte) + 
	  (CoolHeavy.c2975 + CoolHeavy.c1565)*9.132e4*thermal.tsq1;
	CoolAdd("Ne 5",2975,CoolHeavy.c2975);
	CoolAdd("Ne 5",1565,CoolHeavy.c1565);
	CoolAdd("Ne 5",3426,CoolHeavy.c3426);

	/* Ne V 24.2, 14.3 micron
	 * CS from 
	 * >>referold	ne5	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273
	 * A's from 
	 * >>refer	ne5	as	Baluja, K.L. 1985, J.Phys. B, 18, L413 */
	/*cs = (realnum)MIN2(1.84,21.12/(phycon.te10*phycon.te10*phycon.te10/
	  phycon.te003/phycon.te003));*/
	/*>>chng 06 jun 30 Changed-Humeshkar Nemala*/
	/* >>refer	ne5	cs Griffin,D.C., & Badnell,N.R.2000,JPhB,33,4389*/
	/*this line corresponds to the 3P1-3P0 transition - the lowest of the 3P transitions */
	cs = (realnum)(21.917/(phycon.te20 *phycon.te07 *phycon.te02 *phycon.te001 
		*phycon.te0007 *phycon.te0002));
	PutCS(cs,&TauLines[ipTNe24]);
	/*cs = (realnum)MIN2(9.5,261.71/(phycon.te10*phycon.te10*phycon.te10*
	  phycon.te10*phycon.te01*phycon.te003));*/
	/*this line corresponds to the 3P2-3P1 transition*/
	cs = (realnum)(122.49/(phycon.te30 *phycon.te04 *phycon.te005 
		*phycon.te003 *phycon.te0001 ));
	PutCS(cs,&TauLines[ipTNe14]);
	/*we introduce the dummy line 3P2-3P0 transition*/
	/*cs = (realnum)MIN2(3.2,139.86/(phycon.sqrte/phycon.te03*phycon.te003/
	  phycon.te001));*/
	cs = (realnum)(58.788/(phycon.te40*phycon.te001*phycon.te0007));
	PutCS(cs,&TauDummy);
	/* now do the level populations */
	atom_level3(&TauLines[ipTNe24],&TauLines[ipTNe14],&TauDummy);

	/* Ne V 5S - 3P, CS 
	 * >>referold	ne5	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 
	 * A from 
	 * >>referold	ne5	cs	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>referoldcon ed by D.R. Flower, (D. Reidel: Holland), 143 */
	/** \todo	2	transfer these lines */
	/*>>chng 06 july 05 Changed-Humeshkar Nemala*/
	/* >>refer	ne5	cs Griffin,D.C., & Badnell,N.R.2000,JPhB,33,4389*/
	/*this line corresponds to the 2s2 2p2 3P- 2s2p3 5S0 transition*/
	/*The cs was obtained by summing the cs of the transitions
	3P0-5S2,3P1-5S2,3P2-5S2*/
	cs = (realnum)(22.956/((phycon.te30/phycon.te01)*phycon.te001*phycon.te0002));
	/*CoolHeavy.c1134 = atom_pop2(11.9/(phycon.te10*phycon.te10*phycon.te03),
	  9.,5.,4.67e3,1.273e5,dense.xIonDense[ipNEON][4])*1.767e-11;*/
	CoolHeavy.c1134 = atom_pop2(cs,9.,5.,4.67e3,1.273e5,
		dense.xIonDense[ipNEON][4])*1.767e-11;
	CoolAdd("Ne 5",1134,CoolHeavy.c1134);

	/* Ne VI 7.6 micron, A from 
	 * >>refer	ne6	as	Froese Fischer, C. 1983, J.Phys. B, 16, 157
	 * cs from 
	 * >>referold	ne6	cs	Zhang, H.L., Graziani, M., Pradhan, A.K. 1994, A&A, 283, 319
	 * >>chng 96 jul 16 had been just constant 2.0  */
	/*cs = (realnum)MIN2(3.71,23.623/(phycon.te20*phycon.te02/phycon.te003));*/
	/*>>chng 06 jun 30  Changed-Humeshkar Nemala*/
	/*>>refer	ne6	cs	Mitnik,D.M.,Griffin,D.C., & Badnell,N.R. 2001,JPhB,34,4455*/
	cs = (realnum) (35.705/(phycon.te20 *phycon.te07 *phycon.te004 *phycon.te0007 ));
	PutCS(cs,&TauLines[ipxNe0676]);
	atom_level2(&TauLines[ipxNe0676]);
	/* fs76 = atom_pop2(0.37,2.,4.,1.9E-2,1.89E3,dense.xIonDense(10,6))*2.62E-13
	 * dCooldT = dCooldT + fs76*(1890.*tsq1-halfte)
	 * call CoolAdd( 'Ne 6' , 76 , FS76 )
	 *
	 * Ne VII col data from 
	 * >>refer	ne7	cs	Berrington, K.A., Burke, P.G., Dufton, P.L., Kingston, A.E. 1985,
	 * >>refercon At. Data Nucl. Data Tables, 33, 195
	 * low te from 
	 * >>refer	ne7	cs	Dufton, P.L., Doyle, J.G., Kingston, A.E. 1979, A&A, 78, 318
	 * newer fit to 
	 * >>refer	ne7	cs	Ramsbottom, C.A., Berrington, K.A., Bell, K.L. 1995,
	 * >>refercon At. Data Nucl. Data Tables, 61, 105 */
	if( phycon.te < 4e4 )
	{
		cs = (realnum)(0.0352*(phycon.te20/phycon.te03));
	}
	else
	{
		cs = (realnum)(0.736/(phycon.te10*phycon.te02/phycon.te003));
	}
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	PutCS(cs,&TauLines[ipT895]);
	/* AtomSeqBeryllium(CS23,CS24,CS34,tarray,A41)
	 * c895 = AtomSeqBeryllium(.52,.61, 2.0,t895,.0578) * 2.223E-11
	 * A's 
	 * >>refer	ne7	as	Fleming, J., Bell, K.L, Hibbert, A., Vaeck, N., Godefroid, M.R.
	 * >>refercon 1996, MNRAS, 279, 1289 */
	AtomSeqBeryllium(.52,.61,2.0,&TauLines[ipT895],.07066);
	embesq.em895 = (realnum)(atoms.PopLevels[3]*0.0578*2.223e-11);

	/* Ne VIII 774, iso with 1549, extrapolation for omega
	 * >>refer	ne8	??	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(10,&TauLines[ipT770],&TauLines[ipT88],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipT770]);
	PutCS(cs2s2p*0.5,&TauLines[ipT780]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipT780],&TauDummy,&TauLines[ipT770]);

	PutCS(cs2s3p,&TauLines[ipT88]);
	atom_level2(&TauLines[ipT88]);
	return;
}

//neiii_cs calculates the collision strenghts for Ne III
void neiii_cs(double& neiii_cs3P13P0,
		double& neiii_cs3P23P1,
		double& neiii_cs3P23P0,
		double& neiii_cs3P1D2,
		double& neiii_cs3P1S0)
{
    DEBUG_ENTRY( "neiii_cs()" );
    /* Ne III fine structure lines
	 * >>referold	ne3	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1 */
	/*PutCS(0.774,&TauLines[ipTNe16]);
	PutCS(0.244,&TauLines[ipTNe36]);
	PutCS(0.208,&TauDummy);*/
	/* Ne III fine structure lines
	 *>>refer ne3 cs McLaughlin,B.M., & Bell,K.L.2000,JPhB,33,597 */
	/*>>chng 06 jun 30 Changes made by-Humeshkar Nemala*/
	/*Data available over temps 1E3 K to 1E6 K*/
	/* this is the highest of the two transitions, J = 0-1, 36 mm
	 * (3P J levels inverted for this ion) */
	if(phycon.te <6.3E3)
	{
		neiii_cs3P13P0=(realnum)((9.34E-02)*phycon.te10*phycon.te003*phycon.te0002);
    }
	else if(phycon.te < 2.5E4)
	{
		neiii_cs3P13P0=(realnum)((19.8888E-02)*(phycon.te02/(phycon.te003*phycon.te0002)));
	}
	else if(phycon.te < 4E4)
	{
		neiii_cs3P13P0=(realnum)((838.0688E-06)*phycon.sqrte*phycon.te05*(phycon.te007/phycon.te0001));
	}
	else if(phycon.te < 1E5)
	{
		neiii_cs3P13P0=(realnum)((256.2312E-11)*phycon.te*phycon.te70*phycon.te05*phycon.te005*phycon.te0002);
	}
	else if(phycon.te < 2.5E5 )
	{
		neiii_cs3P13P0=(realnum)((238.5789E-06)*phycon.te70*phycon.te05*phycon.te01*
			phycon.te001*phycon.te0004);
	}
	else
	{
		neiii_cs3P13P0=(realnum)(147.59848/(phycon.te30*phycon.te01*phycon.te001*phycon.te0005));
	}

    /* this is the lowest of the two transitions, J = 1-2 16 mm
	 * (3P J levels inverted for this ion) */

	if(phycon.te < 2.5E4)
	{
		neiii_cs3P23P1 = (realnum)(0.3702*phycon.te07*phycon.te005*phycon.te0004);
	}
	else if(phycon.te < 4E4)
	{
		neiii_cs3P23P1 = (realnum)((16.6945E-04)*(phycon.te40*phycon.te20*phycon.te005*phycon.te003*
			phycon.te0005*phycon.te0003));
	}
	else if(phycon.te < 1.6E5)
	{
		neiii_cs3P23P1 = (realnum)((50.5069E-09)*(phycon.te32*phycon.te07*phycon.te02*phycon.te0007
			*phycon.te0001));
	}
	else if(phycon.te < 2.5E5)
	{
		neiii_cs3P23P1 = (realnum)((778.1245E-04)*phycon.te40*(phycon.te002/phycon.te0002));
	}
	else
	{
		neiii_cs3P23P1 = (realnum)(786.6482/(phycon.te30*phycon.te04*phycon.te0001));
	}

	/* this is the transition from highest to lowest J within 3P - it has a tiny A
	 * which we ignore, so use the TauDummy struc */

	if(phycon.te <4E4)
	{
		neiii_cs3P23P0 = (realnum)(0.0999*phycon.te07*phycon.te005*phycon.te001);
	}
	else if( phycon.te < 2.5E5)
	{
		neiii_cs3P23P0=(realnum)(9.02142E-06*
			phycon.te90*phycon.te05*phycon.te004*phycon.te0007*phycon.te0001);
	}
	else
	{
		neiii_cs3P23P0=(realnum)(66.1264/(phycon.te30*phycon.te01*phycon.te007));
	}
	/*cs = 0.207f;*/

	/* Ne III 3869+3968, 3343, A's
	 * >>referold	ne3	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143
	 * CS
	 * >>referold	ne3	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1 */
	/*>>chng 06 jun 30 Changes made-Humeshkar Nemala*/
	/* >>refer	ne3	cs	McLaughlin, B. M., & Bell, K. L. 2000, Journal of Physics B Atomic Molecular Physics, 33, 597
	* (3P J levels inverted for this ion) */
	/*This is the transition between states 1-2 i.e., 3P-1D2*/
	/*The cs for 1-2 transition was obtained by summing the cs of
	transitions 3P2-1D2,3P1-1D2,3P0-1D2*/
	/*The cs are determined over two ranges: as below and above 2.5E5
	 *>>chng 10 feb 24 ML: Updated the coefficients for the cs12 and cs13 values using the McLaughlin and Bell reference.
	 *The coefficients have been chosen to ensure continuity over 2.5E5 K.*/
	if(phycon.te <2.5E5)
	{
		neiii_cs3P1D2 = (realnum)(0.919 *phycon.te03*phycon.te007*
			phycon.te002*phycon.te0007*phycon.te0001);
	}
	else
	{
		neiii_cs3P1D2 = (realnum)((19.7)/(phycon.te20*phycon.te005*
			phycon.te001*phycon.te0007*phycon.te0002));
	}
	/*This is the transition between states 1-3 i.e., 3P-1S0*/
	/*The cs for 1-3 transition was obtained by summing the cs of
	transitions 3P2-1S0,3P1-1S0,3P0-1S0*/
	/*The cs are determined over three ranges: as below 2.5E4,
	above 2.5E4 and below 2.5E5, and above 2.5E5*/
	if(phycon.te < 2.5E4)
	{
		/*The cs remains fairly constant in the range 1E3 to 2.5E4,with an average of 0.152*/
		neiii_cs3P1S0 =  0.152f;
	}
	else if( phycon.te < 2.5E5)
	{
		neiii_cs3P1S0 = (realnum)((798.0776E-05)*(phycon.te30/phycon.te01)*phycon.te001);
	}
	else
	{
		neiii_cs3P1S0 = (realnum)((1026.2621E-02)/(phycon.te20*phycon.te07*phycon.te01*phycon.te005));
	}
	/*cs of 3343 which refers to the 1D2-1S0 cs*/
	/*This is the transition between states 2-3*/
	/*data over temperatures 1E3 to 1E6 remains fairly constant at 0.32
	cs = 0.32f;*/

    return;
}
