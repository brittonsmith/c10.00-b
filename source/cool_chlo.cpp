/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolChlo compute chlorine cooling */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "dense.h"
#include "phycon.h"
#include "thermal.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolChlo(void)
{
	double a21, 
	  a31, 
	  a32, 
	  a41, 
	  a42, 
	  a43, 
	  a51, 
	  a52, 
	  a53, 
	  a54, 
	  cs, 
	  cs01,
	  cs02,
	  cs12, 
	  cs13, 
	  cs14, 
	  cs15, 
	  cs23, 
	  cs24, 
	  cs25, 
	  cs34, 
	  cs35, 
	  cs45, 
	  p[5], 
	  p3,
	  tused;
	realnum
	  p2, 
	  pop2,
	  rate;
	static double gcl3[5]={4.,4.,6.,2.,4.};
	static double excl3[4]={18053.,66.,11693.,95.};

	DEBUG_ENTRY( "CoolChlo()" );

	/* >>chng >>03 nov 09, add this line, from
	 * >>refer	Cl1	cs	Hollenbach, D. & McKee, C.F. 1989, ApJ, 342, 306 */
	/* rates are said to be ok over range 30 - 3000K */
	tused = MAX2( 30. , phycon.te );
	tused = MIN2( 3000. , phycon.te  );
	tused /= 100.;

	/* HM89 give deexcitation rates, must convert into electron collision strength, as expected
	 * by the code's infrastructure */
	/* electron collision strength */
	rate = (realnum)( 4.7e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*8.3e-10*pow( tused , 0.17 ) * dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		8.3e-10*pow( tused , 0.17 ) * dense.xIonDense[ipHYDROGEN][0] );
	/* possible for atomic hydrogen density to be vary small, causing zero rate coef,
	 * which triggers thrown assert - guard against this */
	LineConvRate2CS( &TauLines[ipCl1_11m]  , SDIV(rate) );
	atom_level2(&TauLines[ipCl1_11m]);

	/* [Cl II] 14.3678 mic, 33.281 mic*/
	/* >>refer	cl2	as	Mendoza, C., & Zeippen, C.J. 1983, MNRAS, 202, 981*/
	/* the following cs were about 2x smaller */
	/* >>referold	cl2	cs	Krueger, T.K., & Czyzak, S.J. 1970, Pro Roy Soc Lond, 318, 531 */
	/* >>chgn 03 feb 24, change to following collision strengths */
	/* >>refer	cl2	cs	Wilson, N.J., & Bell, K.L. 2002, MNRAS, 331, 389 */
	/* order of 3P ground term, 2, 1, 0 from lowest to highest */

	/* this is the 3P J=1 -> 2, 14.4 micron line */
	cs12 = 17.5 / (phycon.te10*phycon.te04 );
	/*PutCS(2.17,&TauLines[ipfsCl214]);*/
	PutCS( cs12 , &TauLines[ipfsCl214] );

	/* this is the 3P J=0 -> 1, 33.3 micron line */
	cs01 = 4.85 / (phycon.te10*phycon.te02 );
	PutCS( cs01 , &TauLines[ipfsCl233] );
	/*PutCS(0.93,&TauLines[ipfsCl233]);*/

	/* the 0 - 2 transition */
	cs02 = 4.51 / (phycon.te10*phycon.te04 );
	PutCS( cs02 , &TauDummy );
	/*PutCS(1.00,&TauDummy);*/

	/* atom_level3(  t10,t21,t20) */
	atom_level3(&TauLines[ipfsCl214],&TauLines[ipfsCl233],&TauDummy);

	/* [Cl II] 8578.7, 9123.6, 6161.8, 3677.9 */
	/* >>refer	Cl2	As	Mendoza, C., & Zeippen, C.J. 1983, MNRAS, 202, 981 */
	/* following numbering in terms of level position, 1, 2, then 3 */
	a21 = 0.133;
	a31 = 1.33;
	a32 = 2.06;
	cs13 = 1.01;
	/* this is 10x what is in the paper, as per comment made in 
	 * >>refer	Cl2	cs	Keenan, F.P., Aller, L.H., Exter, K.M., Hyung, S., & 
	 * >>refercon	Pollacco, D.L. 2003, ApJ, 584, 385 */
	cs23 = 1.49;
	cs12 = 8.389;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	p3 = atom_pop3(9.,5.,1.,cs12,cs13,cs23,a21,a31,a32,1.576e4,23344.,
	  &p2,dense.xIonDense[ipCHLORINE][1],0.,0.,0.);
	/*p3 = atom_pop3(9.,5.,1.,3.86,0.456,1.15,a21,a31,a32,1.576e4,23344.,
	  &p2,dense.xIonDense[ipCHLORINE][1],0.);*/

	/* [Cl II] 8578.7, 9123.6 doublet, both together */
	CoolHeavy.c8579 = p2*a21*2.32e-12;
	CoolAdd("Cl 2",8579,CoolHeavy.c8579);

	/* [Cl II] 6161.8 auroral line */
	CoolHeavy.c6164 = p3*a32*3.23e-12;
	CoolAdd("Cl 2",6164,CoolHeavy.c6164);

	/* [Cl II] 3677.9 */
	CoolHeavy.c3679 = p3*a31*5.41e-12;
	CoolAdd("Cl 2",3679,CoolHeavy.c3679);

	/* [Cl III] this is a [SII] - like doublet, vac lam=5519, 5539
	 * all data from 
	 * >>refer	cl3	all	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143 */
	cs12 = 1.26;
	a21 = 4.83e-3;

	cs13 = 1.88;
	a31 = 7.04e-4;

	cs14 = 0.627;
	a41 = 0.305;

	cs15 = 1.26;
	a51 = 0.754;

	cs23 = 3.19;
	a32 = 3.22e-6;

	cs24 = 1.24;
	a42 = 0.303;

	cs25 = 1.91;
	a52 = 0.323;

	cs34 = 1.38;
	a43 = 0.100;

	cs35 = 3.33;
	a53 = 0.316;

	cs45 = 1.34;
	a54 = 7.65e-6;

	double Cooling , CoolingDeriv;
	atom_pop5(gcl3,excl3,cs12,cs13,cs14,cs15,cs23,cs24,cs25,cs34,cs35,
		cs45,a21,a31,a41,a51,a32,a42,a52,a43,a53,a54,p,
		dense.xIonDense[ipCHLORINE][2]	, &Cooling , &CoolingDeriv, 0.,0.,0.,0.);

	CoolHeavy.Cl5539 = p[1]*a21*3.59e-12;
	CoolHeavy.Cl5519 = p[2]*a31*3.61e-12;
	CoolHeavy.Cl3354 = p[3]*a41*5.93e-12;
	CoolHeavy.Cl3344 = p[4]*a51*5.95e-12;
	CoolHeavy.Cl8504 = p[3]*a42*2.34e-12;
	CoolHeavy.Cl8436 = p[4]*a42*2.36e-12;
	CoolHeavy.Cl8552 = p[3]*a43*2.33e-12;
	CoolHeavy.Cl8483 = p[4]*a53*2.35e-12;

	/* following are whole multiplets */
	CoolHeavy.c5525 = CoolHeavy.Cl5539 + CoolHeavy.Cl5519;
	CoolHeavy.c3350 = CoolHeavy.Cl3354 + CoolHeavy.Cl3344;
	CoolHeavy.c8494 = CoolHeavy.Cl8504 + CoolHeavy.Cl8436 + CoolHeavy.Cl8552 + 
	  CoolHeavy.Cl8483;

	// total cooling from 5-level atom
	CoolAdd("Cl 3",5525,Cooling);
	thermal.dCooldT += CoolingDeriv;

	/* [CL IV], like [OIII]
	 * cs from 
	 * >>refer	cl4	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347 */
	a21 = 0.251;
	cs12 = 6.437;
	a32 = 2.80;
	cs23 = MIN2(2.1,0.0450*phycon.te30*phycon.te03*phycon.te03);
	a31 = 2.50;
	cs13 = 1.922;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	p3 = atom_pop3(9.,5.,1.,cs12,cs13,cs23,a21,a31,a32,2.24e4,3.11e4,&pop2,
	  dense.xIonDense[ipCHLORINE][2],0.,0.,0.);
	/* whole 2-1 transition */
	CoolHeavy.c8047 = pop2*a21*2.48e-12;
	CoolAdd("Cl 4",8047,CoolHeavy.c8047);
	thermal.dCooldT += CoolHeavy.c8047*(2.24e4*thermal.tsq1 - thermal.halfte);
	/* 3-1 transition */
	CoolHeavy.c3119 = p3*a31*6.38e-12;
	CoolAdd("Cl 4",3119,CoolHeavy.c3119);
	/* 3-2 transition */
	CoolHeavy.c5324 = p3*a32*3.74e-12;
	CoolAdd("Cl 4",5324,CoolHeavy.c5324);

	/* [Cl IV] fine structure lines, 20.354, 11.741 microns */
	cs = MIN2(2.7,6.637/(phycon.te10*phycon.te03*phycon.te01));
	cs = MAX2(1.6,cs);
	PutCS(cs,&TauLines[ipCl04203]);

	cs = MIN2(8.0,15.65/phycon.te10);
	PutCS(cs,&TauLines[ipCl04117]);

	cs = MIN2(2.0,5.805/(phycon.te10*phycon.te03));
	PutCS(cs,&TauDummy);

	/* atom_level3(  t10,t21,t20) */
	atom_level3(&TauLines[ipCl04203],&TauLines[ipCl04117],&TauDummy);

	/* fixit - add Cl V 6.71 micron using cs from
	 * >>refer	cl5	cs	Saraph, H.E., & Storey, P.J., 1999, A&AS, 134, 369 */

	/* [Cl IX] 7335 A, 
	 * >>referold	cl9	as	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/* >>refer	cl9	as	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 129, 161 */
	/*>>chng 06 jul 18 Changes made-Humeshkar Nemala*/
	if(phycon.te < 2.03E5)
	{
		cs = (realnum)(0.1175*phycon.te07*phycon.te004*phycon.te0001);
	}
	else if(phycon.te < 5.11E5)
	{
		cs = (realnum)((60.7989E-02)/(phycon.te05*phycon.te01*phycon.te0004));
	}
	else if(phycon.te < 1.284E6)
	{
		cs = (realnum)(0.274857); 
	}
	else
	{
		cs = (realnum)(27.327963/(phycon.te30*phycon.te02*phycon.te007));
	}
	/*PutCS(0.28,&TauLines[ipCl973]);*/
	PutCS(cs,&TauLines[ipCl973]);
	atom_level2(&TauLines[ipCl973]);

	/* [Cl X] 3.05m, 9223A */
	cs = MIN2(0.187,6.686e-3*phycon.te30/phycon.te001);
	cs = MAX2(0.096,cs);
	PutCS(cs,&TauLines[ipCl1030]);

	cs = MIN2(0.568,9.432e-3*phycon.te30*phycon.te03*
	  phycon.te03/phycon.te003/phycon.te001);
	cs = MAX2(0.229,cs);
	PutCS(cs,&TauLines[ipCl1092]);

	cs = MIN2(0.148,7.088e-4*phycon.sqrte/phycon.te03/
	  phycon.te003/phycon.te003);
	cs = MAX2(0.048,cs);
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipCl1092],&TauLines[ipCl1030],&TauDummy);
	return;
}
