/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolCalc compute calcium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "doppvel.h"
#include "phycon.h"
#include "ca.h"
#include "dense.h"
#include "thermal.h"
#include "opacity.h"
#include "rfield.h"
#include "ligbar.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolCalc(void)
{
	realnum p2;
	double  a21, 
	  a31, 
	  a41, 
	  a42, 
	  a51, 
	  a52, 
	  a53, 
	  c21, 
	  Ca2pop[5] ,
	  cs, 
	  cs2s2p, 
	  cs2s3p	, 
	  cs01, 
	  cs02, 
	  cs12, 
	  cs14, 
	  cs15, 
	  d41, 
	  d42, 
	  d51, 
	  d52, 
	  d53, 
	  hlgam, 
	  op41, 
	  op51, 
	  opckh, 
	  opcxyz, 
	  PhotoRate2, 
	  p3, 
	  PhotoRate3, 
	  PhotoRate4, 
	  PhotoRate5, 
	  r21, 
	  r31, 
	  r41, 
	  r42, 
	  r51, 
	  r52, 
	  r53;
	static double gCa2[5]={2.,4.,6.,2.,4.};
	static double exCa2[4]={13650.2,60.7,11480.6,222.9};
	static realnum opCax = 0.f;
	static realnum opCay = 0.f;
	static realnum opCaz = 0.f;

	DEBUG_ENTRY( "CoolCalc()" );

	/* Ca I 4228 */
	MakeCS(&TauLines[ipCaI4228]);
	atom_level2(&TauLines[ipCaI4228]);

	/* photoionization of evcited levels by Ly-alpha */
	hlgam = rfield.otslin[ Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].ipCont -1];
	PhotoRate5 = 1.7e-18*hlgam;
	PhotoRate4 = 8.4e-19*hlgam;
	PhotoRate3 = 7.0e-18*hlgam;
	PhotoRate2 = 4.8e-18*hlgam;

	/* spontaneous decays
	 * frist two trans prob from 
	 * >>refer	ca2	as	Zeippen, C.J. 1990, A&A, 229, 248 */
	a21 = 1.02*TauLines[ipT7324].Emis->Pesc;
	a31 = 1.05*TauLines[ipT7291].Emis->Pesc;
	a41 = 1.4e8*TauLines[ipT3969].Emis->Pesc;
	a51 = 1.4e8*TauLines[ipT3934].Emis->Pesc;
	a42 = 7.9e6*TauLines[ipT8662].Emis->Pesc;
	a52 = 8.2e5*TauLines[ipT8498].Emis->Pesc;
	a53 = 7.48e6*TauLines[ipT8542].Emis->Pesc;

	/* destruction of IR triplet by continuous opacities */
	opcxyz = opac.opacity_abs[ TauLines[ipT7324].ipCont -1];

	/* opcxyz = opac(icaxyz) */
	if( opcxyz > 0. )
	{
		d52 = 5.6*opcxyz/(opcxyz + opCax)*(1. - TauLines[ipT8498].Emis->Pesc);
		d53 = 5.6*opcxyz/(opcxyz + opCay)*(1. - TauLines[ipT8542].Emis->Pesc);
		d42 = 5.6*opcxyz/(opcxyz + opCaz)*(1. - TauLines[ipT8662].Emis->Pesc);
	}
	else
	{
		d52 = 0.;
		d53 = 0.;
		d42 = 0.;
	}

	/* near UV dest of KH by background continuum */
	opckh = opac.opacity_abs[ TauLines[ipT3969].ipCont -1];

	/* opckh = opac(icakh) */
	if( opckh > 0. )
	{
		op51 = dense.xIonDense[ipCALCIUM][1]*3.89e-7/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]);
		d51 = 5.6*opckh/(opckh + op51);
		op41 = dense.xIonDense[ipCALCIUM][1]*1.96e-7/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]);
		d41 = 5.6*opckh/(opckh + op41);
	}
	else
	{
		op51 = 0.;
		d51 = 0.;
		op41 = 0.;
		d41 = 0.;
	}
	/* net rates */
	r21 = PhotoRate2 + a21;
	r31 = PhotoRate3 + a31;
	r41 = a41 + PhotoRate4 + d41;
	r51 = a51 + PhotoRate5 + d51;
	r42 = a42 + d42;
	r52 = a52 + d52;
	r53 = a53 + d53;
	cs14 = 0.923*phycon.te10*phycon.te10;
	cs15 = cs14*2.;
	TauLines[ipT3969].Coll.col_str = (realnum)cs14;
	TauLines[ipT3934].Coll.col_str = (realnum)cs15;

	/* following used to correct rec contribution
	 * fcakh = a51 / ( a51 + eden*1.5e-5 / sqrte ) 
	 * cs 1-2 from 
	 * >>refer	ca2	cs	Saraph, H.E. 1970, J.Phys. B, 3, 952
	 * other 
	 * >>refer	ca2	cs	Chidichimo, M.C. 1981, J.Phys. B, 14, 4149 */
	double Cooling , CoolingDeriv;
	atom_pop5(gCa2,exCa2,5.8,8.6,cs14,cs15,20.6,22.9,9.8,3.4,44.4,1.0,
		r21,r31,r41,r51,0.,r42,r52,0.,r53,0.,Ca2pop,
		dense.xIonDense[ipCALCIUM][1],&Cooling , &CoolingDeriv, 0.,0.,0.,0.);

	/* CDSQTE = 8.629E-6*EDEN/SQRTE */
	c21 = 5.8/4.*dense.cdsqte;

	/* remember largest ratio of Ly-al removal to total */
	if( dense.xIonDense[ipCALCIUM][1] > 0. )
		ca.Ca2RmLya = MAX2(ca.Ca2RmLya,(realnum)(PhotoRate2/(PhotoRate2+a21+c21)));

	ca.Cak = (realnum)(Ca2pop[4]*a51*5.06e-12);
	ca.Cah = (realnum)(Ca2pop[3]*a41*5.01e-12);
	ca.Cax = (realnum)(Ca2pop[4]*a52*2.34e-12);
	ca.Cay = (realnum)(Ca2pop[4]*a53*2.33e-12);
	ca.Caz = (realnum)(Ca2pop[3]*a42*2.30e-12);
	ca.Caf1 = (realnum)(Ca2pop[2]*a31*2.73e-12);
	ca.Caf2 = (realnum)(Ca2pop[1]*a21*2.72e-12);
	ca.popca2ex = (realnum)(Ca2pop[1] + Ca2pop[2] + Ca2pop[3] + Ca2pop[4]);

	/* this is the total cooling due to the model atom */
	ca.Cair = ca.Cax + ca.Cay + ca.Caz;
	ca.c7306 = ca.Caf1 + ca.Caf2;
	ca.Cakh = ca.Cak + ca.Cah;

	// total cooling from 5-level atom
	CoolAdd("Ca 2",7306,Cooling);
	thermal.dCooldT += CoolingDeriv;

	/*fprintf(ioQQQ,"DEBUG ca2\t%.2f\t%.5e\t%.4e\t%.4e\n",
		fnzone, phycon.te,ca.Cakh,dense.xIonDense[ipCALCIUM][1]);*/

	/* level populations that will be used for excited state photoionization */
	ca.dstCala = (realnum)(Ca2pop[4]*PhotoRate5 + Ca2pop[3]*PhotoRate4);
	ca.dCakh = (realnum)(ca.dstCala*5.03e-12);
	ca.dCaf12 = (realnum)((Ca2pop[2]*PhotoRate3 + Ca2pop[1]*PhotoRate2)*2.31e-12);
	opCax = (realnum)(Ca2pop[1]*1.13e-8/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]));
	opCay = (realnum)(Ca2pop[2]*6.87e-8/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]));
	opCaz = (realnum)(Ca2pop[1]*5.74e-8/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]));

	/* total rate Lalpha destroys CaII,
	 * this is only used in ioncali to increase ionization rate by
	 * adding it to the ct vector */
	if( dense.xIonDense[ipCALCIUM][1] > 0. )
	{
		ca.dstCala = (realnum)(
			(ca.dstCala + ca.dCaf12/2.31e-12)/dense.xIonDense[ipCALCIUM][1]);
		{
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ," hlgam is %e\n", hlgam);
			}
		}
	}
	else
	{
		ca.dstCala = 0.;
	}
	ca.Ca3d = (realnum)(Ca2pop[1] + Ca2pop[2]);
	ca.Ca4p = (realnum)(Ca2pop[3] + Ca2pop[4]);

	/* incl stimulated emission for Calcium II 5-level atom */
	TauLines[ipT3934].Emis->PopOpc = (Ca2pop[0] - Ca2pop[4]/2.);
	TauLines[ipT3934].Hi->Pop = Ca2pop[4];
	TauLines[ipT3934].Lo->Pop = Ca2pop[0];
	TauLines[ipT3969].Emis->PopOpc = (Ca2pop[0] - Ca2pop[3]);
	TauLines[ipT3969].Hi->Pop = Ca2pop[3];
	TauLines[ipT3969].Lo->Pop = Ca2pop[0];

	TauLines[ipT8498].Emis->PopOpc = (Ca2pop[1] - Ca2pop[4]);
	TauLines[ipT8498].Hi->Pop = Ca2pop[4];
	TauLines[ipT8498].Lo->Pop = Ca2pop[1];
	TauLines[ipT8542].Emis->PopOpc = (Ca2pop[2] - Ca2pop[4]*1.5);
	TauLines[ipT8542].Hi->Pop = Ca2pop[4];
	TauLines[ipT8542].Lo->Pop = Ca2pop[2];
	TauLines[ipT8662].Emis->PopOpc = (Ca2pop[1] - Ca2pop[3]*2.);
	TauLines[ipT8662].Hi->Pop = Ca2pop[3];
	TauLines[ipT8662].Lo->Pop = Ca2pop[1];
	TauLines[ipT7291].Emis->PopOpc = dense.xIonDense[ipCALCIUM][1];
	TauLines[ipT7291].Hi->Pop = 0.;
	TauLines[ipT7291].Lo->Pop = dense.xIonDense[ipCALCIUM][1];
	TauLines[ipT7324].Emis->PopOpc = dense.xIonDense[ipCALCIUM][1];
	TauLines[ipT7324].Hi->Pop = 0.;
	TauLines[ipT7324].Lo->Pop = dense.xIonDense[ipCALCIUM][1];

	/* Ca IV 3.2 micron; data from
	 * >>refer	ca4	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 * Y(ik) from 
	 * >>refer	ca4	cs	Pelan, J., & Berrington, K.A. 1995, A&A Suppl, 110, 209 */
	if( phycon.te <= 1e5 )
	{
		cs = MAX2(1.0,8.854e-3*phycon.sqrte);
	}
	else if( phycon.te < 2.512e5 )
	{
		cs = 2.8;
	}
	else
	{
		cs = 641.1/(phycon.te30*phycon.te10*phycon.te02*phycon.te02/
		  phycon.te003);
	}
	PutCS(cs,&TauLines[ipTCa3]);
	atom_level2(&TauLines[ipTCa3]);

	/* [Ca V] IR 4.16, 11.47 micron; A from
	 * >>refer	ca5	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 * cs from 
	 * >>refer	ca5	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347
	 * >>chng 96 jul 16, big changes in cs */
	cs = MIN2(3.3,0.392*phycon.te20/phycon.te005/phycon.te003);
	cs = MAX2(2.2,cs);
	PutCS(cs,&TauLines[ipTCa4]);

	/* >>chng 96 aug 02, following had error in te dep in ver 90.01 */
	cs = MIN2(0.93,0.162*phycon.te10*phycon.te05*phycon.te003*
	  phycon.te001);
	cs = MAX2(0.67,cs);
	PutCS(cs,&TauLines[ipTCa12]);

	cs = MIN2(0.97,0.0894*phycon.te20*phycon.te01*phycon.te005);
	cs = MAX2(0.60,cs);
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipTCa4],&TauLines[ipTCa12],&TauDummy);

	/* Ca V lines from 1d, 1s; A from 
	 * >>refer	ca5	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 * cs from 
	 * >>refer	ca5	cs	 Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	cs01 = MIN2(4.1,0.533*phycon.te20/phycon.te01);
	cs01 = MAX2(2.8,cs01);
	cs02 = MIN2(0.87,5.22e-03*phycon.sqrte);
	p3 = atom_pop3(9.,5.,1.,cs01,cs02,1.35,2.326,23.2,3.73,2.57e4,3.60e4,
	  &p2,dense.xIonDense[ipCALCIUM][4],0.,0.,0.);

	ca.c3997 = p3*3.73*4.98e-12;
	ca.c2414 = p3*23.1*8.245e-12;
	ca.Ca6087 = p2*0.426*3.268e-12;
	ca.c5311 = p2*1.90*3.747e-12;

	CoolAdd("Ca 5",3997,ca.c3997);
	CoolAdd("Ca 5",2414,ca.c2414);
	CoolAdd("Ca 5",6087,ca.Ca6087);
	CoolAdd("Ca 5",5311,ca.c5311);

	/* Ca VII lines from 1d, 1s
	 * all cs from 
	 * >>refer	ca72	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	cs01 = MIN2(4.4,22.25/(phycon.te20/phycon.te02/phycon.te02));
	cs01 = MAX2(3.5,cs01);

	cs12 = MIN2(1.20,0.303*phycon.te30*phycon.te03);
	cs12 = MAX2(0.62,cs12);

	cs02 = MIN2(0.959,7.889/(phycon.te20*phycon.te05/phycon.te01));
	cs02 = MAX2(0.50,cs02);

	p3 = atom_pop3(9.,5.,1.,cs01,cs02,cs12,3.124,30.4,6.81,2.91e4,3.90e4,
	  &p2,dense.xIonDense[ipCALCIUM][6],0.,0.,0.);

	ca.Ca3688 = p3*6.81*5.40e-12;
	ca.Ca2112 = p3*30.4*9.42e-12;
	ca.Ca5620 = p2*2.15*3.548e-12;
	ca.Ca4941 = p2*0.974*4.037e-12;
	CoolAdd("Ca 7",3688,ca.Ca3688);
	CoolAdd("Ca 7",2112,ca.Ca2112);
	CoolAdd("Ca 7",5620,ca.Ca5620);
	CoolAdd("Ca 7",4941,ca.Ca4941);

	/* all cs from 
	 * >>refer	ca7	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347
	 * [Ca VII] 4.09, 6.15 mic 3P lines */
	cs = MIN2(5.354,0.406*phycon.te20*phycon.te03*phycon.te01);
	cs = MAX2(3.702,cs);
	PutCS(cs,&TauLines[ipCa0741]);

	cs = MIN2(1.59,0.183*phycon.te20);
	cs = MAX2(1.153,cs);
	PutCS(cs,&TauLines[ipCa0761]);

	cs = MIN2(1.497,0.0917*phycon.te20*phycon.te05* phycon.te01);
	cs = MAX2(1.005,cs);
	PutCS(cs,&TauDummy);

	/* atom_level3(  t10,t21,t20) */
	atom_level3(&TauLines[ipCa0761],&TauLines[ipCa0741],&TauDummy);

	/* [Ca VIII]  2.32 microns, cs 
	 * >>refer	ca8	cs	Saraph, H.E., & Storey, P.J. A&AS, 115, 151 */
	cs = MIN2(6.75,22.04/(phycon.te10*phycon.te02*phycon.te005));
	PutCS(cs,&TauLines[ipCa08232]);
	atom_level2(&TauLines[ipCa08232]);

	/* [Ca 12] 3328.78A, cs from 
	 * >>refer	ca12	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	cs = MIN2(0.172,0.0118*phycon.te20*phycon.te01);
	cs = MAX2(0.10,cs);
	PutCS(cs,&TauLines[ipCa12333]);
	atom_level2(&TauLines[ipCa12333]);

	/* Li seq Ca 18 2s2p 2s3p, 2s2p as two separate lines
	 * >>refer	ca18	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(20,&TauLines[ipTCa302],&TauLines[ipTCa19],&cs2s2p,&cs2s3p);

	PutCS(cs2s2p,&TauLines[ipTCa302]);
	atom_level2(&TauLines[ipTCa302]);

	/* funny factor (should have been 0.5) due to energy change */
	PutCS(cs2s2p*0.439,&TauLines[ipTCa345]);
	atom_level2(&TauLines[ipTCa345]);

	PutCS(cs2s3p,&TauLines[ipTCa19]);
	atom_level2(&TauLines[ipTCa19]);
	return;
}
