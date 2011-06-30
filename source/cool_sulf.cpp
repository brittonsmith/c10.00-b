/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSulf compute sulphur cooling */
/*S2cs compute [CoolHeavy] collision strengths
 * compute collision strengths for [SII] transitions 
 * w/in S II ground term. From 
 *>>refer	s2	cs	Ramsbottom, C.A., Bell, K.L., Stafford, R.P. 1996, At. Data Nucl. Data Tables, 63, 57 */
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

void CoolSulf(void)
{
	double cs, 
	  cs2s2p, 
	  cs2s3p, 
	  cs01, 
	  cs02, 
	  cs12, 
	  a21,
	  a31,
	  a32,
	  a41,
	  a42,
	  a43,
	  a51,
	  a52,
	  a53,
	  p[5],
	  tused;
	realnum pop2,
	  rate;
	static double gS2[5]={4.,4.,6.,2.,4.};
	static double exS2[4]={14851.9,31.5,9640.8,48.6};

	long int i;

	DEBUG_ENTRY( "CoolSulf()" );

	/*>>refer	S1	cs	Hollenbach, D. & McKee, C.F. 1989, ApJ, 342, 306 */
	/* >>chng 03 nov 15, add these lines */
	/* rates are said to be ok over range 30 - 3000K */
	tused = MAX2( 30. , phycon.te );
	tused = MIN2( 3000. , phycon.te  );
	tused /= 100.;

	/* the 25.2 micron line */
	rate = (realnum)(3.3e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*7.5e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		7.5e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	LineConvRate2CS( &TauLines[ipS1_25m]  , rate );

	/* the 56.6 micron line */
	rate = (realnum)(1.2e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*4.2e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		4.2e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	LineConvRate2CS( &TauLines[ipS1_56m]  , rate );

	rate = (realnum)(3.3e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*7.1e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		7.1e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	TauDummy.Hi->g = TauLines[ipS1_56m].Hi->g;
	LineConvRate2CS( &TauDummy  , rate );
	/* this says that line is a dummy, not real one */
	TauDummy.Hi->g = 0.;

	atom_level3(&TauLines[ipS1_25m],&TauLines[ipS1_56m],&TauDummy);

	/***********************************************************************************
	**************************************S II*****************************************
	************************************************************************************/
	/* Sulphur II [S II] 6731 + 6716, A 
	 * >>referold	s2	as	trans Mendoza, C., & Zeippen, C.J., 1982, MNRAS, 198, 127
	 * >>referold	s2	as	Froese Fischer, C., Tachiev, G., & Irimia, A. 2006, At. Data Nucl. Data Tables, 92, 607
	 * >>refer	s2	as	Podobedova, L.I., Kelleher, D.E., & Wiese, W.L. 2009, JPCRD, 38, 171
	 * collision strength from 
	 *>>refer	s2	cs	Ramsbottom, C.A., Bell, K.L., Stafford, R.P. 1996,
	 *>>refercon	At. Data Nucl. Data Tables, 63, 57 
	 * this agrees very well with 
	 * >>refer	s2	cs	Tayal, S., 1997, ApJS, 111, 459
	 * */
	double sii_cs4S32D3,
		sii_cs4S32D5,
		sii_cs4S32P1,
		sii_cs4S32P3,
		sii_cs2D32D5,
		sii_cs2D32P1,
		sii_cs2D52P1,
		sii_cs2D32P3,
		sii_cs2D52P3,
		sii_cs2P12P3,
		sii_cs4S34P;

	/*sii_cs compute [CoolHeavy] collision strengths
	 * compute collision strengths for [SII] transitions
	 * w/in S II ground term. From
	 *>>refer	s2	cs	Ramsbottom, C.A., Bell, K.L., Stafford, R.P. 1996, ADNDT, 63, 57 */
	sii_cs(sii_cs4S32D3,
		sii_cs4S32D5,
		sii_cs4S32P1,
		sii_cs4S32P3,
		sii_cs2D32D5,
		sii_cs2D32P1,
		sii_cs2D52P1,
		sii_cs2D32P3,
		sii_cs2D52P3,
		sii_cs2P12P3,
		sii_cs4S34P);
	
	/* cs13 = MIN(4.95 , 3.2640828 + 1.9369099*sexp(te/30337.498) )
	 * cs12 = cs13 / 1.5
	 * cs15 = MIN( 1.59 , 0.10*te10*te10*te10/te01/te003/te001/te001)
	 * cs14 = cs15 * 0.5
	 * FIVEL( G(1-5) , ex(wn,1-5), cs12,cs13,14,15,23,24,25,34,35,45,
	 *  A21,31,41,51,32,42,52,43,53,54, pop(1-5), abund) */

	a21 = 6.84e-4;
	a31 = 2.02e-4;
	a41 = 7.72e-2;
	a42 = 0.135;
	a43 = 6.81e-2;
	a51 = 1.92e-1;
	a52 = 0.115;
	a53 = 0.157;

	double Cooling , CoolingDeriv;
	atom_pop5(gS2,exS2,sii_cs4S32D3,sii_cs4S32D5,sii_cs4S32P1,sii_cs4S32P3,sii_cs2D32D5,
	    sii_cs2D32P1,sii_cs2D32P3,sii_cs2D52P1,sii_cs2D52P3,sii_cs2P12P3,
	    a21,a31,a41,a51,3.35e-7,a42,a52,a43,a53,
	    1.03e-6,p,dense.xIonDense[ipSULPHUR][1],	&Cooling , &CoolingDeriv, 0.,0.,0.,0.);
	CoolHeavy.S6733 = (realnum)(p[1]*a21*2.96e-12);
	CoolHeavy.S6718 = (realnum)(p[2]*a31*2.962e-12);
	CoolHeavy.S4070 = (realnum)(p[4]*a51*4.89e-12);
	CoolHeavy.S4078 = (realnum)(p[3]*a41*4.88e-12);
	CoolHeavy.S10323 = (realnum)(p[4]*a53*1.93e-12);
	CoolHeavy.S10289 = (realnum)(p[4]*a52*1.93e-12);
	CoolHeavy.S10373 = (realnum)(p[3]*a43*1.92e-12);
	CoolHeavy.S10339 = (realnum)(p[3]*a42*1.92e-12);
	CoolHeavy.c6731 = CoolHeavy.S6733 + CoolHeavy.S6718;
	CoolHeavy.c10330 = CoolHeavy.S10323 + CoolHeavy.S10289 + CoolHeavy.S10373 + CoolHeavy.S10339;

	// all cooling from 5-level atom
	CoolAdd("S  2",6731,Cooling);
	thermal.dCooldT += CoolingDeriv;

	
	PutCS(sii_cs4S34P,&TauLines[ipT1256]);
	atom_level2(&TauLines[ipT1256]);

	/***********************************************************************************
	**************************************S III*****************************************
	************************************************************************************/
	double siii_cs3P03P1,
		siii_cs3P03P2,
		siii_cs3P01D2,
		siii_cs3P01S0,
		siii_cs3P13P2,
		siii_cs3P11D2,
		siii_cs3P11S0,
		siii_cs3P21D2,
		siii_cs3P21S0,
		siii_cs1D21S0,
		siii_cs3P3D,
		siii_cs3P5S2;
	// siii_cs calculates collision strengths for transitions of S III
	siii_cs(siii_cs3P03P1,
		siii_cs3P03P2,
		siii_cs3P01D2,
		siii_cs3P01S0,
		siii_cs3P13P2,
		siii_cs3P11D2,
		siii_cs3P11S0,
		siii_cs3P21D2,
		siii_cs3P21S0,
		siii_cs1D21S0,
		siii_cs3P3D,
		siii_cs3P5S2);
	
	PutCS(siii_cs3P03P1,&TauLines[ipTS34]);
	PutCS(siii_cs3P13P2,&TauLines[ipTS19]);
	PutCS(siii_cs3P03P2,&TauDummy);
	atom_level3(&TauLines[ipTS34],&TauLines[ipTS19],&TauDummy);

	/* S III O III-like lines, A from 
	 * >>referold	s3	as	Mendoza, C., & Zeippen, C.J. 1982, MNRAS, 199, 1025
	 * >>referold	s3	as	Froese Fischer, C., Tachiev, G., & Irimia, A. 2006, At. Data Nucl. Data Tables, 92, 607
	 * >>refer	s3	as	Podobedova, L.I., Kelleher, D.E., & Wiese, W.L. 2009, JPCRD, 38, 171
	 * CS from 
	 * >>referold	s3	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347
	 * >>chng 00 Sep 11, cs changed from above to
	 * >>refer	s3	cs	Tayal, S.S., and Gupta, G.P. 1999 ApJ 526, 544 */
	/*cs = MIN2(2.05,0.0821*phycon.te30);*/
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	/*CoolHeavy.c6312 = atom_pop3(9.,5.,1.,7.98,1.14,cs,7.97e-2,0.807,2.22,*/

	a21 = 0.0666;
	a31 = 0.670;
	a32 = 2.08;

	CoolHeavy.c6312 = atom_pop3(9.,5.,1.,siii_cs3P01D2+siii_cs3P11D2+siii_cs3P21D2,
		siii_cs3P01S0+siii_cs3P11S0+siii_cs3P21S0,siii_cs1D21S0,a21,a31,a32,
	  1.55e4,2.28e4,&pop2,dense.xIonDense[ipSULPHUR][2],0.,0.,0.)*2.22*3.15e-12;
	/* folowing is 9532 + 9069 together (OIII-like) */
	CoolHeavy.c9532 = pop2*a21*2.11e-12;
	CoolAdd("S  3",6312,CoolHeavy.c6312);
	CoolAdd("S  3",9532,CoolHeavy.c9532);
	CoolAdd("S  3",3722,CoolHeavy.c6312*0.59);

	PutCS(siii_cs3P3D,&TauLines[ipT1194]);
	atom_level2(&TauLines[ipT1194]);

	PutCS(siii_cs3P5S2,&TauLines[ipTS1720]);
	atom_level2(&TauLines[ipTS1720]);


	
	/***********************************************************************************
	**************************************S IV*****************************************
	************************************************************************************/
	double siv_cs2P12P3;
	// siv_cs calculates collision strengths for transitions of S IV
	// unknown transition
	siv_cs(siv_cs2P12P3);

	PutCS(siv_cs2P12P3,&TauLines[ipTS11]);
	/*atom_level2(&TauLines[ipTS11]);*/

	static vector< pair<transition*,double> > S4Pump;
	S4Pump.reserve(48);

	/* one time initialization if first call */
	if( S4Pump.empty() )
	{
		// set up level 2 pumping lines
		for( i=0; i < nWindLine; ++i )
		{
			/* don't test on nelem==ipIRON since lines on physics, not C, scale */
			if( TauLine2[i].Hi->nelem == 16 && TauLine2[i].Hi->IonStg == 4 )
			{
#				if	0
				DumpLine( &TauLine2[i] );
#				endif
				double branch_ratio;
				// the branching ratios used here ignore cascades via intermediate levels
				// usually the latter are much slower, so this should be reasonable
				if( fp_equal( TauLine2[i].Hi->g, realnum(2.) ) )
					branch_ratio = 2./3.; // 2S upper level
				else if( fp_equal( TauLine2[i].Hi->g, realnum(6.) ) )
					branch_ratio = 1./2.; // 2P upper level
				else if( fp_equal( TauLine2[i].Hi->g, realnum(10.) ) )
					branch_ratio = 1./6.; // 2D upper level
				else
					TotalInsanity();
				pair<transition*,double> pp2( &TauLine2[i], branch_ratio ); 
				S4Pump.push_back( pp2 );
			}
		}
	}

	/* now sum pump rates */
	double pump_rate = 0.;
	vector< pair<transition*,double> >::const_iterator s4p;
	for( s4p=S4Pump.begin(); s4p != S4Pump.end(); ++s4p )
	{
		const transition* t = s4p->first;
		double branch_ratio = s4p->second;
		pump_rate += t->Emis->pump*branch_ratio;
#		if	0
		dprintf( ioQQQ, "S IV %.3e %.3e\n",
			 t->WLAng , t->Emis->pump*branch_ratio );
#		endif
	}

	/* S IV 1404.8, 1398.05, 1423.8, 1416.9, 1406.0 */
	/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
	/* >>refer	s4	cs	Tayal, S.S., 2000, ApJ 530, 1091*/
	AtomSeqBoron(&TauLines[ipTS11], 
	  &TauLines[ipS4_1405], 
	  &TauLines[ipS4_1398], 
	  &TauLines[ipS4_1424], 
	  &TauLines[ipS4_1417], 
	  &TauLines[ipS4_1407], 
	  1.168 , 3.366 , 2.924 , 7.233 , 
	  pump_rate , "S  4");

	/***********************************************************************************
	**************************************S V*****************************************
	************************************************************************************/
	/* S V Be-seq line, A=
	 * >>refer	s5	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 * all cs 
	 * >>refer	s5	as	Dufton, P.L., Hibbert, A., Keenan, F.P, Kingston, A.E., &
	 * >>refercon	Doschek, G.A. 1986, ApJ, 300, 448
	 * and 
	 * >>refer	s5	cs	Dufton, P.L., & Kingston, A.E. 1984, J.Phys. B, 17, 3321 */
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	cs = MIN2(1.58,35.372/(phycon.te10*phycon.te10*phycon.te10));
	PutCS(cs,&TauLines[ipT1198]);
	/* cs01 = MIN(1.98, 29.625/(T**0.25)) */
	cs01 = MIN2(1.98,29.625/(phycon.te20*phycon.te05));
	/* cs02 = MIN(2.26, 13.477/(T**0.165)) */
	cs02 = MIN2(2.26,13.477/(phycon.te10*phycon.te03*phycon.te03*
	  phycon.te005));
	/* cs12 = MIN(7.59, 63.994/(T**0.197)) */
	cs12 = MIN2(7.59,63.994/(phycon.te20/phycon.te003));
	/* AtomSeqBeryllium(cs23,cs24,cs34,tarray,a41)
	 * call AtomSeqBeryllium( 1.98 , 2.26 , 7.59, t1198 ,.066) */
	AtomSeqBeryllium(cs01,cs02,cs12,&TauLines[ipT1198],.066);
	embesq.em1198 = (realnum)(atoms.PopLevels[3]*0.066*1.66e-11);

	/* S V 786
	 * cs from 
	 * >>refer	s5	cs	Dufton, P.L., & Kingston, A.E. 1984, J.Phys. B, 17, 3321 */
	/** \todo	1	upgrade SV to more levels there is intercombination line at 0.7634 ryd
	 * upgrade to atomic data described in
	 *>>refer	s5	cs	Hudson, C.E> & Bell, K.L. 2006, A&A, 452, 1113 */
	PutCS(8.3,&TauLines[ipT786]);
	atom_level2(&TauLines[ipT786]);

	/***********************************************************************************
	**************************************S VI*****************************************
	************************************************************************************/
	/* S VI, Li sequence, A from
	 * >>refer	s6	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 * CS from 
	 * >>refer	s6	cs	Dufton, P.L., & Kingston, A.E. 1987, J.Phys. B, 20, 3899 */
	cs = 11.9;
	PutCS(cs*0.667,&TauLines[ipT933]);
	PutCS(cs*0.333,&TauLines[ipT944]);
	PutCS(1.0,&TauDummy);

	atom_level3(&TauLines[ipT944],&TauDummy,&TauLines[ipT933]);

	/***********************************************************************************
	**************************************S VIII*****************************************
	************************************************************************************/
	double sviii_cs2P32P1;

	// sviii_cs calculates collision strengths for transitions of S VIII
	sviii_cs(sviii_cs2P32P1);

	PutCS(sviii_cs2P32P1,&TauLines[ipfsS810]);

	atom_level2(&TauLines[ipfsS810]);

	/***********************************************************************************
	**************************************S IX*****************************************
	************************************************************************************/
	/* S IX 1715, 1987, 
	 * >>refer	s9	all	Kafatos, M., & Lynch, J.P. 1980, ApJS, 42, 611 */
	CoolHeavy.c1715 = atom_pop2(0.4,9.,5.,80.,8.39e4,dense.xIonDense[ipSULPHUR][8])*
	  1.16e-11;
	CoolAdd("S  9",1715,CoolHeavy.c1715);

	/* S IX 1.25, 3.78 microns
	 * cs 
	 * >>refer	s9	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1 */
	cs = MIN2(3.10,277.2/(phycon.sqrte/phycon.te003));
	PutCS(cs,&TauLines[ipfsS912]);

	cs = MIN2(.50,7.940/(phycon.te10*phycon.te10*phycon.te10));
	PutCS(cs,&TauLines[ipfsS938]);

	cs = MIN2(1.1,195.72/(phycon.sqrte*phycon.te03*phycon.te03*
	  phycon.te01*phycon.te003*phycon.te001));
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipfsS912],&TauLines[ipfsS938],&TauDummy);

	/* S X 1213, 1197, 
	 * >>refer	si10	cs	Bell, K.L., & Ramsbottom, C.A., 2000, ADNDT, 76, 176-190 */
	CoolHeavy.c1213 = atom_pop2(0.131,4.,6.,6.,1.19e5,dense.xIonDense[ipSULPHUR][10-1])*
	  1.64e-11;
	CoolAdd("S 10",1213,CoolHeavy.c1213);

	/* S XI 1.91, 1.40 microns, cs 
	 * >>refer	si11	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	cs = MIN2(0.22,0.0151*phycon.te20*phycon.te03*phycon.te001*
	  phycon.te001);
	PutCS(cs,&TauLines[ipfsS1119]);
	cs = MIN2(0.6,0.0168*phycon.te10*phycon.te10*phycon.te10*
	  phycon.te01);
	PutCS(cs,&TauLines[ipfsS1114]);
	cs = MIN2(0.15,1.5e-3*phycon.te30*phycon.te10*phycon.te005);
	PutCS(cs,&TauDummy);
	atom_level3(&TauLines[ipfsS1119],&TauLines[ipfsS1114],&TauDummy);

	/* S XI 1615, 1826, 
	 * >>refer	si11	cs	Kafatos, M., & Lynch, J.P. 1980, ApJS, 42, 611 */
	CoolHeavy.c1826 = atom_pop2(0.3,9.,5.,100.,8.92e4,dense.xIonDense[ipSULPHUR][11-1])*
	  1.23e-11;
	CoolAdd("S 11",1615,CoolHeavy.c1826);

	/* S XII 7611.2, cs from 
	 * >>refer	s12	cs	Oliva, E., Pasquali, A., & Reconditi, M. 1996, A&A, 305, 210
	 * >>refer	s12	as	Froese Fischer, C. 1983, J.Phys. B, 16, 157
	 * >>refer	s12	as	Saha, H.P., & Trefftz, E. 1983, SoPh, 87, 233 */
	PutCS(0.135,&TauLines[ipfsS1207]);
	atom_level2(&TauLines[ipfsS1207]);

	/* S XII 520, actually group of four intercombination lines
	 * all atomic data from 
	 * >>refer	s12	all	Saha, H.P., & Trefftz, E. 1983, SoPh, 87, 233
	CoolHeavy.c520 = 0.03*1.65e-16*expion(2.8e5,dense.xIonDense[ipSULPHUR][12-1]);
	CoolAdd("S 12",520,CoolHeavy.c520); */
	/* >>chng 03 sep 27, move to simple two level atom */
	PutCS(0.03,&TauLines[ipS12_520]);
	atom_level2(&TauLines[ipS12_520]);

	/* S XIII 488.4, 1909-like */
	CoolHeavy.c488 = atom_pop2(0.11,1.,9.,3.3e5,2.948e5,dense.xIonDense[ipSULPHUR][13-1])*
	  4.07e-11;
	CoolAdd("S 13",488,CoolHeavy.c488);

	/* li seq 2s2p and 2s3p, S XIV 418, 446
	 * >>refer	s14	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(16,&TauLines[ipTSu418],&TauLines[ipTSu30],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipTSu418]);
	PutCS(cs2s2p*0.5,&TauLines[ipTSu446]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipTSu446],&TauDummy,&TauLines[ipTSu418]);

	PutCS(cs2s3p,&TauLines[ipTSu30]);
	atom_level2(&TauLines[ipTSu30]);
	return;
}

/*sii_cs compute [CoolHeavy] collision strengths
 * compute collision strengths for [SII] transitions 
 * w/in S II ground term. From 
 *>>refer	s2	cs	Ramsbottom, C.A., Bell, K.L., Stafford, R.P. 1996, ADNDT, 63, 57 */
void sii_cs(double& sii_cs4S32D3,
	double& sii_cs4S32D5,
	double& sii_cs4S32P1,
	double& sii_cs4S32P3,
	double& sii_cs2D32D5,
	double& sii_cs2D32P1,
	double& sii_cs2D52P1,
	double& sii_cs2D32P3,
	double& sii_cs2D52P3,
	double& sii_cs2P12P3,
	double& sii_cs4S34P)
{
    double telog;
    telog = phycon.alogte;
	/* written by Kirk Korista  */
	double a, 
	  b, 
	  c, 
	  telogn1;

	DEBUG_ENTRY( "sii_cs()" );

	/* limit to stop exceeding bounds */
	telogn1 = MAX2(3.5,telog);
	telogn1 = MIN2(telogn1,5.0);

	/* 2D5/2 - 2P3/2    S II  10320.4   A=0.179 (3-->5) */
	a = 18.335524;
	b = -5.1180248;
	c = 0.44482438;
	sii_cs2D52P3 = MIN2(5.82,a+b*telogn1+c*telogn1*telogn1);
	sii_cs2D52P3 = MAX2(3.87,sii_cs2D52P3);

	/* 2D3/2 - 2P3/2          10286.7   A=0.1335 (2-->5) */
	a = 6.690242;
	b = -1.061514;
	c = 0.034535506;
	sii_cs2D32P3 = MIN2(3.38,a+b*telogn1+c*telogn1*telogn1);
	sii_cs2D32P3 = MAX2(2.24,sii_cs2D32P3);

	/* 2D5/2 - 2P1/2          10373.3   A=0.0779 (3-->4) */
	a = 4.2250081;
	b = -0.46549935;
	c = -0.010172139;
	sii_cs2D52P1 = MIN2(2.46,a+b*telogn1+c*telogn1*telogn1);
	sii_cs2D52P1 = MAX2(1.64,sii_cs2D52P1);

	/* 2D3/2 - 2P1/2          10336.3   A=0.1626 (2-->4) */
	a = 8.274085;
	b = -2.6223732;
	c = 0.2502924;
	sii_cs2D32P1 = MIN2(2.14,a+b*telogn1+c*telogn1*telogn1);
	sii_cs2D32P1 = MAX2(1.42,sii_cs2D32P1);

	/* 2P1/2 - 2P3/2   */
	a = -5.1994665;
	b = 49.334586;
	c = -70.93344;
	sii_cs2P12P3 = MIN2(3.07,a+b/telogn1+c/(telogn1*telogn1));
	sii_cs2P12P3 = MAX2(1.85,sii_cs2P12P3);

	/* 2D3/2 - 2D5/2 */
	a = -27.497273;
	b = 247.27405;
	c = -429.9142;
	sii_cs2D32D5 = MIN2(8.01,a+b/telogn1+c/(telogn1*telogn1));
	sii_cs2D32D5 = MAX2(4.79,sii_cs2D32D5);

	/* 4S3/2 - 2P3/2           4068.6   A=0.220  */
	a = 2.6106784;
	b = -3.2766908e-05;
	c = 6.5105436;
	sii_cs4S32P3 = a+b*pow(telogn1,c);
	sii_cs4S32P3 = MIN2(2.46,sii_cs4S32P3);
	sii_cs4S32P3 = MAX2(1.45,sii_cs4S32P3);

	/* 4S3/2 - 2P1/2           4076.4   A=0.091       */
	sii_cs4S32P1 = 0.5*sii_cs4S32P3;

	/* 4S3/2 - 2D5/2            6716.5   A=2.601e-04   */
	a = 8.1458628;
	b = -0.5389108;
	c = 1.4486586;
	sii_cs4S32D5 = a+b*pow(telogn1,c);
	sii_cs4S32D5 = MIN2(4.77,sii_cs4S32D5);
	sii_cs4S32D5 = MAX2(2.54,sii_cs4S32D5);

	/* 4S3/2 - 2D3/2           6730.8   A=8.82e-04    */
	sii_cs4S32D3 = sii_cs4S32D5/1.5;

	/* SII 1256 */
	sii_cs4S34P = MIN2(8.46,-4.9416304+47.01064/phycon.alogte);
	sii_cs4S34P = MAX2(4.466,sii_cs4S34P);

	return;
}
// siii_cs calculates collision strengths for transitions of S III
void siii_cs(double& siii_cs3P03P1,
		double& siii_cs3P03P2,
		double& siii_cs3P01D2,
		double& siii_cs3P01S0,
		double& siii_cs3P13P2,
		double& siii_cs3P11D2,
		double& siii_cs3P11S0,
		double& siii_cs3P21D2,
		double& siii_cs3P21S0,
		double& siii_cs1D21S0,
		double& siii_cs3P3D,
		double& siii_cs3P5S2)
{
    DEBUG_ENTRY( "siii_cs()" );
    /* S III 18.7M, 33.6M, A
	 * >>referold	s3	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 *
	 * >>refer	s3	cs	Galavis, M.E., Mendoza, C., & Zeippen, C.J. 1995, A&AS, 111, 347
	 * >>chng 99 dec 22, cs changed from above to
	 * >>refer	s3	cs	Tayal, S.S., and Gupta, G.P. 1999 ApJ 526, 544 */
	/* the 1-2 transition */
	if( phycon.te < 5000. )
	{
		siii_cs3P03P1 = 4.44;
		siii_cs3P03P2 = 1.41;
		siii_cs3P01D2 = 0.802;
		siii_cs3P01S0 = 0.129;
		siii_cs3P13P2 = 8.72;
		siii_cs3P11D2 = 2.41;
		siii_cs3P11S0 = 0.388;
		siii_cs3P21D2 = 4.01;
		siii_cs3P21S0 = 0.646;
		siii_cs1D21S0 = 1.31;
	}
	else if( phycon.te > 1e5 )
	{
		siii_cs3P03P1 = 1.9;
		siii_cs3P03P2 = 1.24;
		siii_cs3P01D2 = 0.664;
		siii_cs3P01S0 = 0.136;
		siii_cs3P13P2 = 5.13;
		siii_cs3P11D2 = 1.99;
		siii_cs3P11S0 = 0.407;
		siii_cs3P21D2 = 3.32;
		siii_cs3P21S0 = 0.679;
		siii_cs1D21S0 = 1.84;
	}
	else
	{
		siii_cs3P03P1 = 52.47/(phycon.te30/phycon.te02);
		siii_cs3P03P2 = 1.894/(phycon.te02*phycon.te02);
		siii_cs3P01D2 = 1.34/(phycon.te05*phycon.te01);
		siii_cs3P01S0 = 0.109*phycon.te02;
		siii_cs3P13P2 = 41.3/(phycon.te20/phycon.te02);
		siii_cs3P11D2 = 4.03/(phycon.te05*phycon.te01);
		siii_cs3P11S0 = 0.327*phycon.te02;
		siii_cs3P21D2 = 6.708/(phycon.te05*phycon.te01);
		siii_cs3P21S0 = 0.545*phycon.te02;
		siii_cs1D21S0 = 0.501*phycon.te10*phycon.te01;
	}
	/*cs = MIN2(2.331,7.935*phycon.te/(phycon.te10*phycon.te03*phycon.te003));*/

	/* the 2-3 transition
	if( phycon.te <= 39811. )
	{
		cs = MIN2(5.78,3.114*phycon.te03*phycon.te03);
	}
	else
	{
		cs = 24.93/(phycon.te10*phycon.te03*phycon.te01/phycon.te001/
		  phycon.te001);
	}*/

	/* the 1-3 transition
	cs = MIN2(1.413,0.221*phycon.te*phycon.te20/phycon.te03*phycon.te005);*/

	/* S III 1194, data from
	 * >>refer	s3	cs	Ho, Y.K., & Henry,R.J.W. 1984, ApJ, 282, 816
	 * >>chng 97 may 17, to, about 2x larger than above
	 * >>refer	s3	cs	Tayal, S.S. 1997, ApJ 481, 550 */
	if( phycon.te <= 3e4 )
	{
		siii_cs3P3D = 12.04/(phycon.te02*phycon.te02);
	}
	else if( phycon.te > 3e4 && phycon.te <= 4e4 )
	{
		siii_cs3P3D = 7.97;
	}
	else
	{
		siii_cs3P3D = 55.42/(phycon.te20/phycon.te02*phycon.te003);
	}
	/* S III] 1713.12, 1728.94, cs from
	 * >>refer	s3	cs	Hayes, M.A., 1986, J Phys B 19, 1853.
	 * cs = MIN( 4.0 , 7.794 / (te10/te02/te001/te001)  )
	 * >>chng 97 may 17, about 20% smaller than before
	 * >>refer	s3	cs	Tayal, S.S. 1997, ApJ 481, 550 */
	if( phycon.te <= 3e4 )
	{
		siii_cs3P5S2 = 1.786*phycon.te05*phycon.te01*phycon.te001;
	}
	else
	{
		siii_cs3P5S2 = 9.392/phycon.te10;
	}
	
}

// siv_cs calculates collision strengths for transitions of S IV
// unknown transition
void siv_cs(double& siv_cs2P12P3)
{
    DEBUG_ENTRY( "siv_cs()" );
    /* S IV 1062 */
	/*>>refer	S4	As	Hibbert, A., Brage, T., Fleming, J. 2002, MNRAS 333, 885,
	 * typo noted in CHIANTI data file */
	/*>>refer	S4	cs	Tayal S.S., 2000, ApJ, 530, 1091 */

	/* S IV 10.5MI,
	 * >>refer	s4	as	Johnson, C.T., Kingston, A.E., Dufton, P.L. 1986, 220, 155
	 * >>referold	s4	cs	Johnson, C.T., Kingston, A.E., Dufton, P.L. 1986, MNRAS, 220, 155
	 * >>chng 97 feb 14, error in cs below t = 10,000K
	 * >>chng 96 dec 19, to CS from
	 * >>referold	s4	cs	Saraph, H.E., Storey, P.J., & Tully, J.A. 1995, 5th International
	 * >>referoldcon	Colloquium on Atomic Spectra and Oscillator Strengths, ed. by
	 * >>referoldcon	W.-U L.  Tchang-Brillet, J.-F. Wyart, C.J. Zeippen,
	 * >>referoldcon	(Meudon: Publications de l'Observaroire de Paris), p.110
	 * above said to be A&A in press */
	/* >>refer	s4	cs	Tayal, S.S., 2000, ApJ, 530, 1091
	 */
	if( phycon.te < 1e4 )
	{
		siv_cs2P12P3 = 3.71*phycon.te10/phycon.te01;
	}
	else
	{
		siv_cs2P12P3 = MIN2(8.5,19.472/(phycon.te10/phycon.te01));
	}

	return;
}

// sviii_cs calculates collision strengths for transitions of S VIII
void sviii_cs(double& sviii_cs2P32P1)
{
    DEBUG_ENTRY( "sviii_cs()" );
    /* S VIII 9913
	 * >>referold	s8	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29 */
	/*cs = MIN2(0.291,0.0289*phycon.te20/phycon.te01* phycon.te001);
	cs = MAX2(0.192,cs);*/
	/* >>refer	s8	cs	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 129,161 */
	/*>>chng 06 jul 18 Changes made-Humeshkar Nemala*/
	if(phycon.te < 6.4E5)
	{
		sviii_cs2P32P1 = (realnum)(0.0943*(phycon.te10/(phycon.te01*phycon.te001))*phycon.te0004);
	}
	else
	{
		sviii_cs2P32P1 = (realnum) (8.1555/(phycon.te20*phycon.te04*phycon.te004*phycon.te0002));
	}

    return;
}


    




