/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSili compute silicon cooling */
#include "cddefines.h"
#include "taulines.h"
#include "phycon.h"
#include "dense.h"
#include "ligbar.h"
#include "lines_service.h"
#include "colden.h"
#include "embesq.h"
#include "atoms.h"
#include "sil.h"
#include "cooling.h"

void CoolSili(void)
{
	double cs, 
	  cs2s2p, 
	  cs2s3p, 
	  cs01, 
	  cs02, 
	  cs12,
	  tused,
	  temp;
	realnum
	  p2,
	  rate;
	long int i;

	DEBUG_ENTRY( "CoolSili()" );

	/*>>refer	Si I	cs	Hollenbach, D. & McKee, C.F. 1989, ApJ, 342, 306 */
	/* >>chng 03 nov 15, add these lines */
	/* the Si I 25.2 micron line */
	/* rates are said to be ok over range 30 - 3000K */
	tused = MAX2( 30. , phycon.te );
	tused = MIN2( 3000. , phycon.te  );
	tused /= 100.;

	/* derive their rate, then convert to collision strength */
	rate = (realnum)(7.2e-9 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*3.5e-10*pow(tused, -0.03 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		3.5e-10*pow(tused, -0.03 )*dense.xIonDense[ipHYDROGEN][0] );
	LineConvRate2CS( &TauLines[ipSi1_130m]  , rate );

	/* the Si I 56.6 micron line */
	rate = (realnum)(2.2e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*5.0e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		5.0e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	LineConvRate2CS( &TauLines[ipSi1_68m]  , rate );

	rate = (realnum)(7.2e-9 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*1.7e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		1.7e-10*pow(tused, 0.17 )*dense.xIonDense[ipHYDROGEN][0] );
	TauDummy.Hi->g = TauLines[ipSi1_68m].Hi->g;
	LineConvRate2CS( &TauDummy  , rate );
	/* this says that line is a dummy, not real one */
	TauDummy.Hi->g = 0;

	/* solve model atom for Si I */
	atom_level3(&TauLines[ipSi1_130m],&TauLines[ipSi1_68m],&TauDummy);

	/* Si I 2518 */
	MakeCS(&TauLines[ipSii2518]);
	atom_level2(&TauLines[ipSii2518]);

	/* Si I 2215 */
	MakeCS(&TauLines[ipSii2215]);
	atom_level2(&TauLines[ipSii2215]);

	/* Silicon II 35 micron */
	/* hydrogen collision strength from 
	 * >>referold si2	cs	Tielens, A.G.G., & Hollenbach, D. 1985, ApJ, 291, 722
	 * they give rate de-ex 6.5E-10 cm^3 s^-1, indep of temp */
	/*cs += 6.5e-10/dense.cdsqte*4.*dense.xIonDense[ipHYDROGEN][0];*/
	/* >> chng 05 may 21, GS, rate with hydrogen is updated from 2005, ApJ, 620,537*/
	/* >>refer	si2	cs	Barinovs, G., van Hemert, M., Krems, R. & Dalgarno, A. 2005, ApJ, 620, 537 */
	/* original data only extend up to 2000K 
	 * following is valid up to 2000K */
	temp = MIN2(2e3, phycon.te);
	cs = 1e-10*(3.9436853+ 0.11176758*pow(temp, 0.55762129));
	/* for high temperatures simply extend the power law */
	if( phycon.te>2e3 )
	{
		cs *= pow(phycon.te/2e3, 0.55762129);	
	}
	/* above was rate coef, convert to cs and mult by den of colliders */
	cs *= 4.*dense.xIonDense[ipHYDROGEN][0]/dense.cdsqte;

	/* add on elec cs from 
	 *>>refer	si2	cs	Dufton, P.L., & Kingston, A.E. 1994, At. Data Nucl. Data Tables,
	 *>>refercon 57, 273 */
	cs += 5.77;

	PutCS(cs,&TauLines[ipTSi35]);

	static vector< pair<transition*,double> > Si2Pump;
	Si2Pump.reserve(32);

	/* one time initialization if first call */
	if( Si2Pump.empty() )
	{
		// set up level 1 pumping lines
		pair<transition*,double> ppa( &TauLines[ipT1808], 1./6. );
		Si2Pump.push_back( ppa );
		pair<transition*,double> ppb( &TauLines[ipT1527], 2./3. );
		Si2Pump.push_back( ppb );
		pair<transition*,double> ppc( &TauLines[ipT1305], 2./3. );
		Si2Pump.push_back( ppc );
		pair<transition*,double> ppd( &TauLines[ipT1260], 1./6. );
		Si2Pump.push_back( ppd );
		// set up level 2 pumping lines
		for( i=0; i < nWindLine; ++i )
		{
			/* don't test on nelem==ipIRON since lines on physics, not C, scale */
			if( TauLine2[i].Hi->nelem == 14 && TauLine2[i].Hi->IonStg == 2 )
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
				Si2Pump.push_back( pp2 );
			}
		}
	}

	/* now sum pump rates */
	double pump_rate = 0.;
	vector< pair<transition*,double> >::const_iterator si2p;
	for( si2p=Si2Pump.begin(); si2p != Si2Pump.end(); ++si2p )
	{
		const transition* t = si2p->first;
		double branch_ratio = si2p->second;
		pump_rate += t->Emis->pump*branch_ratio;
#		if	0
		dprintf( ioQQQ, "Si II %.3e %.3e\n",
			 t->WLAng , t->Emis->pump*branch_ratio );
#		endif
	}

	/*atom_level2(&TauLines[ipTSi35]);*/
	/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
	/* >>refer	s4	cs	Tayal, S.S., 2000, ApJ 530, 1091*/
	/*>>refer	si2	cs	Dufton, P.L., & Kingston, A.E., 1991, MNRAS, 248, 827*/
	/*>>refer	si2	as	Dufton, P.L., Keenan, F.P., Hibbert, A., 
	 *>>rerercon Stafford, R.P., Byrne, P.B., & Agnew, D., 1991, MNRAS, 253, 474*/
	AtomSeqBoron(&TauLines[ipTSi35], 
		&TauLines[ipSi2_2334], 
		&TauLines[ipSi2_2329], 
		&TauLines[ipSi2_2350], 
		&TauLines[ipSi2_2344], 
		&TauLines[ipSi2_2336], 
		0.534 , 4.51 , 1.67 , 6.94 , 
		pump_rate ,"Si 2");
	/*fprintf(ioQQQ,"DEBUG Si2\t%.2f\t%.5e\t%.5e\t%.5e\n", 
		fnzone, 
		phycon.te, 
		TauLines[ipTSi35].cool, dense.eden);*/
	for( i=0; i < 5; i++ )
	{
		/* pops and column density for SiII atom */
		colden.Si2Pops[i] = (realnum)atoms.PopLevels[i];
	}

	/* Si II 1808, permitted resonance line,
	 * osc str from 
	 * >>refer	si2	as	morton et al 88 (apj sup); 
	 * all si ii collision data (following 4 lines) are from
	 * >>refer	si2	cs	Dufton, P.L., & Kingston, A.E. 1991, MNRAS, 248, 827
	 * following assumes there is typo in table 1 of dufton and kingston
	 * and that they meant 2s 2p^2 ^2D instead of 3d */

	/* Si II 1814 */
	PutCS(13.01,&TauLines[ipT1808]);
	atom_level2(&TauLines[ipT1808]);

	/* Si II 1531 */
	PutCS(3.61,&TauLines[ipT1527]);
	atom_level2(&TauLines[ipT1527]);

	/* Si II 1307.7 */
	PutCS(2.89,&TauLines[ipT1305]);
	atom_level2(&TauLines[ipT1305]);

	/* Si II 1263.3 */
	PutCS(12.25,&TauLines[ipT1260]);
	atom_level2(&TauLines[ipT1260]);

	/* permitted Si III 1206.5, collision strength from 
	 * >>refer	si3	cs	Callaway, J. 1994, At. Data Nucl. Data Tables, 57, 9 */
	cs = MIN2(7.0,1.442*phycon.te10*phycon.te03*phycon.te03/
	  phycon.te01);
	PutCS(cs,&TauLines[ipT1207]);
	atom_level2(&TauLines[ipT1207]);

	/* Si III] 1895, CS=
	 * >>refer	si3	cs	Dufton, P.L., & Kingston, A.E. 1989, MNRAS, 241, 209
	 * >>refer	si3	cs	Dufton, P.L., & Kingston, A.E. 1994, ADNDT, 57, 273
	 * grnd 3s^2 ^1S, upper lev 3p ^3P^o j=0,1,2 */
	/* >>refer	si3	as	Callegari, F., & Trigueiros, A.G., 1998, ApJS, 119, 181
	 * >>chng 00 nov 01, A about 3x larger than before */
	cs = 106./(phycon.te10*phycon.te10*phycon.te10*phycon.te02);
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	PutCS(cs,&TauLines[ipT1895]);
	AtomSeqBeryllium(1.8,3.6,10.4,&TauLines[ipT1895],.013);
	embesq.em1895 = (realnum)(atoms.PopLevels[3]*0.013*1.05e-11);

	/* Si IV 1394, 1403, data from 
	 * >>refer	si4	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon	ed by D.R. Flower, (D. Reidel: Holland), 143
	 * cs from 
	 * >>refer	si4	cs	Dufton, P.L., & Kingston, A.E. 1987, J.Phys. B, 20, 3899 */
	cs = 6.37*phycon.te10;
	PutCS(cs*0.667,&TauLines[ipT1394]);
	PutCS(cs*0.333,&TauLines[ipT1403]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipT1403],&TauDummy,&TauLines[ipT1394]);

	/* Si VI 1.96 micron
	 * >>referold	si6	cs	Saraph, H.E. & Tully, J.A. 1994, A&AS, 107, 29
	 * >>chng 96 jul 16 had been constant */
	/*cs = MIN2(0.43,0.0448*phycon.te20/phycon.te003/phycon.te003);*/
	/*cs = MAX2(0.3,cs);*/
	/* >>refer	si6	cs	Berrington,K.A.,Saraph, H.E. & Tully, J.A. 1998, A&AS, 124, 161*/
	 /* >>chng 06 jul 11-Humeshkar Nemala*/
	if(phycon.te< 1.43E5)
	{
		cs = (realnum)(0.0207*(phycon.te30/phycon.te04)*phycon.te0001);
	}
	else
	{
		cs = (realnum)(3.9042/((phycon.te20/phycon.te02)*phycon.te001*phycon.te0003));
	}
	PutCS(cs,&TauLines[ipSi619]);
	atom_level2(&TauLines[ipSi619]);

	/* Si VII 2148- OIII like, 
	 * >>refer	si7	cs	Kafatos, M., & Lynch, J.P. 1980, ApJS, 42, 611 */
	sil.c2148 = 
		atom_pop2(0.4,9.,5.,15.,6.7e4,dense.xIonDense[13][6])*9.26e-12;
	CoolAdd("Si 7",2148,sil.c2148);

	/* Si VII ground term, 2.48, 6.51 microns
	 * cs 
	 * >>refer	si7	cs	Butler, K., & Zeippen, C.J. 1994, A&AS, 108, 1 */
	/* more recent paper, for solar case, which does not give thermal averaged
	 * collision strengths, is
	 * >>refer	Si7	data	Bhatia, A.K., & Landi, E. 2003, ApJ, 585, 587-597 */
	/** \todo	2	- update to this reference for As
	 * >>refer	Si7	As	Galavis, M.E., Mendoza, C., * Zeippen, C.J. 1997, A&AS, 123, 159 */
	cs = MIN2(0.217,0.0904*phycon.te05*phycon.te03/phycon.te003/
	  phycon.te001);
	PutCS(cs,&TauLines[ipTSi65]);

	cs = MIN2(0.70,8.79e-2*phycon.te10*phycon.te10/phycon.te02);
	PutCS(cs,&TauLines[ipTSi25]);

	cs = MIN2(0.20,9.751e-3*phycon.te20*phycon.te03*phycon.te03/
	  phycon.te003);
	PutCS(cs,&TauDummy);

	atom_level3(&TauLines[ipTSi25],&TauLines[ipTSi65],&TauDummy);

	/* Si 8 1446, 3727-like, 
	 * >>refer	si8	cs	Kafatos, M., & Lynch, J.P. 1980, ApJS, 42, 611 */
	sil.c1446 = 
		atom_pop2(0.4,4.,10.,1.,9.97e4,dense.xIonDense[13][7])*
	  1.39e-11;
	CoolAdd("Si 8",1446,sil.c1446);

	/* Si 9 1985, 2150
	 * cs, As from
	 * >>refer	si9	cs	Aggarwal, K.M. 1983, J.Phys. B, 16, L59
	 * >>refer	si9	as	Baluja, K.L. 1985, J.Phys. B, 18, L413 */
	sil.c949 = 
		atom_pop3(9.,5.,1.,0.5913,0.0757,0.225,26.3,214.,5.16,
	  7.62e4,7.902e4,&p2,dense.xIonDense[13][8],0.,0.,0.)*214.*2.096e-11;
	sil.c1815 = sil.c949*1.912*0.0516;
	sil.c1985 = p2*26.3*1.0e-11;
	CoolAdd("Si 9",949,sil.c949);
	CoolAdd("Si 9",1815,sil.c1815);
	CoolAdd("Si 9",1985,sil.c1985);

	/* Si 9 3P fine structure lines, A=
	 * >>refer	si9	as	Baluja, K.L. 1985, J.Phys. B, 18, L413
	 * 2.583, 3.9microns
	 * CS=
	 * >>refer	si9	cs	Lennon, D.J. Burke, V.M. 1994, A&AS, 103, 273 */
	cs01 = MIN2(0.98,28.51/(phycon.te10*phycon.te10*phycon.te10*
	  phycon.te10/phycon.te03*phycon.te003*phycon.te001*phycon.te001));

	cs12 = MIN2(2.7,81.21/(phycon.te10*phycon.te10*phycon.te10*
	  phycon.te10/phycon.te01/phycon.te01/phycon.te001/phycon.te001));

	cs02 = MIN2(0.70,19.67/(phycon.te10*phycon.te10*phycon.te10*
	  phycon.te10/phycon.te03*phycon.te001));

	PutCS(cs01,&TauLines[ipTSi4]);
	PutCS(cs12,&TauLines[ipTSi3]);
	PutCS(cs02,&TauDummy);

	atom_level3(&TauLines[ipTSi4],&TauLines[ipTSi3],&TauDummy);

	/* 5S0 - 3P, cs from, A=guess
	 * >>refer	si9	cs	Aggarwal, K.M. 1984, ApJS, 54, 1 */
	sil.c691 = atom_pop2(40.6/phycon.sqrte*phycon.te10,9.,5.,
	  1e4,2.081e5,dense.xIonDense[13][8])*2.88e-11;
	CoolAdd("Si 9",691,sil.c691);

	/* Si 10 606, actually three lines clumped together, ll 621.1, 611.7, 598.6
	 * atomic data 
	 * >>refer	si10	cs	Saha, H.P., & Trefftz, E. 1982, A&A, 116, 224 */
	/* >>chng 03 sep 27, rm expion move to simple two level with rt */
	/*CoolHeavy.c606 = 
		0.10*1.42e-16*expion(2.4e5,dense.xIonDense[13][10-1]);
	CoolAdd("Si10",606,CoolHeavy.c606);*/
	PutCS(0.1,&TauLines[ipSi10_606]);
	atom_level2(&TauLines[ipSi10_606]);

	/* Si 10 1.43m, A from 
	 * >>refer	si10	as	Chandra, S. 1982, SoPh, 75, 133
	 * cs from 
	 * >>refer	si10	cs	Zhang, H.L., Graziani, M., Pradhan, A.K. 1994, A&A, 283, 319 */
	if( phycon.te <= 40500. )
	{
		cs = 0.190*phycon.te20/phycon.te001;
	}
	else
	{
		cs = 24.93/(phycon.te20*phycon.te03*phycon.te01*phycon.te003*
		  phycon.te003);
	}
	PutCS(cs,&TauLines[ipSi10143]);
	atom_level2(&TauLines[ipSi10143]);

	/* SI 11 582.9, 1909-LIKE, CS=
	 * >>refer	si11	cs	Berrington, K.A., Burke, P.G., Dufton, P.L., Kingston, A.E.
	 * >>refercon	1985, At. Data Nucl. Data Tables, 33, 195
	 * A=
	 * >>refer	si11	as	Muhlethaler, H.P., & Nussbaumer, H. 1976, A&A 48, 109 */
	sil.c583 = 
		atom_pop2(0.10,1.,9.,1e5,2.47e5,dense.xIonDense[13][11-1])*
	  3.4e-11;
	CoolAdd("Si11",583,sil.c583);

	/* li seq 2s2p and 2s3p, Si 12 499, 521
	 * >>refer	si12	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(14,&TauLines[ipTSi499],&TauLines[ipTSi41],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipTSi499]);
	PutCS(cs2s2p*0.5,&TauLines[ipTSi521]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipTSi521],&TauDummy,&TauLines[ipTSi499]);

	PutCS(cs2s3p,&TauLines[ipTSi41]);
	atom_level2(&TauLines[ipTSi41]);
	return;
}
