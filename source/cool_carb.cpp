/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolCarb evaluate total cooling due to carbon */
#include "cddefines.h"
#include "physconst.h"
#include "embesq.h"
#include "phycon.h"
#include "taulines.h"
#include "dense.h"
#include "hmi.h"
#include "h2.h"
#include "mole.h"
#include "ligbar.h"
#include "thermal.h"
#include "colden.h"
#include "lines_service.h"
#include "atoms.h"
#include "carb.h"
#include "cooling.h"

void CoolCarb(void)
{
	double SaveAbun, 
	  a21, 
	  a31, 
	  a32, 
	  cs, 
	  cs01, 
	  cs02, 
	  cs12, 
	  cs13, 
	  cs23, 
	  cs2s2p, 
	  cs2s3p ,
	  ortho_frac ,
	  popup,
	  popratio;

	/* added to implement Peter van Hoof additions for new ground term
	 * atomic collision data */
	double cse01,
		cse12,
		cse02,
		csh01,
		csh12,
		csh02,
		csp01,
		csp12,
		csp02,
		csh201,
		csh212,
		csh202 ,
		csh2p01,
		csh2p12,
		csh2p02,
		csh2o01,
		csh2o12,
		csh2o02,
		temp;
	double cs_c2_h12=-1.;
	realnum pciexc ,
		  sum;
	int i;

	DEBUG_ENTRY( "CoolCarb()" );

	TauDummy.Zero();
	TauDummy.Hi->g = 0.;
	TauDummy.Lo->g = 0.;
	TauDummy.Hi->IonStg = 0;
	TauDummy.Lo->IonStg = 0;
	TauDummy.Hi->nelem = 0;
	TauDummy.Lo->nelem = 0;
	TauDummy.Emis->Aul = 0.;
	TauDummy.EnergyErg = 0.;
	TauDummy.EnergyK = 0.;
	TauDummy.EnergyWN = 0.;

	/* subroutine atom_level3( t10,t21,t20)
	 *
	 * Carbon cooling
	 *
	 * C I 1656, collision strength from transition prob */
	/*PutCS(7.3,t1656);
	atom_level2(t1656);*/
	PutCS(7.3, &TauLines[ipT1656] );
	atom_level2(&TauLines[ipT1656]);

	/* C I fine structure lines data originally from
	 * >>refer	c1	cs	Tielens, A.G.G., & Hollenbach, D. 1985, ApJ, 291, 722
	 * >>chng 99 jun 01, to more recent ground term collision data
	 * by Peter van Hoof */

	/* effective collision strength of C I(3P) with e
	 * >>refer	c1	cs	Johnson, C.T., Burke, P.G., Kingston, A.E. 1987, JPhB, 20, 2553
	 * these data are valid for 7.5K <= Te <= 10,000K*/
	if( phycon.te<=3.0e3 )
	{
		/* the first fit is valid for 10K <= Te <= 300K, 
		 * the second 300K <= Te <= 3000K*/
		cse01 = MAX2(4.80E-06*phycon.te*phycon.te20/phycon.te03,
			8.24E-07*phycon.te32/phycon.te01);

		cse12 = MAX2(7.67E-05*phycon.te/phycon.te10/phycon.te03,
			1.47E-06*phycon.te32*phycon.te10/phycon.te03);

		cse02 = MAX2(4.72E-05*phycon.te70*phycon.te03,
			3.05E-07*phycon.te32*phycon.te10);
	}
	else
	{
		/* the first fit is valid for 300K <= Te <= 3000K, 
		 * the second up to 10,000K */
		cse01 = MIN2(8.24E-07*phycon.te32/phycon.te01,
			0.0035*phycon.sqrte*phycon.te01);

		cse12 = MIN2(1.47E-06*phycon.te32*phycon.te10/phycon.te03,
			0.0088*phycon.sqrte*phycon.te01*phycon.te005);

		cse02 = MIN2(3.05E-07*phycon.te32*phycon.te10,
			0.00448*phycon.sqrte/phycon.te10*phycon.te03*phycon.te005);
	}

	/* >>chng 04 nov 24, upper limit of 1000K is too low - for low Z DLA clouds we need
	 * C^0 populations at higher temperature - these are simple power laws - extrapolate them
	 * to 3x too high a temp */
	/* rate coefficients for collisional de-excitation of C I(3P) with neutral H(2S1/2)
	 * >>refer	c1	cs	Launay & Roueff 1977, A&A 56, 289
	 * the first fit is for Te <= 100K, the second for Te >= 250K
	 * these data are valid for 4K <= Te <= 1000K*/
	csh01 = MAX2(1.61e-10,5.66e-11*phycon.te20);

	/* these data are valid for 7K <= Te <= 1000K*/
	csh12 = MAX2(1.93e-10*phycon.te05*phycon.te03,
		5.64e-11*phycon.te30*phycon.te02);

	/* these data are valid for 10K <= Te <= 1000K*/
	csh02 = MAX2(1.08e-10/phycon.te03,
		1.67e-11*phycon.te30*phycon.te02*phycon.te02);

	/* >>chng 05 may 23, collisional de-excitations rate co-efficients (cm3s-1) of C I by proton*/
	/*>>refer	C1	cs	Roueff, E. & Le Bourlot, J. 1990, A&A 236, 515
	 * upward rates are given for 100 to 20,000K*/
	/* csp01 starts increasing below 25 K, but csp02 and csp12 behave properly*/
	if( phycon.te < 25. )
		temp = 25.;
	else if( phycon.te >20000. )
		temp = 20000.;
	else
		temp = phycon.te;
	csp01 = 1e-9*pow2(4.3671821 - 14.39018/log(temp))*(1./3.)*exp(16.4*T1CM/temp);
	csp12 = 1e-9*exp(3.2823932 - 60.99754*(log(temp)/temp))*(1./5.)*exp(37.1*T1CM/temp);
	csp02 = 1e-9/(0.033932579+ (1503.4042/pow(temp,1.5)))*(3./5.)*exp(43.5*T1CM/temp);

	/* >>chng 05 feb 03, this logic had set H2 collisions to zero when
	 * temperature was > 1.2e3.  this is unphysical.  change to use
	 * 1200K collision rate at all higher temperatures.  this is a constant
	 * rate, and the original paper suggested that the rate was fairly
	 * constant at higher tabulated temperatures 
	if( phycon.te<=1.2e3 )*/
	/* >>chng 04 mar 15, use explicit ortho-para densities */
	ortho_frac = h2.ortho_density/SDIV(hmi.H2_total);

	/* rate coefficients for collisional de-excitation of C I(3P) with H2(J=1,0)
	 * >>refer	c1	cs	Schroeder et al. 1991, J.Phys.B 24, 2487
	 * these data are valid for 10K <= Te <= 1200K
	 * the first entry is contribution from ortho H2, the second para H2.*/
	if( phycon.te<=30. )
	{
		csh2p01 = MIN2(8.38E-11*phycon.te05*phycon.te01,
			2.12e-10/phycon.te20/phycon.te05/phycon.te01);

		csh2o01 = MIN2(5.17E-11*phycon.te10*phycon.te05,
			1.07e-10/phycon.te10*phycon.te01);
	}
	else if( phycon.te<=150. )
	{
		csh2p01 = MAX2(6.60e-11,
			2.12e-10/phycon.te20/phycon.te05/phycon.te01);

		csh2o01 = MAX2(7.10e-11,
			1.07e-10/phycon.te10*phycon.te01);
	}
	else
	{
		/* this is high temperature branch - increasing function of T,
		 * so hits cap set by min */
		csh2p01 = MAX2(6.60e-11,3.38e-11*phycon.te10*phycon.te03);
		csh2p01 = MIN2(8.10e-11,csh2p01);

		csh2o01 = MAX2(7.1e-11,3.37e-11*phycon.te10*phycon.te02*phycon.te02);
		csh2o01 = MIN2(8.57e-11,csh2o01);
	}

	/* use computed ortho and para H2 densities to get total collision rate */
	csh201 = ortho_frac*csh2o01 + (1.-ortho_frac)*csh2p01;
	if( phycon.te<=30. )
	{
		csh2p12 = MIN2(1.48E-10*phycon.te05*phycon.te02,
			2.25e-10/phycon.te03/phycon.te03);
	}
	else if( phycon.te <= 100. )
	{
		csh2p12 = MAX2(1.75e-10,
			2.25e-10/phycon.te03/phycon.te03);
	}
	else
	{
		csh2p12 = MAX2(1.75e-10,6.23e-11*phycon.te20*phycon.te01);
		csh2p12 = MIN2(2.61e-10,csh2p12);
	}

	csh2o12 = MIN2(2.83e-10,4.46e-11*phycon.te30/phycon.te03);
	{
		/*csh212 = 0.75*csh2o12 + 0.25*csh2p12;*/
		csh212 = ortho_frac*csh2o12 + (1.-ortho_frac)*csh2p12;
	}

	if( phycon.te<=30 )
	{
		csh2p02 = MIN2(8.67E-11*phycon.te02*phycon.te02,
			1.35e-10/phycon.te10);
	}
	else if( phycon.te<=150. )
	{
		csh2p02 = MAX2(8.40e-11,
			1.35e-10/phycon.te10);
	}
	else
	{
		csh2p02 = MAX2(8.4e-11,4.04e-11*phycon.te10*phycon.te02*phycon.te02);
		csh2p02 = MIN2(1.04e-10,csh2p02);
	}

	csh2o02 = MIN2(1.11e-10,3.16e-11*phycon.te20/phycon.te02);
	/*csh202 = 0.75*csh2o02 + 0.25*csh2p02;*/
	csh202 = ortho_frac*csh2o02 + (1.-ortho_frac)*csh2p02;

	/** \todo	2	add term for protons from Rouef, E., & Le Bourlot, J. 1990, A&A, 236, 515 */
	/** \todo	1	add neutral helium Staemmler, V., & Flower, D. R. 1991, J. Phys. B, 24, 2343 */
	/* assume CS for He^0 is the same as H^0*/
	/*cs01 = cse01 + 3.*(csh01*(dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHELIUM][0]) + csh201*hmi.Hmolec[ipMH2g])/dense.cdsqte;
	cs12 = cse12 + 5.*(csh12*(dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHELIUM][0]) + csh212*hmi.Hmolec[ipMH2g])/dense.cdsqte;
	cs02 = cse02 + 5.*(csh02*(dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHELIUM][0]) + csh202*hmi.Hmolec[ipMH2g])/dense.cdsqte;*/
	cs01 = cse01 + 3.*(csh01*(dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHELIUM][0]) + csp01*dense.xIonDense[ipHYDROGEN][1] + csh201*hmi.H2_total)/dense.cdsqte;
	cs12 = cse12 + 5.*(csh12*(dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHELIUM][0]) + csp12*dense.xIonDense[ipHYDROGEN][1] + csh212*hmi.H2_total)/dense.cdsqte;
	cs02 = cse02 + 5.*(csh02*(dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHELIUM][0]) + csp02*dense.xIonDense[ipHYDROGEN][1] + csh202*hmi.H2_total)/dense.cdsqte;

	PutCS( cs01 , &TauLines[ipT610] );
	PutCS( cs12 , &TauLines[ipT370] );
	PutCS( cs02 , &TauDummy );

	/* ======================================================== */
	/* end changes 99 Jun 01, by Peter van Hoof */
	atom_level3( 
		&TauLines[ipT610],
		&TauLines[ipT370],
		&TauDummy);

	/* now save pops to add col den in radinc */
	for( i=0; i<3; ++i)
	{
		/* >>chng 02 oct 23, bug - had been C1Colden rather than C1Pops */
		colden.C1Pops[i] = (realnum)atoms.PopLevels[i];
	}

	/* C I 9850, 8727, A from 
	 * >>refer	c1	as	Mendoza, C. 1982, in Planetary Nebulae, IAU Symp No. 103,
	 * >>refercon ed by D.R. Flower, (D. Reidel: Holland), 143 */
	if( dense.xIonDense[ipCARBON][0] > 0. && phycon.te < 40000. )
	{
		cs12 = 1.156e-4*phycon.te*(1.09 - 7.5e-6*phycon.te - 2.1e-10*
		  phycon.te*phycon.te);
		cs13 = 2.8e-3*phycon.sqrte;
		cs23 = 2.764e-3*phycon.sqrte;

		a21 = 3.26e-4*TauLines[ipT9830].Emis->Pesc;
		a31 = 2.73e-3;
		a32 = 0.528*TauLines[ipT8727].Emis->Pesc;
		/** \todo	3	change to atom_level3 */
		carb.c8727 = atom_pop3(9.,5.,1.,cs12,cs13,cs23,a21,a31,a32,
		  1.417e4,1.255e4,&pciexc,dense.xIonDense[ipCARBON][0],0.,0.,0.)*a32*
		  2.28e-12;
		TauLines[ipT9830].Emis->PopOpc = dense.xIonDense[ipCARBON][0];
		TauLines[ipT9830].Lo->Pop = dense.xIonDense[ipCARBON][0];
		TauLines[ipT9830].Hi->Pop = 0.;
		TauLines[ipT9830].Coll.col_str = (realnum)cs12;
		TauLines[ipT8727].Emis->PopOpc = (carb.c8727/(a32*2.28e-12));
		TauLines[ipT8727].Lo->Pop = (carb.c8727/(a32*2.28e-12));
		TauLines[ipT8727].Hi->Pop = 0.;
		TauLines[ipT8727].Coll.col_str = (realnum)cs23;

		carb.c9850 = pciexc*a21*2.02e-12;
		thermal.dCooldT += carb.c9850*(1.468e4*thermal.tsq1 + thermal.halfte);
		thermal.dCooldT += carb.c8727*(1.255e4*thermal.tsq1 + thermal.halfte);

		/* C I 9850 correction for deexcitation, needed for rec line */
		carb.r9850 = (realnum)(a21/(a21 + cs12/5.*COLL_CONST/phycon.sqrte*dense.eden));
	}

	else
	{
		carb.c9850 = 0.;
		carb.c8727 = 0.;
		carb.r9850 = 0.;
		TauLines[ipT9830].Emis->PopOpc = 0.;
		TauLines[ipT9830].Lo->Pop = 0.;
		TauLines[ipT9830].Hi->Pop = 0.;
		TauLines[ipT8727].Emis->PopOpc = 0.;
		TauLines[ipT8727].Lo->Pop = 0.;
		TauLines[ipT8727].Hi->Pop = 0.;
	}
	CoolAdd("C  1",8727,carb.c8727);
	CoolAdd("C  1",9850,carb.c9850);

	/* C II 158 micron emission, A=
	 * >>refer	c2	as	Froese Fischer, C. 1983, J.Phys. B, 16, 157
	 * CS From 
	 * >>refer	c2	cs	Blum, R.D., & Pradhan, A.K. 1992, ApJS 80, 425
	 * neutral collision data from 
	 * >>refer	c2	cs	Tielens, A.G.G., & Hollenbach, D. 1985, ApJ, 291, 722
	 * >>chng 96 aug 01, better fit to cs  */
	/* following is a more recent calculation but without extensive tables */
	/* >>refer	c2	cs	Wilson, N.J., & Bell, K.L. 2002, MNRAS, 337, 1027-1034 */
	/* >>chng 03 feb 24, break apart electron and neutral hydrogen for book keeping*/
	/*cs = MIN2(2.20,0.403*phycon.te20/phycon.te02*phycon.te001*phycon.te001) + 
	  5.8e-10*phycon.te02/dense.cdsqte*4.*(dense.xIonDense[ipHYDROGEN][0] + 
	  hmi.Hmolec[ipMH2g]);*/
	/* electron collision strength */
	cse12 = MIN2(2.20,0.403*phycon.te20/phycon.te02*phycon.te001*phycon.te001);

	/* atomic hydrogen collision strength, include H2 with same rate */
	/* >>referold	c2	cs	Tielens, A.G.G., & Hollenbach, D. 1985, ApJ, 291, 722 */
	/*cs_c2_h12 = 5.8e-10*phycon.te02/dense.cdsqte*4.*(dense.xIonDense[ipHYDROGEN][0] + 
	  hmi.Hmolec[ipMH2g]);*/
	/* >> chng 05 may 21, GS, rate with hydrogen is updated from following */
	/* >>refer	c2	cs	Barinovs, G., van Hemert, M., Krems, R. & Dalgarno, A. 2005, ApJ, 620, 537 */
	temp = MIN2(2e3, phycon.te);

	/* evaluate the rate at the temperature set above, if temperature is above 2000 K
	 * then it is evaluated at 2000K - first line is just rate as given in paper */
	cs_c2_h12 = 1e-10*(4.4716028+ 0.69658785*pow(temp, 0.31692387));

	if(phycon.te > 2e3) 
	{
		/* temperature is above 2000 K so extrapolate rate as a power law */
		cs_c2_h12 *= pow(phycon.te/2e3, 0.31692387);
	}

	/* now convert rate into equivalent cs */
	cs_c2_h12 *= 4.*(dense.xIonDense[ipHYDROGEN][0] + hmi.Hmolec[ipMH2g])/dense.cdsqte;

	/* >>chng 05 apr 10, make sure we have good current set of vars */
	ASSERT( fabs(dense.eden + dense.xIonDense[ipHYDROGEN][0]*1.7e-4 * dense.HCorrFac -dense.EdenHCorr )/ 
		dense.EdenHCorr < 1e-8 );

	PutCS( cse12+cs_c2_h12 ,&TauLines[ipT157]);

	/* CII 1335 all collision strengths and A'S from 
	 * >>refer	c2	cs	Lennon, D.J., Dufton, P.L., Hibbert, A., Kingston, A.E. 1985, ApJ, 294, 200
	 * >>refer	c2	cs	Blum, R.D., & Pradhan, A.K. 1992, ApJS 80, 425 */
	cs = MIN2(6.73,2.316*phycon.te10);
	PutCS(cs,&TauLines[ipT1335]);
	atom_level2(&TauLines[ipT1335]);

	static vector< pair<transition*,double> > C2Pump;
	C2Pump.reserve(32);

	/* one time initialization if first call */
	if( C2Pump.empty() )
	{
		// set up level 1 pumping lines
		pair<transition*,double> pp( &TauLines[ipT1335], 1./6. ); 
		C2Pump.push_back( pp );
		// set up level 2 pumping lines
		for( i=0; i < nWindLine; ++i )
		{
			/* don't test on nelem==ipIRON since lines on physics, not C, scale */
			if( TauLine2[i].Hi->nelem == 6 && TauLine2[i].Hi->IonStg == 2 )
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
				C2Pump.push_back( pp2 );
			}
		}
	}

	/* now sum pump rates */
	double pump_rate = 0.;
	vector< pair<transition*,double> >::const_iterator c2p;
	for( c2p=C2Pump.begin(); c2p != C2Pump.end(); ++c2p )
	{
		const transition* t = c2p->first;
		double branch_ratio = c2p->second;
		pump_rate += t->Emis->pump*branch_ratio;
#		if	0
		dprintf( ioQQQ, "C II %.3e %.3e\n",
			 t->WLAng , t->Emis->pump*branch_ratio );
#		endif
	}

	/*atom_level2(&TauLines[ipT157]);*/
	/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
	/* >>refer	c2	cs	Blum, R.D., & Pradhan, A.K., 1992, ApJS 80, 425
	 * >>refer	c2	cs	Lennon, D.J., Dufton, P.L., Hibbert, A., Kingston, A.E. 1985, ApJ, 294, 200*/
	/* >>refer	c2	as	Nahar, S.N. 2003, ADNDT, 80, 205 */
	AtomSeqBoron(&TauLines[ipT157], 
	  &TauLines[ipC2_2325], 
	  &TauLines[ipC2_2324], 
	  &TauLines[ipC2_2329], 
	  &TauLines[ipC2_2328], 
	  &TauLines[ipC2_2327], 
	  0.2349 , 0.8237 , 0.8533 , 1.9818 , pump_rate , "C  2");

	/* following should be set true to print contributors */
	enum {DEBUG_LOC=false};
	if( DEBUG_LOC && nzone > 80 )
	{
		fprintf(ioQQQ,"DEBUG\t%.2f\t%.3e\t%.3e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			fnzone , 
			phycon.te, 
			TauLines[ipT157].Coll.cool , 
			cse12,
			csh12,
			dense.eden,
			(dense.xIonDense[ipHYDROGEN][0] + hmi.Hmolec[ipMH2g])/dense.cdsqte);
	}

	sum = 0.;
	/* now save pops to add col den in radinc */
	for( i=0; i<5; ++i)
	{
		colden.C2Pops[i] = (realnum)atoms.PopLevels[i];
		sum += colden.C2Pops[i];
	}
	ASSERT( fabs(sum-dense.xIonDense[ipCARBON][1])/SDIV(dense.xIonDense[ipCARBON][1]) < 1e-4 );

	/* following used for pumping - cs just made up - no real data */
	PutCS(.1,&TauLines[ipT386]);
	atom_level2(&TauLines[ipT386]);

	PutCS(.1,&TauLines[ipT310]);
	atom_level2(&TauLines[ipT310]);

	PutCS(.1,&TauLines[ipT291]);
	atom_level2(&TauLines[ipT291]);

	PutCS(.1,&TauLines[ipT280]);
	atom_level2(&TauLines[ipT280]);

	PutCS(.1,&TauLines[ipT274]);
	atom_level2(&TauLines[ipT274]);

	PutCS(.1,&TauLines[ipT270]);
	atom_level2(&TauLines[ipT270]);

	/* C III  1909
	 * A for 1909 itself from 
	 * >>refer	c3	as	Kwong, V., Fang, Z., Gibbons, T.T., Parkinson, W.H., Smith, P.L.
	 * >>refercon 1993, ApJ, 411, 431
	 * experimental value of 121 is larger than old NS 96, cs from
	 * >>refer	c3	cs	Berrington, K.A., Burke, P.G., Dufton, P.L., Kingston, A.E. 1985,
	 * >>refercon	At. Data Nucl. Data Tables, 33, 195
	 * AtomSeqBeryllium(CS23,CS24,CS34,tarray,A41) */
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	/* >>refer	c2	as	Nahar, S.N. 2003, ADNDT, 80, 205 */
	cs = MIN2(1.1,2.67/phycon.te10);
	a21 = 5.149e-3;
	PutCS(cs,&TauLines[ipT1909]);
	/* C1909 = AtomSeqBeryllium(.96,.73,2.8 , T1909 ,5.19E-3 )
	 * A's 
	 * >>refer	c3	as	Fleming, J., Bell, K.L, Hibbert, A., Vaeck, N., Godefroid, M.R.
	 * >>refercon 1996, MNRAS, 279, 1289 */
	AtomSeqBeryllium(.96,.73,2.8,&TauLines[ipT1909],a21);
	embesq.em1908 = (realnum)(atoms.PopLevels[3]*a21*1.04e-11);
	/*DumpLine(&TauLines[ipT1909]);*/

	/* >>chng 02 mar 08, add 13C line - this is totally forbidden for 12C
	 * and so provides a mathod of deducing 13C/12C */
	/* >>refer	C3	13C As	Clegg, R.E.S., Storey, P.J., Walsh, J.R., & Neale, L.
	 * >>refercon	1997, MNRAS, 284, 348 */
	a21 = 6.87e-4;
	/* this is the correction for depopulation of the P_0 level due to A21, which is no
	 * present in 12C */
	cs = 2.8*dense.cdsqte/5.*1.667;
	popratio = 	cs/(cs + a21);
	embesq.em13C1910 = (realnum)(a21 * atoms.PopLevels[1]*popratio* 1.04e-11 / co.C12_C13_isotope_ratio);

	/* CIII 1175 excited state line 
	 * following were computed by previous call to AtomSeqBeryllium  */
	/*popup = atoms.PopLevels[1] + atoms.PopLevels[2] + atoms.PopLevels[3];*/
	popup = 0.;
	colden.C3Pops[0] = (realnum)atoms.PopLevels[0];
	for( i=1; i<4; ++i)
	{
		popup += atoms.PopLevels[i];
		colden.C3Pops[i] = (realnum)atoms.PopLevels[i];
	}

	SaveAbun = dense.xIonDense[ipCARBON][2];
	dense.xIonDense[ipCARBON][2] = (realnum)popup;
	/* cs 
	 * >>refer	c3	cs	Berrington, K.A., Burke, P.G., Dufton, P.L., Kingston, A.E. 1985,
	 * >>refercon At. Data Nucl. Data Tables, 33, 195 */
	cs = MIN2(30.,4.806*phycon.te10*phycon.te05/phycon.te01/phycon.te003);
	PutCS(18.45,&TauLines[ipc31175]);
	atom_level2(&TauLines[ipc31175]);
	dense.xIonDense[ipCARBON][2] = (realnum)SaveAbun;

	/* C III 977, cs from 
	 * >>refer	c3	cs	Berrington, K.A. 1985, J.Phys. B, 18, L395 */
	cs = MIN2(7.0,1.556*phycon.te10);
	PutCS(cs,&TauLines[ipT977]);
	atom_level2(&TauLines[ipT977]);

	/* CIV 1548, 1550 doublet
	 * >>refer	c4	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(
		6,
		&TauLines[ipT1548],
		&TauLines[ipT312],
		&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipT1548]);
	PutCS(cs2s2p*0.5,&TauLines[ipT1550]);
	PutCS(1.0,&TauDummy);
	atom_level3(
		&TauLines[ipT1550],
		&TauDummy,
		&TauLines[ipT1548]);

	PutCS(cs2s3p,&TauLines[ipT312]);
	atom_level2(&TauLines[ipT312]);
	return;
}
