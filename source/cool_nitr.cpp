/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolNitr evaluate total cooling due to nitrogen */
#include "cddefines.h"
#include "taulines.h"
#include "dense.h"
#include "phycon.h"
#include "ligbar.h"
#include "thermal.h"
#include "lines_service.h"
#include "embesq.h"
#include "atoms.h"
#include "cooling.h"
#include "nitro.h"
/* xNI_coll_stren: Author: Terry Yun 2006 */
/* This routine interpolates the collision strength of NI 
 * We are only interested in the first 5 levels, therefore there are 10 transitions.
 * Used polynomial fit, eqn: CS = a*t+b*t^1.5+c*t^0.5
 * The range of temp is 0 - 50,000 K */
/* #define PRINT */
static const int N1_SIZE = 10;

static const double aNI[N1_SIZE] = {2.755e-5,4.123e-5,7.536e-6,1.486e-5,4.516e-5,
				    -2.935e-6,4.000e-6,3.751e-6,-2.176e-6,1.024e-5};
static const double bNI[N1_SIZE] = {-8.150e-8,-1.220e-7,-2.226e-8,-4.390e-8,-1.130e-7,
				    8.000e-9,-1.1447e-8,-1.061e-8,5.610e-9,-3.227e-8};
static const double cNI[N1_SIZE] = {2.140e-4,3.272e-4,-4.944e-6,3.473e-6,-8.772e-4,
				    1.654e-3,1.675e-3,1.123e-3,3.867e-3,3.376e-4};

static double NI[6][6][3]; /* coefficients a b c into array */

/*xNI_coll_stren interpolate the collision strength of NI  */
STATIC double xNI_coll_stren(int ipLo, int ipHi)
{
	/* incorporating the N0 collision strengths give in
	 *>>refer	N1	cs	Tayal, S.S., 2006, ApJS, 163, 207 */

	double CS; /* collision strength */

	// Tayal's energy index IS NOT IN INCREASING ENERGY ORDER
	// call routine with proper energy order, swap it around
	// to Tayal order
	int ipLoTayal=-1, ipHiTayal=-1;

	/* Invalid entries returns '-1':the initial indices are smaller than the final indices */
	if(ipLo >= ipHi)
		return( -1. );

	/* Invalid returns '-1': the indices are greater than 5 or smaller than 0 */
	if(ipLo < 1 || ipHi > 5)
		return( -1. );

	int ipTayalOrder[6]={0,1,3,2,4,5};
	ipLoTayal = ipTayalOrder[ipLo];
	ipHiTayal = ipTayalOrder[ipHi];
	// special case of 2-3, must swap order again
	if( ipLo==2 && ipHi==3 )
	{
		ipLoTayal = 2;
		ipHiTayal = 3;
	}

	static bool lgMustInit = true;
	if( lgMustInit )
	{
		lgMustInit = false;
		/* assigning array initially to zero */ 
		for(long i=0; i<6; i++)
		{
			for(long j=0; j<6; j++)
			{
				set_NaN(NI[i][j][0]); /* coeff a */
				set_NaN(NI[i][j][1]); /* coeff b */
				set_NaN(NI[i][j][2]); /* coeff c */
			}
		}

		/* reading coeffs into 3D array */   
		/* The index is in physics scale */
		int index = 0;
		for(int i=1; i<6; i++)
		{
			for(int j=i+1; j<6; j++)
			{
				NI[i][j][0] = aNI[index];
				NI[i][j][1] = bNI[index];
				NI[i][j][2] = cNI[index];
				index++;
			}
		}
	}

#	ifdef PRINT
	/* physical index goes from 1 to 5 */
	for(long i=1; i<6; i++)
	{
		for(long j=i+1; j<6; j++)
		{
			printf("NI %i%i a:%e ", i, j, NI[i][j][0]);
			printf("%i%i b:%e\n", i, j, NI[i][j][1]);
			/* printf("%i%i c:%lf\n", i, j, NI[i][j][2]); */
		}
	}
#	endif

	double temp_max = 50e3;
	double temp_min = 0;
	double temp = phycon.te;
	temp = MAX2(temp, temp_min); /* return lager temp */
	temp = MIN2(temp, temp_max); /* return smaller temp */
	/* eqn: CS = a*t+b*t^1.5+c*t^0.5 */
	CS = NI[ipLoTayal][ipHiTayal][0]*phycon.te + NI[ipLoTayal][ipHiTayal][1]*phycon.te32 + 
		NI[ipLoTayal][ipHiTayal][2]*phycon.sqrte;

	return CS;
}


void CoolNitr(void)
{
	realnum p2;
	double a21, 
	  a31, 
	  a32, 
	  cs, 
	  cs2s2p, 
	  cs2s3p, 
	  cs21, 
	  cs31, 
	  cs32, 
	  p3;
	long int i;

	double a12, 
	  a13, 
	  a14, 
	  a15, 
	  a23, 
	  a24, 
	  a25, 
	  a34, 
	  a35, 
	  a45, 
	  cs12, 
	  cs13, 
	  cs14, 
	  cs15, 
	  cs23, 
	  cs24, 
	  cs25, 
	  cs34, 
	  cs35, 
	  cs45;

	double pop[5];

	DEBUG_ENTRY( "CoolNitr()" );

	static const double gN1[5]={4.,6.,4.,2.,4.};
	static const double exN1[4]={19224.464,8.713,9605.74,0.386};

	/* trans prob from 
	 * >>refer	n1	as	Froese Fischer, C. & Tachiev, G. 2004, ADNDT, 87, 1
	 */
	// 5200, g_up = 6, cs=0.337 @ 1e4K in Tayal06, 
	// NB - his J levels in ^2D are not in energy order
	cs12 = xNI_coll_stren(1, 2);
	a12 = 7.57e-6;
	double a12Tot = a12 + atoms.d5200r;

	// 5198, g_up = 4, cs=0.224 @ 1e4K in Tayal06
	// NB - his J levels in ^2D are not in energy order
	cs13 = xNI_coll_stren(1, 3);
	a13 = 2.03e-5;
	double a13Tot = a13 + atoms.d5200r;

	// 3466.543 0.055 @ 1e4K, g_up = 2
	cs14 = xNI_coll_stren(1, 4);
	a14 = 2.61e-3;

	// 3466.497 0.109 @ 1e4K, g_up = 4
	cs15 = xNI_coll_stren(1,5);
	a15 = 6.50e-3;

	/* this is fraction of 2D excitations that emit a photon in the 5200 multiplet
	 * this is used for collisional quenching in prt lines routine
	 * true line emission is in numerator while all loss terms are in
	 * denominator */
	nitro.quench_5200 = (a12+a13) / ((a12Tot+a13Tot) + (cs12 + cs13) * dense.cdsqte);

	// 0.257 @ 1e4K, g_up = 4
	cs23 = xNI_coll_stren(2, 3);
	a23 = 1.07e-8;

	// 10398.2 0.139 @ 1e4K, g_up = 2
	cs24 = xNI_coll_stren(2, 4);
	a24 = 3.45e-2;

	// 10397.7 0.366 @ 1e4K, g_up = 4
	cs25 = xNI_coll_stren(2, 5);
	a25 = 6.14e-2;

	// 10407.6 0.141 @ 1e4K, g_up = 2
	cs34 = xNI_coll_stren(3, 4);
	a34 = 5.27e-2;

	// 10407.2 0.195 @ 1e4K, g_up = 4
	cs35 = xNI_coll_stren(3, 5);
	a35 = 2.75e-2;

	// 0.123 @ 1e4K, g_up = 4
	cs45 = xNI_coll_stren(4, 5);
	a45 = 9.50e-17;

	// FUV allowed and intercombination lines that pump [NI]
	// the branching ratios used below are based on
	// >>refer	n1	as	Froese Fischer, C. & Tachiev, G. 2004, ADNDT, 87, 1
	//
	// the colission strengths are taken from
	// >>refer	n1	cs	Tayal, S. S., 2006, ApJS, 163, 207
	// the fits are based on data between 500K and 50 kK
	//
	// update rates for the driving transitions
	cs = exp(-11.3423 + 0.8379*min(phycon.alnte,10.82));
	PutCS( cs, &TauLines[ipNI_pumpIndirect] );
	atom_level2( &TauLines[ipNI_pumpIndirect] );

	static const double csN1_a[NI_NDP] = {
		-12.3982, -9.4523, -12.5580, -10.8813, -13.6532, -9.9035, -11.4470, -5.4776, -6.3304
	};
	static const double csN1_b[NI_NDP] = {
		0.7458, 0.3865, 0.7330, 0.6853, 0.7712, 0.3919, 0.6734, 0.1789, 0.1966
	};

  	for( i=0; i < NI_NDP; ++i )
	{
		cs = exp(csN1_a[i] + csN1_b[i]*min(phycon.alnte,10.82));
		PutCS( cs, &TauLines[ipNI_pumpDirect[i]] );
		atom_level2( &TauLines[ipNI_pumpDirect[i]] );
	}
	// units of pump rate s-1
	double RPI = TauLines[ipNI_pumpIndirect].Emis->pump;
	double RPD0 = TauLines[ipNI_pumpDirect[0]].Emis->pump;
	double RPD1 = TauLines[ipNI_pumpDirect[1]].Emis->pump;
	double RPD2 = TauLines[ipNI_pumpDirect[2]].Emis->pump;
	double RPD3 = TauLines[ipNI_pumpDirect[3]].Emis->pump;
	double RPD4 = TauLines[ipNI_pumpDirect[4]].Emis->pump;
	double RPD5 = TauLines[ipNI_pumpDirect[5]].Emis->pump;
	double RPD6 = TauLines[ipNI_pumpDirect[6]].Emis->pump;
	double RPD7 = TauLines[ipNI_pumpDirect[7]].Emis->pump;
	double RPD8 = TauLines[ipNI_pumpDirect[8]].Emis->pump;
	// correct for removal of photons that go straight back to ground
	// beta is the branching ratio directly back to ground
	// for pumping lines not listed below, beta is less than 0.01
	double beta = 0.7955;
	double Premove = TauLines[ipNI_pumpIndirect].Emis->Pesc +
		TauLines[ipNI_pumpIndirect].Emis->Pelec_esc + TauLines[ipNI_pumpIndirect].Emis->Pdest;
	RPI *= (1.-beta)/(1.-beta*(1.-Premove));
	beta = 0.1384;
	Premove = TauLines[ipNI_pumpDirect[0]].Emis->Pesc +
		TauLines[ipNI_pumpDirect[0]].Emis->Pelec_esc + TauLines[ipNI_pumpDirect[0]].Emis->Pdest;
	RPD0 *= (1.-beta)/(1.-beta*(1.-Premove));
	// the branching ratios to the lowest 4 excited levels are calculated in the low-density
	// limit (i.e. no collisions are treated) and exclude transitions straight back to ground
	// from the upper level of the driving line
	double pump14 = 0.0113*RPI + 0.0239*RPD0 + 0.0090*RPD1 + 0.1617*RPD2 + 0.0167*RPD3 + 0.4404*RPD4;
	double pump15 = 0.0112*RPI + 0.0265*RPD0 + 0.6253*RPD1 + 0.5108*RPD2 + 0.0824*RPD3 + 0.2588*RPD4;
	double pump13 = 0.3441*RPI + 0.8621*RPD0 + 0.0233*RPD1 + 0.0895*RPD2 + 0.1068*RPD3 + 0.1644*RPD4;
	double pump12 = 0.0417*RPI + 0.0468*RPD0 + 0.3408*RPD1 + 0.2328*RPD2 + 0.7937*RPD3 + 0.1338*RPD4;

	pump14 += 0.4881*RPD5 + 0.0876*RPD6 + 0.0450*RPD7 + 0.1777*RPD8;
	pump15 += 0.1569*RPD5 + 0.0484*RPD6 + 0.2240*RPD7 + 0.0854*RPD8;
	pump13 += 0.2908*RPD5 + 0.8397*RPD6 + 0.0694*RPD7 + 0.7369*RPD8;
	pump12 += 0.0623*RPD5 + 0.0238*RPD6 + 0.6615*RPD7 + 0.0000*RPD8;

	// estimate of pumping contribution to intensity of [NI] 5199
	nitro.pump5199 = (pump12+pump13+0.9710*pump14+0.9318*pump15) * nitro.quench_5200 * 
		dense.xIonDense[ipNITROGEN][0] * 3.82e-12;

	/**** Terry's addition, recombination from N+ **************/
	/* rate coefficient (cm3 s-1) from Table 3 of
	 *>>refer	NI	rec	Pequignot, D., Petijean, P. & Boisson, C. 1991, A&A, 251, 680 */
	double eff_recrate_2D = 1.108e-13 * pow(phycon.te*1e-4, -0.6085) /
		(1. - 0.0041 * pow(phycon.te*1e-4, -0.3975) );
	double eff_recrate_2P = 0.659e-13 * pow(phycon.te*1e-4, -0.6158);
	double fac_n1;
	if( dense.xIonDense[ipNITROGEN][0] > 0. )
		fac_n1 = dense.eden * dense.xIonDense[ipNITROGEN][1] / dense.xIonDense[ipNITROGEN][0];
	else
		fac_n1 = 0.;

	// assume levels are populated according to statistical weight
	double rec14 = eff_recrate_2P * fac_n1 * 2./6.;
	double rec15 = eff_recrate_2P * fac_n1 * 4./6.;
	double rec13 = eff_recrate_2D * fac_n1 * 4./10.;
	double rec12 = eff_recrate_2D * fac_n1 * 6./10.;

	// estimate of recombination contribution to intensity of [NI] 5199
	nitro.rec5199 = (rec12+rec13+0.9710*rec14+0.9318*rec15) * nitro.quench_5200 * 
		dense.xIonDense[ipNITROGEN][0] * 3.82e-12;

	pump14 += rec14;
	pump15 += rec15;
	pump13 += rec13;
	pump12 += rec12;

	double Cooling , CoolingDeriv;
	atom_pop5(gN1,exN1,cs12,cs13,cs14,cs15,cs23,cs24,cs25,cs34,cs35,
		cs45,a12Tot,a13Tot,a14,a15,a23,a24,a25,a34,a35,a45,pop,
		dense.xIonDense[ipNITROGEN][0],	&Cooling , &CoolingDeriv ,
		pump12 , pump13 , pump14 , pump15 );

	nitro.xN5200 = pop[1]*a12*3.818e-12;
	nitro.xN5198 = pop[2]*a13*3.821e-12;
	nitro.xN3467 = pop[3]*a14*5.729e-12;
	nitro.xN3466 = pop[4]*a15*5.729e-12;
	nitro.xN10398 = pop[3]*a24*1.910e-12;
	nitro.xN10397 = pop[4]*a25*1.910e-12;
	nitro.xN10408 = pop[3]*a34*1.908e-12;
	nitro.xN10407 = pop[4]*a35*1.908e-12;

	// population of excited NI - used in ionization balance routines
	atoms.p2nit = (realnum)(pop[1]+pop[2])/SDIV(dense.xIonDense[ipNITROGEN][0]);

	// total cooling from 5-level system
	thermal.dCooldT += CoolingDeriv;
	CoolAdd("N  1",5199, Cooling );

	/* N I 1200, cs from trans probability */
	PutCS(4.1,&TauLines[ipT1200]);
	atom_level2(&TauLines[ipT1200]);

	/* N II 1084, cs from trans prob */
	PutCS(5.5,&TauLines[ipT1085]);
	atom_level2(&TauLines[ipT1085]);

	/* [N II] coll data from 
	 * >>referold	n2	cs	Stafford, R.P., Bell, K.L, Hibbert, A. & Wijesundera, W.P.,
	 * >>rereroldcon 1994, MNRAS 268, 816,
	 * at 10,000K (v weak T dep)
	 * >>chng 00 dec 11, to Lennon & Burke
	 * >>referold	n2	cs	Lennon, D.J., & Burke, V.M., 1994, A&AS, 103, 273-277
	 * transit prob from 
	 * >>referold	n2	as	Nussbaumer, H., & Rusca, C. 1979, A&A, 72, 129
	 * >>chng 00 dec 11, to Lennon & Burke cs
	 *>>refer	n2	as	Froese Fischer, C., & Tachiev, G. 2004, At. Data Nucl. Data Tables, 87, 1
	 */
	a21 = 3.889e-3;
	a31 = 0.03140;
	a32 = 1.136;
	/* >>chng 02 may 02, put in option to switch between Lennon Burke and Stafford et al. ,
	 * default state of this var is true */
	/** \todo	1	update cs these to following reference:
	>>refer	n2	cs	Hudson, C.E. & Bell, K.L. 2004, MNRAS, 348, 1275 and A&A, 430, 725 
	 * they agree with Lennon & Burke
	 * >>chng 10 feb 24 ML: cs values updated with Hudson & Bell. Not sure why they were wrong.*/
	cs21 = 2.722;
	cs31 = 0.3162;
	cs32 = 0.806;

	/* POP3(   G1,G2,G3,   O12,  O13,  O23,  A21, A31, A32,*/
	p3 = atom_pop3(  9.,5.,1.,  cs21, cs31, cs32, a21, a31, a32 ,
	  /*   E12,E23,P2,ABUND,GAM2) */
	  21955.,24982.,&p2,dense.xIonDense[ipNITROGEN][1],0.,0.,0.);

	nitro.c5755 = p3*a32*3.46e-12;

	nitro.c6584 = p2*a21*3.03e-12;
	thermal.dCooldT += nitro.c6584*(2.2e4*thermal.tsq1 - thermal.halfte);
	CoolAdd("N  2",5755,nitro.c5755);
	CoolAdd("N  2",6584,nitro.c6584);

	/* nitro.xN2_A3_tot is fraction of excitations that produce a photon 
	 * and represents the correction for collisiona deexcitation */
	nitro.xN2_A3_tot = (a31+a32) /(a31+a32 + (cs31+cs32)/1.*dense.cdsqte );
	ASSERT( nitro.xN2_A3_tot <= 1. );

	/* N II fine structure lines, */
	/* >>chng 00 dec 11, to Lennon & Burke CS
	* >>refer	n2	cs	Hudson, C.E. & Bell, K.L. 2004, MNRAS, 348, 1275 and A&A, 430, 725
	 */
	cs21 = 0.431;
	cs32 = 1.15;
	cs31 = 0.273;

	PutCS(cs21,&TauLines[ipT205]);
	PutCS(cs32,&TauLines[ipT122]);
	PutCS(cs31,&TauDummy);
	atom_level3(&TauLines[ipT205],&TauLines[ipT122],&TauDummy);

	/* N II 2140, data 
	 * >>referold	n2	cs	Stafford, R.P., Bell, K.L, Hibbert, A. & Wijesundera, W.P. 1994, MNRAS, 268, 816
	 * >>referold	n2	cs	Lennon, D.J., & Burke, V.M., 1994, A&AS, 103, 273-277
	 * A from 
	 * >>referold	n2	as	Brage, T., Hibbert, A., Leckrone, D.S. 1997, ApJ, 478, 423
	 * >>refer	n2	as	Froese Fischer, C., & Tachiev, G. 2004, At. Data Nucl. Data Tables, 87, 1
	 * >>refer	n2	cs	Hudson, C.E. & Bell, K.L. 2004, MNRAS, 348, 1275 and A&A, 430, 725
	 */
	cs21 = 1.125;

	PutCS(cs21,&TauLines[ipT2140]);
	atom_level2(&TauLines[ipT2140]);

	/* N III 989.8, cs from 
	 * >>refer	N3	cs	Blum, R.D., & Pradhan, A.K. 1992, ApJS 80, 425 */
	PutCS(7.12,&TauLines[ipT990]);
	atom_level2(&TauLines[ipT990]);

	/* 57 micron N III, A=
	 * >>refer	n3	as	Froese Fischer, C. 1983, J.Phys. B, 16, 157
	 * collision strength from 
	 * >>refer	n3	cs	Blum, R.D., & Pradhan, A.K. 1992, ApJS 80, 425 */
	cs = MIN2(1.90,0.2291*phycon.te10*phycon.te10);
	PutCS(cs,&TauLines[ipT57]);

	static vector< pair<transition*,double> > N3Pump;
	N3Pump.reserve(32);

	/* one time initialization if first call */
	if( N3Pump.empty() )
	{
		// set up level 1 pumping lines
		pair<transition*,double> ppa( &TauLines[ipT990], 1./6. ); 
		N3Pump.push_back( ppa );
		pair<transition*,double> ppb( &TauLines[ipT374g], 1./6. ); 
		N3Pump.push_back( ppb );
		pair<transition*,double> ppc( &TauLines[ipT315], 1./6. ); 
		N3Pump.push_back( ppc );
		pair<transition*,double> ppd( &TauLines[ipT324], 1./2. ); 
		N3Pump.push_back( ppd );
		pair<transition*,double> ppe( &TauLines[ipT333], 2./3. ); 
		N3Pump.push_back( ppe );
		// set up level 2 pumping lines
		for( i=0; i < nWindLine; ++i )
		{
			/* don't test on nelem==ipIRON since lines on physics, not C, scale */
			if( TauLine2[i].Hi->nelem == 7 && TauLine2[i].Hi->IonStg == 3 )
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
				N3Pump.push_back( pp2 );
			}
		}
	}

	/* now sum pump rates */
	double pump_rate = 0.;
	vector< pair<transition*,double> >::const_iterator n3p;
	for( n3p=N3Pump.begin(); n3p != N3Pump.end(); ++n3p )
	{
		const transition* t = n3p->first;
		double branch_ratio = n3p->second;
		pump_rate += t->Emis->pump*branch_ratio;
#		if	0
		dprintf( ioQQQ, "N III %.3e %.3e\n",
			 t->WLAng , t->Emis->pump*branch_ratio );
#		endif
	}

	/* N III] N 3, N  3, 1765 multiplet */
	/*atom_level2(&TauLines[ipT57]);*/
	/*AtomSeqBoron compute cooling from 5-level boron sequence model atom */
	AtomSeqBoron(&TauLines[ipT57], 
	  &TauLines[ipN3_1749], 
	  &TauLines[ipN3_1747], 
	  &TauLines[ipN3_1754], 
	  &TauLines[ipN3_1752], 
	  &TauLines[ipN3_1751], 
	  0.201 , 1.088 , 0.668 , 2.044 , pump_rate,"N  3");
	/*fprintf(ioQQQ," n4 %.3e\n",  (
		TauLines[ipN3_1749].xIntensity + 
		TauLines[ipN3_1747].xIntensity + 
		TauLines[ipN3_1754].xIntensity+ 
		TauLines[ipN3_1752].xIntensity+ 
		TauLines[ipN3_1751].xIntensity) / dense.xIonDense[ipNITROGEN][2] );*/

	/* N IV 1486, N 4, N  4,collisions within 3P just guess
	 * cs to ground from 
	 * >>refer	n4	cs	Ramsbottom, C.A., Berrington, K.A., Hibbert, A., Bell, K.L. 1994,
	 * >>refercon Physica Scripta, 50, 246 */
	if( phycon.te > 1.584e4 )
	{
		cs = 21.346/(phycon.te10*phycon.te10*phycon.te10*phycon.te02);
	}
	else
	{
		cs = 75.221/(phycon.sqrte/phycon.te03/phycon.te02);
	}
	/* >>chng 01 sep 09, AtomSeqBeryllium will reset this to 1/3 so critical density correct */
	cs = MAX2(0.01,cs);
	PutCS(cs,&TauLines[ipT1486]);
	AtomSeqBeryllium(.9,.9,3.0,&TauLines[ipT1486],.0115);
	embesq.em1486 = (realnum)(atoms.PopLevels[3]*0.0115*1.34e-11);
	/*fprintf(ioQQQ," n4 %.3e\n", (TauLines[ipT1486].xIntensity + embesq.em1486 ) / dense.xIonDense[ipNITROGEN][3] );*/

	/* N IV 765, cs from 
	 * >>refer	n4	cs	Ramsbottom, C.A., Berrington, K.A., Hibbert, A., Bell, K.L. 1994,
	 * >>refercon Physica Scripta, 50, 246 */
	/* >>refer	n4	as	Flemming, J., Brage, T., Bell, K.L., Vaeck, N., Hibbert, A., 
	 * >>refercon	Godefroid, M., & Froese Fischer, C., 1995, ApJ, 455, 758*/
	cs = MIN2(4.0,1.864*phycon.te03*phycon.te03);
	PutCS(cs,&TauLines[ipT765]);
	atom_level2(&TauLines[ipT765]);

	/* N V 1240
	 * >>refer	n5	cs	Cochrane, D.M., & McWhirter, R.W.P. 1983, PhyS, 28, 25 */
	ligbar(7,&TauLines[ipT1239],&TauLines[ipT209],&cs2s2p,&cs2s3p);
	PutCS(cs2s2p,&TauLines[ipT1239]);
	PutCS(cs2s2p*0.5,&TauLines[ipT1243]);
	PutCS(1.0,&TauDummy);
	atom_level3(&TauLines[ipT1243],&TauDummy,&TauLines[ipT1239]);
	/*fprintf(ioQQQ," n5 %.3e\n", (TauLines[ipT1243].xIntensity + TauLines[ipT1239].xIntensity ) / dense.xIonDense[ipNITROGEN][4] );*/

	/* N V 209 */
	PutCS(cs2s3p,&TauLines[ipT209]);
	atom_level2(&TauLines[ipT209]);
	return;
}
