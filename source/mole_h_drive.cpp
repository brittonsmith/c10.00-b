/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*hmole determine populations of hydrogen molecules */
/*hmole_old determine populations of hydrogen molecules */
/*hmole_init - initialize some hmole vars */
/*hmole_reactions update hmole reactions */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "called.h"
#include "thermal.h"
#include "gammas.h"
#include "colden.h"
#include "thermal.h"
#include "secondaries.h"
#include "h2.h"
#include "mole.h"
#include "radius.h"
#include "doppvel.h"
#include "rfield.h"
#include "ionbal.h"
#include "rt.h"
#include "opacity.h"
#include "iso.h"
#include "conv.h"
#include "phycon.h"
#include "hmi.h"
#include "cosmology.h"

/* Define to verify chemistry solution */
#if 0
#if !defined(NDEBUG)
/* >>chng 02 dec 21, line 14 changed from 
	if(fabs(total) > 1e-20 && fabs(total) > 1e-14*mtotal) { 
	to
	if(fabs(total) > 1e-20 && fabs(total) > 1e-8*mtotal) {  */
#define AUDIT(a)	{ \
		double total, mtotal; \
		for( i=0;i<N_H_MOLEC;i++) { \
			total = 0.; \
			for( j=0;j<N_H_MOLEC;j++) { \
				total += c[i][j]*nprot[j]; \
			} \
			if( fabs(total) > 1e-6*fabs(c[i][i]*nprot[i])) { \
					fprintf(ioQQQ,"PROBLEM Subtotal1 %c %.2e\n",a,fabs(total)/fabs(c[i][i]*nprot[i])); \
					fprintf(ioQQQ,"Species %li Total %g Diag %g\n",i,total,c[i][i]*nprot[i]); \
			} \
		} \
    total = mtotal = 0.;for( j=0;j<N_H_MOLEC;j++) { total += bvec[j]*nprot[j]; mtotal += fabs(bvec[j]*nprot[j]); }\
			if( fabs(total) > 1e-30 && fabs(total) > 1e-10*rtot) { \
					fprintf(ioQQQ,"PROBLEM Subtotal2 %c %.2e\n",a,fabs(total)/mtotal); \
					fprintf(ioQQQ,"RHS Total %g cf %g\n",total,mtotal); \
			} else if( a == '.' && fabs(total) > 1e-7*mtotal)  { \
					fprintf(ioQQQ,"PROBLEM zone %li Hmole RHS conservation error %.2e of %.2e\n",nzone,total,mtotal); \
					fprintf(ioQQQ,"(may be due to high rate equilibrium reactions)\n"); \
			} \
	}
#else
#define AUDIT /* nothing */
#endif

#endif

void hmole( void )
{
	int nFixup, i;
	double error;
	realnum oxy_error=0.;
	static realnum abund0=-BIGFLOAT , abund1=-BIGFLOAT;
	realnum save1=dense.xIonDense[ipOXYGEN][1], 
		save0=dense.xIonDense[ipOXYGEN][0];

	DEBUG_ENTRY( "hmole()" );

	/* will be used to keep track of neg solns */
	nFixup = 0;
	error = 1.;

	/* >>chng 04 apr 29, use PDR solution, scaled to density of H0, as first guess*/
	/* if very first call on this sim, set H mol abundances to scaled TH85 PDR values */
	if( conv.nTotalIoniz==0 && iteration==0 )
	{
		realnum pdr_mole_h2[N_H_MOLEC] = {1.00E+00f,
			3.18E-05f,
			1.95E-11f,
			4.00E-08f,
			1.08E-14f,
			1.08E-20f,
			3.85E-07f,
			8.04E-14f};

		/* we should have an H0 soln at this point */
		ASSERT( dense.xIonDense[ipHYDROGEN][0]>SMALLFLOAT );
		/* make sure order of molecules has not changed */
		ASSERT( ipMH==0&&
			ipMHp  == 1&&
			ipMHm  == 2&&
			ipMH2g == 3&&
			ipMH2p == 4&&
			ipMH3p == 5&&
			ipMH2s == 6&&
			ipMHeHp== 7&&
			N_H_MOLEC==8 );

		for( i=0; i<N_H_MOLEC; ++i )
		{
			hmi.Hmolec[i] = dense.xIonDense[ipHYDROGEN][0]*pdr_mole_h2[i];
		}
	}

	/* update hmole reactions that depend on temperature */
	hmole_reactions();

	/* will test against this error for quitting */
#	define BIGERROR		1e-4
	/* >>chng 04 jul 20, upper limit had been 6, why?  change to 20
	 * to get closer to soln */
#	define LIM_LOOP	20
	/* loop until run to limit, or no fix ups needed and error is small */
	for(i=0; i < LIM_LOOP && ((error > BIGERROR||nFixup || oxy_error>conv.EdenErrorAllowed));i++) 
	{
		{
			/* option to print deduced abundances */
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC && (nzone>140) )
			{
				fprintf(ioQQQ,"DEBUG hmole in \t%.2f\t%.5e\t%.5e",
					fnzone,
					phycon.te,
					dense.eden);
				for( i=0; i<N_H_MOLEC; i++ )
					fprintf(ioQQQ,"\t%.2e", hmi.Hmolec[i] );
				fprintf(ioQQQ,"\n" );
			}
		}
		/* nFixup is number of times negative abundances were fixed, should be zero 
		 * at return for valid soln */
		nFixup = 0;
		/* >>chng 04 jul 16, bring oxygen into the H solver */
		if( !conv.lgSearch )
		{
			IonOxyge();
		}

		save1 = dense.xIonDense[ipOXYGEN][1];
		save0 = dense.xIonDense[ipOXYGEN][0];
		/* >>chng 04 jul 29, necessary to damp out small changes in the O^0, O^+ density
		 * when calling the hydrogen solver, since this can cause changes in H+/H0, which feed
		 * back into the O+/O0 - small oscillations develop - shown in finely converged
		 * dynamics models */
		if( nzone )
		{
#			define OLD	0.75f
			abund0 = abund0*OLD+dense.xIonDense[ipOXYGEN][0]*(1.f-OLD);
			abund1 = abund1*OLD+dense.xIonDense[ipOXYGEN][1]*(1.f-OLD);
		}
		else
		{
			abund0 = dense.xIonDense[ipOXYGEN][0];
			abund1 = dense.xIonDense[ipOXYGEN][1];
		}
		dense.xIonDense[ipOXYGEN][0] = abund0;
		/* conserve sum of O0 + O+ */
		dense.xIonDense[ipOXYGEN][1] = abund1;
		/* error in this averaging */
		oxy_error = (realnum)fabs(save1-abund1)/SDIV(dense.gas_phase[ipOXYGEN]);
		/*fprintf(ioQQQ,"DEBUG hmole\t%.2f\t%.5e\t%.5e\t%.5e\n", fnzone,save0 , save1 , oxy_error);*/
		/*dense.xIonDense[ipOXYGEN][1] = abund1;*/
		/* call the hydrogen solver */
		hmole_step(&nFixup, &error);
		dense.xIonDense[ipOXYGEN][0] = save0;
		dense.xIonDense[ipOXYGEN][1] = save1;
		{
			/* option to print deduced abundances */
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC /*&& (nzone>68)*/ )
			{
				fprintf(ioQQQ,"DEBUG hmole out\t%i\t%.2f\t%.5e\t%.5e",
					i,
					fnzone,
					phycon.te,
					dense.eden);
				fprintf(ioQQQ,
					"\terror\t%.4e\tH0\t%.4e\tH+\t%.4e\tsink\t%.4e\t%.4e\tsour\t%.4e\t%.4e\tion\t%.4e\trec\t%.4e", 
					error,
					dense.xIonDense[ipHYDROGEN][0],
					dense.xIonDense[ipHYDROGEN][1],
					mole.sink[ipHYDROGEN][0],
					mole.sink[ipHYDROGEN][1],
					mole.source[ipHYDROGEN][0] , 
					mole.source[ipHYDROGEN][1] ,
					ionbal.RateIonizTot(ipHYDROGEN,0),
					ionbal.RateRecomTot[ipHYDROGEN][0]);

				/*for( j=0; j<N_H_MOLEC; j++ )
					fprintf(ioQQQ,"\t%.4e", hmi.Hmolec[j] );*/
				fprintf(ioQQQ,"\n" );
			}
		}
	}

	if( (i == LIM_LOOP && error > BIGERROR)  || nFixup != 0 )
	{
		/* most of these failures occur just one time during search phase -
		 * that is not a serious problem */
		conv.lgConvPops = false;

		if( !conv.lgSearch && called.lgTalk )
		{
			fprintf(ioQQQ," PROBLEM  hmole, zone %li: %d iters, %d bad; final error: %g lgSearch %i\n",
			nzone,
			i,
			nFixup,
			error,
			conv.lgSearch);
		}
		ConvFail( "pops" , "Hmole");
	}

	/* check that density is still correct - CDEN is constant density */
	if( strcmp( dense.chDenseLaw, "CDEN" )==0 && !cosmology.lgDo )
		/* check that we are conserving hydrogen nuclei */
		ASSERT( fabs( dense.gas_phase[ipHYDROGEN] - dense.den0 )/
		dense.gas_phase[ipHYDROGEN]<1e-4 );


	/* total number of H per unit vol in molecules,
	 * of course not including H0/H+ */
	dense.xMolecules[ipHYDROGEN] = 0.;
	for(i=0;i<N_H_MOLEC;i++) 
	{
		dense.xMolecules[ipHYDROGEN] += hmi.Hmolec[i]*hmi.nProton[i];
	}
	/* remove the atom/ion which was just counted */
	dense.xMolecules[ipHYDROGEN] -= (hmi.Hmolec[ipMH]+hmi.Hmolec[ipMHp]);

	/* now add on all H in heavy element molecules */
	for( i=0; i < mole.num_comole_calc; i++ )
	{
		dense.xMolecules[ipHYDROGEN] += COmole[i]->hevmol*COmole[i]->nElem[ipHYDROGEN];
	}

	/*fprintf(ioQQQ," hmole return value is %.3e\n",timesc.time_H2_Dest_here);*/
	return;
}

/*hmirat compute radiative association rate for H- */
double hmirat(double te)
{
	double hmirat_v;

	DEBUG_ENTRY( "hmirat()" );

	/* fits to radiative association rate coefficients */
	if( te < 31.62 )
	{
		hmirat_v = 8.934e-18*phycon.sqrte*phycon.te003*phycon.te001*
		  phycon.te001;
	}
	else if( te < 90. )
	{
		hmirat_v = 5.159e-18*phycon.sqrte*phycon.te10*phycon.te03*
		  phycon.te03*phycon.te003*phycon.te001;
	}
	else if( te < 1200. )
	{
		hmirat_v = 2.042e-18*te/phycon.te10/phycon.te03;
	}
	else if( te < 3800. )
	{
		hmirat_v = 8.861e-18*phycon.te70/phycon.te03/phycon.te01*
		  phycon.te003;
	}
	else if( te <= 1.4e4 )
	{
		/* following really only optimized up to 1e4 */
		hmirat_v = 8.204e-17*phycon.sqrte/phycon.te10/phycon.te01*
		  phycon.te003;
	}
	else
	{
		/* >>chng 00 sep 28, add this branch */
		hmirat_v = 5.44e-16*phycon.te20/phycon.te01;
	}
	return( hmirat_v );
}


/* hmole_init - one-time initialize of some hmole vars, called by cdinit  */
void hmole_init(void)
{

	DEBUG_ENTRY( "hmole_init()" );

	/* do we need to fill in the molecular labels? */
	strcpy( hmi.chLab[ipMH], "H0  " );
	strcpy( hmi.chLab[ipMHp], "H+  " );
	strcpy( hmi.chLab[ipMHm], "H-  " );
	strcpy( hmi.chLab[ipMH2g], "H2g " );
	strcpy( hmi.chLab[ipMH2p], "H2+ " );
	strcpy( hmi.chLab[ipMH3p], "H3+ " );
	strcpy( hmi.chLab[ipMH2s], "H2* " );	
	strcpy( hmi.chLab[ipMHeHp], "HeH+" );	

	/* rate for He+ + H2 -> */
	hmi.rheph2hpheh = 0.;
	return;

}

/*hmole_reactions update hmole reactions */
void hmole_reactions( void )
{
	static double teused=-1;
	double exph2,
		exph2s,
		exphp,
	  ex3hp;
	long i;
	double h2esc, 
	  th2,
	  cr_H2s ,
	  cr_H2dis,
	  cr_H2dis_ELWERT_H2g,
	  cr_H2dis_ELWERT_H2s;

	DEBUG_ENTRY( "hmole_reactions()" );

	/* everything here depends only on temperature - don't do anything if we don't
	 * need to */
	if( ! fp_equal( phycon.te, teused ) )
	{
		teused = phycon.te;

		/* get LTE populations */
		/* related to population of H- in LTE
		* IP is 0.754 eV */
		hmi.exphmi = sexp(8.745e3/phycon.te);
		if( hmi.exphmi > 0. )
		{
			/* these are ratio n*(H-)/[  n*(ne) n*(Ho)  ] */
			hmi.rel_pop_LTE_Hmin = SAHA/(phycon.te32*hmi.exphmi)*(1./(2.*2.));
		}
		else
		{
			hmi.rel_pop_LTE_Hmin = 0.;
		}

		/* population of H2 in LTE
		 * dissociation energy of H2g is 4.477eV, for TH85 model */
		/* >>chng 05 oct 17, GS, averaged in big H2 in order to consider correct statistical weight*/
		if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			/* the terms on the right were computed in the large molecule */
			hmi.rel_pop_LTE_H2g = hmi.H2g_LTE_bigH2;
			hmi.rel_pop_LTE_H2s = hmi.H2s_LTE_bigH2;
#			if 0
			exph2 = sexp((5.195e4-hmi.H2_BigH2_H2g_av*T1CM)/phycon.te);
			/* >>chng 05 oct 21, GS, protect against zero Boltzmann factor */
			if( exph2 > 0. )
				hmi.rel_pop_LTE_H2g =  SAHA/SDIV(phycon.te32*exph2)*(1./(2.*2.))*3.634e-5;
			else
				hmi.rel_pop_LTE_H2g =  0.;
			/*fprintf(ioQQQ,"test\t%.2e\t%.1e\t%.1e\n", 
					hmi.rel_pop_LTE_H2g, 
					exph2,
					SAHA/(phycon.te32*exph2)*(1./(2.*2.))*3.634e-5);*/

			exph2s = sexp((5.195e4-hmi.H2_BigH2_H2s_av * T1CM)/phycon.te);
			/* >>chng 05 oct 21, GS, protect against zero Boltzmann factor */
			if( exph2s > 0. )
				hmi.rel_pop_LTE_H2s = SAHA/SDIV(phycon.te32*exph2s)*(1./(2.*2.))*3.634e-5; 
			else
				hmi.rel_pop_LTE_H2s = 0.; 
			/*fprintf(ioQQQ,"testh2s\t%.2e\t%.1e\t%.1e\n", 
					hmi.rel_pop_LTE_H2s, 
					exph2s,
					SAHA/(phycon.te32*exph2s)*(1./(2.*2.))*3.634e-5);*/
#			endif
		}
		else
		{

			/* H2 ground */
			exph2 = sexp((5.195e4)/phycon.te);
			/* >>chng 05 oct 17, GS, note that statical weight of H2g is assumed to be 1 if big H2 is not turned on*/

			if( exph2 > 0. ) 
			{
				/* >>chng 05 oct 17, GS, note that statical weight of H2g is assumed to be 1 if big H2 is not turned on*/
				hmi.rel_pop_LTE_H2g = SAHA/(phycon.te32*exph2)*(1./(2.*2.))*3.634e-5;
			}
			else
			{
				hmi.rel_pop_LTE_H2g = 0.;
			}

			/* H2 star */
			/* population of H2s in LTE
			 * dissociation energy is 1.877eV, if h2s = 2.6eV, assumed for TH85 model */
			/* >>chng 05 oct 17, GS, averaged in big H2 in order to consider correct statistical weight*/
			exph2s = sexp(2.178e4/phycon.te);

			if( exph2s > 0. ) 
			{
			/* >>chng 05 oct 17, GS, note that statical weight of H2s is assumed to be 1 if big H2 is not turned on*/
				hmi.rel_pop_LTE_H2s = SAHA/(phycon.te32*exph2s)*(1./(2.*2.))*3.634e-5;
			}
			else
			{
				hmi.rel_pop_LTE_H2s = 0.;
			}
		}
		{
			/*@-redef@*/
			/* often the H- route is the most efficient formation mechanism for H2,
			* will be through rate called Hmolec_old[ipMH]*hmi.assoc_detach
			* this debug print statement is to trace h2 oscillations */
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC && nzone>187&& iteration > 1/**/)
			{
				/* rapid increase in H2 density caused by rapid increase in hmi.rel_pop_LTE_H2g */
				fprintf(ioQQQ,"ph2lteee\t%.2e\t%.1e\t%.1e\n", 
					hmi.rel_pop_LTE_H2g, 
					sexp(2.178e4/phycon.te),
					phycon.te);
			}
		}

		/* population of H2+ in LTE, hmi.rel_pop_LTE_H2p is H_2+/H / H+
		* dissociation energy is 2.647 */
		exphp = sexp(3.072e4/phycon.te);
		if( exphp > 0. )
		{
			/* stat weight of H2+ is 4
			* last factor was put in ver 85.23, missing before */
			hmi.rel_pop_LTE_H2p = SAHA/(phycon.te32*exphp)*(4./(2.*1.))*3.634e-5;
		}
		else
		{
			hmi.rel_pop_LTE_H2p = 0.;
		}

		/* related to population of H3+ in LTE, hmi.rel_pop_LTE_H3p is H_3+/( H2+ H+ )
		* dissociation energy is 2.647 */
		ex3hp = sexp(1.882e4/phycon.te);
		if( ex3hp > 0. )
		{
			/* stat weight of H2+ is 4
			* last factor was put in ver 85.23, missing before */
			hmi.rel_pop_LTE_H3p = SAHA/(phycon.te32*ex3hp)*(4./(2.*1.))*3.634e-5;
		}
		else
		{
			hmi.rel_pop_LTE_H3p = 0.;
		}
	}
	/* end constant temperature - */


	hmi.hminus_rad_attach = hmirat(phycon.te);
	/* cooling due to radiative attachment */
	hmi.hmicol = hmi.hminus_rad_attach*EN1RYD*phycon.te*1.15e-5;

	/*fprintf(ioQQQ,"%.2e %.2e %.2e %.2e\n", phycon.te, hmi.hminus_rad_attach , hmi.hmicol,
		hmi.hmicol/(hmi.hminus_rad_attach*EN1RYD*phycon.te*1.15e-5) );*/

	/* get per unit vol */
	hmi.hminus_rad_attach *= dense.eden;
	hmi.hmicol *= dense.eden*hmi.Hmolec[ipMH]; /* was dense.xIonDense[ipHYDROGEN][0]; */

	/* ================================================================= */
	/* evaluate H- photodissociation rate, induced rec and rec cooling rates */
	/* >>chng 00 dec 24, add test so that photo rate only reevaluated two times per zone.
	 * in grain-free models this was sometimes dominated by Lya and so oscillated.  
	 * especially bad in primal.in - change 2 to 4 and primal.in will stop due to Lya oscil */
	/** \todo	2	following always true, why?  either remove test or use it - 
	 * it is here to save time - this step routine is called very often */
	/* >>chng 02 feb 16, add damper on H- photo rate, wild oscillations in Lya photo rate in 
	 * grain free models */

	t_phoHeat photoHeat;

	hmi.HMinus_photo_rate = GammaBn( hmi.iphmin-1 , iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0] , opac.iphmop ,
		0.055502 , &hmi.HMinus_induc_rec_rate , &hmi.HMinus_induc_rec_cooling, &photoHeat );

	/* save H- photodissociation heating */
	hmi.HMinus_photo_heat = photoHeat.HeatNet;

	{
		/* following should be set true to print populations */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC)
		{
			fprintf(ioQQQ,"hminphoto\t%li\t%li\t%.2e\n", nzone, conv.nPres2Ioniz , hmi.HMinus_photo_rate );
		}
	}

	/* induced recombination */
	hmi.HMinus_induc_rec_rate *= hmi.rel_pop_LTE_Hmin*dense.eden;

	/* induced recombination cooling per unit volume */
	/** \todo	2	this should be done with new populations after converged soln */
	hmi.HMinus_induc_rec_cooling *= hmi.rel_pop_LTE_Hmin*dense.eden*hmi.Hmolec[ipMH]; /* dense.xIonDense[ipHYDROGEN][0]; */

	{
		/* following should be set true to debug H- photoionization rates */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && nzone>400/*&& iteration > 1*/)
		{
			fprintf(ioQQQ,"hmoledebugg %.2f ",fnzone);
			GammaPrt(
				hmi.iphmin-1 , iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0] , opac.iphmop ,
				/* io unit we will write to */
				ioQQQ, 
				/* total photo rate from previous call to GammaK */
				hmi.HMinus_photo_rate, 
				/* we will print contributors that are more than this rate */
				hmi.HMinus_photo_rate*0.05);
		}
	}
	/* add on high energy ionization, assume hydrogen cross section
	 * n.b.; HGAMNC with secondaries */
	/* >>chng 00 dec 24, above goes to HeI edge, no need for this, and was not important*/
	/*hmi.HMinus_photo_rate += iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s];*/

	/* ================================================================= */
	/* photodissociation by Lyman band absorption: esc prob treatment,
	* treatment based on 
	* >>refer	HI	abs	Tielens & Hollenbach 1985 ApJ 291, 722. */
	/* do up to carbon photo edge if carbon is turned on */
	/* >>>chng 00 apr 07, add test for whether element is turned on */
	hmi.UV_Cont_rel2_Habing_TH85_depth = 0.;
	/* >>chng 00 apr 07 from explicit ipHeavy to ipLo */
	/* find total intensity over carbon-ionizing continuum */
	/* >>chng 03 jun 09, use exact bounds rather than CI photo threshold for lower bound */
	/*for( i=ipLo; i < iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]; i++ )*/
	/* the integral is from 6eV to 13.6, or 2060 - 912 Ang */
	for( i=rfield.ipG0_TH85_lo; i < rfield.ipG0_TH85_hi; ++i )
	{
		hmi.UV_Cont_rel2_Habing_TH85_depth += (rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
			rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1])*rfield.anu[i-1];
	}

	/* now convert to Habing ISM units
	 * UV_Cont_rel2_Habing_TH85_face is FUV continuum relative to Habing value 
	 * 1.6e-3 ergs em-2 s-1 is the Habing 1968 field, quoted on page 723, end of first
	 * full paragraph on left */
	hmi.UV_Cont_rel2_Habing_TH85_depth = (realnum)(hmi.UV_Cont_rel2_Habing_TH85_depth*EN1RYD/1.6e-3);
	/* if start of calculation remember G0 at illuminated face */
	if( nzone<=1 )
	{
		hmi.UV_Cont_rel2_Habing_TH85_face = hmi.UV_Cont_rel2_Habing_TH85_depth;
	}


	/* >>chng 05 jan 09, add special ver of G0 */
	hmi.UV_Cont_rel2_Habing_spec_depth = 0.; 
	for( i=rfield.ipG0_spec_lo; i < rfield.ipG0_spec_hi; ++i )
	{
		hmi.UV_Cont_rel2_Habing_spec_depth += (rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
			rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1])*rfield.anu[i-1];
	}
	hmi.UV_Cont_rel2_Habing_spec_depth = (realnum)(hmi.UV_Cont_rel2_Habing_spec_depth*EN1RYD/1.6e-3);

	/* the Draine & Bertoldi version of the same thing, defined over their energy range */
	/* >>chng 04 feb 07, only evaluate at the illuminated face */
	if( hmi.UV_Cont_rel2_Draine_DB96_face ==0 )
	{
		/* this is sum of photon number between 912A and 1110, as per BD96 */
		for( i=rfield.ipG0_DB96_lo; i < rfield.ipG0_DB96_hi; ++i )
		{
			hmi.UV_Cont_rel2_Draine_DB96_face += (rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
				rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1]);
		}
		/* convert into scaled ISM background field, total number of photons over value for 1 ISM field,
		 * the coefficient 1.232e7 is the number of photons over this wavelength range for 1H and is
		 * given in BD96, page 225, 4th line from start of section 4, also page 272, table 1, 2nd line
		 * from bottom */
		/* >>chng 04 mar 16, introduce the 1.71 */
		/* equation 20 of DB96 gives chi as flux over 1.21e7, to produce one Habing field.
		 * to get the Draine field we need to further divide by 1.71 as stated on the first
		 * line after equation 23 */
		hmi.UV_Cont_rel2_Draine_DB96_face = hmi.UV_Cont_rel2_Draine_DB96_face/(1.232e7f * 1.71f);
	}

	/* escape prob takes into account line shielding, 
	 * next is opacity then optical depth in H2 UV lines, using eqn A7 of TH85,
	 * LIMELM+2 points to H2 */
	hmi.H2Opacity = (realnum)(1.2e-14*(1e5/GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])));
	/* the typical Lyman -Werner H2 line optical depth eq A7 of TH85a */
	th2 = (colden.colden[ipCOL_H2g]+ colden.colden[ipCOL_H2s])*hmi.H2Opacity;
	/* the escape probability - chance that continuum photon will penetrate to
	 * this depth to pump the Lyman Werner bands */
	h2esc = esc_PRD_1side(th2,1e-4);

	/* cross section is eqn A8 of 
	 * >>refer	H2	dissoc	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722
	 * branching ratio of 10% is already included, so 10x smaller than their number
	 * 10% lead to dissociation through H_2 + h nu => 2H */
	/* >>chng 05 mar 10, by TE, break into 2g and 2s 
	 * note use of same shielding column in below - can do far better */
	hmi.H2_Solomon_dissoc_rate_TH85_H2g = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_Solomon_dissoc_rate_TH85_H2s = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_H2g_to_H2s_rate_TH85 = hmi.H2_Solomon_dissoc_rate_TH85_H2g*9.;

	/* these are Burton et al. 1990 rates */
	hmi.H2_Solomon_dissoc_rate_BHT90_H2g = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_Solomon_dissoc_rate_BHT90_H2s = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_H2g_to_H2s_rate_BHT90 = hmi.H2_Solomon_dissoc_rate_BHT90_H2g*9.;

	{ 
		/* the following implements Drain & Bertoldi 1996, equation 37 from
		 * >>refer	H2	destruction	Draine, B.T., & Bertoldi, F., 1996, ApJ, 468, 269-289
		 * but the constant 4.6e-11 comes from Bertoldi & Draine equation 23,
		 * this is the normalized H2 column density */
		double x = (colden.colden[ipCOL_H2g]+colden.colden[ipCOL_H2s]) / 5e14;
		double sqrtx = sqrt(1.+x);
		/* Doppler with of H2 in km/s */
		double b5 = GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])/1e5;
		/* the molecular hydrogen line self-shielding factor */
		double fshield = 0.965/POW2(1.+x/b5) + 0.035/sqrtx *
			sexp(8.5e-4*sqrtx);

		/*double fshield = pow( MAX2(1.,colden.colden[ipCOLH2]/1e14) , -0.75 );*/
		/* this is the Bertoldi & Draine version of the Habing field,
		 * with dilution of radiation and extinction due to grains */
		/* >>chng 04 apr 18, moved fshield, the line shielding factor, from this defn to
		 * the following defn of dissociation rate, since following should measure continuum */
		hmi.UV_Cont_rel2_Draine_DB96_depth = hmi.UV_Cont_rel2_Draine_DB96_face * 
			(realnum)(sexp( opac.TauAbsFace[rfield.ip1000A-1] )/radius.r1r0sq);

		/* the following comes from Bertoldi & Draine 1996, equation 23,
		 * hmi.UV_Cont_rel2_Draine_DB96_depth already includes a factor of 1.71 to correct back from TH85 */
		/* >>chng 05 mar 10, TE, break into 2s and 2s */
		if( co.lgUMISTrates )
		{
			/* this is default, when set Leiden hack UMIST rates not entered */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g = 4.6e-11 * hmi.UV_Cont_rel2_Draine_DB96_depth * fshield;
			hmi.H2_Solomon_dissoc_rate_BD96_H2s = 4.6e-11 * hmi.UV_Cont_rel2_Draine_DB96_depth * fshield;
		}
		else
		{
			/* done when set Leiden hack UMIST rates command entered */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g = 5.18e-11* (hmi.UV_Cont_rel2_Habing_TH85_face/1.66f)
				*sexp(3.02*rfield.extin_mag_V_point)* fshield;
			hmi.H2_Solomon_dissoc_rate_BD96_H2s = 5.18e-11* (hmi.UV_Cont_rel2_Habing_TH85_face/1.66f)
				*sexp(3.02*rfield.extin_mag_V_point)* fshield;
		}

		/* BD do not give an excitation rate, so used 9 times the dissociation
		 * rate by analogy with 90% going back into X, as per TH85 */
		/*>>chng 04 feb 07, had used 90% relax into X from TH85,
		 * BD say that 15% dissociate, so 85/15 = 5.67 is factor */
		hmi.H2_H2g_to_H2s_rate_BD96 = 5.67* hmi.H2_Solomon_dissoc_rate_BD96_H2g;
	}

	/* do Elwert approximation to the dissociation rate */
	if( hmi.UV_Cont_rel2_Draine_DB96_face > SMALLFLOAT )
	{
		/* this evaluates the new H2 dissociation rate by Torsten Elwert */
		/* chng 05 jun 23, TE
		 * >>chng 05 sep 13, update master source with now approximation */

		/* Doppler with of H2 in km/s */
		double b5 = GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])/1e5;			

		/* split the Solomon rates in H2g and H2s */
		/* >>chng 05 oct 19, 
		 * >>chng 05 dec 05, TE, define new approximation for the heating due to the destruction of H2
		 *	use this approximation for the specified cloud parameters, otherwise
		 * use limiting cases for 1 <= G0, G0 >= 1e7, n >= 1e7, n <= 1 */

		double x_H2g, x_H2s,
				fshield_H2g, fshield_H2s,
				f_H2s;
		static double a_H2g, a_H2s,
						e1_H2g, e1_H2s,
						e2_H2g,
						b_H2g,
						sl_H2g, sl_H2s,
						k_f_H2s,
						k_H2g_to_H2s, 
						log_G0_face = -1;

		/* define parameter range for the new approximation
		 * test for G0 
		 *BUGFIX - this tested on lg_G0_face < 0 for initialization needed so did not work
		 * in grids - change to evaluate in zone 0 */
		/* >>chng 07 feb 24, BUGFIX crash when G0=0 at start and radiation
		 * field builds up due to diffuse fields - soln is to always reevaluate */
		/*if( !nzone )*/
		{
			if(hmi.UV_Cont_rel2_Draine_DB96_face <= 1.) 
			{ 
				log_G0_face = 0.;
			}
			else if(hmi.UV_Cont_rel2_Draine_DB96_face >= 1e7) 
			{ 
				log_G0_face = 7.;
			}
			else 
			{ 
				log_G0_face = log10(hmi.UV_Cont_rel2_Draine_DB96_face); 
			}

			/* terms only dependent on G0_face */

			/* coefficients and exponents */
			a_H2g = 0.06 * log_G0_face + 1.32;
			a_H2s = 0.375 * log_G0_face + 2.125;

			e1_H2g = -0.05 * log_G0_face + 2.25;
			e1_H2s = -0.125 * log_G0_face + 2.625;

			e2_H2g = -0.005 * log_G0_face + 0.625;

			b_H2g = -4.0e-3  * log_G0_face + 3.2e-2;

			/* scalelength for H2g and H2s */
			sl_H2g = 4.0e14;
			sl_H2s = 9.0e15;

			/* coefficient for 2nd term of Solomon H2s */
			k_f_H2s = MAX2(0.1,2.375 * log_G0_face - 1.875 );

			/* coefficient for branching ratio */
			k_H2g_to_H2s =  MAX2(1.,-1.75 * log_G0_face + 11.25);

			/*fprintf( ioQQQ, "e1_H2g%.2e, e1_H2s%.2e, e2_H2g%.2e, b_H2g%.2e, a_H2g%.2e, a_H2s%.2e,sl_H2g: %.2e,sl_H2s: %.2e\n",
						e1_H2g, e1_H2s, e2_H2g, b_H2g, a_H2g, a_H2s, sl_H2g, sl_H2s);
			*/
		}

		/* Solomon H2s ~G0^0.2 at large depth*/
		f_H2s = k_f_H2s * pow((double)hmi.UV_Cont_rel2_Draine_DB96_depth, 0.2 );

		/* scale length for absorption of UV lines */
		x_H2g = (colden.colden[ipCOL_H2g]) / sl_H2g;
		x_H2s = (colden.colden[ipCOL_H2s]) / sl_H2s;

		/* the molecular hydrogen line self-shielding factor */
		fshield_H2g = 0.965/pow(1.+x_H2g/b5,e1_H2g) + b_H2g/pow(1.+x_H2g/b5,e2_H2g);
		fshield_H2s = 0.965/pow(1.+x_H2s/b5,e1_H2s);

		/* the Elwert Solomon rates for H2g and H2s  hmi.chH2_small_model_type == 'E' */
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2g = a_H2g * 4.6e-11 * fshield_H2g * hmi.UV_Cont_rel2_Draine_DB96_depth;
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2s = a_H2s * 4.6e-11 * fshield_H2s * (hmi.UV_Cont_rel2_Draine_DB96_depth + f_H2s);

		/* assume branching ratio dependent on G0*/
		hmi.H2_H2g_to_H2s_rate_ELWERT = k_H2g_to_H2s * hmi.H2_Solomon_dissoc_rate_ELWERT_H2g;

		/* use G0_BD96 as this definition declines faster with depth which is physical as
		 * the longer wavelengths in the definition of G0_TH85 cannot dissociate
		 * H2s directly */
		hmi.H2_photodissoc_ELWERT_H2s = hmi.UV_Cont_rel2_Draine_DB96_depth*1e-11;
		hmi.H2_photodissoc_ELWERT_H2g = hmi.H2_photodissoc_ELWERT_H2s * 1.0e-10;
	}
	else
	{
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2g = 0.;
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2s = 0.;
		hmi.H2_photodissoc_ELWERT_H2s = 0.;
		hmi.H2_photodissoc_ELWERT_H2g = 0.;
	}

	/* this is rate of photodissociation of H2*, A12 of TH85 */
	hmi.H2_photodissoc_TH85 = hmi.UV_Cont_rel2_Habing_TH85_depth*1e-11;

	/* dissociation rate from Burton et al. 1990 */
	hmi.H2_photodissoc_BHT90 = hmi.UV_Cont_rel2_Habing_TH85_depth*1e-11;

	/* rates for cosmic ray excitation of singlet bound electronic bound excited states 
	 * only add this to small molecule since automatically included in large 
	 *>>refer	H2	cr excit	Dalgarno, A., Yan, Min, & Liu, Weihong 1999, ApJS, 125, 237
	 * this is excitation of H2* */
	/* >>chng 05 sep 13, do not include this process when Leiden hacks are in place */
	cr_H2s = secondaries.x12tot*0.9 / 2. * hmi.lgLeidenCRHack;
	/* this is the fraction that dissociate */
	/* >>chng 05 sep 13, do not include this process when Leiden hacks are in place */
	cr_H2dis = secondaries.x12tot*0.1 / 2. * hmi.lgLeidenCRHack;

	/* >>chng 05 sep 13, TE, synchronize treatment of CR */
	/* cosmic ray rates for dissociation of ground and H2s 
	 * two factors done to agree with large H2 deep in the cloud where
	 * cosmic rays are important */
	cr_H2dis_ELWERT_H2g = secondaries.x12tot*5e-8 * hmi.lgLeidenCRHack; 
	cr_H2dis_ELWERT_H2s = secondaries.x12tot*4e-2 * hmi.lgLeidenCRHack;

	/* at this point there are two or three independent estimates of the H2 dissociation rate.
	 * if the large H2 molecule is on, then H2 Solomon rates has been defined in the last
	 * call to the large molecule.  Just above we have defined hmi.H2_Solomon_dissoc_rate_TH85,
	 * the dissociation rate from Tielens & Hollenbach 1985, and hmi.H2_Solomon_dissoc_rate_BD96,
	 * the rate from Bertoldi & Draine 1996.  We can use any defined rate.  If the big H2
	 * molecule is on, use its rate.  If not, for now use the TH85 rate, since that is the
	 * rate we always have used in the past.
	 * The actual rate we will use is given by hmi.H2_Solomon_dissoc_rate_used
	 */
	/* this is the Solomon process dissociation rate actually used */
	if( h2.lgH2ON  && hmi.lgBigH2_evaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		/* only update after big H2 molecule has been evaluated,
		 * when very little H2 and big molecule not used, leave at previous (probably TH85) value,
		 * since that value is always known */

		/* Solomon process rate from X into the X continuum with units s-1
		 * rates are total rate, and rates from H2g and H2s */ 
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_BigH2_H2g;
		if( hmi.H2_Solomon_dissoc_rate_used_H2g <= 0. )
		{
			/*fprintf(ioQQQ,"DEBUG his sol dis <0\n");*/
			hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_TH85_H2g;
		}

		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_BigH2_H2s;
		if( hmi.H2_Solomon_dissoc_rate_used_H2s <= 0. )
			hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_TH85_H2g;

		/* photoexcitation from H2g to H2s */
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_BigH2;
		if( hmi.H2_H2g_to_H2s_rate_used <= 0. )
			hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_TH85;

		/* add up H2s + hnu (continuum) => 2H + KE, continuum photodissociation,
		 * this is not the Solomon process, true continuum, units s-1 */
		/* only include rates from H2s since this is only open channel, this process is well
		 * shielded against Lyman continuum destruction by atomic hydrogen */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_BigH2_H2s;
		if( hmi.H2_photodissoc_used_H2s <= 0. )
			hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_TH85;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_BigH2_H2g;
		if( hmi.H2_photodissoc_used_H2g <= 0. )
			hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_TH85*1.0e-10f;
	}
	else if( hmi.chH2_small_model_type == 'T' )
	{
		/* the TH85 rate  */
		/*>>chng 05 jun 23, add cosmic rays */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_TH85_H2g + cr_H2dis;
		/* >>chng 05 sep 13, cr_H2dis was not included */
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_TH85_H2s + cr_H2dis;
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_TH85 + cr_H2s;

		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_TH85;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_TH85*1.0e-10f;
	}

	else if( hmi.chH2_small_model_type == 'H' )
	{
		/* the improved H2 formalism given by 
		*>>refer	H2	dissoc	Burton, M.G., Hollenbach, D.J. & Tielens, A.G.G.M 
		>>refcon	1990, ApJ, 365, 620 */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_BHT90_H2g + cr_H2dis;
		/* >>chng 05 sep 13, cr_H2dis was not included */
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_BHT90_H2s + cr_H2dis;
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_BHT90 + cr_H2s;

		/* continuum photodissociation H2s + hnu => 2H, ,
		* this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_BHT90;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		* for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_BHT90*1.0e-10f;
	}

	else if( hmi.chH2_small_model_type == 'B' )
	{
		/* the Bertoldi & Draine rate - this is the default */
		/*>>chng 05 jun 23, add cosmic rays */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_BD96_H2g + cr_H2dis;
		/* >>chng 05 sep 13, cr_H2dis was not included */
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_BD96_H2s + cr_H2dis;
		/* they did not do the excitation or dissoc rate, so use TH85 */
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_BD96 + cr_H2s;


		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_TH85;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_TH85*1.0e-10f;
	}
	else if( hmi.chH2_small_model_type == 'E' )
	{
		/* the Elwert et al. rate 
		 *>>chng 05 jun 23, add cosmic rays */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_ELWERT_H2g + cr_H2dis_ELWERT_H2g;
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_ELWERT_H2s + cr_H2dis_ELWERT_H2s;
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_ELWERT + cr_H2s;


		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_ELWERT_H2s;
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_ELWERT_H2g;
	}
	else
		TotalInsanity();

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && h2.lgH2ON )
		{
			fprintf(ioQQQ," Solomon H2 dest rates: TH85 %.2e BD96 %.2e Big %.2e excit rates: TH85 %.2e Big %.2e\n",
				hmi.H2_Solomon_dissoc_rate_TH85_H2g,
				hmi.H2_Solomon_dissoc_rate_BD96_H2g,
				hmi.H2_Solomon_dissoc_rate_BigH2_H2g ,
				hmi.H2_H2g_to_H2s_rate_TH85 ,
				hmi.H2_H2g_to_H2s_rate_BigH2);
		}
	}
	return;
}

