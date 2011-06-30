/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolEvaluate main routine to call others, to evaluate total cooling */
#include "cddefines.h"
#include "physconst.h"
#include "hydrogenic.h"
#include "taulines.h"
#include "wind.h"
#include "coolheavy.h"
#include "radius.h"
#include "conv.h"
#include "h2.h"
#include "rt.h"
#include "doppvel.h"
#include "opacity.h"
#include "ionbal.h"
#include "dense.h"
#include "trace.h"
#include "dynamics.h"
#include "rfield.h"
#include "grainvar.h"
#include "atmdat.h"
#include "atoms.h"
#include "called.h"
#include "hmi.h"
#include "numderiv.h"
#include "magnetic.h"
#include "phycon.h"
#include "lines_service.h"
#include "hyperfine.h"
#include "iso.h"
#include "thermal.h"
#include "cooling.h"
#include "pressure.h"
/*fndneg search cooling array to find negative values */
STATIC void fndneg(void);
/*fndstr search cooling stack to find strongest values */
STATIC void fndstr(double tot, 
  double dc);

/* set true to debug derivative of heating and cooling */
static const bool PRT_DERIV = false;

void CoolEvaluate(double *tot)
{
	static long int nhit = 0, 
	  nzSave=0;

	static double TeEvalCS = 0., TeEvalCS_21cm=0.;
	static double TeUsedBrems=-1.f;
	static int nzoneUsedBrems=-1;

	static double 	electron_rate_21cm,
		atomic_rate_21cm,
		proton_rate_21cm;

	double 
	  cs ,
	  deriv, 
	  factor, 
	  qn, 
	  rothi=-SMALLFLOAT, 
	  rotlow=-SMALLFLOAT, 
	  x;

	static double oltcool=0., 
	  oldtemp=0.;

	DEBUG_ENTRY( "CoolEvaluate()" );

	/* returns tot, the total cooling,
	 * and dc, the derivative of the cooling */

	/* routine atom_level2( t10 )
	 * routine atom_level3( abund , t10,t21,t20)
	 * tsq1 = 1. / (te**2)
	 * POPEXC( O12,g1,g2,A21,excit,abund); result already*a21
	 * POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2)
	 * AtomSeqBeryllium(cs23,cs24,cs34,tarray,a41)
	 * FIVEL( G(1-5) , ex(wn,1-5), cs12,cs13,14,15,23,24,25,34,35,45,
	 *  A21,31,41,51,32,42,52,43,53,54, pop(1-5), abund) */

	if( trace.lgTrace )
		fprintf( ioQQQ, "   COOLR TE:%.4e zone %li %li Cool:%.4e Heat:%.4e eden:%.4e edenTrue:%.4e\n", 
		phycon.te, 
		nzone, conv.nPres2Ioniz ,
		thermal.ctot , thermal.htot,dense.eden,dense.EdenTrue );

	/* must call TempChange since ionization has changed, there are some
	 * terms that affect collision rates (H0 term in electron collision) */
	TempChange(phycon.te , false);

	/* now zero out the cooling stack */
	CoolZero();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  0 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);
	if( gv.lgGrainPhysicsOn )
	{
		/* grain heating and cooling */
		/* grain recombination cooling, evaluated elsewhere
		* can either heat or cool the gas, do cooling here */
		CoolAdd("dust",0,MAX2(0.,gv.GasCoolColl));

		/* grain cooling proportional to temperature ^3/2 */
		thermal.dCooldT += MAX2(0.,gv.GasCoolColl)*3./(2.*phycon.te);

		/* these are the various heat agents from grains */
		/* options to force gas heating or cooling by grains to zero - for tests only ! */
		if( gv.lgDustOn() && gv.lgDHetOn )
		{
			/* rate dust heats gas by photoelectric effect */
			thermal.heating[0][13] = gv.GasHeatPhotoEl;

			/* if grains hotter than gas then collisions with gas act
			* to heat the gas, add this in here
			* a symmetric statement appears in COOLR, where cooling is added on */
			thermal.heating[0][14] = MAX2(0.,-gv.GasCoolColl);

			/* this is gas heating due to thermionic emissions */
			thermal.heating[0][25] = gv.GasHeatTherm;
		}
		else
		{
			thermal.heating[0][13] = 0.;
			thermal.heating[0][14] = 0.;
			thermal.heating[0][25] = 0.;
		}
	}
	else if( gv.lgBakesPAH_heat )
	{
		/* >>chng 06 jul 21, option to include Bakes PAH hack with grain physics off,
		 * needed to test dynamics models */
		thermal.heating[0][13] = gv.GasHeatPhotoEl;
	}

	/* find net heating - cooling due to large H2 molecule */
	H2_Cooling("CoolEvaluate");
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  1 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* molecular molecules molecule cooling */
	if( hmi.lgNoH2Mole )
	{
		/* this branch - do not include molecules */
		hmi.hmicol = 0.;
		CoolHeavy.brems_cool_hminus = 0.;
		/* line cooling within simple H2 molecule - zero when big used */
		CoolHeavy.h2line = 0.;
		/*  H + H+ => H2+ cooling */
		CoolHeavy.H2PlsCool = 0.;
		CoolHeavy.HD = 0.;

		/* thermal.heating[0][8] is heating due to collisions within X of H2 */
		thermal.heating[0][8] = 0.;
		/* thermal.heating[0][15] is H minus heating*/
		thermal.heating[0][15] = 0.;
		/* thermal.heating[0][16] is H2+ heating */
		thermal.heating[0][16] = 0.;
		hmi.HeatH2Dish_used = 0.;
		hmi.HeatH2Dexc_used = 0.;
	}

	else
	{
		/* save various molecular heating/cooling agent */
		thermal.heating[0][15] = hmi.hmihet;
		thermal.heating[0][16] = hmi.h2plus_heat;
		/* now get heating from H2 molecule, either simple or from big one */
		if( h2.lgH2ON  && hmi.lgH2_Thermal_BigH2 )
		{
			if( hmi.lgBigH2_evaluated )
			{
				/* these are explicitly from big H2 molecule,
				 * first is heating due to radiative pump of excited states, followed by
				 * radiative decay into continuum of X, followed by dissociation of molecule
				 * with kinetic energy, typically 0.25 - 0.5 eV per event */
				hmi.HeatH2Dish_used = hmi.HeatH2Dish_BigH2;
				hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BigH2;
				/*fprintf(ioQQQ,"DEBUG big %.2f\t%.5e\t%.2e\t%.2e\t%.2e\n", 
					fnzone , phycon.te, hmi.HeatH2Dexc_used,
					hmi.H2_total , hmi.H2_star_BigH2);*/
				/* negative sign because right term is really deriv of heating,
				 * but will be used below as deriv of cooling */
				hmi.deriv_HeatH2Dexc_used = -hmi.deriv_HeatH2Dexc_BigH2;
			}
			else
			{
				hmi.HeatH2Dish_used = 0;
				hmi.HeatH2Dexc_used = 0;
				hmi.deriv_HeatH2Dexc_used = 0;
			}
		}

		else if( hmi.chH2_small_model_type == 'T' )
		{
			/* TH85 dissociation heating */
			/* these come from approximations in TH85, see comments above */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_TH85;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_TH85;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_TH85;
		}
		else if( hmi.chH2_small_model_type == 'H' )
		{
			/* Burton et al. 1990 */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_BHT90;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BHT90;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_BHT90;
		}
		else if( hmi.chH2_small_model_type == 'B')
		{
			/* Bertoldi & Draine */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_BD96;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BD96;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_BD96;
		}
		else if(hmi.chH2_small_model_type == 'E')
		{
			/* this is the default when small H2 used */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_ELWERT;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_ELWERT;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_ELWERT;
		}
		else
			TotalInsanity();

		/* heating due to photodissociation heating */
		thermal.heating[0][17] = hmi.HeatH2Dish_used;

		/* heating (usually cooling in big H2) due to collisions within X */
		/* add to heating is net heating is positive */
		thermal.heating[0][8] = MAX2(0.,hmi.HeatH2Dexc_used);

		/* add to cooling if net heating is negative */
		CoolAdd("H2cX",0,MAX2(0.,-hmi.HeatH2Dexc_used));
		/*fprintf(ioQQQ,"DEBUG coolh2\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
			fnzone, phycon.te, dense.eden, hmi.H2_total, thermal.ctot, -hmi.HeatH2Dexc_used );*/
		/* add to net derivative */
		/*thermal.dCooldT += MAX2(0.,-hmi.HeatH2Dexc_used)* ( 30172. * thermal.tsq1 - thermal.halfte );*/
		/* >>chng 04 jan 25, check sign to prevent cooling from entering here,
		 * also enter neg sign since going into cooling stack (bug), in heatsum
		 * same term adds to deriv of heating */
		if( hmi.HeatH2Dexc_used < 0. )
			thermal.dCooldT -= hmi.deriv_HeatH2Dexc_used;

		/*  H + H+ => H2+ cooling */
		CoolHeavy.H2PlsCool = (realnum)(MAX2((2.325*phycon.te-1875.)*1e-20,0.)*
		  dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipHYDROGEN][1]*1.66e-11);

		if( h2.lgH2ON )
		{
			/* this is simplified approximation to H2 rotation cooling,
			 * big molecule goes this far better */
			CoolHeavy.h2line = 0.;
		}
		else
		{
			/*  rate for rotation lines from 
			*  >>refer	h2	cool	Lepp, S., & Shull, J.M. 1983, ApJ, 270, 578 */
			x = phycon.alogte - 4.;
			if( phycon.te > 1087. )
			{
				rothi = 3.90e-19*sexp(6118./phycon.te);
			}
			else
			{
				rothi = pow(10.,-19.24 + 0.474*x - 1.247*x*x);
			}

			/*  low density rotation cooling */
			/*&qn = pow(MAX2(hmi.Hmolec[ipMH2g],1e-37),0.77) + 1.2*pow(MAX2(dense.xIonDense[ipHYDROGEN][0],1e-37),0.77);*/
			qn = pow(MAX2(hmi.H2_total,1e-37),0.77) + 1.2*pow(MAX2(dense.xIonDense[ipHYDROGEN][0],1e-37),0.77);
			/* these are equations 11 from LS83 */
			if( phycon.te > 4031. )
			{
				rotlow = 1.38e-22*sexp(9243./phycon.te)*qn;
			}
			else
			{
				rotlow = pow(10.,-22.90 - 0.553*x - 1.148*x*x)*qn;
			}

			CoolHeavy.h2line = 0.;
			if( rotlow > 0. )
				CoolHeavy.h2line += hmi.H2_total*rothi/(1. + rothi/rotlow);
			/* \todo 1 add this from LS83 or (better yet) update to another form.  See Galli & Palla 1998, A5-7. */
			//if( viblow > 0. )
			//	CoolHeavy.h2line += hmi.H2_total*vibhi/(1. + vibhi/viblow);
		}

		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && nzone>187&& iteration > 1/**/)
			{
				fprintf(ioQQQ,"h2coolbug\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
					phycon.te, 
					CoolHeavy.h2line, 
					hmi.H2_total, 
					hmi.Hmolec[ipMHm], 
					hmi.HMinus_photo_rate,
					rothi,
					rotlow );
			}
		}

		/* >>chng 02 mar 07, add DH cooling using rates (eqn 6) from
		 * >>refer	HD	cooling	Puy, D., Grenacher, L, & Jetzer, P., 1999, A&A, 345, 723 */
		factor = sexp(128.6/phycon.te);
		/*CoolHeavy.HD = 2.66e-21 * hydro.D2H_ratio * POW2(hmi.Hmolec[ipMH2g]) * phycon.sqrte *
			factor/(1416.+phycon.sqrte*hmi.Hmolec[ipMH2g] * (1. + 3.*factor));*/
		CoolHeavy.HD = 2.66e-21 * hydro.D2H_ratio * POW2((double)hmi.H2_total) * phycon.sqrte *
			/*factor/(1416.+phycon.sqrte*hmi.Hmolec[ipMH2g] * (1. + 3.*factor));*/
			factor/(1416.+phycon.sqrte*hmi.H2_total * (1. + 3.*factor));
	}

	/* cooling due to charge transfer ionization / recombination */
	CoolAdd("CT C" , 0. , thermal.char_tran_cool );

	/*  H- FB; H + e -> H- + hnu */
	/*  H- FF is in with H ff */
	CoolAdd("H-fb",0,hmi.hmicol);

	/* >>chng 96 nov 15, fac of 2 in deriv to help convergence in very dense
	 * models where H- is important, this takes change in eden into
	 * partial account */
	thermal.dCooldT += 2.*hmi.hmicol*phycon.teinv;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  2 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	CoolAdd("H2ln",0,CoolHeavy.h2line);
	/* >>chng 00 oct 21, added coef of 3.5, sign had been wrong */
	/*thermal.dCooldT += CoolHeavy.h2line*phycon.teinv;*/
	/* >>chng 03 mar 17, change 3.5 to 15 as per behavior in primal.in */
	/*thermal.dCooldT += 3.5*CoolHeavy.h2line*phycon.teinv;*/
	/* >>chng 03 may 18, from 15 to 30 as per behavior in primal.in - overshoots happen */
	/*thermal.dCooldT += 15.*CoolHeavy.h2line*phycon.teinv;*/
	/*>>chng 03 oct 03, from 30 to 3, no more overshoots in primalin */
	/*thermal.dCooldT += 30.*CoolHeavy.h2line*phycon.teinv;*/
	thermal.dCooldT += 3.0*CoolHeavy.h2line*phycon.teinv;

	{
		/* problems with H2 cooling */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /*&& nzone>300 && iteration > 1*/)
		{
			fprintf(ioQQQ,"CoolEvaluate debuggg\t%.2f\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",
			fnzone, 
			phycon.te, 
			hmi.H2_total , 
			CoolHeavy.h2line,
			hmi.Hmolec[ipMHm] , 
			dense.eden);
		}
	}

	CoolAdd("HDro",0,CoolHeavy.HD);
	thermal.dCooldT += CoolHeavy.HD*phycon.teinv;

	CoolAdd("H2+ ",0,CoolHeavy.H2PlsCool);
	thermal.dCooldT += CoolHeavy.H2PlsCool*phycon.teinv;

	/* heating due to three-body, will be incremented in iso_cool*/
	thermal.heating[0][3] = 0.;
	/* heating due to hydrogen lines */
	thermal.heating[0][23] = 0.;
	/* heating due to photoionization of all excited states of hydrogen species */
	thermal.heating[0][1] = 0.;

	/* isoelectronic species cooling, mainly lines, and level ionization */
	for( long int ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long int nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always call iso_cool since we must zero those variables
			* that would have been set had the species been present */
			iso_cool( ipISO , nelem );
		}
	}

	/* >>chng 02 jun 18, don't reevaluate needlessly */
	/* >>chng 03 nov 28, even faster - special logic for when ff is pretty
	 * negligible - eval of ff is pretty slow */
	/* >>chng 04 feb 19, must not test on temp not changing, since ionization
	 * can change at constant temperature 
	 * >>chng 04 sep 14, above introduced bug since brems never reevaluated
	 * now test is zone or temp has changed */
	if( fabs(CoolHeavy.brems_cool_net)/SDIV(thermal.ctot) > conv.HeatCoolRelErrorAllowed/5. || 
	    conv.lgSearch || 
	    !fp_equal(phycon.te,TeUsedBrems) ||
	    nzone != nzoneUsedBrems )
	{
		double BremsThisEnergy;
		/*double OpacityThisIon;*/
		long int limit;
		/* free-free free free brems emission for all ions */

		TeUsedBrems = phycon.te;
		nzoneUsedBrems = nzone;
		/* highest frequency where we have non-zero Boltzmann factors */
		limit = MIN2( rfield.ipMaxBolt , rfield.nflux );

		CoolHeavy.brems_cool_hminus = 0.;
		CoolHeavy.brems_cool_h = 0.;
		CoolHeavy.brems_cool_metals = 0.;
		CoolHeavy.brems_cool_he = 0.;
		CoolHeavy.brems_heat_total = 0.;

		{
			double bhfac, bhMinusfac;
			realnum sumion[LIMELM+1];
			long int ion_lo , ion_hi;

			ASSERT(rfield.ipEnergyBremsThin < rfield.nupper);
			ASSERT(limit < rfield.nupper);

			/* ipEnergyBremsThin is index to energy where gas becomes optically thin to brems,
			 * so this loop is over optically thin frequencies 
			 * do not include optically thick part as net emission since self absorbed */

			/* do hydrogen first, before main loop since want to break out as separate
			 * coolant, and what to add on H- brems */
			CoolHeavy.brems_cool_h = 0.;
			CoolHeavy.brems_cool_hminus = 0.;
			/* this is all done in opacity_addtotal - why do here too? */
			for( long int i=rfield.ipEnergyBremsThin; i < limit; i++ )
			{
				long int ion = 1;

				/* in all following CoolHeavy.lgFreeOn is flag set with 'no free-free' to
				 * turn off brems heating and cooling */
				BremsThisEnergy = rfield.gff[ion][i]*rfield.widflx[i]*rfield.ContBoltz[i];
				/*ASSERT( BremsThisEnergy >= 0. );*/
				CoolHeavy.brems_cool_h += BremsThisEnergy;

				/* for case of hydrogen, add H- brems - OpacStack contains the ratio
				 * of the H- to H brems cross section - multiply by this and H(1s) population */
				CoolHeavy.brems_cool_hminus += BremsThisEnergy * opac.OpacStack[i-1+opac.iphmra];
			}
			bhfac = (dense.xIonDense[ipHYDROGEN][1]+hmi.Hmolec[ipMH2p]+hmi.Hmolec[ipMH3p])*CoolHeavy.lgFreeOn* dense.eden*1.032e-11/phycon.sqrte*EN1RYD;
			bhMinusfac = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop*CoolHeavy.lgFreeOn* dense.eden*1.032e-11/phycon.sqrte*EN1RYD;
			CoolHeavy.brems_cool_h *= bhfac;
			CoolHeavy.brems_cool_hminus *= bhMinusfac;

			/* now do helium, both He+ and He++ */
			CoolHeavy.brems_cool_he = 0.;
			for( long int i=rfield.ipEnergyBremsThin; i < limit; i++ )
			{
				long int nelem = ipHELIUM;
				/* eff. charge is ion, so first rfield.gff argument must be "ion".	*/
				BremsThisEnergy = 
					(dense.xIonDense[nelem][1]*rfield.gff[1][i] + 4.*dense.xIonDense[nelem][2]*rfield.gff[2][i])*
					rfield.widflx[i]*rfield.ContBoltz[i];
				CoolHeavy.brems_cool_he += BremsThisEnergy;
			}
			CoolHeavy.brems_cool_he *=
				CoolHeavy.lgFreeOn * dense.eden*1.032e-11/phycon.sqrte*EN1RYD;

			/* >>chng 05 jul 13, rewrite this for speed */
			/* gaunt factors depend only on photon energy and ion charge, so do
			* sum of ions here before entering into loop over photon energy */
			sumion[0] = 0.;
			for( long int ion=1; ion<=LIMELM; ++ion )
			{
				sumion[ion] = 0.;
				for( long int nelem=ipLITHIUM; nelem < LIMELM; ++nelem )
				{
					if( dense.lgElmtOn[nelem] && ion<=nelem+1 )
					{
						sumion[ion] += dense.xIonDense[nelem][ion];
					}
				}
				/* now include the charge, density, and temperature */
				sumion[ion] *= POW2((realnum)ion);
			}

			// HeH+ is most important ion missing above, add here
			sumion[1] += hmi.Hmolec[ipMHeHp] * POW2((realnum)1.f);

			/* now find lowest and highest ion we need to consider - following loop is over
			 * full continuum and eats time
			 * >>chng 05 oct 19, bounds check had been on ion, rather than ion_lo and ion_hi, so
			 * array bounds were exceeded */
			ion_lo = 1;
			while( sumion[ion_lo]==0 && ion_lo<LIMELM-1 )
				++ion_lo;
			ion_hi = LIMELM;
			while( sumion[ion_hi]==0 && ion_hi>0 )
				--ion_hi;

			/* heavy elements */
			CoolHeavy.brems_cool_metals = 0.;
			CoolHeavy.brems_heat_total = 0.;
			for( long int i=rfield.ipEnergyBremsThin; i < limit; i++ )
			{
				BremsThisEnergy = 0.;
				for(long int ion=ion_lo; ion<=ion_hi; ++ion )
				{
					BremsThisEnergy += sumion[ion]*rfield.gff[ion][i];
				}
				CoolHeavy.brems_cool_metals += BremsThisEnergy*rfield.widflx[i]*rfield.ContBoltz[i];
				/* the total heating due to bremsstrahlung */
				CoolHeavy.brems_heat_total += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu[i]*EN1RYD*CoolHeavy.lgFreeOn;
			}
			CoolHeavy.brems_cool_metals *=
				CoolHeavy.lgFreeOn * dense.eden*1.032e-11/phycon.sqrte*EN1RYD;
			{
				enum {DEBUG_LOC=false};
				if( DEBUG_LOC && nzone>60 /*&& iteration > 1*/)
				{
					double sumfield = 0., sumtot=0., sum1=0., sum2=0.;
					for( long int i=rfield.ipEnergyBremsThin; i<limit;  i++ )
					{
						sumtot += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu[i];
						sumfield += rfield.flux[0][i]*rfield.anu[i];
						sum1 += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu[i];
						sum2 += opac.FreeFreeOpacity[i]*rfield.flux[0][i];
					}
					fprintf(ioQQQ,"DEBUG brems heat\t%.2f\t%.3e\t%.3e\t%.3e\t%e\t%.3e\t%.3e\n",
						fnzone,
						CoolHeavy.brems_heat_total,
						sumtot/SDIV(sumfield) ,
						sum1/SDIV(sum2),
						phycon.te , 
						rfield.gff[1][1218],
						opac.FreeFreeOpacity[1218]);
				}
			}
		}
	}

	/* these two terms are both large, nearly canceling, near lte */
	CoolHeavy.brems_cool_net = 
		CoolHeavy.brems_cool_h + 
		CoolHeavy.brems_cool_he + 
		CoolHeavy.brems_cool_hminus + 
		CoolHeavy.brems_cool_metals - 
		CoolHeavy.brems_heat_total;
	/*fprintf(ioQQQ,"DEBUG brems\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
		fnzone,
		phycon.te,
		CoolHeavy.brems_cool_net,
		CoolHeavy.brems_cool_h ,
		CoolHeavy.brems_cool_he ,
		CoolHeavy.brems_cool_hminus,
		CoolHeavy.brems_cool_metals ,
		CoolHeavy.brems_heat_total);*/

	/* net free free brems cooling, count as cooling if positive */
	CoolAdd( "FF c" , 0, MAX2(0.,CoolHeavy.brems_cool_net) );

	/* now stuff into heating array if negative */
	thermal.heating[0][11] = MAX2(0.,-CoolHeavy.brems_cool_net);

	/* >>chng 96 oct 30, from HFFNet to just FreeFreeCool,
	 * since HeatSum picks up CoolHeavy.brems_heat_total */
	thermal.dCooldT += CoolHeavy.brems_cool_h*thermal.halfte;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  3 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* >>chng 02 jun 21, net cooling already includes this */
	/* end of brems cooling */

	/* heavy element recombination cooling, do not count hydrogenic since
	 * already done above, also helium singlets have been done */
	/* >>chng 02 jul 21, put in charge dependent rec term */
	CoolHeavy.heavfb = 0.;
	for( long int nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* note that detailed iso seq atoms are done in iso_cool */
			long limit_lo = MAX2( 1 , dense.IonLow[nelem] );
			/* there is a -1 on the NISO since loop will be up to <= not < */
			long limit_hi = MIN2( nelem-NISO+1-1 , dense.IonHigh[nelem] );
			for( long int ion=limit_lo; ion<=limit_hi; ++ion )
			{
				/* factor of 0.9 is roughly correct for nebular conditions, see
				 * >>refer	H	rec cooling	LaMothe, J., & Ferland, G.J., 2001, PASP, 113, 165 */
				/* note that ionbal.RR_rate_coef_used is the rate coef, cm3 s-1, needs eden */
				/* >>chng 02 nov 07, move rec arrays around, this now has ONLY rad rec,
				 * previously had included grain rec and three body */
				/* recombination cooling for iso-seq atoms are done in iso_cool */
				double one = dense.xIonDense[nelem][ion] * ionbal.RR_rate_coef_used[nelem][ion-1]*
					dense.eden * phycon.te * BOLTZMANN;
				/*fprintf(ioQQQ,"debugggfb\t%li\t%li\t%.3e\t%.3e\t%.3e\n", nelem, ion, one, 
					dense.xIonDense[nelem][ion] , ionbal.RR_rate_coef_used[nelem][ion]);*/
				CoolHeavy.heavfb += one;
			}
		}
	}

	/*fprintf(ioQQQ,"debuggg hvFB\t%i\t%.2f\t%.2e\t%.2e\n",iteration, fnzone,CoolHeavy.heavfb, dense.eden);*/

	CoolAdd("hvFB",0,CoolHeavy.heavfb);
	thermal.dCooldT += CoolHeavy.heavfb*.113*phycon.teinv;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  4 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* electron-electron brems, approx form from 
	 * >>refer	ee	brems	Stepney and Guilbert, MNRAS 204, 1269 (1983)
	 * ok for T<10**9 */
	CoolHeavy.eebrm = POW2(dense.eden*phycon.te*1.84e-21);

	/* >>chng 97 mar 12, added deriv */
	thermal.dCooldT += CoolHeavy.eebrm*thermal.halfte;
	CoolAdd("eeff",0,CoolHeavy.eebrm);

	/* add advective heating and cooling */
	/* this is cooling due to loss of matter from this region */
	CoolAdd("adve",0,dynamics.Cool() );
	/* >>chng02 dec 04, rm factor of 8 in front of dCooldT */
	thermal.dCooldT += dynamics.dCooldT();
	/* local heating due to matter moving into this location */
	thermal.heating[1][5] = dynamics.Heat();
	thermal.dHeatdT += dynamics.dHeatdT;

	/* total Compton cooling */
	CoolHeavy.tccool = rfield.cmcool*phycon.te;
	CoolAdd("Comp",0,CoolHeavy.tccool);
	thermal.dCooldT += rfield.cmcool;

	/* option for "extra" cooling, expressed as power-law in temperature, these
	 * are set with the CEXTRA command */
	if( thermal.lgCExtraOn )
	{
		CoolHeavy.cextxx = 
			(realnum)(thermal.CoolExtra*pow(phycon.te/1e4,(double)thermal.cextpw));
	}
	else
	{
		CoolHeavy.cextxx = 0.;
	}
	CoolAdd("Extr",0,CoolHeavy.cextxx);

	realnum dDensityDT;

	/* cooling due to wind expansion, only for winds expansion cooling */
	if( wind.lgBallistic() )
	{
		dDensityDT = -(realnum)(wind.AccelTotalOutward/wind.windv + 2.*wind.windv/
		  radius.Radius);
		CoolHeavy.expans = -2.5*pressure.PresGasCurr*dDensityDT;
	}
	else if( dynamics.lgTimeDependentStatic && 
				iteration > dynamics.n_initial_relax)
	{
		dDensityDT = 
			(realnum)((dense.gas_phase[ipHYDROGEN]-dynamics.Upstream_hden)/
			(dynamics.timestep*0.5*(dense.gas_phase[ipHYDROGEN]+dynamics.Upstream_hden)));
		// pdV work term
		CoolHeavy.expans = -pressure.PresGasCurr*dDensityDT;
	}
	else
	{
		dDensityDT = 0.;
		CoolHeavy.expans = 0.;
	}
	CoolAdd("Expn",0,CoolHeavy.expans);
	thermal.dCooldT += CoolHeavy.expans/phycon.te;

	/* cyclotron cooling */
	/* coef is 4/3 /8pi /c * vtherm(elec) */
	CoolHeavy.cyntrn = 4.5433e-25f*magnetic.pressure*PI8*dense.eden*phycon.te;
	CoolAdd("Cycl",0,CoolHeavy.cyntrn);
	thermal.dCooldT += CoolHeavy.cyntrn/phycon.te;

	/* heavy element collisional ionization
	 * derivative should be zero since increased coll ion rate
	 * decreases neutral fraction by proportional amount */
	CoolAdd("Hvin",0,CoolHeavy.colmet);

	/* evaluate H 21 cm spin changing collisions */
	if( !fp_equal(phycon.te,TeEvalCS_21cm) )
	{
		{
			/* this prints table of rates at points given in original data paper */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
#				define N21CM_TE	16
				int n;
				double teval[N21CM_TE]={2.,5.,10.,20.,50.,100.,200.,500.,1000.,
					2000.,3000.,5000.,7000.,10000.,15000.,20000.};
				for( n = 0; n<N21CM_TE; ++n )
				{
					fprintf(
						ioQQQ,"DEBUG 21 cm deex Te=\t%.2e\tH0=\t%.2e\tp=\t%.2e\te=\t%.2e\n",
						teval[n], 
						H21cm_H_atom( teval[n] ),
						H21cm_proton( teval[n] ),
						H21cm_electron( teval[n] ) );
				}
				cdEXIT(EXIT_FAILURE);
#				undef N21CM_TE
			}
		}
		/*only evaluate T dependent part when Te changes, but update
		 * density part below since densities may constantly change */
		atomic_rate_21cm = H21cm_H_atom( phycon.te );
		proton_rate_21cm = H21cm_proton( phycon.te );
		electron_rate_21cm = H21cm_electron( phycon.te );
		TeEvalCS_21cm = phycon.te;
	}
	/* H 21 cm emission/population,
	* cs will be sum of e cs and H cs converted from rate */
	cs = (electron_rate_21cm * dense.eden + 
		    atomic_rate_21cm * dense.xIonDense[ipHYDROGEN][0] +
		    proton_rate_21cm * dense.xIonDense[ipHYDROGEN][1] ) *
			3./	dense.cdsqte;
	PutCS(  cs , &HFLines[0] );

	/* fine structure lines */
	if( !fp_equal(phycon.te,TeEvalCS) )
	{
		/* H 21 cm done above, now loop over remaining lines to get collision strengths */
		for( long int i=1; i < nHFLines; i++ )
		{
			cs = HyperfineCS( i );
			/* now generate the collision strength and put into the line array */
			PutCS(  cs , &HFLines[i] );
		}
		TeEvalCS = phycon.te;
	}

	/* do level pops for H 21 cm which is a special case since Lya pumping in included */
	RT_line_one( &HFLines[0], true,0.f, GetDopplerWidth(dense.AtomicWeight[HFLines[0].Hi->nelem-1]) );
	H21_cm_pops();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  5 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* find total cooling due to hyperfine structure lines */
	hyperfine.cooling_total = HFLines[0].Coll.cool;

	/* now do level pops for all except 21 cm  */
	for( long int i=1; i < nHFLines; i++ )
	{
		/* remember current gas-phase abundance of this isotope */
		realnum save = dense.xIonDense[HFLines[i].Hi->nelem-1][HFLines[i].Hi->IonStg-1];

		/* bail if no abundance */
		if( save<=0. ) 
			continue;

		RT_line_one( &HFLines[i], true,0.f, GetDopplerWidth(dense.AtomicWeight[HFLines[i].Hi->nelem-1]) );

		/* set gas-phase abundance to total times isotope ratio */
		dense.xIonDense[HFLines[i].Hi->nelem-1][HFLines[i].Hi->IonStg-1] *= 
			hyperfine.HFLabundance[i];

		/* use the collision strength generated above and find pops and cooling */
		atom_level2( &HFLines[i] );

		/* put the correct gas-phase abundance back in the array */
		dense.xIonDense[HFLines[i].Hi->nelem-1][HFLines[i].Hi->IonStg-1] = save;

		/* find total cooling due to hyperfine structure lines */
		hyperfine.cooling_total += HFLines[i].Coll.cool;
	}

	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  6 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	realnum xIonDenseSave[LIMELM][LIMELM+1];
	if( atmdat.lgChiantiOn )
	{
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{		
			for( int ion=0; ion<=nelem+1; ++ion )
			{
				xIonDenseSave[nelem][ion] = dense.xIonDense[nelem][ion];
				// zero abundance of species if we are using Chianti for this ion
				if( dense.lgIonChiantiOn[nelem][ion] )
					dense.xIonDense[nelem][ion] = 0.;
			}
		}
	}

	/* Carbon cooling */
	CoolCarb();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  C %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Nitrogen cooling */
	CoolNitr();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  N %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Oxygen cooling */
	CoolOxyg();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  7 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Fluorine cooling */
	CoolFluo();

	/* Neon cooling */
	CoolNeon();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT Ne %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Magnesium cooling */
	CoolMagn();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  8 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Sodium cooling */
	CoolSodi();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT Na %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Aluminum cooling */
	CoolAlum();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT Al %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Silicon cooling */
	CoolSili();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  9 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Phosphorus */
	CoolPhos();

	/* Sulphur cooling */
	CoolSulf();

	/* Chlorine cooling */
	CoolChlo();

	/* Argon cooling */
	CoolArgo();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 10 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Potasium cooling */
	CoolPota();

	/* Calcium cooling */
	CoolCalc();

	/* Scandium cooling */
	CoolScan();

	/* Titanium cooling */
	CoolTita();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 11 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Vanadium cooling */
	CoolVana();

	/* Chromium cooling */
	CoolChro();

	/* Iron cooling */
	CoolIron();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 12 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Manganese cooling */
	CoolMang();

	/* Cobalt cooling */
	CoolCoba();

	/* Nickel cooling */
	CoolNick();

	/* Zinc cooling */
	CoolZinc();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 13 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* do all the thousands of lines Dima added with g-bar approximation */
	CoolDima();
	
	if( atmdat.lgChiantiOn )
	{
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{		
			for( int ion=0; ion<=nelem+1; ++ion )
			{
				dense.xIonDense[nelem][ion] = xIonDenseSave[nelem][ion];
			}
		}
	}

	/* do external database lines */
 	dBase_solve();
	/* Print number of levels for each species */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			static bool lgMustPrintHeader=true;
			if( lgMustPrintHeader )
			{
				lgMustPrintHeader = false;
				printf("DEBUG Levels\t%s",Species[0].chLabel );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					printf("\t%s",Species[ipSpecies].chLabel );
				}
				printf("\n" );
				printf("DEBUG Max\t%li" ,Species[0].numLevels_max );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					printf( "\t%li" ,Species[ipSpecies].numLevels_max );
				}
				printf("\n");
			}
			printf("DEBUG Local\t%li" ,Species[0].numLevels_local );
			for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
			{
				printf("\t%li" ,Species[ipSpecies].numLevels_local );
			}
			printf("\n");
		}
	}

	/* now add up all the coolants */
	CoolSum(tot);
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 14 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* negative cooling */
	if( *tot <= 0. )
	{
		fprintf( ioQQQ, " COOLR; cooling is <=0, this is impossible.\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* bad derivative */
	if( thermal.dCooldT == 0. )
	{
		fprintf( ioQQQ, " COOLR; cooling slope <=0, this is impossible.\n" );
		if( *tot > 0. && dense.gas_phase[ipHYDROGEN] < 1e-4 )
		{
			fprintf( ioQQQ, " Probably due to very low density.\n" );
		}
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace )
	{
		fndstr(*tot,thermal.dCooldT);
	}

	/* lgTSetOn true for constant temperature model */
	if( (((((!thermal.lgTemperatureConstant) && *tot < 0.) && called.lgTalk) && 
	  !conv.lgSearch) && thermal.lgCNegChk) && nzone > 0 )
	{
		fprintf( ioQQQ, 
			" NOTE Negative cooling, zone %4ld, =%10.2e coola=%10.2e CHION=%10.2e Te=%10.2e\n", 
		  nzone, 
		  *tot, 
		  iso.cLya_cool[ipH_LIKE][ipHYDROGEN], 
		  iso.coll_ion[ipH_LIKE][ipHYDROGEN], 
		  phycon.te );
		fndneg();
	}

	/* possibility of getting empirical cooling derivative
	 * normally false, set true with 'set numerical derivatives' command */
	if( NumDeriv.lgNumDeriv )
	{
		if( ((nzone > 2 && nzone == nzSave) && ! fp_equal( oldtemp, phycon.te )) && nhit > 4 )
		{
			/* hnit is number of tries on this zone - use to stop numerical problems
			 * do not evaluate numerical deriv until well into solution */
			deriv = (oltcool - *tot)/(oldtemp - phycon.te);
			thermal.dCooldT = deriv;
		}
		else
		{
			deriv = thermal.dCooldT;
		}
		if( nzone != nzSave )
			nhit = 0;

		nzSave = nzone;
		nhit += 1;
		oltcool = *tot;
		oldtemp = phycon.te;
	}
	return;
}

/*  */
#ifdef EPS
#	undef EPS
#endif
#define	EPS	0.01

/*fndneg search cooling array to find negative values */
STATIC void fndneg(void)
{
	long int i;
	double trig;

	DEBUG_ENTRY( "fndneg()" );

	trig = fabs(thermal.htot*EPS);
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( thermal.cooling[i] < 0. && fabs(thermal.cooling[i]) > trig )
		{
			fprintf( ioQQQ, " negative line=%s %.2f fraction of heating=%.3f\n", 
			  thermal.chClntLab[i], thermal.collam[i], thermal.cooling[i]/
			  thermal.htot );
		}

		if( thermal.heatnt[i] > trig )
		{
			fprintf( ioQQQ, " heating line=%s %.2f fraction of heating=%.3f\n", 
			  thermal.chClntLab[i], thermal.collam[i], thermal.heatnt[i]/
			  thermal.htot );
		}
	}
	return;
}

/*fndstr search cooling stack to find strongest values */
STATIC void fndstr(double tot, 
  double dc)
{
	char chStrngLab[NCOLNT_LAB_LEN+1];
	long int i;
	realnum wl;
	double str, 
	  strong;

	DEBUG_ENTRY( "fndstr()" );

	strong = 0.;
	wl = -FLT_MAX;
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( fabs(thermal.cooling[i]) > strong )
		{
			/* this is the wavelength of the coolant, 0 for a continuum*/
			wl = thermal.collam[i];
			/* make sure labels are all valid*/
			/*>>chng 06 jun 06, bug fix, assert length was ==4, should be <=NCOLNT_LAB_LEN */
			ASSERT( strlen( thermal.chClntLab[i] ) <= NCOLNT_LAB_LEN );
			strcpy( chStrngLab, thermal.chClntLab[i] );
			strong = fabs(thermal.cooling[i]);
		}
	}

	str = strong;

	fprintf( ioQQQ, 
		"   fndstr cool: TE=%10.4e Ne=%10.4e C=%10.3e dC/dT=%10.2e ABS(%s %.1f)=%.2e nz=%ld\n", 
	  phycon.te, dense.eden, tot, dc, chStrngLab
	  , wl, str, nzone );

	/* option for extensive printout of lines */
	if( trace.lgCoolTr )
	{
		realnum ratio;

		/* flag all significant coolants, first zero out the array */
		coolpr(ioQQQ,(char*)thermal.chClntLab[0],1,0.,"ZERO");

		/* push all coolants onto the stack */
		for( i=0; i < thermal.ncltot; i++ )
		{
			/* usually positive, although can be neg for coolants that heats, 
			 * only do positive here */
			ratio = (realnum)(thermal.cooling[i]/thermal.ctot);
			if( ratio >= EPS )
			{
				/*>>chng 99 jan 27, only cal when ratio is significant */
				coolpr(ioQQQ,(char*)thermal.chClntLab[i],thermal.collam[i], ratio,"DOIT");
			}
		}

		/* complete the printout for positive coolants */
		coolpr(ioQQQ,"DONE",1,0.,"DONE");

		/* now do cooling that was counted as a heat source if significant */
		if( thermal.heating[0][22]/thermal.ctot > 0.05 )
		{
			fprintf( ioQQQ, 
				"     All coolant heat greater than%6.2f%% of the total will be printed.\n", 
			  EPS*100. );

			coolpr(ioQQQ,"ZERO",1,0.,"ZERO");
			for( i=0; i < thermal.ncltot; i++ )
			{
				ratio = (realnum)(thermal.heatnt[i]/thermal.ctot);
				if( fabs(ratio) >=EPS )
				{
					coolpr(ioQQQ,(char*)thermal.chClntLab[i],thermal.collam[i],
					  ratio,"DOIT");
				}
			}
			coolpr(ioQQQ,"DONE",1,0.,"DONE");
		}
	}
	return;
}
