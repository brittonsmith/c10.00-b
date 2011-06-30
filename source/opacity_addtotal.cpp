/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityAddTotal derive total opacity for this position,
 * called by ConvBase */
#include "cddefines.h"
#include "physconst.h"
#include "iso.h"
#include "grainvar.h"
#include "ca.h"
#include "rfield.h"
#include "oxy.h"
#include "mole.h"
#include "dense.h"
#include "atoms.h"
#include "conv.h"
#include "ionbal.h"
#include "trace.h"
#include "hmi.h"
#include "phycon.h"
#include "opacity.h"
#include "taulines.h"

void OpacityAddTotal(void)
{
	long int i, 
	  ion,
	  ipHi,
	  ipISO,
	  ipop,
	  limit,
	  low,
	  nelem,
	  n; 
	double DepartCoefInv ,
	  fac, 
	  sum;
	realnum SaveOxygen1 , 
		SaveCarbon1;

	DEBUG_ENTRY( "OpacityAddTotal()" );

	/* OpacityZero will zero out scattering and absorption opacities,
	 * and set OldOpacSave to opac to save it */
	OpacityZero();

	/* free electron scattering opacity, Compton recoil energy loss */
	for( i=0; i < rfield.nflux; i++ )
	{
		/* scattering part of total opacity */
		opac.opacity_sct[i] += opac.OpacStack[i-1+opac.iopcom]*
		  dense.eden;
	}

	/* opacity due to Compton bound recoil ionization */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{ 
			for( ion=0; ion<nelem+1; ++ion )
			{
				realnum factor = dense.xIonDense[nelem][ion];
				/*>>chng 05 nov 26, add molecular hydrogen assuming same
				 * as two free hydrogen atoms - hence mult density by two
				 *>>KEYWORD	H2 bound Compton ionization */
				if( nelem==ipHYDROGEN )
					factor += hmi.H2_total*2.f;
				if( factor > 0. )
				{
					// loop_min and loop_max are needed to work around a bug in icc 10.0
					long loop_min = ionbal.ipCompRecoil[nelem][ion]-1;
					long loop_max = rfield.nflux;
					/* ionbal.nCompRecoilElec number of electrons in valence shell
					 * that can compton recoil ionize */
					factor *= ionbal.nCompRecoilElec[nelem-ion];
					for( i=loop_min; i < loop_max; i++ )
					{
						/* add in bound hydrogen electron scattering, treated as absorption */
						opac.opacity_abs[i] += opac.OpacStack[i-1+opac.iopcom]*factor;
					}
				}
			}
		}
	}

	/* opacity due to pair production - does not matter what form these
	 * elements are in */
	/** \todo	2	add charged heavy elements */
	sum = dense.gas_phase[ipHYDROGEN] + 4.*dense.gas_phase[ipHELIUM];
	OpacityAdd1Subshell(opac.ioppr,opac.ippr,rfield.nflux,(realnum)sum,'s');

	/* free-free free free brems emission for all ions */

	/* >>chng 02 jul 21, use full set of ions and gaunt factor */
	/* ipEnergyBremsThin is index to energy where gas goes optically thin to brems,
	 * so this loop is over optically thick frequencies */
	static double *TotBremsAllIons;
	static bool lgFirstTime=true;
	double BremsThisEner,bfac, sumion[LIMELM+1];
	long int ion_lo , ion_hi;

	if( lgFirstTime )
	{
		/* rfield.nupper will not change in one coreload, so just malloc this once */
		TotBremsAllIons = (double *)MALLOC((unsigned long)rfield.nupper*sizeof(double));
		lgFirstTime = false;
	}
	bfac = (dense.eden/1e20)/phycon.sqrte/1e10;
	/* gaunt factors depend only on photon energy and ion charge, so do
	 * sum of ions here before entering into loop over photon energy */
	sumion[0] = 0.;
	for(ion=1; ion<=LIMELM; ++ion )
	{
		sumion[ion] = 0.;
		for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] && ion<=nelem+1 )
			{
				sumion[ion] += dense.xIonDense[nelem][ion];
			}
		}
		/* now include the charge, density, and temperature */
		sumion[ion] *= POW2((double)ion)*bfac;
	}

	// HeH+ is most important ion missing below, add here
	sumion[1] += hmi.Hmolec[ipMHeHp] * POW2((realnum)1.f)*bfac;

	/* now find lowest and highest ion we need to consider - following loop is over
	 * full continuum and eats time 
	 * >>chng 05 oct 20, following had tests on ion being within bounds, bug,
	 * changed to ion_lo and ion_hi */
	ion_lo = 1;
	while( sumion[ion_lo]==0 && ion_lo<LIMELM-1 )
		++ion_lo;
	ion_hi = LIMELM;
	while( sumion[ion_hi]==0 && ion_hi>0 )
		--ion_hi;

	for( i=0; i < rfield.nflux; i++ )
	{
		/* do hydrogen first, before main loop since want to add on H- brems */
		nelem = ipHYDROGEN;
		ion = 1;/* >>chng 02 nov 07, add charged ions as H+ */

		BremsThisEner =  bfac * rfield.gff[ion][i] * (dense.xIonDense[ipHYDROGEN][1]+hmi.Hmolec[ipMH2p]+hmi.Hmolec[ipMH3p]);
		/* for case of hydrogen, add H- brems - OpacStack contains the ratio
		 * of the H- to H brems cross section - multiply by this and H(1s) population */
		BremsThisEner += bfac * rfield.gff[ion][i] * opac.OpacStack[i-1+opac.iphmra] * StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;

		/* there is at least one standard test (grainlte) which has ZERO ionization -
		 * this assert will pass that test (the ==0) and also the usual case */
		/* ASSERT( (dense.xIonDense[nelem][ion]==0.) || (TotBremsAllIons[i] > 0.) );*/

		/* chng 02 may 16, by Ryan...do all brems for all ions in one fell swoop,
		* using gaunt factors from rfield.gff.	*/
		/* >>chng 05 jul 11 reorganize loop for speed */
		for(ion=ion_lo; ion<=ion_hi; ++ion )
		{
			BremsThisEner += sumion[ion]*rfield.gff[ion][i];
		}
		ASSERT( !isnan( BremsThisEner ) );
		TotBremsAllIons[i] = BremsThisEner;

		/* >>>chng 05 jul 12, bring these two loops together */
		if( rfield.ContBoltz[i] < 0.995 )
		{
			TotBremsAllIons[i] *= opac.OpacStack[i-1+opac.ipBrems]*
				(1. - rfield.ContBoltz[i]);
		}
		else
		{
			TotBremsAllIons[i] *= opac.OpacStack[i-1+opac.ipBrems]*
				rfield.anu[i]*TE1RYD/phycon.te;
		}
		opac.FreeFreeOpacity[i] = TotBremsAllIons[i];
		/* >>chng 02 jul 23, from >0 to >=0, some models have no ionization,
		 * like grainlte.in */
		/*ASSERT( (opac.FreeFreeOpacity[i] > 0.) || (dense.xIonDense[ipHYDROGEN][1] == 0.) );*/
		opac.opacity_abs[i] += opac.FreeFreeOpacity[i];
	}


	/* H minus absorption, with correction for stimulated emission */
	if( hmi.hmidep > SMALLFLOAT )
	{
		DepartCoefInv = 1./hmi.hmidep;
	}
	else
	{
		/* the hmidep departure coef can become vastly small in totally
		 * neutral gas (no electrons) */
		DepartCoefInv = 1.;
	}

	limit = iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0];
	if(limit > rfield.nflux)
		limit = rfield.nflux;

	for( i=hmi.iphmin-1; i < limit; i++ )
	{
		double factor;
		factor = 1. - rfield.ContBoltz[i]*DepartCoefInv;
		if(factor > 0)
			opac.opacity_abs[i] += opac.OpacStack[i-hmi.iphmin+opac.iphmop]*
											 hmi.Hmolec[ipMHm]*factor;
	}

	/* H2P h2plus bound free opacity */
	limit = opac.ih2pnt[1]; 
	if(limit > rfield.nflux)
		limit = rfield.nflux;
	for( i=opac.ih2pnt[0]-1; i < limit; i++ )
	{
		opac.opacity_abs[i] += hmi.Hmolec[ipMH2p]*opac.OpacStack[i-opac.ih2pnt[0]+opac.ih2pof];
		ASSERT( !isnan( opac.opacity_abs[i] ) );
	}

	/* get total population of hydrogen ground to do Rayleigh scattering */
	if( dense.xIonDense[ipHYDROGEN][1] <= 0. )
	{
		fac = dense.xIonDense[ipHYDROGEN][0];
	}
	else
	{
		fac = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;
	}

	/* Ly a damp wing opac (Rayleigh scattering) */
	limit = iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s];
	if(limit > rfield.nflux)
		limit = rfield.nflux;
	for( i=0; i < limit; i++ )
	{
		opac.opacity_sct[i] += (fac*opac.OpacStack[i-1+opac.ipRayScat]);
	}

	 /* remember largest correction for stimulated emission */
	if( iso.DepartCoef[ipH_LIKE][ipHYDROGEN][ipH1s] > 1e-30 && !conv.lgSearch )
	{
		realnum factor;
		factor = (realnum)(rfield.ContBoltz[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1]/iso.DepartCoef[ipH_LIKE][ipHYDROGEN][ipH1s]);
		if(opac.stimax[0] < factor)
			opac.stimax[0] = factor;
	}

	if( iso.DepartCoef[ipH_LIKE][ipHYDROGEN][ipH2p] > 1e-30 && !conv.lgSearch )
	{
		realnum factor;
		factor = (realnum)(rfield.ContBoltz[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH2p]-1]/iso.DepartCoef[ipH_LIKE][ipHYDROGEN][ipH2p]);
		if(opac.stimax[1] < factor)
			opac.stimax[1] = factor;
	}

#	if 0
	/* check whether hydrogen or Helium singlets mased, if not in search mode */
	if( !conv.lgSearch )
	{
		if( Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].PopOpc < 0. )
		{
			hydro.lgHLyaMased = true;
		}
	}
#	endif

	/* >>chng 05 nov 25, use Yan et al. H2 photo cs
	 * following reference gives cross section for all energies
	 * >>refer	H2	photo cs	Yan, M., Sadeghpour, H.R., & Dalgarno, A., 1998, ApJ, 496, 1044 
	 * Wilms, J., Allen, A., & McCray, R. 2000, ApJ, 542, 914 */
	/* >>chng 02 jan 16, approximate inclusion of H_2 photoelectric opacity 
	 * include H_2 in total photoelectric opacity as twice H0 cs */
	/* set lower and upper limits to this range */
	/*>>KEYWORD	H2	photoionization opacity */
	low = opac.ipH2_photo_thresh;
	ipHi = rfield.nupper;
	ipop = opac.ipH2_photo_opac_offset;
	/* OpacityAdd1Subshell just returns for static opacities if opac.lgRedoStatic not set*/
	/* >>chng 05 nov 27, change on nov 25 had left 2*density from H0, so
	 * twice the H2 density was used - 	 * also changed to static opacity 
	 * this assumes that all v,J levels contribute the same opacity */
	OpacityAdd1Subshell( ipop , low , ipHi , hmi.H2_total , 's' );

	/*>>KEYWORD	CO photoionization opacity */
	/* include photoionization of CO - assume C and O in CO each have 
	 * same photo cs as atom - this should only be significant in highly
	 * shielded regions where only very hard photons penetrate 
	 * also H2O condensed onto grain surfaces - very important deep in cloud */
	SaveOxygen1 = dense.xIonDense[ipOXYGEN][0];
	SaveCarbon1 = dense.xIonDense[ipCARBON][0];
	dense.xIonDense[ipOXYGEN][0] += findspecies("CO")->hevmol + findspecies("H2Ogrn")->hevmol;
	dense.xIonDense[ipCARBON][0] += findspecies("CO")->hevmol;

	/* following loop adds standard opacities for first 30 elements 
	 * most heavy element opacity added here */
	for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		/* this element may be turned off */
		if( dense.lgElmtOn[nelem] )
		{ 
			OpacityAdd1Element(nelem);
		}
	}

	/* now reset the abundances */
	dense.xIonDense[ipOXYGEN][0] = SaveOxygen1;
	dense.xIonDense[ipCARBON][0] = SaveCarbon1;

	/* following are opacities due to specific excited levels */

	/* nitrogen opacity
	 * excited level of N+ */
	OpacityAdd1Subshell(opac.in1[2],opac.in1[0],opac.in1[1],
		dense.xIonDense[ipNITROGEN][0]*atoms.p2nit , 'v' );

	/* oxygen opacity
	 * excited level of Oo */
	OpacityAdd1Subshell(opac.ipo1exc[2],opac.ipo1exc[0],opac.ipo1exc[1],
	  dense.xIonDense[ipOXYGEN][0]*oxy.poiexc,'v');

	/* O2+ excited states */
	OpacityAdd1Subshell(opac.ipo3exc[2],opac.ipo3exc[0],opac.ipo3exc[1],
	  dense.xIonDense[ipOXYGEN][2]*oxy.poiii2,'v');

	OpacityAdd1Subshell(opac.ipo3exc3[2],opac.ipo3exc3[0],opac.ipo3exc3[1],
	  dense.xIonDense[ipOXYGEN][2]*oxy.poiii3,'v');

	/* magnesium opacity
	 * excited level of Mg+ */
	OpacityAdd1Subshell(opac.ipOpMgEx,opac.ipmgex,iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s],
		dense.xIonDense[ipMAGNESIUM][1]* atoms.popmg2,'v');

	/* calcium opacity
	 * photoionization of excited levels of Ca+ */
	OpacityAdd1Subshell(opac.ica2op,opac.ica2ex[0],opac.ica2ex[1],
	  ca.popca2ex,'v');

	/*******************************************************************
	 *
	 * complete evaluation of total opacity by adding in the static part and grains
	 *
	 *******************************************************************/

	/* this loop defines the variable iso.ConOpacRatio[ipH_LIKE][nelem][n],
	 * the ratio of not H to Hydrogen opacity.  for grain free environments
	 * at low densities this is nearly zero.  The correction includes 
	 * stimulated emission correction */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* this element may be turned off */
			if( dense.lgElmtOn[nelem] )
			{ 
				/* this branch is for startup only */
				if( nzone < 1 )
				{
					/* >>chng 06 aug 17, should go to numLevels_local instead of _max */
					for( n=0; n < iso.numLevels_local[ipISO][nelem]; n++ )
					{
						if(iso.ipIsoLevNIonCon[ipISO][nelem][n] < rfield.nflux )
						{
							/*>>chng 04 dec 12, had tested against 1e-30, set to zero if not
							 * greater than this - caused oscillations as opacity fell below
							 * and around this value - change to greater than 0 */
							/*if( opac.opacity_abs[iso.ipIsoLevNIonCon[ipISO][nelem][n]-1] > 1e-30 )*/
							if( opac.opacity_abs[iso.ipIsoLevNIonCon[ipISO][nelem][n]-1] > 0. )
							{
								/* >>chng 02 may 06, use general form of threshold cs */
								/*double t1 = atmdat_sth(n)/(POW2(nelem+1.-ipISO));*/
								long int ip = iso.ipIsoLevNIonCon[ipISO][nelem][n];
								double t2 = csphot(
									ip ,
									ip ,
									iso.ipOpac[ipISO][nelem][n] );

								iso.ConOpacRatio[ipISO][nelem][n] = 1.-
									(StatesElemNEW[nelem][nelem-ipISO][n].Pop*t2)/
									opac.opacity_abs[ip-1];
							}
							else
							{
								iso.ConOpacRatio[ipISO][nelem][n] = 0.;
							}
						}
					}
				}
				/* end branch for startup only, start branch for all zones including startup */
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max */
				for( n=0; n < iso.numLevels_local[ipISO][nelem]; n++ )
				{
					/* ratios of other to total opacity for continua of all atoms done with iso model */
					if(iso.ipIsoLevNIonCon[ipISO][nelem][n] < rfield.nflux )
					{
						/*>>chng 04 dec 12, had tested against 1e-30, set to zero if not
						 * greater than this - caused oscillations as opacity fell below
						 * and around this value - change to greater than 0 */
						/*if( opac.opacity_abs[iso.ipIsoLevNIonCon[ipISO][nelem][n]-1] > 1e-30 )*/
						if( opac.opacity_abs[iso.ipIsoLevNIonCon[ipISO][nelem][n]-1] > 0. )
						{
							/* first get departure coef */
							if( iso.DepartCoef[ipISO][nelem][n] > 1e-30 && (!conv.lgSearch ) )
							{
								/* this is the usual case, use inverse of departure coef */
								fac = 1./iso.DepartCoef[ipISO][nelem][n];
							}
							else if( conv.lgSearch )
							{
								/* do not make correction for stim emission during search
								* for initial temperature solution, since trys are very non-equil */
								fac = 0.;
							}
							else
							{
								fac = 1.;
							}

							/** \todo	1	stupid - why this test on opacity_abs ? - we only get here
							 * if we already passed above test on this very thing */
							/* now get opaicty ratio with correction for stimulated emission */
							/*>>chng 04 dec 12, had tested against 1e-30, set to zero if not
							 * greater than this - caused oscillations as opacity fell below
							 * and around this value - change to greater than 0 */
							/*if( opac.opacity_abs[iso.ipIsoLevNIonCon[ipISO][nelem][n]-1] > 1e-30 )*/
							if( opac.opacity_abs[iso.ipIsoLevNIonCon[ipISO][nelem][n]-1] > 0. )
							{
								/* >>chng 02 may 06, use general form of threshold cs */
								long int ip = iso.ipIsoLevNIonCon[ipISO][nelem][n];

								double t2 = csphot(
									ip ,
									ip ,
									iso.ipOpac[ipISO][nelem][n] );

								double opacity_this_species = 
									StatesElemNEW[nelem][nelem-ipISO][n].Pop*t2*
									(1. - fac*rfield.ContBoltz[ip-1]);

								double opacity_fraction =  1. - opacity_this_species / opac.opacity_abs[ip-1];
								if(opacity_fraction < 0)
									opacity_fraction = 0.;

								/* use mean of old and new ratios */
								iso.ConOpacRatio[ipISO][nelem][n] = 
									iso.ConOpacRatio[ipISO][nelem][n]* 0.75 + 0.25*opacity_fraction;

								if(iso.ConOpacRatio[ipISO][nelem][n] < 0.)
									iso.ConOpacRatio[ipISO][nelem][n] = 0.;
							}
							else
							{
								iso.ConOpacRatio[ipISO][nelem][n] = 0.;
							}
						}
						else
						{
							iso.ConOpacRatio[ipISO][nelem][n] = 0.;
						}
					}
					else
					{
						iso.ConOpacRatio[ipISO][nelem][n] = 0.;
					}
				}
			}
		}
	}

	/* add dust grain opacity if dust present */
	if( gv.lgDustOn() )
	{
		/* generate current grain opacities since may be function of depth */
		/* >>chng 01 may 11, removed code to update grain opacities, already done by GrainChargeTemp */
		for( i=0; i < rfield.nflux; i++ )
		{
			/* units cm-1 */
			opac.opacity_sct[i] += gv.dstsc[i]*dense.gas_phase[ipHYDROGEN];
			opac.opacity_abs[i] += gv.dstab[i]*dense.gas_phase[ipHYDROGEN];
		}
	}

	/* check that opacity is sane */
	for( i=0; i < rfield.nflux; i++ )
	{
		/* OpacStatic was zeroed in OpacityZero, incremented in opacityadd1subshell */
		opac.opacity_abs[i] += opac.OpacStatic[i];
		/* make sure that opacity is positive */
		/*ASSERT( opac.opacity_abs[i] > 0. );*/
	}

	/* compute gas albedo here */
	for( i=0; i < rfield.nflux; i++ )
	{
		opac.albedo[i] = opac.opacity_sct[i]/
			(opac.opacity_sct[i] + opac.opacity_abs[i]);
	}

	/* during search phase set old opacity array to current value */
	if( conv.lgSearch )
		OpacityZeroOld();

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "     OpacityAddTotal returns; grd rec eff (opac) for Hn=1,4%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e HeI,II:%10.2e%10.2e\n", 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH1s], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH2s], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH2p], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH3s], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH3p], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH3d], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH4s], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH4p], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH4d], 
		  iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][ipH4f],
		  iso.ConOpacRatio[ipHE_LIKE][ipHELIUM][ipHe1s1S],
		  iso.ConOpacRatio[ipH_LIKE][ipHELIUM][ipH1s] );
	}

	{
		/* following should be set true to print strongest ots contributors */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && (nzone>=378)/**/ )
		{
			if( nzone > 380 ) 
				cdEXIT( EXIT_FAILURE );
			for( i=0; i<rfield.nflux; ++i )
			{
				fprintf(ioQQQ,"rtotsbugggg\t%li\t%.3e\t%.3e\t%.3e\t%.3e\n",
					conv.nPres2Ioniz,
					rfield.anu[i],
					opac.opacity_abs[i],
					rfield.otscon[i],
					rfield.otslin[i]);
			}
		}
	}
	return;
}
