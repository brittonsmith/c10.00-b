/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PresTotCurrent determine the gas and line radiation pressures for current conditions,
 * this sets values of pressure.PresTotlCurr, also calls tfidle */
#include "cddefines.h"
#include "physconst.h"
#include "taulines.h"
#include "opacity.h"
#include "hextra.h"
#include "elementnames.h"
#include "hydrogenic.h"
#include "conv.h"
#include "iso.h"
#include "wind.h"
#include "magnetic.h"
#include "doppvel.h"
#include "rfield.h"
#include "phycon.h"
#include "thermal.h"
#include "hmi.h"
#include "h2.h"
#include "dense.h"
#include "atomfeii.h"
#include "mole.h"
#include "dynamics.h"
#include "trace.h"
#include "rt.h"
#include "atmdat.h"
#include "lines_service.h"
#include "pressure.h"
#include "radius.h"

/* this sets values of pressure.PresTotlCurr, also calls tfidle */
void PresTotCurrent(void)
{
	static long int
	  /* used in debug print statement to see which line adds most pressure */
	  ipLineTypePradMax=-1 , 
	  /* points to line causing radiation pressure */
	  ipLinePradMax=-LONG_MAX,
	  ip2=-LONG_MAX,
	  ip3=-LONG_MAX,
	  ip4=-LONG_MAX;

	/* the line radiation pressure variables that must be preserved since
	 * a particular line radiation pressure may not be evaluated if it is
	 * not important */
	static double Piso_seq[NISO]={0.,0.},
		PLevel1=0.,
		PLevel2=0.,
		PHfs=0.,
		PCO=0.,
		P_H2=0.,
		P_FeII=0.,
		P_dBase=0.;

	double 
		RadPres1, 
		TotalPressure_v, 
		pmx;
	realnum den_hmole ,
		den_comole ,
		DenseAtomsIons,
		DensElements;

	DEBUG_ENTRY( "PresTotCurrent()" );

	if( !conv.nTotalIoniz )
	{
		/* resetting ipLinePradMax, necessary for 
		 * multiple cloudy calls in single coreload. */ 
		ipLinePradMax=-LONG_MAX;
		//pressure.PresTotlCurr = 0.;
	}

	/* PresTotCurrent - evaluate total pressure, dyne cm^-2
	 * and radiative acceleration */

	/* several loops over the emission lines for radiation pressure and
	 * radiative acceleration are expensive due to the number of lines and
	 * the memory they occupy.  Many equations of state do not include
	 * radiation pressure or radiative acceleration.  Only update these
	 * when included in EOS - when not included only evaluate once per zone.
	 * do it once per zone since we will still report these quantities.
	 * this flag says whether we must update all terms */
	bool lgMustEvaluate = false;

	/* this says we already have an evaluation in this zone so do not 
	 * evaluate terms known to be trivial, even when reevaluating major terms */
	bool lgMustEvaluateTrivial = false;
	/* if an individual agent is larger than this fraction of the total current
	 * radiation pressure then it is not trivial 
	 * conv.PressureErrorAllowed is relative error allowed in pressure */
	double TrivialLineRadiationPressure = conv.PressureErrorAllowed/10. * 
		pressure.pres_radiation_lines_curr;

	/* reevaluate during search phase when conditions are changing dramatically */
	if( conv.lgSearch )
	{
		lgMustEvaluate = true;
		lgMustEvaluateTrivial = true;
	}

	/* reevaluate if zone or iteration has changed */
	static long int nzoneEvaluated=0, iterationEvaluated=0;
	if( nzone!=nzoneEvaluated || iteration!=iterationEvaluated )
	{
		lgMustEvaluate = true;
		lgMustEvaluateTrivial = true;
		/* this flag says which source of radiation pressure was strongest */
		ipLineTypePradMax = -1;
	}

	/* constant pressure or dynamical sim - we must reevaluate terms 
	 * but do not need to reevaluate trivial contributors */
	if( (strcmp(dense.chDenseLaw,"WIND") == 0 ) ||
		(strcmp(dense.chDenseLaw,"CPRE") == 0 ) )
		lgMustEvaluate = true;

	if( lgMustEvaluate )
	{
		/* we are reevaluating quantities in this zone and iteration,
		 * remember that we did it */
		nzoneEvaluated = nzone;
		iterationEvaluated = iteration;
	}

	/* update density - temperature variables which may have changed */
	/* must call TempChange since ionization has changed, there are some
	* terms that affect collision rates (H0 term in electron collision) */
	TempChange(phycon.te , false);

	/* find total number of atoms and ions */
	DenseAtomsIons = 0.;
	DensElements = 0.;
	realnum factor = 1.1f;
	if( conv.lgSearch )
		factor = 2.f;
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		/* only check element solution if it is on, and ionization 
		 * has not been set with TABLE ELEMENT command */
		if( dense.lgElmtOn[nelem]  )
		{
			/* gas_phase is density of all atoms and ions, but not molecules */
			DensElements += dense.gas_phase[nelem];
			realnum DenseOne = 0;
			for( long ion=0; ion<=nelem+1; ++ion )
				DenseOne += dense.xIonDense[nelem][ion];

			/* during search phase sums may not add up correctly, and during
			 * early parts of search phase may be badly off.  This test 
			 * was introduced as result off fpe due to another
			 * reason - Te had changed too much during initial search
			 * for sim in which chemistry was important - fix was to not
			 * let Te change.  During resulting insanity caused by large
			 * change, linearization did not work, CO and ionization solvers
			 * could diverge leading to molecule or ionization population
			 * larger than 1e38 limit to a realnum, fpe followed.  this is to
			 * catch that in its early stages */
			if( conv.nTotalIoniz>5 && !dense.lgSetIoniz[nelem] &&
				DenseOne > dense.gas_phase[nelem]*factor ) 
			{
				/* species is not conserved */
				fprintf(ioQQQ,"\n\n PROBLEM DISASTER PressureTotal: the chemical species "
					"are not conserved.\n");
				fprintf(ioQQQ," The sum of the densities of the atoms and ions for "
					"%s is %.3e which is greater than the gas-phase abundance "
					"of the element, %.3e.\n",
					elementnames.chElementName[nelem],DenseOne,
					dense.gas_phase[nelem]);
				/* not including cosmic rays is the most likely cause of this */ 
				if( hextra.cryden == 0. )
				{
					fprintf(ioQQQ," The chemistry network is known to fail this way"
					" in cold molecular environments when cosmic rays are not"
					" included.  They were not - try including them with the"
					" COSMIC RAY BACKBROUND command.\n");
					fprintf(ioQQQ," The calculation is aborting. conv.nTotalIoniz=%li\n"
						" Sorry.\n\n", conv.nTotalIoniz );
					lgAbort = true;
				}
				/* if they were included then something seriously wrong has happened*/
				else
					TotalInsanity(); 
			}
			DenseAtomsIons += DenseOne;
		}
	}
	/* DensElements is sum of all gas-phase densities of all elements,
	 * DenseAtomsIons is sum of density of atoms and ions, does not 
	 * include molecules */
	ASSERT( conv.nTotalIoniz<5 || DenseAtomsIons <= DensElements*factor );
	ASSERT( DenseAtomsIons > 0. );

	/* evaluate the total ionization energy on a scale where neutral atoms
	 * correspond to an energy of zero, so the ions have a positive energy */
	phycon.EnergyIonization = 0.;
#if 0
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
		{
			/* lgMETALS is true, set false with "no advection metals" command */
			int kadvec = dynamics.lgMETALS;
			/* this gives the iso sequence for this ion - should it be included in the
			 * advected energy equation? lgISO is true, set false with 
			 * "no advection H-like" and "he-like" - for testing*/
			ipISO = nelem-ion;
			fixit(); /* should this be kadvec = kadvec && dynamics.lgISO[ipISO]; ? */
			if( ipISO >= 0 && ipISO<NISO )
				kadvec = dynamics.lgISO[ipISO];
			for( long i=1; i<=ion; ++i )
			{
				long int nelec = nelem+1 - i;
				/* this is the sum of all the energy needed to bring the atom up
				 * to the ion+1 level of ionization - at this point a positive number */
				phycon.EnergyIonization += dense.xIonDense[nelem][ion] * 
					t_ADfA::Inst().ph1(Heavy.nsShells[nelem][i-1]-1,nelec,nelem,0)/EVRYD* 0.9998787*kadvec;
			}
		}
	}
	/* convert to ergs/cm^3 */
	phycon.EnergyIonization *= EN1RYD;
#endif

	/** \todo	2	this is the total binding energy of the molecules, and 
	 * is negative, the energy need to get back to free atoms 
	 * never set and only appears in print statements */
	phycon.EnergyBinding = 0.;

	/* find number of molecules of the heavy elements in the gas phase. 
	 * Atoms and ions were counted above.  Do not include ices, solids
	 * on grains */
	den_comole = 0.;
	for( long i=0; i < mole.num_comole_calc; i++ )
	{
		/* term COmole[i]->lgGas_Phase is 0 or 1 depending on whether molecule 
		 * is on grain or in gas phase 
		 * n_nuclei is number of nuclei in molecule, this tests whether
		 * actually an atom or molecule - atoms were already counted above */
		if(COmole[i]->n_nuclei > 1)
			den_comole += COmole[i]->hevmol * COmole[i]->lgGas_Phase;
	}

	/* number of hydrogen molecules, loop over all H molecular species,
	 * do not include H0, H+ */
	den_hmole = 0.;
	for( long i=0; i<N_H_MOLEC; ++i )
	{
		if( i!=ipMH && i!=ipMHp )
			den_hmole += hmi.Hmolec[i];
	}

	/* total number of atoms, ions, and molecules in gas phase per unit vol */
	dense.xNucleiTotal = den_hmole + den_comole + DenseAtomsIons;
	if( dense.xNucleiTotal > BIGFLOAT )
	{
		fprintf(ioQQQ, "PROBLEM DISASTER pressure_total has found "
			"dense.xNucleiTotal with an insane density.\n");
		fprintf(ioQQQ, "The density was %.2e\n",
			dense.xNucleiTotal);
		TotalInsanity();
	}
	ASSERT( dense.xNucleiTotal > 0. );

	/* particle density that enters into the pressure includes electrons
	 * cm-3 */
	dense.pden = (realnum)(dense.eden + dense.xNucleiTotal);

	/* the current gas pressure */
	pressure.PresGasCurr = dense.pden*phycon.te*BOLTZMANN;
	/*fprintf(ioQQQ,"DEBUG gassss %.2f %.4e %.4e %.4e \n", 
		fnzone, pressure.PresGasCurr , dense.pden , pressure.PresInteg );*/

	/* dense.wmole is molecular weight, AtomicMassUnit per particle */
	dense.wmole = 0.;
	for( long i=0; i < LIMELM; i++ )
	{
		dense.wmole += dense.gas_phase[i]*(realnum)dense.AtomicWeight[i];
	}

	/* dense.wmole is now molecular weight, AtomicMassUnit per particle */
	dense.wmole /= dense.pden;
	ASSERT( dense.wmole > 0. && dense.pden > 0. );

	/* xMassDensity is density in grams per cc */
	/** \todo	2	- should this include mass in grains? */
	/** \todo	2	- should this include mass in grain mantle ice deposits? */
	dense.xMassDensity = (realnum)(dense.wmole*ATOMIC_MASS_UNIT*dense.pden);

	/*>>chng 04 may 25, introduce following test on xMassDensity0 
	 * WJH 21 may 04, this is the mass density that corresponds to the hden 
	 * specified in the init file. It is used to calculate the mass flux in the 
	 * dynamics models. It may not necessarily be the same as struc.DenMass[0],
	 * since the pressure solver may have jumped to a different density at the
	 * illuminated face from that specified.*/   
	if( dense.xMassDensity0 < 0.0 )
		 dense.xMassDensity0 = dense.xMassDensity;

	/* >>chng 01 nov 02, move to here from dynamics routines */
	/* >>chng 02 jun 18 (rjrw), fix to use local values */
	/* WJH 21 may 04, now recalculate wind v for the first zone too.
	 * This is necessary when we are forcing the sub or supersonic branch */
	if(!(wind.lgBallistic() || wind.lgStatic()))
	{
		wind.windv = DynaFlux(radius.depth)/(dense.xMassDensity);
	}

	/* this is the current ram pressure, at this location */
	pressure.PresRamCurr = dense.xMassDensity*POW2( wind.windv );

	/* this is the current turbulent pressure, not now added to equation of state 
	 * >>chng 06 mar 29, add Heiles_Troland_F / 6. term*/
	pressure.PresTurbCurr = dense.xMassDensity*POW2( DoppVel.TurbVel ) *
		DoppVel.Heiles_Troland_F / 6.;

	/** \todo	0	add this press term due to cosmic rays - hextra.cr_energydensity */

	/* radiative acceleration, lines and continuum */
	if( lgMustEvaluate )
	{
		/* continuous radiative acceleration */
		double rforce = 0.;
		double relec = 0.;
		for( long i=0; i < rfield.nflux; i++ )
		{
			rforce += (rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i]+ rfield.ConInterOut[i])*
				rfield.anu[i]*(opac.opacity_abs[i] + opac.opacity_sct[i]);

			/* electron scattering acceleration - used only to derive force multiplier */
			relec += 
				(rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i]+ 
				rfield.ConInterOut[i]) * 
				opac.OpacStack[i-1+opac.iopcom]*dense.eden*rfield.anu[i];
		}
		/* total continuum radiative acceleration */
		wind.AccelCont = (realnum)(rforce*EN1RYD/SPEEDLIGHT/dense.xMassDensity);

		wind.AccelElectron = (realnum)(relec*EN1RYD/SPEEDLIGHT/dense.xMassDensity);

		/* line acceleration; xMassDensity is gm per cc */
		wind.AccelLine = (realnum)(RT_line_driving()/SPEEDLIGHT/dense.xMassDensity);

		/* total acceleration cm s^-2 */
		wind.AccelTotalOutward = wind.AccelCont + wind.AccelLine;

		/* find wind.fmul - the force multiplier; wind.AccelElectron can be zero for very low 
		 * densities - fmul is of order unity - wind.AccelLine and wind.AccelCont 
		 * are both floats to will underflow long before wind.AccelElectron will - fmul is only used
		 * in output, not in any physics */
		if( wind.AccelElectron > SMALLFLOAT )
			wind.fmul = (realnum)( (wind.AccelLine + wind.AccelCont) / wind.AccelElectron);
		else
			wind.fmul = 0.;

		double reff = radius.Radius-radius.drad/2.;
		/* inward acceleration of gravity cm s^-2 */
		wind.AccelGravity = (realnum)(
			/* wind.comass - mass of star in solar masses */
			GRAV_CONST*wind.comass*SOLAR_MASS/POW2(reff)*
			/* wind.DiskRadius normally zero, set with disk option on wind command */
			(1.-wind.DiskRadius/reff) );
#		if 0
		if( fudge(-1) )
			fprintf(ioQQQ,"DEBUG pressure_total updates AccelTotalOutward to %.2e grav %.2e\n",
			wind.AccelTotalOutward , wind.AccelGravity );
#		endif
	}

	/* must always evaluate H La width since used elsewhere */
	if( Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->PopOpc > 0. )
	{
		hydro.HLineWidth = (realnum)RT_LineWidth(&Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s],GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]));
#		if 0
		fprintf(ioQQQ,"DEBUG widLya\t%li\t%.2f\t%.3e",
			iteration,
			fnzone,
			hydro.HLineWidth);
		hydro.HLineWidth = (realnum)(
			RT_LyaWidth(
			Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].TauIn,
			Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].TauTot,
			4.72e-2/phycon.sqrte,
			GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]) ) );
		fprintf(ioQQQ,"\t%.3e\n",
			hydro.HLineWidth);
#		endif
	}
	else
		hydro.HLineWidth = 0.;


	/* find max rad pressure even if capped
	 * lgLineRadPresOn is turned off by NO RADIATION PRESSURE command */
	if( trace.lgTrace )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent nzone %li iteration %li lgMustEvaluate:%c "
			"lgMustEvaluateTrivial:%c " 
			"pressure.lgLineRadPresOn:%c " 
			"rfield.lgDoLineTrans:%c \n", 
			nzone , iteration , TorF(lgMustEvaluate) , TorF(lgMustEvaluateTrivial),
			TorF(pressure.lgLineRadPresOn), TorF(rfield.lgDoLineTrans) );
	}

	if( lgMustEvaluate && pressure.lgLineRadPresOn && rfield.lgDoLineTrans )
	{
		/* RadPres is pressure due to lines, lgPres_radiation_ON turns off or on */
		pressure.pres_radiation_lines_curr = 0.;
		/* used to remember largest radiation pressure source */
		pmx = 0.;
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			if( lgMustEvaluateTrivial || Piso_seq[ipISO] > TrivialLineRadiationPressure )
			{
				Piso_seq[ipISO] = 0.;
				for( long nelem=ipISO; nelem < LIMELM; nelem++ )
				{
					/* does this ion stage exist? */
					if( dense.IonHigh[nelem] >= nelem + 1 - ipISO  )
					{
						/* do not include highest levels since maser can occur due to topoff,
						 * and pops are set to small number in this case */
						for( long ipHi=1; ipHi <iso.numLevels_local[ipISO][nelem] - iso.nCollapsed_local[ipISO][nelem]; ipHi++ )
						{
							for( long ipLo=0; ipLo < ipHi; ipLo++ )
							{
								if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont <= 0 )
									continue;

								ASSERT( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul > iso.SmallA );

								if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc > SMALLFLOAT &&
									/* test that have not overrun optical depth scale */
									( (Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauTot - Transitions[ipISO][nelem][ipHi][ipLo].Emis->TauIn) > SMALLFLOAT ) &&
									Transitions[ipISO][nelem][ipHi][ipLo].Emis->Pesc > FLT_EPSILON*100. )
								{
									RadPres1 =  PressureRadiationLine( &Transitions[ipISO][nelem][ipHi][ipLo], GetDopplerWidth(dense.AtomicWeight[nelem]) );
									
									if( RadPres1 > pmx )
									{
										ipLineTypePradMax = 2;
										pmx = RadPres1;
										ip4 = ipISO;
										ip3 = nelem;
										ip2 = ipHi;
										ipLinePradMax = ipLo;
									}
									Piso_seq[ipISO] += RadPres1;
									{
										/* option to print particulars of some line when called */
										enum {DEBUG_LOC=false};
										if( DEBUG_LOC && ipISO==ipH_LIKE && ipLo==3 && ipHi==5 && nzone > 260 )
										{
											fprintf(ioQQQ,
												"DEBUG prad1 \t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t\n",
												fnzone,
												RadPres1,
												Transitions[ipISO][nelem][ipHi][ipLo].Emis->PopOpc,
												StatesElemNEW[nelem][nelem-ipISO][ipLo].Pop,
												StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop,
												Transitions[ipISO][nelem][ipHi][ipLo].Emis->Pesc);
										}
									}
								}
							}
						}
					}
				}
				ASSERT( Piso_seq[ipISO] >= 0. );
			}
			pressure.pres_radiation_lines_curr += Piso_seq[ipISO];
		}
		{
			/* option to print particulars of some line when called */
			enum {DEBUG_LOC=false};
#			if 0
			if( DEBUG_LOC /*&& iteration > 1*/ && nzone > 260 )
			{
				fprintf(ioQQQ,
					"DEBUG prad2 \t%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
					nzone,
					pmx,
					Transitions[ipISO][ip3][ipLinePradMax][ip2].Emis->PopOpc,
					StatesElemNEW[ip3][ip3-ipISO][0].Pop,
					StatesElemNEW[ip3][ip3-ipISO][2].Pop,
					StatesElemNEW[ip3][ip3-ipISO][3].Pop,
					StatesElemNEW[ip3][ip3-ipISO][4].Pop,
					StatesElemNEW[ip3][ip3-ipISO][5].Pop,
					StatesElemNEW[ip3][ip3-ipISO][6].Pop);
			}
#			endif
			if( DEBUG_LOC /*&& iteration > 1 && nzone > 150 */)
			{
				fprintf(ioQQQ,
					"DEBUG prad3\t%.2e\t%li\t%li\t%li\t%li\t%.2e\t%.2e\t%.2e\n",
					pmx,
					ip4,
					ip3,
					ip2,
					ipLinePradMax,
					Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->PopOpc,
					StatesElemNEW[ip3][ip3-ip4][ip2].Pop,
					1.-Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->Pesc );
			}
		}

		if( lgMustEvaluateTrivial || PLevel1 > TrivialLineRadiationPressure )
		{
			PLevel1 = 0.;
			/* line radiation pressure from large set of level 1 lines */
			for( long i=1; i <= nLevel1; i++ )
			{
				if( TauLines[i].Hi->Pop > 1e-30 )
				{
					RadPres1 = PressureRadiationLine( &TauLines[i], GetDopplerWidth(dense.AtomicWeight[TauLines[i].Hi->nelem-1]) );

					if( RadPres1 > pmx )
					{
						ipLineTypePradMax = 3;
						pmx = RadPres1;
						ipLinePradMax = i;
					}
					PLevel1 += RadPres1;
				}
			}
			ASSERT( PLevel1 >= 0. );
		}
		pressure.pres_radiation_lines_curr += PLevel1;

		if( lgMustEvaluateTrivial || PLevel2 > TrivialLineRadiationPressure )
		{
			/* level 2 lines */
			PLevel2 = 0.;
			for( long i=0; i < nWindLine; i++ )
			{
				if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
				{
					if( TauLine2[i].Hi->Pop > 1e-30 )
					{
						RadPres1 = PressureRadiationLine( &TauLine2[i], GetDopplerWidth(dense.AtomicWeight[TauLine2[i].Hi->nelem-1]) );

						PLevel2 += RadPres1;
						if( RadPres1 > pmx )
						{
							ipLineTypePradMax = 4;
							pmx = RadPres1;
							ipLinePradMax = i;
						}
					}
				}
			}
			ASSERT( PLevel2 >= 0. );
		}
		pressure.pres_radiation_lines_curr += PLevel2;

		/* fine structure lines */
		if( lgMustEvaluateTrivial || PHfs > TrivialLineRadiationPressure )
		{
			PHfs = 0.;
			for( long i=0; i < nHFLines; i++ )
			{
				if( HFLines[i].Hi->Pop > 1e-30 )
				{
					RadPres1 = PressureRadiationLine( &HFLines[i], GetDopplerWidth(dense.AtomicWeight[HFLines[i].Hi->nelem-1]) );

					PHfs += RadPres1;
					if( RadPres1 > pmx )
					{
						ipLineTypePradMax = 8;
						pmx = RadPres1;
						ipLinePradMax = i;
					}
				}
			}
			ASSERT( PHfs >= 0. );
		}
		pressure.pres_radiation_lines_curr += PHfs;

		/* radiation pressure due to H2 lines */
		if( lgMustEvaluateTrivial || P_H2 > TrivialLineRadiationPressure )
		{
			P_H2 = H2_RadPress();
			/* flag to remember H2 radiation pressure */
			if( P_H2 > pmx )
			{
				pmx = P_H2;
				ipLineTypePradMax = 9;
			}
			ASSERT( P_H2 >= 0. );
		}
		pressure.pres_radiation_lines_curr += P_H2;

		/* radiation pressure due to FeII lines and large atom */
		if( lgMustEvaluateTrivial || P_FeII > TrivialLineRadiationPressure )
		{
			P_FeII = FeIIRadPress();
			if( P_FeII > pmx )
			{
				pmx = P_FeII;
				ipLineTypePradMax = 7;
			}
			ASSERT( P_FeII >= 0. );
		}
		pressure.pres_radiation_lines_curr += P_FeII;

		/* do lines from third-party databases */
		if( lgMustEvaluateTrivial || P_dBase > TrivialLineRadiationPressure )
		{
			P_dBase = 0.;
			for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
			{
				if( Species[ipSpecies].lgActive )
				{
					for( long ipHi=1; ipHi<Species[ipSpecies].numLevels_local; ipHi++ )
					{
						transition *t = &dBaseTrans[ipSpecies][ipHi][0];
						for(long ipLo = 0; ipLo < ipHi; ipLo++ )
						{
							if( t->ipCont > 0 && t->Hi->Pop > 1e-30 )								{
								RadPres1 =  PressureRadiationLine( t, GetDopplerWidth( t->Hi->sp->fmolweight ) );

								if( RadPres1 > pmx )
								{
									ipLineTypePradMax = 10;
									pmx = RadPres1;
									ip3 = ipSpecies;
									ip2 = ipHi;
									ipLinePradMax = ipLo;
								}
								P_dBase += RadPres1;
							}
							++t;
						}
					}
				}
			}
			ASSERT( P_dBase >= 0. );
		}
		pressure.pres_radiation_lines_curr += P_dBase;

	}
	else if( !pressure.lgLineRadPresOn || !rfield.lgDoLineTrans )
	{
		/* case where radiation pressure is not turned on */
		ipLinePradMax = -1;
		ipLineTypePradMax = 0;
	}

	fixit();  // all other line stacks need to be included here.
	// can we just sweep over line stack?  Is that ready yet?

	/* the ratio of radiation to total (gas + continuum + etc) pressure */
	pressure.pbeta = (realnum)(pressure.pres_radiation_lines_curr/SDIV(pressure.PresTotlCurr));

	/* this would be a major logical error */
	if( pressure.pres_radiation_lines_curr < 0. )
	{
		fprintf(ioQQQ,
			"PresTotCurrent: negative pressure, constituents follow %e %e %e %e %e \n",
		Piso_seq[ipH_LIKE],
		Piso_seq[ipHE_LIKE],
		PLevel1,
		PLevel2,
		PCO);
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* following test will never succeed, here to trick lint, ipLineTypePradMax is only used
	 * when needed for debug */
	if( trace.lgTrace && ipLineTypePradMax <0 )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent, pressure.pbeta = %e, ipLineTypePradMax%li ipLinePradMax=%li \n", 
			pressure.pbeta,ipLineTypePradMax,ipLinePradMax );
	}

	/* this is the total internal energy of the gas, the amount of energy needed
	 * to produce the current level populations, relative to ground,
	 * only used for energy terms in advection */
	phycon.EnergyExcitation = 0.;
#if 0
	fixit(); /* the quantities phycon.EnergyExcitation, phycon.EnergyIonization, phycon.EnergyBinding
		  * are not used anywhere, except in print statements... */
	broken(); /* the code below contains serious bugs! It is supposed to implement loops
		   * over quantum states in order to evaluate the potential energy stored in
		   * excited states of all atoms, ions, and molecules. The code below however
		   * often implements loops over all combinations of upper and lower levels!
		   * This code needs to be rewritten once quantumstates are fully implemented. */
	ipLo = 0;
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.IonHigh[nelem] == nelem + 1 - ipISO )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( long ipHi=1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
				{
					phycon.EnergyExcitation += 
						StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop * 
						Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg*
						/* last term is option to turn off energy term, no advection hlike, he-like */
						dynamics.lgISO[ipISO];
				}
			}
		}
	}

	if( dynamics.lgISO[ipH_LIKE] )
		/* internal energy of H2 */
		phycon.EnergyExcitation += H2_InterEnergy();

	/* this is option to turn off energy effects of advection, no advection metals */
	if( dynamics.lgMETALS )
	{
		for( long i=1; i <= nLevel1; i++ )
		{
			phycon.EnergyExcitation += 
				TauLines[i].Hi->Pop* TauLines[i].EnergyErg;
		}
		for( long i=0; i < nWindLine; i++ )
		{
			if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
			{
				phycon.EnergyExcitation += 
					TauLine2[i].Hi->Pop* TauLine2[i].EnergyErg;
			}
		}
		for( long i=0; i < nHFLines; i++ )
		{
			phycon.EnergyExcitation += 
				HFLines[i].Hi->Pop* HFLines[i].EnergyErg;
		}

		/* internal energy of large FeII atom */
		phycon.EnergyExcitation += FeII_InterEnergy();
	}
#endif

	/* ==================================================================
	 * end internal energy of atoms and molecules */

	/* evaluate some parameters to do with magnetic field */
	Magnetic_evaluate();

	/*lint -e644 Piso_seq not init */
	if( trace.lgTrace && (pressure.PresTotlCurr > SMALLFLOAT) )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent zn %.2f Ptot:%.5e Pgas:%.3e Prad:%.3e Pmag:%.3e Pram:%.3e "
			"gas parts; H:%.2e He:%.2e L1:%.2e L2:%.2e CO:%.2e fs%.2e h2%.2e turb%.2e\n",
			fnzone,
			pressure.PresTotlCurr, 
			pressure.PresGasCurr/pressure.PresTotlCurr,
			pressure.pres_radiation_lines_curr*pressure.lgPres_radiation_ON/pressure.PresTotlCurr,
			magnetic.pressure/pressure.PresTotlCurr,
			pressure.PresRamCurr*pressure.lgPres_ram_ON/pressure.PresTotlCurr,
			Piso_seq[ipH_LIKE]/pressure.PresTotlCurr,
			Piso_seq[ipHE_LIKE]/pressure.PresTotlCurr,
			PLevel1/pressure.PresTotlCurr,
			PLevel2/pressure.PresTotlCurr,
			PCO/pressure.PresTotlCurr,
			PHfs/pressure.PresTotlCurr,
			P_H2/pressure.PresTotlCurr,
			pressure.PresTurbCurr*DoppVel.lgTurb_pressure/pressure.PresTotlCurr);
		/*lint +e644 Piso_seq not initialized */
	}

	/* Conserved quantity in steady-state energy equation for dynamics:
	 * thermal energy, since recombination is treated as cooling
	 * (i.e. is loss of electron k.e. to emitted photon), so don't
	 * include
	 * ...phycon.EnergyExcitation + phycon.EnergyIonization + phycon.EnergyBinding
	 * */

	/* constant density case is also hypersonic case */
	if( dynamics.lgTimeDependentStatic )
	{
		/* this branch is time dependent single-zone */
		/* \todo	1	this has 3/2 on the PresGasCurr while true dynamics case below
		 * has 5/2 - so this is not really the enthalpy density - better
		 * would be to always use this term and add the extra PresGasCurr
		 * when the enthalpy is actually needed */
		phycon.EnthalpyDensity =  
			0.5*POW2(wind.windv)*dense.xMassDensity +	/* KE */
			3./2.*pressure.PresGasCurr +				/* thermal plus PdV work */
			magnetic.EnthalpyDensity;					/* magnetic terms */
			//pressure.RhoGravity * (radius.Radius/(radius.Radius+radius.drad)); /* gravity */
	}
	else
	{
		/* this branch either advective or constant pressure */
		/*fprintf(ioQQQ,"DEBUG enthalpy HIT2\n");*/
		/* this is usual dynamics case, or time-varying case where pressure
		 * is kept constant and PdV work is added to the cell */
		phycon.EnthalpyDensity =  
			0.5*POW2(wind.windv)*dense.xMassDensity +	/* KE */
			5./2.*pressure.PresGasCurr +				/* thermal plus PdV work */
			magnetic.EnthalpyDensity;					/* magnetic terms */
			//pressure.RhoGravity * (radius.Radius/(radius.Radius+radius.drad)); /* gravity */
	}

	/* stop model from running away on first iteration, when line optical
	 * depths are not defined correctly anyway.
	 * if( iter.le.itermx .and. RadPres.ge.GasPres ) then
	 * >>chng 97 jul 23, only cap radiation pressure on first iteration */
	if( iteration <= 1 && pressure.pres_radiation_lines_curr >= pressure.PresGasCurr )
	{
		/* stop RadPres from exceeding the gas pressure on first iteration */
		pressure.pres_radiation_lines_curr = 
			MIN2(pressure.pres_radiation_lines_curr,pressure.PresGasCurr);
		pressure.lgPradCap = true;
	}

	/* remember the globally most important line, in the entire model 
	 * test on nzone so we only do test after solution is stable */
	if( pressure.pbeta > pressure.RadBetaMax && nzone )
	{
		pressure.RadBetaMax = pressure.pbeta;
		pressure.ipPradMax_line = ipLinePradMax;
		pressure.ipPradMax_nzone = nzone;
		if( ipLineTypePradMax == 2 )
		{
			/* hydrogenic */
			strcpy( pressure.chLineRadPres, "ISO   ");
			ASSERT( ip4 < NISO && ip3<LIMELM );
			ASSERT( ipLinePradMax>=0 && ip2>=0 && ip3>=0 && ip4>=0 );
			strcat( pressure.chLineRadPres, chLineLbl(&Transitions[ip4][ip3][ip2][ipLinePradMax]) );
			{
				/* option to print particulars of some line when called */
				enum {DEBUG_LOC=false};
				/*lint -e644 Piso_seq not initialized */
				/* trace serious constituents in radiation pressure */
				if( DEBUG_LOC  )
				{
					fprintf(ioQQQ,"DEBUG iso prad\t%li\t%li\tISO,nelem=\t%li\t%li\tnu, no=\t%li\t%li\t%.4e\t%.4e\t%e\t%e\t%e\n",
					iteration, nzone, 
					ip4,ip3,ip2,ipLinePradMax,
					Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->TauIn,
					Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->TauTot,
					Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->Pesc,
					Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->Pelec_esc,
					Transitions[ip4][ip3][ip2][ipLinePradMax].Emis->Pdest);
					if( ip2==5 && ipLinePradMax==4 )
					{
						double width;
						fprintf(ioQQQ,"hit it\n");
						width = RT_LineWidth(&Transitions[ip4][ip3][ip2][ipLinePradMax],GetDopplerWidth(dense.AtomicWeight[ip3]));
						fprintf(ioQQQ,"DEBUG width %e\n", width);
					}
				}
			}
		}
		else if( ipLineTypePradMax == 3 )
		{
			/* level 1 lines */
			ASSERT( ipLinePradMax>=0  );
			strcpy( pressure.chLineRadPres, "Level1 ");
			strcat( pressure.chLineRadPres, chLineLbl(&TauLines[ipLinePradMax]) );
		}
		else if( ipLineTypePradMax == 4 )
		{
			/* level 2 lines */
			ASSERT( ipLinePradMax>=0  );
			strcpy( pressure.chLineRadPres, "Level2 ");
			strcat( pressure.chLineRadPres, chLineLbl(&TauLine2[ipLinePradMax]) );
		}
		else if( ipLineTypePradMax == 5 )
		{
			cdEXIT( EXIT_FAILURE );
		}
		else if( ipLineTypePradMax == 6 )
		{
			cdEXIT( EXIT_FAILURE );
		}
		else if( ipLineTypePradMax == 7 )
		{
			/* FeII lines */
			strcpy( pressure.chLineRadPres, "Fe II ");
		}
		else if( ipLineTypePradMax == 8 )
		{
			/* hyperfine struct lines */
			strcpy( pressure.chLineRadPres, "hyperf ");
			ASSERT( ipLinePradMax>=0  );
			strcat( pressure.chLineRadPres, chLineLbl(&HFLines[ipLinePradMax]) );
		}
		else if( ipLineTypePradMax == 9 )
		{
			/* large H2 lines */
			strcpy( pressure.chLineRadPres, "large H2 ");
		}
		else if( ipLineTypePradMax == 10 )
		{
			/* database lines */
			strcpy( pressure.chLineRadPres, "dBaseLin " );
			strcat( pressure.chLineRadPres, dBaseTrans[ip3][ip2][ipLinePradMax].Hi->sp->chLabel );
		}
		else
		{
			fprintf(ioQQQ," PresTotCurrent ipLineTypePradMax set to %li, this is impossible.\n", ipLineTypePradMax);
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( trace.lgTrace && pressure.pbeta > 0.5 )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent Pbeta:%.2f due to %s\n",
			pressure.pbeta ,
			pressure.chLineRadPres
			);
	}

	/* >>chng 02 jun 27 - add in magnetic pressure */
	/* start to bring total pressure together */
	TotalPressure_v = pressure.PresGasCurr;

	/* these flags are set false by default since constant density is default,
	 * set true for constant pressure or dynamics */
	TotalPressure_v += pressure.PresRamCurr * pressure.lgPres_ram_ON;

	/* magnetic pressure, evaluated in magnetic.c - this can be negative for an ordered field! 
	 * option is on by default, turned off with constant density, or constant gas pressure, cases */
	/** \todo	0	code has variable magnetic energydensity and pressure, which are equal,
	 * as they must be - del one or the other */
	TotalPressure_v += magnetic.pressure * pressure.lgPres_magnetic_ON;

	/* turbulent pressure
	 * >>chng 06 mar 24, include this by default */
	TotalPressure_v += pressure.PresTurbCurr * DoppVel.lgTurb_pressure;

	/* radiation pressure only included in total eqn of state when this flag is
	 * true, set with constant pressure command */
	/* option to add in internal line radiation pressure */
	TotalPressure_v += pressure.pres_radiation_lines_curr * pressure.lgPres_radiation_ON;

	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC && pressure.PresTotlCurr > SMALLFLOAT /*&& iteration > 1*/ )
		{
			fprintf(ioQQQ,"pressureee%li\t%.4e\t%.4e\t%.4e\t%.3f\t%.3f\t%.3f\n", 
				nzone,
				pressure.PresTotlCorrect,
			    pressure.PresTotlCurr, 
				TotalPressure_v ,
				pressure.PresGasCurr/pressure.PresTotlCurr,
				pressure.pres_radiation_lines_curr/pressure.PresTotlCurr ,
				pressure.PresRamCurr/pressure.PresTotlCurr
				);
		}
	}

	if( TotalPressure_v< 0. )
	{
		ASSERT( magnetic.pressure < 0. );

		/* negative pressure due to ordered field overwhelms total pressure - punt */
		fprintf(ioQQQ," The negative pressure due to ordered magnetic field overwhelms the total outward pressure - please reconsider the geometry & field.\n");
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( TotalPressure_v > 0. );

	/* remember highest pressure anywhere */
	pressure.PresMax = MAX2(pressure.PresMax,(realnum)TotalPressure_v);

	/* this is what we came for - set the current pressure */
	pressure.PresTotlCurr = TotalPressure_v;

	return;
}
