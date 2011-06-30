/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_diffuse evaluate local diffuse emission for this zone,
 * fill in ConEmitLocal[depth][energy] with diffuse emission,
 * called by Cloudy, this routine adds energy to the outward beam 
 * OTS rates for this zone were set in RT_OTS - not here */
#include "cddefines.h"
#include "physconst.h"
#include "taulines.h"
#include "grains.h"
#include "grainvar.h"
#include "iso.h"
#include "dense.h"
#include "opacity.h"
#include "trace.h"
#include "coolheavy.h"
#include "rfield.h"
#include "phycon.h"
#include "hmi.h"
#include "radius.h"
#include "atmdat.h"
#include "heavy.h"
#include "atomfeii.h"
#include "lines_service.h"
#include "h2.h"
#include "ipoint.h"
#include "rt.h"

#define TwoS    (1+ipISO)

#if defined (__ICC) && defined(__ia64) && __INTEL_COMPILER < 910
#pragma optimization_level 0
#endif
void RT_diffuse(void)
{
	/* set this true to print a table of 2 photon emission coefficients
	 * comparable to Brown & Mathews (1971) */
	static long lgPrt2PhtEmisCoef = false;

	/* arrays used in this routine
	 * rfield.ConEmitLocal[depth][energy] local emission per unit vol 
	 * rfield.ConOTS_local_photons - local ots two-photon rates
	 * rfield.rfield.DiffuseEscape is the spectrum of diffuse emission that escapes this zone,
	 * at end of this routine part is thrown into the outward beam
	 * by adding to rfield.ConInterOut
	 * units are photons s-1 cm-3 
	 * one-time init done on first call */

	/* rfield.rfield.DiffuseEscape and rfield.ConEmitLocal are same except that 
	 * rfield.ConEmitLocal is local emission, would be source function if div by opac
	 * rfield.rfield.DiffuseEscape is part that escapes so has RT built into it 
	 * rfield.rfield.DiffuseEscape is used to define rfield.ConInterOut below as per this statement
	 * rfield.ConInterOut[nu] += rfield.DiffuseEscape[nu]*(realnum)radius.dVolOutwrd;
	 */
	/* \todo	0	define only rfield.ConEmitLocal as it is now done, 
	 * do not define rfield.DiffuseEscape at all
	 * at bottom of this routine use inward and outward optical depths to define
	 * local and escaping parts 
	 * this routine only defines 
	 * rfield.ConInterOut - set to rfield.DiffuseEscape times vol element 
	 * rfield.ConOTS_local_photons (two photon only)
	 * so this is only var that
	 * needs to be set 
	 */

	long int i=-100000, 
	  ip=-100000, 
	  ipISO=-100000,
	  ipHi=-100000, 
	  ipLo=-100000, 
	  ipla=-100000,
	  nelem=-100000,
	  limit=-100000, 
	  n=-100000,
	  nu=-10000;

	double EnergyLimit;

	double EdenAbund, 
	  difflya, 
	  arg, 
	  fac, 
	  factor, 
	  gamma, 
	  gion, 
	  gn, 
	  photon, 
	  sum,
	  Sum1level,
	  SumCaseB;

	realnum two_photon_emiss;

	DEBUG_ENTRY( "RT_diffuse()" );

	/* many arrays were malloced to nupper, and we will add unit flux to [nflux] -
	 8 this must be true to work */
	ASSERT(rfield.nflux < rfield.nupper );

	/* this routine evaluates the local diffuse fields
	 * it fills in all of the following vectors */
	memset(rfield.DiffuseEscape               , 0 , (unsigned)rfield.nupper*sizeof(realnum) );
	memset(rfield.ConEmitLocal[nzone]  , 0 , (unsigned)rfield.nupper*sizeof(realnum) );
	memset(rfield.ConOTS_local_photons , 0 , (unsigned)rfield.nupper*sizeof(realnum) );
	memset(rfield.TotDiff2Pht          , 0 , (unsigned)rfield.nupper*sizeof(realnum) );
	memset(rfield.DiffuseLineEmission  , 0 , (unsigned)rfield.nupper*sizeof(realnum) );

	/* must abort after setting all of above to zero because some may be
	 * used in various ways before abort is complete */
	if( lgAbort )
	{
		/* quit if we are aborting */
		return;
	}

	/* loop over iso-sequences of all elements 
	 * to add all recombination continua and lines*/
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* >>chng 01 sep 23, rewrote for iso sequences */
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* this will be the sum of recombinations to all excited levels */
			SumCaseB = 0.;

			/* the product of the densities of the parent ion and electrons */
			EdenAbund = dense.eden*dense.xIonDense[nelem][nelem+1-ipISO];

			/* recombination continua for all iso seq - 
			 * if this stage of ionization exists */
			if( dense.IonHigh[nelem] >= nelem+1-ipISO  )
			{
				/* loop over all levels to include recombination diffuse continua,
				 * pick highest energy continuum point that opacities extend to 
				 * for ground continuum, go up to highest defined Boltzmann factor,
				 * at bottom of loop will be reset to ground state photo edge */
				ipHi = rfield.nflux;
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( n=0; n < iso.numLevels_local[ipISO][nelem]; n++ )
				{
					Sum1level = 0.;
					iso.RadRecCon[ipISO][nelem][n] = 0.;
					/* >>chng 02 nov 20 - pull the plug on old he treatment */
					/* the number is (2 pi me k/h^2) ^ -3/2 * 8 pi/c^2 / ge - it includes
					 * the stat weight of the free electron in the demominator */
					/** \todo	2	This doesn't really seem to be the expression above!!!	*/
					/*gamma = 2.0618684e11*StatesElemNEW[nelem][nelem-ipISO][n].g/iso.stat_ion[ipISO]/phycon.te/phycon.sqrte;*/
					gamma = 0.5*MILNE_CONST*StatesElemNEW[nelem][nelem-ipISO][n].g/iso.stat_ion[ipISO]/phycon.te/phycon.sqrte;

					/* loop over all recombination continua 
					 * escaping part of recombinations are added to rfield.ConEmitLocal 
					 * added to ConInterOut at end of routine */
					for( nu=iso.ipIsoLevNIonCon[ipISO][nelem][n]-1; nu < ipHi; nu++ )
					{
						/* dwid used to adjust where within WIDFLX exp is evaluated -
						 * weighted to lower energy due to exp(-energy/T) */
						double dwid = 0.2;

						/* this is the term in the negative exponential Boltzmann factor
						 * for continuum emission */
						arg = (rfield.anu[nu]-iso.xIsoLevNIonRyd[ipISO][nelem][n]+
							rfield.widflx[nu]*dwid)/phycon.te_ryd;
						arg = MAX2(0.,arg);
						/* don't bother evaluating for this or higher energies if 
						 * Boltzmann factor is tiny. */
						if( arg > SEXP_LIMIT ) 
							break;

						/* photon is in photons cm^3 s^-1 per cell */
						photon = gamma*exp(-arg)*rfield.widflx[nu]*
							opac.OpacStack[ nu-iso.ipIsoLevNIonCon[ipISO][nelem][n] + iso.ipOpac[ipISO][nelem][n] ] *
							rfield.anu2[nu];

						Sum1level += photon;

						/* total local diffuse emission units photons cm-3 s-1 cell-1,*/
						rfield.ConEmitLocal[nzone][nu] += (realnum)(photon*EdenAbund);

						/* iso.RadRecomb[ipISO][nelem][n][ipRecEsc] is escape probability
						 * rfield.DiffuseEscape is local emission that escapes this zone */
						rfield.DiffuseEscape[nu] += 
							(realnum)(photon*EdenAbund*iso.RadRecomb[ipISO][nelem][n][ipRecEsc]);

						// total RRC radiative recombination continuum
						iso.RadRecCon[ipISO][nelem][n] += rfield.anu[nu] *
							emergent_line( photon*EdenAbund/2. , photon*EdenAbund/2. , 
							// energy on fortran scale
							nu+1 );
					}

					// convert to erg cm-3 s-1
					iso.RadRecCon[ipISO][nelem][n] *= EN1RYD;
					/* this will be used below to confirm case B sum */
					if( n > 0 )
					{
						/* SumCaseB will be sum to all excited */
						SumCaseB += Sum1level;
					}

					/* on entry to this loop ipHi was set to the upper limit of the code,
					 * and ground state recombination continua for all energies was added.  For
					 * excited states (which are the ones that will be done after
					 * first pass through) only go to ground state threshold,
					 * since that will be so much stronger than excited state recom*/
					ipHi = iso.ipIsoLevNIonCon[ipISO][nelem][0]-1;
				}

				// no RRC emission from levels that do not exist
				for( n=iso.numLevels_local[ipISO][nelem];n<iso.numLevels_max[ipISO][nelem]; n++ )
					iso.RadRecCon[ipISO][nelem][n] = 0.;

				/* this is check on self-consistency */
				iso.CaseBCheck[ipISO][nelem] = MAX2(iso.CaseBCheck[ipISO][nelem],
				  (realnum)(SumCaseB/iso.RadRec_caseB[ipISO][nelem]));

				/* this add line emission from the model atoms,
				 * do not include very highest level since disturbed by topoff */
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( ipLo=0; ipLo < (iso.numLevels_local[ipISO][nelem] - 2); ipLo++ )
				{
					/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
					for( ipHi=ipLo+1; ipHi < iso.numLevels_local[ipISO][nelem] - 1; ipHi++ )
					{
						/* must not include fake transitions (2s-1s, two photon, or truly
						 * forbidden transitions */
						if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont < 1 )
							continue;

						/* number of photons in the line has not been defined up until now,
						 * do so now.  this is redone in lines.  */
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->phots = 
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul*
							StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop*
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->Pesc;
						
						// Would be better to enable checks (and remove argument) --
						// present state is to ensure backwards compatibility with previous
						// unchecked code.
						// First argument is fraction of line not emitted by scattering --
						// would be better to do this on the basis of line physics rather than
						// fiat...
						const bool lgDoChecks = false;
						Transitions[ipISO][nelem][ipHi][ipLo].outline(1.0, lgDoChecks );
					}
				}

				/*Iso treatment of two photon emission.  */
				/* NISO could in the future be increased, but we want this assert to blow
				 * so that it is understood this may not be correct for other iso sequences,
				 * probably should break since will not be present */
				ASSERT( ipISO <= ipHE_LIKE );

				/* upper limit to 2-phot is energy of 2s to ground */
				EnergyLimit = Transitions[ipISO][nelem][TwoS][0].EnergyWN/RYD_INF;

				/* this could fail when pops very low and Pop2Ion is zero */
				ASSERT( iso.ipTwoPhoE[ipISO][nelem]>0 && EnergyLimit>0. );

				sum = 0.;
				/* two photon emission, iso.ipTwoPhoE[ipISO][nelem] is 
				 * continuum array index for Lya energy */
				for( nu=0; nu<iso.ipTwoPhoE[ipISO][nelem]; nu++ )
				{
					/* We do not assert rfield.anu[nu]<=EnergyLimit because the maximum 
					 * index could be set to point to the next highest bin.	*/
					ASSERT( rfield.anu[nu]/EnergyLimit < 1.01 || rfield.anu[nu-1]<EnergyLimit);

					/* iso.As2nu[ipISO][nelem][nu] is transition probability A per bin
					 * So sum is the total transition probability - this sum should
					 * add up to the A two photon */
					sum += iso.As2nu[ipISO][nelem][nu];

					/* flag saying whether to also do he triplets */
					const bool lgDo2TripSToo = false;

					if( ipISO == ipHE_LIKE && lgDo2TripSToo )
					{
						two_photon_emiss = 2.f*iso.As2nu[ipISO][nelem][nu]*
							((realnum)StatesElemNEW[nelem][nelem-ipISO][TwoS].Pop+0.0001f*(realnum)StatesElemNEW[nelem][nelem-ipISO][ipHe2s3S].Pop);
					}
					else
					{
						/* StatesElemNEW[nelem][nelem-ipISO][TwoS].Pop has dimension cm^-3.  The factor of 2 is for two photons per 
						 * transition. Thus two_photon_emiss has dimension photons cm-3 s-1.	*/
						two_photon_emiss = 2.f*(realnum)StatesElemNEW[nelem][nelem-ipISO][TwoS].Pop*iso.As2nu[ipISO][nelem][nu];
					}

					/* >>chng 03 mar 26, only do if induced processes turned on,
					 * otherwise inconsistent with rate solver treatment.	*/
					/* >>chng 02 aug 14, add induced two photon emission */
					/* product of occupation numbers for induced emission */
					if( iso.lgInd2nu_On)
					{
						/* this is the total rate */
						two_photon_emiss *= (1.f + rfield.SummedOcc[nu]) *
							(1.f+rfield.SummedOcc[iso.ipSym2nu[ipISO][nelem][nu]-1]);
					}

					/* information - only used in save output */
					rfield.TotDiff2Pht[nu] += two_photon_emiss;

					/* total local diffuse emission */
					rfield.ConEmitLocal[nzone][nu] += two_photon_emiss;

					/* this is escaping part of two-photon emission, 
					 * as determined from optical depth to illuminated face */
					rfield.DiffuseEscape[nu] += two_photon_emiss*opac.ExpmTau[nu];

					/* save locally destroyed OTS two-photon continuum - this is only
					 * place this is set in this routine */
					rfield.ConOTS_local_photons[nu] += two_photon_emiss*(1.f - opac.ExpmTau[nu]);
				}

				/* a sanity check on the code, see Spitzer and Greenstein (1951), eqn 4.	*/
				/* >>refer	HI	2nu	Spitzer, L., & Greenstein, J., 1951, ApJ, 114, 407 */
				ASSERT( fabs( 1.f - (realnum)sum/Transitions[ipISO][nelem][TwoS][0].Emis->Aul ) < 0.01f );
			}

			/* option to print hydrogen and helium two-photon emission coefficients?	*/
			if( lgPrt2PhtEmisCoef )
			{
				long yTimes20;
				double y, E2nu;

				lgPrt2PhtEmisCoef = false;

				fprintf( ioQQQ, "\ny\tGammaNot(2q)\n");

				for( yTimes20=1; yTimes20<=10; yTimes20++ )
				{
					y = yTimes20/20.;

					fprintf( ioQQQ, "%.3e\t", (double)y );	

					E2nu = Transitions[0][0][1][0].EnergyWN/RYD_INF;
					i = ipoint(y*E2nu);
					fprintf( ioQQQ, "%.3e\t", 
						8./3.*HPLANCK*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][1].Pop/dense.eden/dense.xIonDense[0][1]
						*y*iso.As2nu[0][0][i]*E2nu/rfield.widflx[i] );

					E2nu = Transitions[1][1][2][0].EnergyWN/RYD_INF;
					i = ipoint(y*E2nu);
					fprintf( ioQQQ, "%.3e\n", 
						8./3.*HPLANCK*StatesElemNEW[ipHELIUM][ipHELIUM-ipHE_LIKE][2].Pop/dense.eden/dense.xIonDense[1][1]
						*y*iso.As2nu[1][1][i]*E2nu/rfield.widflx[i] );
				}
				fprintf( ioQQQ, "eden%.3e\n", dense.eden );
			}
		}
	}

	/* add recombination continua for elements heavier than those done with iso seq */
	for( nelem=NISO; nelem < LIMELM; nelem++ )
	{
		/* do not include species with iso-sequence in following */
		/* >>chng 03 sep 09, upper bound was wrong, did not include NISO */
		for( long ion=dense.IonLow[nelem]; ion < nelem-NISO+1; ion++ )
		{
			if( dense.xIonDense[nelem][ion+1] > 0. )
			{
				long int ns, nshell,igRec , igIon,
					iplow , iphi , ipop;

				ip = Heavy.ipHeavy[nelem][ion]-1;
				ASSERT( ip >= 0 );

				/* nflux was reset upward in ConvInitSolution to encompass all
				 * possible line and continuum emission.  this test should not
				 * possibly fail.  It could if the ionization were to increase with depth
				 * although the continuum mesh is designed to deal with this.
				 * This test is important because the nflux cell in ConInterOut 
				 * is used to carry out the unit integration, and if it gets 
				 * clobbered by diffuse emission the code will declare 
				 * insanity in PrtComment */
				if( ip >= rfield.nflux )
					continue;

				/* get shell number, stat weights for this species */
				atmdat_outer_shell( nelem+1 , nelem+1-ion , &nshell, &igRec , &igIon );
				gn = (double)igRec;
				gion = (double)igIon;

				/* shell number */
				ns = Heavy.nsShells[nelem][ion]-1;
				ASSERT( ns == (nshell-1) );

				/* lower and upper energies, and offset for opacity stack */
				iplow = opac.ipElement[nelem][ion][ns][0]-1;
				iphi = opac.ipElement[nelem][ion][ns][1];
				iphi = MIN2( iphi , rfield.nflux );
				ipop = opac.ipElement[nelem][ion][ns][2];

				/* now convert ipop to the offset in the opacity stack from threshold */
				ipop = ipop - iplow;

				gamma = 0.5*MILNE_CONST*gn/gion/phycon.te/phycon.sqrte*dense.eden*dense.xIonDense[nelem][ion+1];

				/* this is ground state continuum from stored opacities */
				Heavy.RadRecCon[nelem][ion] = 0;
				if( rfield.ContBoltz[iplow] > SMALLFLOAT )
				{
					for( nu=iplow; nu < iphi; ++nu )
					{
						photon = gamma*rfield.ContBoltz[nu]/rfield.ContBoltz[iplow]*
							rfield.widflx[nu]*opac.OpacStack[nu+ipop]*rfield.anu2[nu];
						/* add heavy rec to ground in active beam,*/
						/** \todo	2	should use ConEmitLocal for all continua, but not followed
						 * by rfield.DiffuseEscape - put that at the end.  Once continua all
						 * bundled this way, it will be easy to save them as a function
						 * of depth and then do exact rt */
						rfield.ConEmitLocal[nzone][nu] += (realnum)photon;
						rfield.DiffuseEscape[nu] += (realnum)photon*opac.ExpmTau[nu];

						// escaping RRC
						Heavy.RadRecCon[nelem][ion] += rfield.anu[nu] *
						  emergent_line( photon/2. , photon/2. , 
						  // energy on fortran scale
						  nu+1 );
					}
				}
				// units erg cm-3 s-1
				Heavy.RadRecCon[nelem][ion] *= EN1RYD;

				/* now do the recombination Lya */
				ipla = Heavy.ipLyHeavy[nelem][ion]-1;
				ASSERT( ipla >= 0 );
				/* xLyaHeavy is set to a fraction of the total rad rec in ion_recomb, includes eden */
				difflya = Heavy.xLyaHeavy[nelem][ion]*dense.xIonDense[nelem][ion+1];
				rfield.DiffuseLineEmission[ipla] += (realnum)difflya;

				/* >>chng 03 jul 10, here and below, use outlin_noplot */
				rfield.outlin_noplot[ipla] += (realnum)(difflya*radius.dVolOutwrd*opac.tmn[ipla]*opac.ExpmTau[ipla]);

				/* now do the recombination Balmer photons */
				ipla = Heavy.ipBalHeavy[nelem][ion]-1;
				ASSERT( ipla >= 0 );
				/* xLyaHeavy is set to fraction of total rad rec in ion_recomb, includes eden */
				difflya = Heavy.xLyaHeavy[nelem][ion]*dense.xIonDense[nelem][ion+1];
				rfield.outlin_noplot[ipla] += (realnum)(difflya*radius.dVolOutwrd*opac.tmn[ipla]*opac.ExpmTau[ipla]);
			}
		}
	}

	/* free-free free free brems emission for all ions */
	limit = MIN2( rfield.ipMaxBolt , rfield.nflux );
	for( nu=0; nu < limit; nu++ )
	{
		double TotBremsAllIons = 0., BremsThisIon;

		/* First add H- brems.  Reaction is H(1s) + e -> H(1s) + e + hnu.
		 * OpacStack contains the ratio of the H- to H brems cross section.
		 * Multiply H brems by this and the population of H(1s). */
		TotBremsAllIons += rfield.gff[1][nu] * opac.OpacStack[nu-1+opac.iphmra] * StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;

		/* chng 02 may 16, by Ryan...do all brems for all ions in one fell swoop,
		 * using gaunt factors from rfield.gff.	*/
		for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
		{
			/* MAX2 occurs because we want to start at first ion (or above)
			 * and do not want atom */
			for( long ion=MAX2(1,dense.IonLow[nelem]); ion<=dense.IonHigh[nelem]; ++ion )
			{
				/* eff. charge is ion, so first rfield.gff argument must be "ion".	*/
				BremsThisIon = POW2( (realnum)ion )*dense.xIonDense[nelem][ion]*rfield.gff[ion][nu];
				TotBremsAllIons += BremsThisIon;
			}
		}

		/* add molecular ions */
		long ion = 1;
		BremsThisIon = POW2( (double)ion )*rfield.gff[ion][nu]*(hmi.Hmolec[ipMH2p] + hmi.Hmolec[ipMH3p] + hmi.Hmolec[ipMHeHp]);
		TotBremsAllIons += BremsThisIon;

		/** \todo	2	Replace this constant with the appropriate macro, if any */
		/* >>chng 06 apr 05, no free free also turns off emission */
		TotBremsAllIons *= dense.eden*1.032e-11*rfield.widflx[nu]*rfield.ContBoltz[nu]/rfield.anu[nu]/phycon.sqrte *
			CoolHeavy.lgFreeOn;
		ASSERT( TotBremsAllIons >= 0.);

		/* >>chng 01 jul 01, move thick brems back to ConEmitLocal but do not add
		 * to outward beam - ConLocNoInter array removed as result
		 * if problems develop with very dense blr clouds, this may be reason */
		/*rfield.ConLocNoInter[nu] += (realnum)fac;*/
		/*rfield.ConEmitLocal[nzone][nu] += (realnum)TotBremsAllIons;*/

		if( nu >= rfield.ipEnergyBremsThin )
		{
			/* >>chng 05 feb 20, move into this test on brems opacity - should not be
			 * needed since would use expmtau to limit outward beam */
			/* >>chng 01 jul 01, move thick brems back to ConEmitLocal but do not add
			* to outward beam - ConLocNoInter array removed as result
			* if problems develop with very dense BLR clouds, this may be reason */
			/*rfield.ConLocNoInter[nu] += (realnum)fac;*/
			rfield.ConEmitLocal[nzone][nu] += (realnum)TotBremsAllIons;

			/* do not add optically thick part to outward beam since self absorbed
			* >>chng 96 feb 27, put back into outward beam since do not integrate
			* over it anyway. */
			/* >>chng 99 may 28, take back out of beam since DO integrate over it
			* in very dense BLR clouds */
			/* >>chng 01 jul 10, add here, in only one loop, where optically thin */
			rfield.DiffuseEscape[nu] += (realnum)TotBremsAllIons;
		}
	}

	/* grain dust emission */
	/* >>chng 01 nov 22, moved calculation of grain flux to qheat.c, PvH */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		/* this calculates diffuse emission from grains,
		 * and stores the result in gv.GrainEmission */
		GrainMakeDiffuse();

		for( nu=0; nu < rfield.nflux; nu++ )
		{
			rfield.ConEmitLocal[nzone][nu] += gv.GrainEmission[nu];
			rfield.DiffuseEscape[nu] += gv.GrainEmission[nu];
		}
	}

	/* hminus emission */
	fac = dense.eden*(double)dense.xIonDense[ipHYDROGEN][0];
	gn = 1.;
	gion = 2.;
	gamma = 0.5*MILNE_CONST*gn/gion/phycon.te/phycon.sqrte;
	/* >>chng 00 dec 15 change limit to -1 of H edge */
	limit = MIN2(iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1,rfield.nflux);

	if( rfield.ContBoltz[hmi.iphmin-1] > 0. )
	{
		for( nu=hmi.iphmin-1; nu < limit; nu++ )
		{
			/* H- flux photons cm-3 s-1 
			 * ContBoltz is ratio of Boltzmann factor for each freq */
			factor = gamma*rfield.ContBoltz[nu]/rfield.ContBoltz[hmi.iphmin-1]*rfield.widflx[nu]*
			  opac.OpacStack[nu-hmi.iphmin+opac.iphmop]*
			  rfield.anu2[nu]*fac;
			rfield.ConEmitLocal[nzone][nu] += (realnum)factor;
			rfield.DiffuseEscape[nu] += (realnum)factor;
		}
	}
	else
	{
		for( nu=hmi.iphmin-1; nu < limit; nu++ )
		{
			arg = MAX2(0.,TE1RYD*(rfield.anu[nu]-0.05544)/phycon.te);
			/* this is the limit sexp normally uses */
			if( arg > SEXP_LIMIT ) 
				break;
			/* H- flux photons cm-3 s-1 
			 * flux is in photons per sec per ryd */
			factor = gamma*exp(-arg)*rfield.widflx[nu]*
				opac.OpacStack[nu-hmi.iphmin+opac.iphmop]*
			  rfield.anu2[nu]*fac;
			rfield.ConEmitLocal[nzone][nu] += (realnum)factor;
			rfield.DiffuseEscape[nu] += (realnum)factor;
		}
	}

	/* outward level 1 line photons, 0 is dummy line */
	for( i=1; i <= nLevel1; i++ )
		TauLines[i].outline_resonance();

	for( ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && iso.lgDielRecom[ipISO] ) 
			{
				for( i=0; i<iso.numLevels_max[ipISO][nelem]; i++ )
				{

					SatelliteLines[ipISO][nelem][i].Emis->phots = 
						SatelliteLines[ipISO][nelem][i].Emis->Aul*
						SatelliteLines[ipISO][nelem][i].Hi->Pop*
						(SatelliteLines[ipISO][nelem][i].Emis->Pesc+
						SatelliteLines[ipISO][nelem][i].Emis->Pelec_esc);

					SatelliteLines[ipISO][nelem][i].outline_resonance();
				}
			}
		}
	}

	/* outward level 2 line photons */
	for( i=0; i < nWindLine; i++ )
	{
		/* must not also do lines that were already done as part
		 * of the isoelectronic sequences */
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			{
				enum {DEBUG_LOC=false};
				if( DEBUG_LOC /*&& nzone > 10*/ && i==4821 )
				{
					/* set up to dump the Fe 9 169A line */
					fprintf(ioQQQ,"DEBUG dump lev2 line %li\n", i );
					DumpLine( &TauLine2[i] );/**/
					fprintf(ioQQQ,"DEBUG dump %.3e %.3e %.3e\n",
						rfield.outlin[0][TauLine2[i].ipCont-1],
						TauLine2[i].Emis->phots*TauLine2[i].Emis->FracInwd*radius.BeamInOut*opac.tmn[i]*TauLine2[i].Emis->ColOvTot,
						TauLine2[i].Emis->phots*(1. - TauLine2[i].Emis->FracInwd)*radius.BeamOutOut* TauLine2[i].Emis->ColOvTot );
				}
			}
			TauLine2[i].outline_resonance();
			/*if( i==2576 ) fprintf(ioQQQ,"DEBUG dump %.3e %.3e \n",
				rfield.outlin[0][TauLine2[i].ipCont-1] , rfield.outlin_noplot[TauLine2[i].ipCont-1]);*/
		}
	}

	/* outward hyperfine structure line photons */
	for( i=0; i < nHFLines; i++ )
	{
		HFLines[i].outline_resonance();
	}

	/* external database lines */
	for( i=0; i < linesAdded2; i++ )
	{
		dBaseLines[i].tran->outline_resonance();
	}

	/* H2 emission */
	H2_RT_diffuse();

	/* do outward parts of FeII lines, if large atom is turned on */
	FeII_RT_Out();
	/** \todo	2	add fegrain to outward beams, but within main formalism by including grains
	 * in all x-ray processes */

	if( trace.lgTrace )
		fprintf( ioQQQ, " RT_diffuse returns.\n" );

	/* >>chng 02 jul 25, zero out all light below plasma freq */
	for( nu=0; nu < rfield.ipPlasma-1; nu++ )
	{
		rfield.flux_beam_const[nu] = 0.;
		rfield.flux_beam_time[nu] = 0.;
		rfield.flux_isotropic[nu] = 0.;
		rfield.flux[0][nu] = 0.;
		rfield.ConEmitLocal[nzone][nu] = 0.;
		rfield.otscon[nu] = 0.;
		rfield.otslin[nu] = 0.;
		rfield.outlin[0][nu] = 0.;
		rfield.outlin_noplot[nu] = 0.;
		rfield.reflin[0][nu] = 0.;
		rfield.TotDiff2Pht[nu] = 0.;
		rfield.ConOTS_local_photons[nu] = 0.;
		rfield.ConInterOut[nu] = 0.;
	}

	/* find occupation number, also assert that no continua are negative */
	for( nu=0; nu < rfield.nflux; nu++ )
	{
		/* >>chng 00 oct 03, add diffuse continua */
		/* local diffuse continua */
		rfield.OccNumbDiffCont[nu] = 
			rfield.ConEmitLocal[nzone][nu]*rfield.convoc[nu];

		/* units are photons cell-1 cm-2 s-1  */
		rfield.ConSourceFcnLocal[nzone][nu] = 
			/* units photons cm-3 s-1 cell-1, */
			(realnum)safe_div( (double)rfield.ConEmitLocal[nzone][nu],
			/* units cm-1 */
			opac.opacity_abs[nu] );

		/* confirm that all are non-negative */
		ASSERT( rfield.flux_beam_const[nu] >= 0.);
		ASSERT( rfield.flux_beam_time[nu] >= 0.);
		ASSERT( rfield.flux_isotropic[nu] >= 0.);
		ASSERT( rfield.flux[0][nu] >= 0.);
		ASSERT( rfield.ConEmitLocal[nzone][nu] >= 0.);
		ASSERT( rfield.otscon[nu] >= 0.);
		ASSERT( rfield.otslin[nu] >= 0.);
		ASSERT( rfield.outlin[0][nu] >= 0.);
		ASSERT( rfield.outlin_noplot[nu] >= 0.);
		ASSERT( rfield.reflin[0][nu] >= 0.);
		ASSERT( rfield.TotDiff2Pht[nu] >= 0.);
		ASSERT( rfield.ConOTS_local_photons[nu] >= 0.);
		ASSERT( rfield.ConInterOut[nu] >= 0.);
	}

	/* option to kill outward lines with no outward lines command*/
	if( rfield.lgKillOutLine )
	{
		for( nu=0; nu < rfield.nflux; nu++ )
		{
			rfield.outlin[0][nu] = 0.;
			rfield.outlin_noplot[nu] = 0.;
		}
	}

	/* option to kill outward continua with no outward continua command*/
	if( rfield.lgKillOutCont )
	{
		for( nu=0; nu < rfield.nflux; nu++ )
		{
			rfield.ConInterOut[nu] = 0.;
		}
	}
	return;
}
