/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_create create data for hydrogen and helium, 1 per coreload, called by ContCreatePointers 
 * in turn called after commands parsed */
/*iso_zero zero data for hydrogen and helium */
#include "cddefines.h"
#include "atmdat.h"
#include "dense.h"
#include "elementnames.h"
#include "helike.h"
#include "helike_einsta.h"
#include "hydro_bauman.h"
#include "hydrogenic.h"
#include "hydroeinsta.h"
#include "iso.h"
#include "lines_service.h"
#include "opacity.h"
#include "phycon.h"
#include "physconst.h"
#include "secondaries.h"
#include "taulines.h"
#include "thirdparty.h"

/*iso_zero zero data for hydrogen and helium */
STATIC void iso_zero(void);

/* allocate memory for iso sequence structures */
STATIC void iso_allocate(void);

/* define levels of iso sequences and assign quantum numbers to those levels */
STATIC void iso_assign_quantum_numbers(void);

STATIC void FillExtraLymanLine( transition *t, long ipISO, long nelem, long nHi );

STATIC void iso_satellite( void );

char chL[21]={'S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V','W','X','Y','Z'};

void iso_create(void)
{
	long int ipHi, 
		ipLo,
		ipISO,
		nelem;

	static int nCalled = 0;

	double HIonPoten;

	DEBUG_ENTRY( "iso_create()" );

	/* > 1 if not first call, then just zero arrays out */
	if( nCalled > 0 )
	{
		iso_zero();
		return;
	}

	/* this is first call, increment the nCalled counterso never do this again */
	++nCalled;

	/* these are the statistical weights of the ions */
	iso.stat_ion[ipH_LIKE] = 1.f;
	iso.stat_ion[ipHE_LIKE] = 2.f;

	/* this routine allocates all the memory
	 * needed for the iso sequence structures */
	iso_allocate();

	/* loop over iso sequences and assign quantum numbers to all levels */
	iso_assign_quantum_numbers();

	/* this is a dummy line, junk it too. */
	TauDummy.Junk();
	TauDummy.Hi = AddState2Stack();
	TauDummy.Lo = AddState2Stack();
	TauDummy.Emis = AddLine2Stack( true );

	/********************************************/
	/**********  Line and level energies ********/
	/********************************************/
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* main hydrogenic arrays, fill with sane values */
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				/* Dima's array has ionization potentials in eV, but not on same
				 * scale as cloudy itself*/
				/* extra factor accounts for this */
				HIonPoten = t_ADfA::Inst().ph1(0,0,nelem,0)/EVRYD* 0.9998787;
				ASSERT(HIonPoten > 0.);

				/* go from ground to the highest level */
				for( ipHi=0; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					double EnergyWN, EnergyRyd;

					if( ipISO == ipH_LIKE )
					{
						EnergyRyd = HIonPoten/POW2((double)N_(ipHi));
					}
					else if( ipISO == ipHE_LIKE )
					{
						EnergyRyd = helike_energy( nelem, ipHi ) * WAVNRYD;
					}
					else
					{
						/* Other iso sequences don't exist yet. */
						TotalInsanity();
					}

					/* >>chng 02 feb 09, change test to >= 0 since we now use 0 for 2s-2p */
					ASSERT(EnergyRyd >= 0.);

					iso.xIsoLevNIonRyd[ipISO][nelem][ipHi] = EnergyRyd;

					/* now loop from ground to level below ipHi */
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						EnergyWN = RYD_INF * (iso.xIsoLevNIonRyd[ipISO][nelem][ipLo] -
							iso.xIsoLevNIonRyd[ipISO][nelem][ipHi]);

						/* This is the minimum line energy we will allow.  */
						/* \todo 2 wire this to lowest energy of code. */
						if( EnergyWN==0 && ipISO==ipHE_LIKE )
							EnergyWN = 0.0001;

						if( EnergyWN < 0. )
							EnergyWN = -1.0 * EnergyWN;

						/* transition energy in various units: */
						Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN = (realnum)EnergyWN;
						Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg = (realnum)(EnergyWN*WAVNRYD*EN1RYD);
						Transitions[ipISO][nelem][ipHi][ipLo].EnergyK = (realnum)(EnergyWN*WAVNRYD*TE1RYD );

						ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN >= 0.);
						ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg >= 0.);
						ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].EnergyK >= 0.);

						/** \todo 2 this should be changed for j-resolved levels */
						if( N_(ipLo)==N_(ipHi) && ipISO==ipH_LIKE )
						{
							Transitions[ipISO][nelem][ipHi][ipLo].WLAng = 0.;
						}
						else
						{
							/* make following an air wavelength */
							Transitions[ipISO][nelem][ipHi][ipLo].WLAng = 
								(realnum)(1.0e8/
							  Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN/
							  RefIndex( Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN));
							ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].WLAng > 0.);
						}

					}
				}

				/* fill the extra Lyman lines */
				for( ipHi=2; ipHi < iso.nLyman_malloc[ipISO]; ipHi++ )
				{
					FillExtraLymanLine( &ExtraLymanLines[ipISO][nelem][ipHi], ipISO, nelem, ipHi );
				}
			}
		}
	}

	/***************************************************************/
	/***** Set up recombination tables for later interpolation *****/
	/***************************************************************/
	/* NB - the above is all we need if we are compiling recombination tables. */
	iso_recomb_malloc();
	iso_recomb_setup( ipH_LIKE );
	iso_recomb_setup( ipHE_LIKE );
	iso_recomb_auxiliary_free();

	/* set up helium collision tables */
	HeCollidSetup();

	/***********************************************************************************/
	/**********  Transition Probabilities, Redistribution Functions, Opacitites ********/
	/***********************************************************************************/
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		if( ipISO == ipH_LIKE )
		{
			/* do nothing here */
		}
		else if( ipISO == ipHE_LIKE )
		{
			/* This routine reads in transition probabilities from a file. */ 
			HelikeTransProbSetup();
		}
		else
		{
			TotalInsanity();
		}

		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipLo=ipH1s; ipLo < (iso.numLevels_max[ipISO][nelem] - 1); ipLo++ )
				{
					for( ipHi=ipLo + 1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
					{
						realnum Aul;

						/* transition prob, EinstA uses current H atom indices */
						if( ipISO == ipH_LIKE )
						{
							Aul = hydro_transprob( nelem, ipHi, ipLo );
						}
						else if( ipISO == ipHE_LIKE )
						{
							Aul = helike_transprob(nelem, ipHi, ipLo);
						}
						else
						{
							TotalInsanity();
						}

						if( Aul <= iso.SmallA )
							Transitions[ipISO][nelem][ipHi][ipLo].Emis = AddLine2Stack( false );
						else
							Transitions[ipISO][nelem][ipHi][ipLo].Emis = AddLine2Stack( true );

						Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul = Aul;

						ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul > 0.);

						if( ipLo == 0 && ipHi == iso.nLyaLevel[ipISO] )
						{
							/* this is La, special redistribution, default is ipLY_A */
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->iRedisFun = iso.ipLyaRedist[ipISO];
						}
						else if( ipLo == 0 )
						{
							/* these are rest of Lyman lines, 
							 * complete redistribution, doppler core only, K2 core, default ipCRD */
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->iRedisFun = iso.ipResoRedist[ipISO];
						}
						else
						{
							/* all lines coming from excited states, default is complete
							 * redis with wings, ipCRDW*/
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->iRedisFun = iso.ipSubRedist[ipISO];
						}

						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA ||
							Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN <= 0.)
						{
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->gf = 0.;
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->opacity = 0.;
						}
						else
						{
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->gf = 
								(realnum)(GetGF(Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul,
								Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN,
								Transitions[ipISO][nelem][ipHi][ipLo].Hi->g));
							ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].Emis->gf > 0.);

							/* derive the abs coef, call to function is gf, wl (A), g_low */
							Transitions[ipISO][nelem][ipHi][ipLo].Emis->opacity = 
								(realnum)(abscf(Transitions[ipISO][nelem][ipHi][ipLo].Emis->gf,
								Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN,
								Transitions[ipISO][nelem][ipHi][ipLo].Lo->g));
							ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].Emis->opacity > 0.);
						}
					}
				}
			}
		}
	}

	/************************************************/
	/**********  Fine Structure Mixing - FSM ********/
	/************************************************/
	if( iso.lgFSM[ipHE_LIKE] )
	{
		/* set some special optical depth values */
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipHi=ipHe2s3S; ipHi<iso.numLevels_max[ipHE_LIKE][nelem]; ipHi++ )
				{
					for( ipLo=ipHe1s1S; ipLo<ipHi; ipLo++ )
					{
						DoFSMixing( nelem, ipLo, ipHi );
					}
				}
			}
		}
	}

	/* following comes out very slightly off, correct here */
	Transitions[ipH_LIKE][ipHELIUM][ipH3s][ipH2s].WLAng = 1640.f;
	Transitions[ipH_LIKE][ipHELIUM][ipH3s][ipH2p].WLAng = 1640.f;
	if( iso.n_HighestResolved_max[ipH_LIKE][ipHELIUM] >=3 )
	{
		Transitions[ipH_LIKE][ipHELIUM][ipH3p][ipH2s].WLAng = 1640.f;
		Transitions[ipH_LIKE][ipHELIUM][ipH3p][ipH2p].WLAng = 1640.f;
		Transitions[ipH_LIKE][ipHELIUM][ipH3d][ipH2s].WLAng = 1640.f;
		Transitions[ipH_LIKE][ipHELIUM][ipH3d][ipH2p].WLAng = 1640.f;
	}

	/****************************************************/
	/**********  lifetimes and damping constants ********/
	/****************************************************/
	for( ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* define these for H and He always */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				/* these are not defined and must never be used */
				StatesElemNEW[nelem][nelem-ipISO][0].lifetime = -FLT_MAX;

				for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime = SMALLFLOAT;
					
					for( ipLo=0; ipLo < ipHi; ipLo++ )  
					{
						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
							continue;

						StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime += Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul;
					}

					/* sum of A's was just stuffed, now invert for lifetime. */
					StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime = 1./StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime;

					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						if( Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN <= 0. )
							continue;

						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
							continue;

						Transitions[ipISO][nelem][ipHi][ipLo].Emis->dampXvel = (realnum)( 
							(1.f/StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime)/
							PI4/Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN);

						ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].Emis->dampXvel> 0.);
					}
				}
			}
		}
	}

	/********************************************/
	/**********  Fix some 2-photon stuff ********/
	/********************************************/
	for( ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		/* set some special optical depth values */
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				/* Must think out two-photon treatment for other sequences.  */
				ASSERT( ipISO <= ipHE_LIKE );

				/* total optical depth in 1s - 2s */
				Transitions[ipISO][nelem][1+ipISO][0].Emis->TauTot = 0.;

				/* opacity in two-photon transition is incorrect - actually a continuum,
				 * so divide by typical width*/
				Transitions[ipISO][nelem][1+ipISO][0].Emis->opacity /= 1e4f;

				/* wavelength of 0 is sentinel for two-photon emission */
				Transitions[ipISO][nelem][1+ipISO][0].WLAng = 0.;
			}
		}
	}


	/* zero out some line information */
	iso_zero();

	/* loop over iso sequences */
	for( ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
			{
				/* calculate cascade probabilities, branching ratios, and associated errors. */
				iso_cascade( ipISO, nelem);
			}
		}
	}

	iso_satellite();

	iso_satellite_update();

	/***************************************/
	/**********  Stark Broadening **********/
	/***************************************/
	/* fill in iso.strkar array */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipLo=ipH1s; ipLo < (iso.numLevels_max[ipISO][nelem] - 1); ipLo++ )
				{
					for( ipHi= ipLo + 1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
					{
						long nHi = StatesElemNEW[nelem][nelem-ipISO][ipHi].n;
						long nLo = StatesElemNEW[nelem][nelem-ipISO][ipLo].n;

						iso.strkar[ipISO][nelem][ipLo][ipHi] = (realnum)pow((realnum)( nLo * nHi ),(realnum)1.2f);
						iso.pestrk[ipISO][nelem][ipLo][ipHi] = 0.;
					}
				}
			}
		}
	}

	return;
}

/* ============================================================================== */
STATIC void iso_zero(void)
{
	long int i, 
	  ipHi, 
	  ipISO,
	  ipLo, 
	  nelem;

	DEBUG_ENTRY( "iso_zero()" );

	fixit(); /* this routine appears to be completely unnecessary because all of 
			  * this is done in RT_tau_init shortly later. */

	hydro.HLineWidth = 0.;

	/****************************************************/
	/**********  initialize some variables **********/
	/****************************************************/
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipHi=0; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop = 0.;

					iso.PopLTE[ipISO][nelem][ipHi] = 0.;
					iso.ColIoniz[ipISO][nelem][ipHi] = 0.;
					iso.gamnc[ipISO][nelem][ipHi] = -DBL_MAX;
					iso.RecomInducRate[ipISO][nelem][ipHi] = -DBL_MAX;
					iso.DepartCoef[ipISO][nelem][ipHi] = -DBL_MAX;
					iso.RateLevel2Cont[ipISO][nelem][ipHi] = 0.;
					iso.RateCont2Level[ipISO][nelem][ipHi] = 0.;
					iso.ConOpacRatio[ipISO][nelem][ipHi] = 1.;
					iso.RadRecCon[ipISO][nelem][ipHi] = 0.;
					iso.RadRecomb[ipISO][nelem][ipHi][ipRecRad] = 0.;
					iso.RadRecomb[ipISO][nelem][ipHi][ipRecNetEsc] = 1.;
					iso.RadRecomb[ipISO][nelem][ipHi][ipRecEsc] = 1.;
					iso.DielecRecomb[ipISO][nelem][ipHi] = 0.;
				}
			}
		}
	}

	/* ground state of H and He is different since totally determine
	 * their own opacities */
	iso.ConOpacRatio[ipH_LIKE][ipHYDROGEN][0] = 1e-5;
	iso.ConOpacRatio[ipH_LIKE][ipHELIUM][0] = 1e-5;
	iso.ConOpacRatio[ipHE_LIKE][ipHELIUM][0] = 1e-5;


	/* main isoelectronic arrays, fill with sane values */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* >>chng 05 dec 30, move nelem.numLevel to numLevels_local and int here */
			iso.numLevels_local[ipISO][nelem] = iso.numLevels_max[ipISO][nelem];

			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				SatelliteLines[ipISO][nelem][0].Zero();
				for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					SatelliteLines[ipISO][nelem][ipHi].Zero();

					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						Transitions[ipISO][nelem][ipHi][ipLo].Zero();
					}
				}

				for( ipHi=2; ipHi < iso.nLyman[ipISO]; ipHi++ )
				{
					ExtraLymanLines[ipISO][nelem][ipHi].Zero();
				}
			}
		}
	}

	/* initialize heavy element line array */
	for( i=0; i <= nLevel1; i++ )
	{
		TauLines[i].Zero();
	}

	/* this is a dummy optical depth array for non-existant lines 
	 * when this goes over to struc, make sure all are set to zero here since
	 * init in grid may depend on it */
	TauDummy.Zero();

	for( i=0; i < nUTA; i++ )
	{
		/* heat is special for this array - it is heat per pump */
		double hsave = UTALines[i].Coll.heat;
		UTALines[i].Zero();
		UTALines[i].Coll.heat = hsave;
	}

	for( i=0; i < nWindLine; i++ )
	{
		TauLine2[i].Zero();
	}

	for( i=0; i < nHFLines; i++ )
	{
		HFLines[i].Zero();
	}

	/* database lines - DB lines */
	for( i=0; i <linesAdded2; i++)
	{
		dBaseLines[i].tran->Zero();
	}

	return;
}

STATIC void iso_allocate(void)
{

	DEBUG_ENTRY( "iso_allocate()" );

	iso.strkar.reserve( NISO );
	iso.ipOpac.reserve( NISO );
	iso.RadRecomb.reserve( NISO );
	iso.RadRecCon.reserve( NISO );
	iso.DielecRecombVsTemp.reserve( NISO );
	iso.Boltzmann.reserve( NISO );
	iso.QuantumNumbers2Index.reserve( NISO );
	iso.BranchRatio.reserve( NISO );
	iso.CascadeProb.reserve( NISO );
	iso.CachedAs.reserve( NISO );

	/* the hydrogen and helium like iso-sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		iso.strkar.reserve( ipISO, LIMELM );
		iso.ipOpac.reserve( ipISO, LIMELM );
		iso.RadRecomb.reserve( ipISO, LIMELM );
		iso.RadRecCon.reserve( ipISO, LIMELM );
		iso.DielecRecombVsTemp.reserve( ipISO, LIMELM );
		iso.Boltzmann.reserve( ipISO, LIMELM );
		iso.QuantumNumbers2Index.reserve( ipISO, LIMELM );
		iso.BranchRatio.reserve( ipISO, LIMELM );
		iso.CascadeProb.reserve( ipISO, LIMELM );
		iso.CachedAs.reserve( ipISO, LIMELM );

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				iso.numLevels_malloc[ipISO][nelem] = iso.numLevels_max[ipISO][nelem];

				/*fprintf(ioQQQ,"assert failll %li\t%li\t%li\n", ipISO, nelem , iso.numLevels_max[ipISO][nelem] );*/
				ASSERT( iso.numLevels_max[ipISO][nelem] > 0 );

				iso.nLyman_malloc[ipISO] = iso.nLyman[ipISO];

				iso.strkar.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				iso.ipOpac.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				iso.RadRecomb.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				iso.RadRecCon.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				iso.DielecRecombVsTemp.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				iso.Boltzmann.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );

				iso.CachedAs.reserve( ipISO, nelem, MAX2(1, iso.nCollapsed_max[ipISO][nelem]) );

				for( long i = 0; i < iso.nCollapsed_max[ipISO][nelem]; ++i )
				{
					iso.CachedAs.reserve( ipISO, nelem, i, iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem] );
					for( long i1 = 0; i1 < iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem]; ++i1 )
					{
						/* allocate two spaces delta L +/- 1 */
						iso.CachedAs.reserve( ipISO, nelem, i, i1, 2 );
					}
				}

				iso.QuantumNumbers2Index.reserve( ipISO, nelem, iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] + 1 );


				for( long i = 1; i <= (iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem]); ++i )
				{
					/* Allocate proper number of angular momentum quantum number.	*/
					iso.QuantumNumbers2Index.reserve( ipISO, nelem, i, i );

					for( long i1 = 0; i1 < i; ++i1 )
					{
						/* This may have to change for other iso sequences. */
						ASSERT( NISO == 2 );
						/* Allocate 4 spaces for multiplicity.	H-like will be accessed with "2" for doublet,
						 * He-like will be accessed via "1" for singlet or "3" for triplet.  "0" will not be used. */
						iso.QuantumNumbers2Index.reserve( ipISO, nelem, i, i1, 4 );
					}
				}

				for( long n=0; n < iso.numLevels_max[ipISO][nelem]; ++n )
				{
					iso.RadRecomb.reserve( ipISO, nelem, n, 3 );

					/* sec to last dim is upper level n,
					 * last dim of this array is lower level, will go from 0 to n-1 */
					if( n > 0 )
						iso.Boltzmann.reserve( ipISO, nelem, n, n );

					/* malloc space for number of temperature points in Badnell grids */
					iso.DielecRecombVsTemp.reserve( ipISO, nelem, n, NUM_DR_TEMPS );
				}

				/* In this array are stored the C values described in Robbins 68. */
				iso.CascadeProb.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				iso.BranchRatio.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem]+1 );

				for( long i = 1; i < iso.numLevels_max[ipISO][nelem]; i++ )
					iso.BranchRatio.reserve( ipISO, nelem, i, i+1 );

				/* reserve final dimension of cascade probabilities. */
				for( long i = 0; i < iso.numLevels_max[ipISO][nelem]; ++i )
				{
					iso.strkar.reserve( ipISO, nelem, i, iso.numLevels_max[ipISO][nelem] );
					iso.CascadeProb.reserve( ipISO, nelem, i, i+1 );
				}
			}
		}
	}

	iso.strkar.alloc();
	iso.pestrk.alloc( iso.strkar.clone() );
	iso.ipOpac.alloc();
	secondaries.Hx12.alloc( iso.ipOpac.clone() );
	iso.ipIsoLevNIonCon.alloc( iso.ipOpac.clone() );
	iso.xIsoLevNIonRyd.alloc( iso.ipOpac.clone() );
	iso.DielecRecomb.alloc( iso.ipOpac.clone() );
	iso.DepartCoef.alloc( iso.ipOpac.clone() );
	iso.RateLevel2Cont.alloc( iso.ipOpac.clone() );
	iso.RateCont2Level.alloc( iso.ipOpac.clone() );
	iso.ConOpacRatio.alloc( iso.ipOpac.clone() );
	iso.gamnc.alloc( iso.ipOpac.clone() );
	iso.RecomInducRate.alloc( iso.ipOpac.clone() );
	iso.RecomInducCool_Coef.alloc( iso.ipOpac.clone() );
	iso.PhotoHeat.alloc( iso.ipOpac.clone() );
	iso.PopLTE.alloc( iso.ipOpac.clone() );
	iso.ColIoniz.alloc( iso.ipOpac.clone() );
	iso.RadEffec.alloc( iso.ipOpac.clone() );
	iso.RadRecCon.alloc();
	iso.RadRecomb.alloc();
	iso.Boltzmann.alloc();
	iso.DielecRecombVsTemp.alloc();
	iso.CachedAs.alloc();
	iso.QuantumNumbers2Index.alloc();
	iso.BranchRatio.alloc();
	iso.CascadeProb.alloc();

	iso.bnl_effective.alloc( iso.QuantumNumbers2Index.clone() );

	iso.DielecRecombVsTemp.zero();
	iso.QuantumNumbers2Index.invalidate();
	iso.bnl_effective.invalidate();
	iso.BranchRatio.invalidate();
	iso.CascadeProb.invalidate();

	SatelliteLines.reserve( NISO );
	Transitions.reserve( NISO );
	ExtraLymanLines.reserve( NISO );

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		SatelliteLines.reserve( ipISO, LIMELM );
		Transitions.reserve( ipISO, LIMELM );
		ExtraLymanLines.reserve( ipISO, LIMELM );

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				ASSERT( iso.numLevels_max[ipISO][nelem] > 0 );

				SatelliteLines.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				Transitions.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
				ExtraLymanLines.reserve( ipISO, nelem, iso.nLyman_malloc[ipISO] );

				for( long ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
					Transitions.reserve( ipISO, nelem, ipHi, ipHi );
			}
		}
	}

	SatelliteLines.alloc();
	Transitions.alloc();
	ExtraLymanLines.alloc();

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( long ipLo=0; ipLo<iso.numLevels_max[ipISO][nelem]; ipLo++ )
				{
					/* Upper level is continuum, use a generic state
					 * lower level is the same as the index. */
					SatelliteLines[ipISO][nelem][ipLo].Junk();
					SatelliteLines[ipISO][nelem][ipLo].Hi = AddState2Stack();
					SatelliteLines[ipISO][nelem][ipLo].Lo = &StatesElemNEW[nelem][nelem-ipISO][ipLo];
					SatelliteLines[ipISO][nelem][ipLo].Emis = AddLine2Stack( true );
				}

				for( long ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{
						/* set ENTIRE array to impossible values, in case of bad pointer */
						Transitions[ipISO][nelem][ipHi][ipLo].Junk();
						Transitions[ipISO][nelem][ipHi][ipLo].Hi = &StatesElemNEW[nelem][nelem-ipISO][ipHi];
						Transitions[ipISO][nelem][ipHi][ipLo].Lo = &StatesElemNEW[nelem][nelem-ipISO][ipLo];
					}
				}

				/* junk the extra Lyman lines */
				for( long ipHi=2; ipHi < iso.nLyman_malloc[ipISO]; ipHi++ )
				{
					ExtraLymanLines[ipISO][nelem][ipHi].Junk();
					ExtraLymanLines[ipISO][nelem][ipHi].Hi = AddState2Stack();
					/* lower level is just ground state of the ion */
					ExtraLymanLines[ipISO][nelem][ipHi].Lo = &StatesElemNEW[nelem][nelem-ipISO][0];
					ExtraLymanLines[ipISO][nelem][ipHi].Emis = AddLine2Stack( true );
				}
			}
		}
	}

	/* malloc space for random error generation, if appropriate flags set. */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		static bool lgNotSetYet = true;

		if( iso.lgRandErrGen[ipISO] && lgNotSetYet )
		{
			iso.Error.reserve( NISO );
			iso.SigmaCascadeProb.reserve( NISO );
			lgNotSetYet = false;
		}

		if( iso.lgRandErrGen[ipISO] )
		{
			ASSERT( !lgNotSetYet );
			iso.Error.reserve( ipISO, LIMELM );
			iso.SigmaCascadeProb.reserve( ipISO, LIMELM );

			for( long nelem=ipISO; nelem < LIMELM; ++nelem )
			{
				if( nelem == ipHELIUM || dense.lgElmtOn[nelem] )
				{
					/* Here we add one slot for rates involving the continuum.  */
					iso.Error.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem]+1 );

					for( long i = 1; i< iso.numLevels_max[ipISO][nelem] + 1; i++ )
					{
						iso.Error.reserve( ipISO, nelem, i, i+1 );

						for( long j = 0; j<i+1; j++ )
							iso.Error.reserve( ipISO, nelem, i, j, 3 );
					}

					/* Uncertainty in cascade probability and effective recombination,
					 * only when error generation is turned on. */
					iso.SigmaCascadeProb.reserve( ipISO, nelem, iso.numLevels_max[ipISO][nelem] );
					for( long i=0; i < iso.numLevels_max[ipISO][nelem]; i++ )
						/* error in cascade probabilities, for all levels */
						iso.SigmaCascadeProb.reserve( ipISO, nelem, i, i+1 );
				}
			}
			iso.Error.alloc();
			iso.ErrorFactor.alloc( iso.Error.clone() );
			iso.SigmaAtot.alloc( iso.ipOpac.clone() );
			iso.SigmaRadEffec.alloc( iso.RadEffec.clone() );
			iso.SigmaCascadeProb.alloc();

			iso.Error.invalidate();
			iso.ErrorFactor.invalidate();		
		}
	}

	return;
}

STATIC void iso_assign_quantum_numbers(void)
{
	long int
	  ipISO,
	  ipLo,
	  level,
	  nelem,
	  i,
	  in,
	  il,
	  is,
	  ij;

	DEBUG_ENTRY( "iso_assign_quantum_numbers()" );

	ipISO = ipH_LIKE;
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		/* only check elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem] )
		{
			i = 0;

			/* 2 for doublet */
			is = 2;

			/* fill in StatesElemNEW and iso.QuantumNumbers2Index arrays. */
			/* this loop is over quantum number n */
			for( in = 1L; in <= iso.n_HighestResolved_max[ipISO][nelem]; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					StatesElemNEW[nelem][nelem-ipISO][i].n = in;
					StatesElemNEW[nelem][nelem-ipISO][i].S = is;
					StatesElemNEW[nelem][nelem-ipISO][i].l = il;
					StatesElemNEW[nelem][nelem-ipISO][i].j = -1;
					iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] = i;
					++i;
				}
			}
			/* now do the collapsed levels */
			in = iso.n_HighestResolved_max[ipISO][nelem] + 1;
			for( level = i; level< iso.numLevels_max[ipISO][nelem]; ++level)
			{
				StatesElemNEW[nelem][nelem-ipISO][level].n = in;
				StatesElemNEW[nelem][nelem-ipISO][level].S = -LONG_MAX;
				StatesElemNEW[nelem][nelem-ipISO][level].l = -LONG_MAX;
				StatesElemNEW[nelem][nelem-ipISO][level].j = -1;
				/* Point every l to same index for collapsed levels.	*/
				for( il = 0; il < in; ++il )
				{
					iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] = level;
				}
				++in;
			}
			--in;

			/* confirm that we did not overrun the array */
			ASSERT( i <= iso.numLevels_max[ipISO][nelem] );

			/* confirm that n is positive and not greater than the max n. */
			ASSERT( (in > 0) && (in < (iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] + 1) ) );

			/* Verify StatesElemNEW and iso.QuantumNumbers2Index[ipISO] agree in all cases	*/
			for( in = 2L; in <= iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem]; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					ASSERT( StatesElemNEW[nelem][nelem-ipISO][ iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] ].n == in );
					if( in <= iso.n_HighestResolved_max[ipISO][nelem] )
					{
						/* Must only check these for resolved levels...
						 * collapsed levels in StatesElemNEW have pointers for l and s that will blow if used.	*/
						ASSERT( StatesElemNEW[nelem][nelem-ipISO][ iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] ].l == il );
						ASSERT( StatesElemNEW[nelem][nelem-ipISO][ iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] ].S == is );
					}
				}
			}
		}
	}

	/* then do he-like */
	ipISO = ipHE_LIKE;
	for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{
		/* only check elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem] )
		{
			i = 0;

			/* fill in StatesElemNEW and iso.QuantumNumbers2Index arrays. */
			/* this loop is over quantum number n */
			for( in = 1L; in <= iso.n_HighestResolved_max[ipISO][nelem]; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					for( is = 3L; is >= 1L; is -= 2 )
					{
						/* All levels except singlet P follow the ordering scheme:	*/
						/*	lower l's have lower energy	*/
						/* 	triplets have lower energy	*/
						if( (il == 1L) && (is == 1L) )
							continue;
						/* n = 1 has no triplet, of course.	*/
						if( (in == 1L) && (is == 3L) )
							continue;

						/* singlets */
						if( is == 1 )
						{
							StatesElemNEW[nelem][nelem-ipISO][i].n = in;
							StatesElemNEW[nelem][nelem-ipISO][i].S = is;
							StatesElemNEW[nelem][nelem-ipISO][i].l = il;
							/* this is not a typo, J=L for singlets.  */
							StatesElemNEW[nelem][nelem-ipISO][i].j = il;
							iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] = i;
							++i;
						}
						/* 2 triplet P is j-resolved */
						else if( (in == 2) && (il == 1) && (is == 3) )
						{
							ij = 0;
							do 
							{
								StatesElemNEW[nelem][nelem-ipISO][i].n = in;
								StatesElemNEW[nelem][nelem-ipISO][i].S = is;
								StatesElemNEW[nelem][nelem-ipISO][i].l = il;
								StatesElemNEW[nelem][nelem-ipISO][i].j = ij;
								iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] = i;
								++i;
								++ij;
								/* repeat this for the separate j-levels within 2^3P. */
							}	while ( ij < 3 );
						}
						else
						{
							StatesElemNEW[nelem][nelem-ipISO][i].n = in;
							StatesElemNEW[nelem][nelem-ipISO][i].S = is;
							StatesElemNEW[nelem][nelem-ipISO][i].l = il;
							StatesElemNEW[nelem][nelem-ipISO][i].j = -1L;
							iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] = i;
							++i;
						}
					}
				}
				/*	Insert singlet P at the end of every sequence for a given n.	*/
				if( in > 1L )
				{
					StatesElemNEW[nelem][nelem-ipISO][i].n = in;
					StatesElemNEW[nelem][nelem-ipISO][i].S = 1L;
					StatesElemNEW[nelem][nelem-ipISO][i].l = 1L;
					StatesElemNEW[nelem][nelem-ipISO][i].j = 1L;
					iso.QuantumNumbers2Index[ipISO][nelem][in][1][1] = i;
					++i;
				}
			}
			/* now do the collapsed levels */
			in = iso.n_HighestResolved_max[ipISO][nelem] + 1;
			for( level = i; level< iso.numLevels_max[ipISO][nelem]; ++level)
			{
				StatesElemNEW[nelem][nelem-ipISO][level].n = in;
				StatesElemNEW[nelem][nelem-ipISO][level].S = -LONG_MAX;
				StatesElemNEW[nelem][nelem-ipISO][level].l = -LONG_MAX;
				StatesElemNEW[nelem][nelem-ipISO][level].j = -1;
				/* Point every l and s to same index for collapsed levels.	*/
				for( il = 0; il < in; ++il )
				{
					for( is = 1; is <= 3; is += 2 )
					{
						iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] = level;
					}
				}
				++in;
			}
			--in;

			/* confirm that we did not overrun the array */
			ASSERT( i <= iso.numLevels_max[ipISO][nelem] );

			/* confirm that n is positive and not greater than the max n. */
			ASSERT( (in > 0) && (in < (iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] + 1) ) );

			/* Verify StatesElemNEW and iso.QuantumNumbers2Index[ipISO] agree in all cases	*/
			for( in = 2L; in <= iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem]; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					for( is = 3L; is >= 1; is -= 2 )
					{
						/* Ground state is not triplicate.	*/
						if( (in == 1L) && (is == 3L) )
							continue;

						ASSERT( StatesElemNEW[nelem][nelem-ipISO][ iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] ].n == in );
						if( in <= iso.n_HighestResolved_max[ipISO][nelem] )
						{
							/* Must only check these for resolved levels...
							 * collapsed levels in StatesElemNEW have pointers for l and s that will blow if used.	*/
							ASSERT( StatesElemNEW[nelem][nelem-ipISO][ iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] ].l == il );
							ASSERT( StatesElemNEW[nelem][nelem-ipISO][ iso.QuantumNumbers2Index[ipISO][nelem][in][il][is] ].S == is );
						}
					}
				}
			}
		}
	}

	for( ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipLo=ipH1s; ipLo < iso.numLevels_max[ipISO][nelem]; ipLo++ )
				{
					StatesElemNEW[nelem][nelem-ipISO][ipLo].nelem = (int)(nelem+1);
					StatesElemNEW[nelem][nelem-ipISO][ipLo].IonStg = (int)(nelem+1-ipISO);

					if( StatesElemNEW[nelem][nelem-ipISO][ipLo].j >= 0 )
					{
						StatesElemNEW[nelem][nelem-ipISO][ipLo].g = 2.f*StatesElemNEW[nelem][nelem-ipISO][ipLo].j+1.f;
					}
					else if( StatesElemNEW[nelem][nelem-ipISO][ipLo].l >= 0 )
					{
						StatesElemNEW[nelem][nelem-ipISO][ipLo].g = (2.f*StatesElemNEW[nelem][nelem-ipISO][ipLo].l+1.f) *
							StatesElemNEW[nelem][nelem-ipISO][ipLo].S;
					}
					else
					{
						if( ipISO == ipH_LIKE )
							StatesElemNEW[nelem][nelem-ipISO][ipLo].g = 2.f*(realnum)POW2( StatesElemNEW[nelem][nelem-ipISO][ipLo].n );
						else if( ipISO == ipHE_LIKE )
							StatesElemNEW[nelem][nelem-ipISO][ipLo].g = 4.f*(realnum)POW2( StatesElemNEW[nelem][nelem-ipISO][ipLo].n );
						else
						{
							/* replace this with correct thing if more sequences are added. */
							TotalInsanity();
						}
					}
					char chConfiguration[11] = "          ";
					long nCharactersWritten = 0;

					ASSERT( StatesElemNEW[nelem][nelem-ipISO][ipLo].n < 1000 );

					/* include j only if defined.  */
					if( StatesElemNEW[nelem][nelem-ipISO][ipLo].n > iso.n_HighestResolved_max[ipISO][nelem] )
					{
						nCharactersWritten = sprintf( chConfiguration, "n=%3li", 
							StatesElemNEW[nelem][nelem-ipISO][ipLo].n );
					}
					else if( StatesElemNEW[nelem][nelem-ipISO][ipLo].j > 0 )
					{
						nCharactersWritten = sprintf( chConfiguration, "%3li^%li%c_%li", 
							StatesElemNEW[nelem][nelem-ipISO][ipLo].n, 
							StatesElemNEW[nelem][nelem-ipISO][ipLo].S,
							chL[ MIN2( 20, StatesElemNEW[nelem][nelem-ipISO][ipLo].l ) ],
							StatesElemNEW[nelem][nelem-ipISO][ipLo].j );
					}
					else
					{
						nCharactersWritten = sprintf( chConfiguration, "%3li^%li%c", 
							StatesElemNEW[nelem][nelem-ipISO][ipLo].n, 
							StatesElemNEW[nelem][nelem-ipISO][ipLo].S,
							chL[ MIN2( 20, StatesElemNEW[nelem][nelem-ipISO][ipLo].l) ] );
					}

					ASSERT( nCharactersWritten <= 10 );
					chConfiguration[10] = '\0';

					strncpy( StatesElemNEW[nelem][nelem-ipISO][ipLo].chConfig, chConfiguration, 10 );
				}
			}
		}
	}
	return;
}

#if defined(__ICC) && defined(__i386)
#pragma optimization_level 1
#endif
STATIC void FillExtraLymanLine( transition *t, long ipISO, long nelem, long nHi )
{
	double Enerwn, Aul;

	DEBUG_ENTRY( "FillExtraLymanLine()" );

	/* atomic number or charge and stage: */
	t->Hi->nelem = (int)(nelem+1);
	t->Hi->IonStg = (int)(nelem+1-ipISO);
	
	/* statistical weight is same as statistical weight of corresponding LyA. */
	t->Hi->g = StatesElemNEW[nelem][nelem-ipISO][iso.nLyaLevel[ipISO]].g;

	/* energies */
	Enerwn = iso.xIsoLevNIonRyd[ipISO][nelem][0] * RYD_INF * (  1. - 1./POW2((double)nHi) );

	/* transition energy in various units:*/
	t->EnergyWN = (realnum)(Enerwn);
	t->EnergyErg = (realnum)(Enerwn * WAVNRYD * EN1RYD);
	t->EnergyK = (realnum)(Enerwn * WAVNRYD * TE1RYD);
	t->WLAng = (realnum)(1.0e8/ Enerwn/ RefIndex(Enerwn));

	if(	ipISO == ipH_LIKE )
	{
		Aul = H_Einstein_A( nHi, 1, 1, 0, nelem+1 );
	}
	else
	{
		if( nelem == ipHELIUM )
		{
			/* A simple fit for the calculation of Helium lyman Aul's.	*/
			Aul = (1.508e10) / pow((double)nHi,2.975);
		}
		else 
		{
			/* Fit to values given in 
			 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
			 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
			/* originally astro.ph. 0201454  */
			Aul = 1.375E10 * pow((double)nelem, 3.9) / pow((double)nHi,3.1);
		}
	}

	t->Emis->Aul = (realnum)Aul;

	t->Hi->lifetime = iso_state_lifetime( ipISO, nelem, nHi, 1 );

	t->Emis->dampXvel = (realnum)( 1.f / t->Hi->lifetime / PI4 / t->EnergyWN );

	t->Emis->iRedisFun = iso.ipResoRedist[ipISO];

	t->Emis->gf = (realnum)(GetGF(t->Emis->Aul,	t->EnergyWN, t->Hi->g));

	/* derive the abs coef, call to function is Emis.gf, wl (A), g_low */
	t->Emis->opacity = (realnum)(abscf(t->Emis->gf, t->EnergyWN, t->Lo->g));

	/* create array indices that will blow up */
	t->ipCont = INT_MIN;
	t->Emis->ipFine = INT_MIN;

	{
		/* option to print particulars of some line when called
			* a prettier print statement is near where chSpin is defined below
			* search for "pretty print" */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"%li\t%li\t%.2e\t%.2e\n",
				nelem+1,
				nHi,
				t->Emis->Aul , 	
				t->Emis->opacity
				);
		}
	}
	return;
}

/* calculate radiative lifetime of an individual iso state */
double iso_state_lifetime( long ipISO, long nelem, long n, long l )
{
	/* >>refer hydro lifetimes	Horbatsch, M. W., Horbatsch, M. and Hessels, E. A. 2005, JPhysB, 38, 1765 */

	double tau, t0, eps2;
	/* mass of electron */
	double m = ELECTRON_MASS;
	/* nuclear mass */
	double M = (double)dense.AtomicWeight[nelem] * ATOMIC_MASS_UNIT;
	double mu = (m*M)/(M+m);
	long z = 1;
	long Z = nelem + 1 - ipISO;
	
	DEBUG_ENTRY( "iso_state_lifetime()" );

	/* this should not be used for l=0 per the Horbatsch et al. paper */
	ASSERT( l > 0 );

	eps2 = 1. - ( l*l + l + 8./47. - (l+1.)/69./n ) / POW2( (double)n );

	t0 = 3. * H_BAR * pow( (double)n, 5.) / 
		( 2. * POW4( (double)( z * Z ) ) * pow( FINE_STRUCTURE, 5. ) * mu * POW2( SPEEDLIGHT ) ) *
		POW2( (m + M)/(Z*m + z*M) );

	tau = t0 * ( 1. - eps2 ) / 
		( 1. + 19./88.*( (1./eps2 - 1.) * log( 1. - eps2 ) + 1. - 
		0.5 * eps2 - 0.025 * eps2 * eps2 ) );

	if( ipISO == ipHE_LIKE )
	{
		/* iso_state_lifetime is not spin specific, must exclude helike triplet here. */
		tau /= 3.;
		/* this is also necessary to correct the helike lifetimes */
		tau *= 1.1722 * pow( (double)nelem, 0.1 );
	}

	/* would probably need a new lifetime algorithm for any other iso sequences. */
	ASSERT( ipISO <= ipHE_LIKE );
	ASSERT( tau > 0. );

	return tau;
}

/* calculate cascade probabilities, branching ratios, and associated errors. */
void iso_cascade( long ipISO, long nelem )
{
	/* The sum of all A's coming out of a given n,
	 * Below we assert a monotonic trend. */
	double *SumAPerN;

	long int i, j, ipLo, ipHi;

	DEBUG_ENTRY( "iso_cascade()" );

	SumAPerN = ((double*)MALLOC( (size_t)(iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] + 1 )*sizeof(double )));
	memset( SumAPerN, 0, (iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] + 1 )*sizeof(double ) );

	/* Initialize some ground state stuff, easier here than in loops.	*/
	iso.CascadeProb[ipISO][nelem][0][0] = 1.;
	if( iso.lgRandErrGen[ipISO] )
	{
		iso.SigmaAtot[ipISO][nelem][0] = 0.;
		iso.SigmaCascadeProb[ipISO][nelem][0][0] = 0.;
	}

	/***************************************************************************/
	/****** Cascade probabilities, Branching ratios, and associated errors *****/
	/***************************************************************************/
	for( ipHi=1; ipHi<iso.numLevels_max[ipISO][nelem]; ipHi++ )
	{
		double SumAs = 0.;

		/** Cascade probabilities are as defined in Robbins 68,
		 * generalized here for cascade probability for any iso sequence.	
		 * >>refer He triplets	Robbins, R.R. 1968, ApJ 151, 497R	
		 * >>refer He triplets	Robbins, R.R. 1968a, ApJ 151, 511R	*/

		/* initialize variables. */
		iso.CascadeProb[ipISO][nelem][ipHi][ipHi] = 1.;
		iso.CascadeProb[ipISO][nelem][ipHi][0] = 0.;
		iso.BranchRatio[ipISO][nelem][ipHi][0] = 0.;

		if( iso.lgRandErrGen[ipISO] )
		{
			iso.SigmaAtot[ipISO][nelem][ipHi] = 0.;
			iso.SigmaCascadeProb[ipISO][nelem][ipHi][ipHi] = 0.;
		}

		long ipLoStart = 0;
		if( opac.lgCaseB && L_(ipHi)==1 && (ipISO==ipH_LIKE || S_(ipHi)==1) )
			ipLoStart = 1;

		for( ipLo=ipLoStart; ipLo<ipHi; ipLo++ )
		{
			SumAs += Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul;
		}

		for( ipLo=ipLoStart; ipLo<ipHi; ipLo++ )
		{
			if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
			{
				iso.CascadeProb[ipISO][nelem][ipHi][ipLo] = 0.;
				iso.BranchRatio[ipISO][nelem][ipHi][ipLo] = 0.;
				continue;
			}

			iso.CascadeProb[ipISO][nelem][ipHi][ipLo] = 0.;
			iso.BranchRatio[ipISO][nelem][ipHi][ipLo] = 
				Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul / SumAs;

			ASSERT( iso.BranchRatio[ipISO][nelem][ipHi][ipLo] <= 1.0000001 );

			SumAPerN[N_(ipHi)] += Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul;

			/* there are some negative energy transitions, where the order
			 * has changed, but these are not optically allowed, these are
			 * same n, different L, forbidden transitions */
			ASSERT( Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN > 0. ||
				Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA );

			if( iso.lgRandErrGen[ipISO] )
			{
				ASSERT( iso.Error[ipISO][nelem][ipHi][ipLo][IPRAD] >= 0. );
				/* Uncertainties in A's are added in quadrature, square root is taken below. */
				iso.SigmaAtot[ipISO][nelem][ipHi] += 
					pow( iso.Error[ipISO][nelem][ipHi][ipLo][IPRAD] * 
					(double)Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul, 2. );
			}
		}

		if( iso.lgRandErrGen[ipISO] )
		{
			/* Uncertainties in A's are added in quadrature above, square root taken here. */
			iso.SigmaAtot[ipISO][nelem][ipHi] = sqrt( iso.SigmaAtot[ipISO][nelem][ipHi] );
		}

		/* cascade probabilities */
		for( ipLo=0; ipLo<ipHi; ipLo++ )
		{
			double SigmaA, SigmaCul = 0.;

			for( i=ipLo; i<ipHi; i++ )
			{
				iso.CascadeProb[ipISO][nelem][ipHi][ipLo] += iso.BranchRatio[ipISO][nelem][ipHi][i] * iso.CascadeProb[ipISO][nelem][i][ipLo];
				if( iso.lgRandErrGen[ipISO] )
				{
					if( Transitions[ipISO][nelem][ipHi][i].Emis->Aul > iso.SmallA )
					{
						/* Uncertainties in A's and cascade probabilities */
						SigmaA = iso.Error[ipISO][nelem][ipHi][i][IPRAD] * 
							Transitions[ipISO][nelem][ipHi][i].Emis->Aul;
						SigmaCul += 
							pow(SigmaA*iso.CascadeProb[ipISO][nelem][i][ipLo]*StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime, 2.) +
							pow(iso.SigmaAtot[ipISO][nelem][ipHi]*iso.BranchRatio[ipISO][nelem][ipHi][i]*iso.CascadeProb[ipISO][nelem][i][ipLo]*StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime, 2.) +
							pow(iso.SigmaCascadeProb[ipISO][nelem][i][ipLo]*iso.BranchRatio[ipISO][nelem][ipHi][i], 2.);
					}
				}
			}
			if( iso.lgRandErrGen[ipISO] )
			{
				SigmaCul = sqrt(SigmaCul);
				iso.SigmaCascadeProb[ipISO][nelem][ipHi][ipLo] = SigmaCul;
			}
		}
	}

	/************************************************************************/
	/*** Allowed decay conversion probabilities. See Robbins68b, Table 1. ***/
	/************************************************************************/
	{
		enum {DEBUG_LOC=false};

		if( DEBUG_LOC && (nelem == ipHELIUM) && (ipISO==ipHE_LIKE) )
		{
			/* To output Bm(n,l; ipLo), set ipLo, hi_l, and hi_s accordingly.	*/
			long int hi_l,hi_s;
			double Bm;

			/* these must be set for following output to make sense
			 * as is, a dangerous bit of code - set NaN for safety */
			hi_s = -100000;
			hi_l = -100000;
			ipLo = -100000;
			/* tripS to 2^3P	*/
			//hi_l = 0, hi_s = 3, ipLo = ipHe2p3P0;

			/* tripD to 2^3P	*/
			//hi_l = 2, hi_s = 3, ipLo = ipHe2p3P0;

			/* tripP to 2^3S	*/
			//hi_l = 1, hi_s = 3, ipLo = ipHe2s3S;	

			ASSERT( hi_l != StatesElemNEW[nelem][nelem-ipISO][ipLo].l );

			fprintf(ioQQQ,"Bm(n,%ld,%ld;%ld)\n",hi_l,hi_s,ipLo);
			fprintf(ioQQQ,"m\t2\t\t3\t\t4\t\t5\t\t6\n");

			for( ipHi=ipHe2p3P2; ipHi<iso.numLevels_max[ipISO][nelem]-iso.nCollapsed_max[ipISO][nelem]; ipHi++ )
			{
				/* Pick out excitations from metastable 2tripS to ntripP.	*/
				if( (StatesElemNEW[nelem][nelem-ipISO][ipHi].l == 1) && (StatesElemNEW[nelem][nelem-ipISO][ipHi].S == 3) )
				{
					fprintf(ioQQQ,"\n%ld\t",StatesElemNEW[nelem][nelem-ipISO][ipHi].n);
					j = 0;
					Bm = 0;
					for( i = ipLo; i<=ipHi; i++)
					{
						if( (StatesElemNEW[nelem][nelem-ipISO][i].l == hi_l) && (StatesElemNEW[nelem][nelem-ipISO][i].S == hi_s)  )
						{
							if( (ipLo == ipHe2p3P0) && (i > ipHe2p3P2) )
							{
								Bm += iso.CascadeProb[ipISO][nelem][ipHi][i] * ( iso.BranchRatio[ipISO][nelem][i][ipHe2p3P0] + 
									iso.BranchRatio[ipISO][nelem][i][ipHe2p3P1] + iso.BranchRatio[ipISO][nelem][i][ipHe2p3P2] );
							}
							else
								Bm += iso.CascadeProb[ipISO][nelem][ipHi][i] * iso.BranchRatio[ipISO][nelem][i][ipLo];

							if( (i == ipHe2p3P0) || (i == ipHe2p3P1) || (i == ipHe2p3P2) )
							{
								j++;
								if(j == 3)
								{
									fprintf(ioQQQ,"%2.4e\t",Bm);
									Bm = 0;
								}
							}
							else
							{
								fprintf(ioQQQ,"%2.4e\t",Bm);
								Bm = 0;
							}
						}
					}
				}
			}
			fprintf(ioQQQ,"\n\n");
		}
	}

	/******************************************************/
	/***  Lifetimes should increase monotonically with  ***/
	/***  increasing n...Make sure the A's decrease.    ***/
	/******************************************************/
	for( i=2; i < iso.n_HighestResolved_max[ipISO][nelem]; ++i)
	{
		ASSERT( (SumAPerN[i] > SumAPerN[i+1]) || opac.lgCaseB );
	}

	// Now check collapsed.  Don't bother with first two since 
	// these are so strongly coupled to the resolved levels	
	for( i=iso.n_HighestResolved_max[ipISO][nelem]+3; i < iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem] - 1; ++i)
	{
		ASSERT( (SumAPerN[i] > SumAPerN[i+1]) || opac.lgCaseB );
	}

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /* && (ipISO == ipH_LIKE) && (nelem == ipHYDROGEN) */)
		{
			for( i = 2; i<= (iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem]); ++i)
			{
				fprintf(ioQQQ,"n %ld\t lifetime %.4e\n", i, 1./SumAPerN[i]);
			}
		}
	}

	free( SumAPerN );

	return;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void iso_satellite( void )
{
	long i, ipISO, nelem;
	
	DEBUG_ENTRY( "iso_satellite()" );

	for( ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
			{
				for( i=0; i<iso.numLevels_max[ipISO][nelem]; i++ )
				{
					char chLab[5]="    ";

					SatelliteLines[ipISO][nelem][i].Zero();
			
					/* Make approximation that all levels have energy of H-like 2s level */
					/* Lines to 1s2s have roughly energy of parent Ly-alpha.  So lines to 1snL will have an energy
					 * smaller by the difference between nL and 2s energies.  Therefore, the following has
					 * energy of parent Ly-alpha MINUS the difference between daughter level and daughter n=2 level. */
					SatelliteLines[ipISO][nelem][i].WLAng = (realnum)(RYDLAM/
						((iso.xIsoLevNIonRyd[ipISO-1][nelem][0] - iso.xIsoLevNIonRyd[ipISO-1][nelem][1]) -
						 (iso.xIsoLevNIonRyd[ipISO][nelem][1]- iso.xIsoLevNIonRyd[ipISO][nelem][i])) );

					SatelliteLines[ipISO][nelem][i].EnergyWN = 1.e8f / 
						SatelliteLines[ipISO][nelem][i].WLAng;

					SatelliteLines[ipISO][nelem][i].EnergyErg = (realnum)ERG1CM * 
						SatelliteLines[ipISO][nelem][i].EnergyWN;

					SatelliteLines[ipISO][nelem][i].EnergyK = (realnum)T1CM * 
						SatelliteLines[ipISO][nelem][i].EnergyWN;

					/* generate label for this ion */
					sprintf( chLab, "%2s%2ld",elementnames.chElementSym[nelem], nelem+1-ipISO );

					SatelliteLines[ipISO][nelem][i].Emis->iRedisFun = ipCRDW;
					/* this is not the usual nelem, is it atomic not C scale. */
					SatelliteLines[ipISO][nelem][i].Hi->nelem = nelem + 1;
					SatelliteLines[ipISO][nelem][i].Hi->IonStg = nelem + 1 - ipISO;
					fixit(); /* what should the stat weight of the upper level be? For now say 2. */
					SatelliteLines[ipISO][nelem][i].Hi->g = 2.f;
					SatelliteLines[ipISO][nelem][i].Lo->g = StatesElemNEW[nelem][nelem-ipISO][i].g;
					SatelliteLines[ipISO][nelem][i].Emis->PopOpc = 
						SatelliteLines[ipISO][nelem][i].Lo->Pop;

					SatelliteLines[ipISO][nelem][i].Emis->pump = 0.;

				}
			}
		}
	}

	return;
}

void iso_satellite_update( void )
{
	double ConBoltz, LTE_pop=SMALLFLOAT+FLT_EPSILON, factor1, ConvLTEPOP;
	long i, ipISO, nelem;
	double dr_rate;
	
	DEBUG_ENTRY( "iso_satellite()" );

	for( ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
			{
				for( i=0; i<iso.numLevels_max[ipISO][nelem]; i++ )
				{
					dr_rate = MAX2( iso.SmallA, iso.DielecRecomb[ipISO][nelem][i] * iso.lgDielRecom[ipISO] );

					SatelliteLines[ipISO][nelem][i].Emis->phots = 
						dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO];
					
					SatelliteLines[ipISO][nelem][i].Emis->xIntensity = 
						SatelliteLines[ipISO][nelem][i].Emis->phots * 
						ERG1CM * SatelliteLines[ipISO][nelem][i].EnergyWN;

					/* We set line intensity above using a rate, but here we need a transition probability.  
					 * We can obtain this by dividing dr_rate by the population of the autoionizing level.
					 * We assume this level is in statistical equilibrium. */
					factor1 = HION_LTE_POP*dense.AtomicWeight[nelem]/
						(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);

					/* term in () is stat weight of electron * ion */
					ConvLTEPOP = pow(factor1,1.5)/(2.*iso.stat_ion[ipISO])/phycon.te32;

					/* This Boltzmann factor is exp( +ioniz energy / Te ).  For simplicity, we make
					 * the fair approximation that all of the autoionizing levels have an energy
					 * equal to the parents n=2. */
					ConBoltz = dsexp(iso.xIsoLevNIonRyd[ipISO-1][nelem][1]/phycon.te_ryd);

					if( ConBoltz >= SMALLDOUBLE )
					{
						/* The energy used to calculate ConBoltz above
						 * should be negative since this is above the continuum, but 
						 * to be safe we calculate ConBoltz with a positive energy above
						 * and multiply by it here instead of dividing.  */
						LTE_pop = SatelliteLines[ipISO][nelem][i].Hi->g * ConBoltz * ConvLTEPOP;
					}
					
					LTE_pop = max( LTE_pop, 1e-30f );

					/* Now the transition probability is simply dr_rate/LTE_pop. */
					SatelliteLines[ipISO][nelem][i].Emis->Aul = (realnum)(dr_rate/LTE_pop);
					SatelliteLines[ipISO][nelem][i].Emis->Aul = 
						max( iso.SmallA, SatelliteLines[ipISO][nelem][i].Emis->Aul );

					SatelliteLines[ipISO][nelem][i].Emis->gf = (realnum)GetGF(
						SatelliteLines[ipISO][nelem][i].Emis->Aul,
						SatelliteLines[ipISO][nelem][i].EnergyWN,
						SatelliteLines[ipISO][nelem][i].Hi->g);

					SatelliteLines[ipISO][nelem][i].Emis->gf = 
						max( 1e-20f, SatelliteLines[ipISO][nelem][i].Emis->gf );

					SatelliteLines[ipISO][nelem][i].Hi->Pop = LTE_pop * dense.xIonDense[nelem][nelem+1-ipISO] * dense.eden;

					SatelliteLines[ipISO][nelem][i].Emis->PopOpc = 
						SatelliteLines[ipISO][nelem][i].Lo->Pop - 
						SatelliteLines[ipISO][nelem][i].Hi->Pop * 
						SatelliteLines[ipISO][nelem][i].Lo->g/SatelliteLines[ipISO][nelem][i].Hi->g;

					SatelliteLines[ipISO][nelem][i].Emis->opacity = 
						(realnum)(abscf(SatelliteLines[ipISO][nelem][i].Emis->gf,
						SatelliteLines[ipISO][nelem][i].EnergyWN,
						SatelliteLines[ipISO][nelem][i].Lo->g));

					/* a typical transition probability is of order 1e10 s-1 */
					double lifetime = 1e-10;

					SatelliteLines[ipISO][nelem][i].Emis->dampXvel = (realnum)( 
							(1.f/lifetime)/PI4/SatelliteLines[ipISO][nelem][i].EnergyWN);
				}
			}
		}
	}

	return;
}

long iso_get_total_num_levels( long ipISO, long nmaxResolved, long numCollapsed )
{
	DEBUG_ENTRY( "iso_get_total_num_levels()" );

	long tot_num_levels;

	/* return the number of levels up to and including nmaxResolved PLUS 
	 * the number (numCollapsed) of collapsed n-levels */		

	if( ipISO == ipH_LIKE )
	{
		tot_num_levels = (long)( nmaxResolved * 0.5 *( nmaxResolved + 1 ) ) + numCollapsed;
	}
	else if( ipISO == ipHE_LIKE )
	{
		tot_num_levels = nmaxResolved*nmaxResolved + nmaxResolved + 1 + numCollapsed;
	}
	else
		TotalInsanity();

	return tot_num_levels;
}

void iso_update_num_levels( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_update_num_levels()" );

	/* This is the minimum resolved nmax. */
	ASSERT( iso.n_HighestResolved_max[ipISO][nelem] >= 3 );

	iso.numLevels_max[ipISO][nelem] = 
		iso_get_total_num_levels( ipISO, iso.n_HighestResolved_max[ipISO][nelem], iso.nCollapsed_max[ipISO][nelem] );

	if( iso.numLevels_max[ipISO][nelem] > iso.numLevels_malloc[ipISO][nelem] )
	{
		fprintf( ioQQQ, "The number of levels for ipISO %li, nelem %li, has been increased since the initial coreload.\n",
			ipISO, nelem );
		fprintf( ioQQQ, "This cannot be done.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* set local copies to the max values */
	iso.numLevels_local[ipISO][nelem] = iso.numLevels_max[ipISO][nelem];
	iso.nCollapsed_local[ipISO][nelem] = iso.nCollapsed_max[ipISO][nelem];
	iso.n_HighestResolved_local[ipISO][nelem] = iso.n_HighestResolved_max[ipISO][nelem];

	/* find the largest number of levels in any element in all iso sequences
	 * we will allocate one matrix for ionization solver, and just use a piece of that memory
	 * for smaller models. */
	max_num_levels = MAX2( max_num_levels, iso.numLevels_max[ipISO][nelem]);

	return;
}

void iso_collapsed_bnl_set( long ipISO, long nelem )
{

	DEBUG_ENTRY( "iso_collapsed_bnl_set()" );

	double bnl_array[4][3][4][10] = {
		{
			/* H */
			{
				{6.13E-01,	2.56E-01,	1.51E-01,	2.74E-01,	3.98E-01,	4.98E-01,	5.71E-01,	6.33E-01,	7.28E-01,	9.59E-01},
				{1.31E+00,	5.17E-01,	2.76E-01,	4.47E-01,	5.87E-01,	6.82E-01,	7.44E-01,	8.05E-01,	9.30E-01,	1.27E+00},
				{1.94E+00,	7.32E-01,	3.63E-01,	5.48E-01,	6.83E-01,	7.66E-01,	8.19E-01,	8.80E-01,	1.02E+00,	1.43E+00},
				{2.53E+00,	9.15E-01,	4.28E-01,	6.16E-01,	7.42E-01,	8.13E-01,	8.60E-01,	9.22E-01,	1.08E+00,	1.56E+00}
			},
			{
				{5.63E-01,	2.65E-01,	1.55E-01,	2.76E-01,	3.91E-01,	4.75E-01,	5.24E-01,	5.45E-01,	5.51E-01,	5.53E-01},
				{1.21E+00,	5.30E-01,	2.81E-01,	4.48E-01,	5.80E-01,	6.62E-01,	7.05E-01,	7.24E-01,	7.36E-01,	7.46E-01},
				{1.81E+00,	7.46E-01,	3.68E-01,	5.49E-01,	6.78E-01,	7.51E-01,	7.88E-01,	8.09E-01,	8.26E-01,	8.43E-01},
				{2.38E+00,	9.27E-01,	4.33E-01,	6.17E-01,	7.38E-01,	8.05E-01,	8.40E-01,	8.65E-01,	8.92E-01,	9.22E-01}
			},
			{
				{2.97E-01,	2.76E-01,	2.41E-01,	3.04E-01,	3.66E-01,	4.10E-01,	4.35E-01,	4.48E-01,	4.52E-01,	4.53E-01},
				{5.63E-01,	5.04E-01,	3.92E-01,	4.67E-01,	5.39E-01,	5.85E-01,	6.10E-01,	6.20E-01,	6.23E-01,	6.23E-01},
				{7.93E-01,	6.90E-01,	4.94E-01,	5.65E-01,	6.36E-01,	6.79E-01,	7.00E-01,	7.09E-01,	7.11E-01,	7.11E-01},
				{1.04E+00,	8.66E-01,	5.62E-01,	6.31E-01,	7.01E-01,	7.43E-01,	7.63E-01,	7.71E-01,	7.73E-01,	7.73E-01}
			}
		},
		{
			/* He+ */
			{
				{6.70E-02,	2.93E-02,	1.94E-02,	4.20E-02,	7.40E-02,	1.12E-01,	1.51E-01,	1.86E-01,	2.26E-01,	3.84E-01},
				{2.39E-01,	1.03E-01,	6.52E-02,	1.31E-01,	2.11E-01,	2.91E-01,	3.61E-01,	4.17E-01,	4.85E-01,	8.00E-01},
				{4.26E-01,	1.80E-01,	1.10E-01,	2.09E-01,	3.18E-01,	4.15E-01,	4.93E-01,	5.54E-01,	6.34E-01,	1.04E+00},
				{6.11E-01,	2.55E-01,	1.51E-01,	2.74E-01,	3.99E-01,	5.02E-01,	5.80E-01,	6.41E-01,	7.30E-01,	1.21E+00}
			},
			{
				{6.79E-02,	3.00E-02,	2.00E-02,	4.30E-02,	7.48E-02,	1.11E-01,	1.44E-01,	1.70E-01,	1.87E-01,	1.96E-01},
				{2.40E-01,	1.04E-01,	6.62E-02,	1.32E-01,	2.11E-01,	2.87E-01,	3.51E-01,	3.98E-01,	4.32E-01,	4.58E-01},
				{4.26E-01,	1.81E-01,	1.11E-01,	2.10E-01,	3.17E-01,	4.12E-01,	4.89E-01,	5.53E-01,	6.14E-01,	6.84E-01},
				{6.12E-01,	2.55E-01,	1.51E-01,	2.73E-01,	3.97E-01,	4.98E-01,	5.77E-01,	6.51E-01,	7.82E-01,	1.18E+00}
			},
			{
				{4.98E-02,	3.47E-02,	2.31E-02,	4.54E-02,	7.14E-02,	9.37E-02,	1.08E-01,	1.13E-01,	1.13E-01,	1.11E-01},
				{1.75E-01,	1.16E-01,	7.36E-02,	1.36E-01,	2.01E-01,	2.50E-01,	2.76E-01,	2.84E-01,	2.81E-01,	2.77E-01},
				{3.38E-01,	1.97E-01,	1.18E-01,	2.13E-01,	3.06E-01,	3.72E-01,	4.06E-01,	4.15E-01,	4.10E-01,	4.04E-01},
				{6.01E-01,	2.60E-01,	1.53E-01,	2.76E-01,	3.95E-01,	4.87E-01,	5.45E-01,	5.76E-01,	5.93E-01,	6.05E-01}
			}
		},
		{
			/* He singlets */
			{
				{1.77E-01,	3.59E-01,	1.54E-01,	2.75E-01,	3.98E-01,	4.94E-01,	5.51E-01,	5.68E-01,	5.46E-01,	4.97E-01},
				{4.09E-01,	7.23E-01,	2.83E-01,	4.48E-01,	5.89E-01,	6.78E-01,	7.22E-01,	7.30E-01,	7.07E-01,	6.65E-01},
				{6.40E-01,	1.02E+00,	3.74E-01,	5.49E-01,	6.85E-01,	7.63E-01,	7.98E-01,	8.03E-01,	7.84E-01,	7.53E-01},
				{8.70E-01,	1.28E+00,	4.42E-01,	6.17E-01,	7.44E-01,	8.13E-01,	8.42E-01,	8.46E-01,	8.34E-01,	8.13E-01}
			},
			{
				{1.78E-01,	3.62E-01,	1.55E-01,	2.73E-01,	3.91E-01,	4.73E-01,	5.10E-01,	5.04E-01,	4.70E-01,	4.32E-01},
				{4.08E-01,	7.26E-01,	2.83E-01,	4.45E-01,	5.79E-01,	6.54E-01,	6.78E-01,	6.64E-01,	6.30E-01,	5.98E-01},
				{6.37E-01,	1.03E+00,	3.73E-01,	5.46E-01,	6.75E-01,	7.40E-01,	7.57E-01,	7.43E-01,	7.15E-01,	6.92E-01},
				{8.65E-01,	1.28E+00,	4.41E-01,	6.14E-01,	7.35E-01,	7.92E-01,	8.05E-01,	7.95E-01,	7.74E-01,	7.59E-01}
			},
			{
				{2.07E-01,	3.73E-01,	1.73E-01,	2.85E-01,	4.03E-01,	4.76E-01,	5.06E-01,	5.03E-01,	4.84E-01,	4.63E-01},
				{4.32E-01,	7.13E-01,	3.06E-01,	4.54E-01,	5.81E-01,	6.44E-01,	6.59E-01,	6.49E-01,	6.28E-01,	6.11E-01},
				{6.40E-01,	9.85E-01,	3.98E-01,	5.53E-01,	6.74E-01,	7.27E-01,	7.36E-01,	7.26E-01,	7.10E-01,	6.98E-01},
				{8.38E-01,	1.21E+00,	4.67E-01,	6.20E-01,	7.34E-01,	7.79E-01,	7.87E-01,	7.79E-01,	7.69E-01,	7.63E-01}
			}
		},
		{
			/* He triplets */
			{
				{9.31E-02,	3.96E-01,	1.36E-01,	2.74E-01,	3.99E-01,	4.95E-01,	5.52E-01,	5.70E-01,	5.48E-01,	4.96E-01},
				{2.25E-01,	8.46E-01,	2.49E-01,	4.46E-01,	5.89E-01,	6.79E-01,	7.23E-01,	7.31E-01,	7.08E-01,	6.64E-01},
				{3.59E-01,	1.24E+00,	3.30E-01,	5.47E-01,	6.85E-01,	7.63E-01,	7.98E-01,	8.04E-01,	7.85E-01,	7.53E-01},
				{4.93E-01,	1.60E+00,	3.91E-01,	6.15E-01,	7.44E-01,	8.13E-01,	8.42E-01,	8.47E-01,	8.35E-01,	8.12E-01}
			},
			{
				{9.32E-02,	3.99E-01,	1.35E-01,	2.72E-01,	3.91E-01,	4.75E-01,	5.14E-01,	5.09E-01,	4.73E-01,	4.31E-01},
				{2.25E-01,	8.49E-01,	2.49E-01,	4.44E-01,	5.79E-01,	6.56E-01,	6.81E-01,	6.68E-01,	6.31E-01,	5.96E-01},
				{3.58E-01,	1.25E+00,	3.29E-01,	5.44E-01,	6.76E-01,	7.42E-01,	7.60E-01,	7.46E-01,	7.16E-01,	6.91E-01},
				{4.92E-01,	1.60E+00,	3.90E-01,	6.12E-01,	7.36E-01,	7.93E-01,	8.07E-01,	7.97E-01,	7.74E-01,	7.58E-01}
			},
			{
				{1.13E-01,	4.21E-01,	1.47E-01,	2.83E-01,	4.04E-01,	4.80E-01,	5.13E-01,	5.12E-01,	4.93E-01,	4.71E-01},
				{2.52E-01,	8.56E-01,	2.61E-01,	4.50E-01,	5.82E-01,	6.48E-01,	6.66E-01,	6.56E-01,	6.35E-01,	6.16E-01},
				{3.85E-01,	1.23E+00,	3.41E-01,	5.49E-01,	6.75E-01,	7.30E-01,	7.41E-01,	7.31E-01,	7.15E-01,	7.02E-01},
				{5.14E-01,	1.56E+00,	4.01E-01,	6.15E-01,	7.34E-01,	7.82E-01,	7.90E-01,	7.83E-01,	7.72E-01,	7.65E-01}
			}
		}
	};

	double temps[4] = {5000., 10000., 15000., 20000. };
	double log_dens[3] = {2., 4., 6.};
	long ipTe, ipDens;

	ASSERT( nelem <= 1 );

	/* find temperature in tabulated values.  */
	ipTe = hunt_bisect( temps, 4, phycon.te );			
	ipDens = hunt_bisect( log_dens, 3, log10(dense.eden) );			

	ASSERT( (ipTe >=0) && (ipTe < 3)  );
	ASSERT( (ipDens >=0) && (ipDens < 2)  );

	for( long nHi=iso.n_HighestResolved_max[ipISO][nelem]+1; nHi<=iso.n_HighestResolved_max[ipISO][nelem]+iso.nCollapsed_max[ipISO][nelem]; nHi++ )
	{
		for( long lHi=0; lHi<nHi; lHi++ )
		{
			for( long sHi=1; sHi<4; sHi++ )
			{
				if( ipISO == ipH_LIKE && sHi != 2 )
					continue;
				else if( ipISO == ipHE_LIKE && sHi != 1 && sHi != 3 )
					continue;

				double bnl_at_lo_den, bnl_at_hi_den, bnl;
				double bnl_max, bnl_min, temp, dens;

				long ipL = MIN2(9,lHi);
				long ip1;

				if( nelem==ipHYDROGEN )
					ip1 = 0;
				else if( nelem==ipHELIUM )
				{
					if( ipISO==ipH_LIKE )
						ip1 = 1;
					else if( ipISO==ipHE_LIKE )
					{
						if( sHi==1 )
							ip1 = 2;
						else if( sHi==3 )
							ip1 = 3;
						else
							TotalInsanity();
					}
					else
						TotalInsanity();
				}
				else
					TotalInsanity();

				temp = MAX2( temps[0], phycon.te );
				temp = MIN2( temps[3], temp );

				dens = MAX2( log_dens[0], log10(dense.eden) );
				dens = MIN2( log_dens[2], dens );

				/* Calculate the answer...must interpolate on two variables.
				 * First interpolate on T, at both the lower and upper densities.
				 * Then interpolate between these results for the right density.	*/
			
				if( temp < temps[0] && dens < log_dens[0] )
					bnl = bnl_array[ip1][0][0][ipL];
				else if( temp < temps[0] && dens >= log_dens[2] )
					bnl = bnl_array[ip1][2][0][ipL]; 
				else if( temp >= temps[3] && dens < log_dens[0] )
					bnl = bnl_array[ip1][0][3][ipL]; 
				else if( temp >= temps[3] && dens >= log_dens[2] )
					bnl = bnl_array[ip1][2][3][ipL];
				else
				{
					bnl_at_lo_den = ( temp - temps[ipTe]) / (temps[ipTe+1] - temps[ipTe]) *
						(bnl_array[ip1][ipDens][ipTe+1][ipL] - bnl_array[ip1][ipDens][ipTe][ipL]) + bnl_array[ip1][ipDens][ipTe][ipL];

					bnl_at_hi_den = ( temp - temps[ipTe]) / (temps[ipTe+1] - temps[ipTe]) *
						(bnl_array[ip1][ipDens+1][ipTe+1][ipL] - bnl_array[ip1][ipDens+1][ipTe][ipL]) + bnl_array[ip1][ipDens+1][ipTe][ipL];

					bnl = ( dens - log_dens[ipDens]) / (log_dens[ipDens+1] - log_dens[ipDens]) * 
						(bnl_at_hi_den - bnl_at_lo_den) + bnl_at_lo_den;
				}

				/** these are just sanity checks, the interpolated value should be between values at interpolation points */
				{
					bnl_max = MAX4( bnl_array[ip1][ipDens][ipTe+1][ipL], bnl_array[ip1][ipDens+1][ipTe+1][ipL],
						bnl_array[ip1][ipDens][ipTe][ipL], bnl_array[ip1][ipDens+1][ipTe][ipL] );
					ASSERT( bnl <= bnl_max );

					bnl_min = MIN4( bnl_array[ip1][ipDens][ipTe+1][ipL], bnl_array[ip1][ipDens+1][ipTe+1][ipL],
						bnl_array[ip1][ipDens][ipTe][ipL], bnl_array[ip1][ipDens+1][ipTe][ipL] );
					ASSERT( bnl >= bnl_min );
				}

				iso.bnl_effective[ipISO][nelem][nHi][lHi][sHi] = bnl;

				ASSERT( iso.bnl_effective[ipISO][nelem][nHi][lHi][sHi]  > 0. );
			}
		}
	}

	return;
}


void iso_collapsed_bnl_print( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collapsed_bnl_print()" );

	for( long is = 1; is<=3; ++is)
	{
		if( ipISO == ipH_LIKE && is != 2 )
			continue;
		else if( ipISO == ipHE_LIKE && is != 1 && is != 3 )
			continue;

		char chSpin[3][9]= {"singlets", "doublets", "triplets"};

		/* give element number and spin */
		fprintf(ioQQQ," %s %s  %s bnl\n",
			iso.chISO[ipISO],
			elementnames.chElementSym[nelem],
			chSpin[is-1]);

		/* header with the l states */
		fprintf(ioQQQ," n\\l=>    ");
		for( long i =0; i < iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem]; ++i)
		{
			fprintf(ioQQQ,"%2ld         ",i);
		}
		fprintf(ioQQQ,"\n");

		/* loop over prin quant numbers, one per line, with l across */
		for( long in = 1; in <= iso.n_HighestResolved_max[ipISO][nelem] + iso.nCollapsed_max[ipISO][nelem]; ++in)
		{
			if( is==3 && in==1 )
				continue;

			fprintf(ioQQQ," %2ld      ",in);

			for( long il = 0; il < in; ++il)
			{
				fprintf( ioQQQ, "%9.3e ", iso.bnl_effective[ipISO][nelem][in][il][is] );
			}
			fprintf(ioQQQ,"\n");
		}
	}

	return;
}

void iso_collapsed_Aul_update( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collapsed_Aul_update()" );

	long ipFirstCollapsed = iso.numLevels_max[ipISO][nelem] - iso.nCollapsed_max[ipISO][nelem];

	for( long ipLo=0; ipLo<ipFirstCollapsed; ipLo++ )
	{
		long spin = StatesElemNEW[nelem][nelem-ipISO][ipLo].S;

		/* calculate effective Aul's from collapsed levels */
		for( long nHi=iso.n_HighestResolved_max[ipISO][nelem]+1; nHi<=iso.n_HighestResolved_max[ipISO][nelem]+iso.nCollapsed_max[ipISO][nelem]; nHi++ )
		{
			realnum Auls[2] = {
				iso.CachedAs[ipISO][nelem][ nHi-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][0],
				iso.CachedAs[ipISO][nelem][ nHi-iso.n_HighestResolved_max[ipISO][nelem]-1 ][ ipLo ][1] };

			realnum EffectiveAul = 
				Auls[0]*spin*(2.f*(L_(ipLo)+1.f)+1.f)*(realnum)iso.bnl_effective[ipISO][nelem][nHi][ L_(ipLo)+1 ][spin];
			
			/* this is for n,L-1 -> n',L
			 * make sure L-1 exists. */
			if( L_(ipLo) > 0 )
			{
				EffectiveAul += 
					Auls[1]*spin*(2.f*(L_(ipLo)-1.f)+1.f)*(realnum)iso.bnl_effective[ipISO][nelem][nHi][ L_(ipLo)-1 ][spin];
			}

			if( ipISO==ipH_LIKE )
				EffectiveAul /= (2.f*nHi*nHi);
			else if( ipISO==ipHE_LIKE )
				EffectiveAul /= (4.f*nHi*nHi);
			else
				TotalInsanity();

			long ipHi = iso.QuantumNumbers2Index[ipISO][nelem][nHi][ L_(ipLo)+1 ][spin];

			/* FINALLY, put the effective A in the proper Emis structure. */
			Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul = EffectiveAul;

			ASSERT( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul > 0. );
		}
	}

	return;
}

void iso_collapsed_lifetimes_update( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collapsed_Aul_update()" );

	for( long ipHi=iso.numLevels_max[ipISO][nelem]- iso.nCollapsed_max[ipISO][nelem]; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
	{
		StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime = SMALLFLOAT;

		for( long ipLo=0; ipLo < ipHi; ipLo++ )  
		{
			if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
				continue;

			StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime += Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul;
		}

		/* sum of A's was just stuffed, now invert for lifetime. */
		StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime = 1./StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime;

		for( long ipLo=0; ipLo < ipHi; ipLo++ )
		{
			if( Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN <= 0. )
				continue;

			if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
				continue;

			Transitions[ipISO][nelem][ipHi][ipLo].Emis->dampXvel = (realnum)( 
				(1.f/StatesElemNEW[nelem][nelem-ipISO][ipHi].lifetime)/
				PI4/Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN);

			ASSERT(Transitions[ipISO][nelem][ipHi][ipLo].Emis->dampXvel> 0.);
		}
	}

	return;
}

#if	0
STATIC void Prt_AGN_table( void )
{
	/* the designation of the levels, chLevel[n][string] */
	multi_arr<char,2> chLevel(max_num_levels,10);

	/* create spectroscopic designation of labels */
	for( long ipLo=0; ipLo < iso.numLevels_max[ipISO][ipISO]-iso.nCollapsed_max[ipISO][ipISO]; ++ipLo )
	{
		long nelem = ipISO;
		sprintf( &chLevel.front(ipLo) , "%li %li%c", N_(ipLo), S_(ipLo), chL[MIN2(20,L_(ipLo))] );
	}

	/* option to print cs data for AGN */
	/* create spectroscopic designation of labels */
	{
		/* option to print particulars of some line when called */
		enum {AGN=false};
		if( AGN )
		{
#			define NTEMP 6
			double te[NTEMP]={6000.,8000.,10000.,15000.,20000.,25000. };
			double telog[NTEMP] ,
				cs ,
				ratecoef;
			long nelem = ipHELIUM;
			fprintf( ioQQQ,"trans");
			for( long i=0; i < NTEMP; ++i )
			{
				telog[i] = log10( te[i] );
				fprintf( ioQQQ,"\t%.3e",te[i]);
			}
			for( long i=0; i < NTEMP; ++i )
			{
				fprintf( ioQQQ,"\t%.3e",te[i]);
			}
			fprintf(ioQQQ,"\n");

			for( long ipHi=ipHe2s3S; ipHi< iso.numLevels_max[ipHE_LIKE][ipHELIUM]; ++ipHi )
			{
				for( long ipLo=ipHe1s1S; ipLo < ipHi; ++ipLo )
				{

					/* deltaN = 0 transitions may be wrong because 
					 * COLL_CONST below is only correct for electron colliders */
					if( N_(ipHi) == N_(ipLo) )
						continue;

					/* print the designations of the lower and upper levels */
					fprintf( ioQQQ,"%s - %s",
						 &chLevel.front(ipLo) , &chLevel.front(ipHi) );

					/* print the interpolated collision strengths */
					for( long i=0; i < NTEMP; ++i )
					{
						phycon.alogte = telog[i];
						/* print cs */
						cs = HeCSInterp( nelem , ipHi , ipLo, ipELECTRON );
						fprintf(ioQQQ,"\t%.2e", cs );
					}

					/* print the rate coefficients */
					for( long i=0; i < NTEMP; ++i )
					{
						phycon.alogte = telog[i];
						phycon.te = pow(10.,telog[i] );
						tfidle(false);
						cs = HeCSInterp( nelem , ipHi , ipLo, ipELECTRON );
						/* collisional deexcitation rate */
						ratecoef = cs/sqrt(phycon.te)*COLL_CONST/StatesElemNEW[nelem][nelem-ipHE_LIKE][ipLo].g *
							sexp( Transitions[ipHE_LIKE][nelem][ipHi][ipLo].EnergyK / phycon.te );
						fprintf(ioQQQ,"\t%.2e", ratecoef );
					}
					fprintf(ioQQQ,"\n");
				}
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	return;
}
#endif
