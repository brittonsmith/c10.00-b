/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ContCreatePointers set up pointers for lines and continua called by cloudy after input read in 
 * and continuum mesh has been set */
/*fiddle adjust energy bounds of certain cells so that match ionization edges exactly */
/*ipShells assign continuum energy pointers to shells for all atoms,
 * called by ContCreatePointers */
/*LimitSh sets upper energy limit to subshell integrations */
/*ContBandsCreate - read set of continuum bands to enter total emission into line stack*/
#include "cddefines.h"
#include "physconst.h"
#include "lines_service.h"
#include "iso.h"
#include "secondaries.h"
#include "taulines.h"
#include "elementnames.h"
#include "ionbal.h"
#include "rt.h"
#include "opacity.h"
#include "yield.h"
#include "dense.h"
#include "he.h"
#include "fe.h"
#include "rfield.h"
#include "oxy.h"
#include "atomfeii.h"
#include "atoms.h"
#include "trace.h"
#include "hmi.h"
#include "heavy.h"
#include "atmdat.h"
#include "ipoint.h"
#include "h2.h"
#include "continuum.h"

/*LimitSh sets upper energy limit to subshell integrations */
STATIC long LimitSh(long int ion, 
  long int nshell, 
  long int nelem);

STATIC void ipShells(
	/* nelem is the atomic number on the C scale, Li is 2 */
	long int nelem);

/*ContBandsCreate - read set of continuum bands to enter total emission into line*/
STATIC void ContBandsCreate(
	/* chFile is optional filename, if void then use default bands,
	 * if not void then use file specified,
	 * return value is 0 for success, 1 for failure */
	 const char chFile[] );

/* upper level for two-photon emission in H and He iso sequences */
#define TwoS	(1+ipISO)

/*fiddle adjust energy bounds of certain cells so that match ionization edges exactly */
STATIC void fiddle(long int ipnt, 
  double exact);

void ContCreatePointers(void)
{
	char chLab[5];
	long int 
	  i, 
	  ion, 
	  ipHi, 
	  ipLo, 
	  ipISO,
	  iWL_Ang,
	  j, 
	  nelem,
	  nshells;
	double energy,
		xnew;
	/* counter to say whether pointers have ever been evaluated */
	static int nCalled = 0;

	DEBUG_ENTRY( "ContCreatePointers()" );

	/* create the hydrogen atom for this core load, routine creates space then zeros it out
	 * on first call, on second and later calls it only zeros things out */
	iso_create();

	/* create internal static variables needed to do the H2 molecule */
	H2_Create();

	/* nCalled is local static variable defined 0 when defined. 
	 * it is incremented below, so that space only allocated one time per coreload. */
	if( nCalled > 0 )
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, " ContCreatePointers called, not evaluating.\n" );
		return;
	}
	else
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, " ContCreatePointers called first time.\n" );
		++nCalled;
	}

	for( i=0; i < rfield.nupper; i++ )
	{
		/* this is array of labels for lines and continua, set to blanks at first */
		strcpy( rfield.chContLabel[i], "    ");
		strcpy( rfield.chLineLabel[i], "    ");
	}

	/* we will generate a set of array indices to ionization edges for
	 * the first thirty elements.  First set all array indices to
	 * totally bogus values so we will crash if misused */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( ion=0; ion<LIMELM; ++ion )
			{
				for( nshells=0; nshells<7; ++nshells )
				{
					for( j=0; j<3; ++j )
					{
						opac.ipElement[nelem][ion][nshells][j] = INT_MIN;
					}
				}
			}
		}
	}

	/* pointer to excited state of O+2 */
	opac.ipo3exc[0] = ipContEnergy(3.855,"O3ex");

	/* main hydrogenic arrays - THIS OCCURS TWICE!! H and He here, then the
	 * remaining hydrogenic species near the bottom.  This is so that H and HeII get
	 * their labels stuffed into the arrays, and the rest of the hydrogenic series 
	 * get whatever is left over after the level 1 lines.
	 * to find second block, search on "ipZ=2" */
	/* NB note that no test for H or He being on exists here - we will always
	 * define the H and He arrays even when He is off, since we need to
	 * know where the he edges are for the bookkeeping that occurs in continuum
	 * binning routines */
	/* this loop is over H, He-like only - DO NOT CHANGE TO NISO */
	for( ipISO=ipH_LIKE; ipISO<=ipHE_LIKE; ++ipISO )
	{
		/* this will be over HI, HeII, then HeI only */
		for( nelem=ipISO; nelem < 2; nelem++ )
		{
			/* generate label for this ion */
			sprintf( chLab, "%2s%2ld",elementnames.chElementSym[nelem], nelem+1-ipISO );

			/* array index for continuum edges for ground */
			iso.ipIsoLevNIonCon[ipISO][nelem][0] = ipContEnergy(iso.xIsoLevNIonRyd[ipISO][nelem][0],chLab);
			for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				/* array index for continuum edges for excited levels */
				iso.ipIsoLevNIonCon[ipISO][nelem][ipHi] = ipContEnergy(iso.xIsoLevNIonRyd[ipISO][nelem][ipHi],chLab);

				/* define all line array indices */
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
						continue;

					/* some lines have negative or zero energy */
					/* >>chng 03 apr 22, this check was if less than or equal to zero,
					 * changed to lowest energy point so that ultra low energy transitions are
					 * not considered.	*/
					if( Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN * WAVNRYD < continuum.filbnd[0] )
						continue;

					/* some energies are negative for inverted levels */
					Transitions[ipISO][nelem][ipHi][ipLo].ipCont = 
						ipLineEnergy(Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN * WAVNRYD , chLab ,
						iso.ipIsoLevNIonCon[ipISO][nelem][ipLo]);
					Transitions[ipISO][nelem][ipHi][ipLo].Emis->ipFine = 
						ipFineCont(Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN * WAVNRYD );
					/* check that energy scales are the same, to within energy resolution of arrays */
#					ifndef NDEBUG
					if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->ipFine > 0 )
					{
						realnum anuCoarse = rfield.anu[Transitions[ipISO][nelem][ipHi][ipLo].ipCont-1];
						realnum anuFine = rfield.fine_anu[Transitions[ipISO][nelem][ipHi][ipLo].Emis->ipFine];
						realnum widCoarse = rfield.widflx[Transitions[ipISO][nelem][ipHi][ipLo].ipCont-1];
						realnum widFine = anuFine - rfield.fine_anu[Transitions[ipISO][nelem][ipHi][ipLo].Emis->ipFine-1];
						realnum width = MAX2( widFine , widCoarse );
						/* NB - can only assert more than width of coarse cell 
						 * since close to ionization limit, coarse index is 
						 * kept below next ionization edge 
						 * >>chng 05 mar 02, pre factor below had been 1.5, chng to 2
						 * tripped near H grnd photo threshold */
						ASSERT( fabs(anuCoarse - anuFine) / anuCoarse < 
							2.*width/anuCoarse);
					}
#					endif
				}
			}/* ipHi loop */
		}/* nelem loop */
	}/* ipISO */

	/* need to increase the cell for the HeII balmer continuum by one, so that
	 * hydrogen Lyman continuum will pick it up. */
	nelem = ipHELIUM;
	/* copy label over to next higher cell */
	if( strcmp( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s]-1] , "He 2" ) == 0)
	{
		strcpy( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s]], 
				 rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s]-1] );
		/* set previous spot to blank */
		strcpy( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s]-1] , "    ");
	}
	else if( strcmp( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s]-1] , "H  1" ) == 0)
	{
		/* set previous spot to blank */
		strcpy( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s]] , "He 2");
	}
	else
	{
		fprintf(ioQQQ," insanity heii pointer fix contcreatepointers\n");
	}
	/* finally increment the two HeII pointers so that they are above the Lyman continuum */
	++iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2s];
	++iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH2p];

	/* we will define an array of either 1 or 0 to show whether photooelectron
	 * energy is large enough to produce secondary ionizations
	 * below 100eV no secondary ionization - this is the threshold*/
	secondaries.ipSecIon = ipoint(7.353);

	/* this is highest energy where k-shell opacities are counted
	 * can be adjusted with "set kshell" command */
	continuum.KshellLimit = ipoint(continuum.EnergyKshell);

	/* pointers for molecules
	 * H2+ dissociation energy 2.647 eV but cs small below 0.638 Ryd */
	opac.ih2pnt[0] = ipContEnergy(0.502,"H2+ ");
	opac.ih2pnt[1] = ipoint(1.03);

	/* pointers for most prominent PAH features
	 * energies given to ipContEnergy are only to put lave in the right place
	 * wavelengths are rough observed values of blends
	 * 7.6 microns */
	i = ipContEnergy(0.0117, "PAH " );

	/* feature near 6.2 microns */
	i = ipContEnergy(0.0147, "PAH " );

	/* 3.3 microns */
	i = ipContEnergy(0.028, "PAH " );

	/* 11.2 microns */
	i = ipContEnergy(0.0080, "PAH " );

	/* 12.3 microns */
	i = ipContEnergy(0.0077, "PAH " );

	/* 13.5 microns */
	i = ipContEnergy(0.0069, "PAH " );


	/* fix pointers for hydrogen and helium */

	/* pointer to Bowen 374A resonance line */
	he.ip374 = ipLineEnergy(1.92,"He 2",0);
	he.ip660 = ipLineEnergy(1.38,"He 2",0);

	/* pointer to energy defining effect x-ray column */
	rt.ipxry = ipoint(73.5);

	/* pointer to Hminus edge at 0.754eV */
	hmi.iphmin = ipContEnergy(0.05544,"H-  ");

	/* pointer to threshold for H2 photoionization at 15.4 eV */
	opac.ipH2_photo_thresh = ipContEnergy( 15.4/EVRYD , "H2  ");

	hmi.iheh1 = ipoint(1.6);
	hmi.iheh2 = ipoint(2.3);

	/* pointer to carbon k-shell ionization */
	opac.ipCKshell = ipoint(20.6);

	/* pointer to threshold for pair production */
	opac.ippr = ipoint(7.51155e4) + 1;

	/* pointer to x-ray - gamma ray bound; 100 kev */
	rfield.ipEnerGammaRay = ipoint(rfield.EnerGammaRay);

	t_fe2ovr_la::Inst().init_pointers();

	/* make low energy edges of some cells exact,
	 * index is on c scale
	 * 0.99946 is correct energy of hydrogen in inf mass ryd */
	fiddle(iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1,0.99946);
	fiddle(iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH2p]-1,0.24986);
	/* confirm that labels are in correct location */
	ASSERT( strcmp( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1], "H  1" ) ==0 );
	ASSERT( strcmp( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH2p]-1], "H  1" ) ==0 );

	fiddle(iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s]-1,4.00000);
	ASSERT( strcmp( rfield.chContLabel[iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s]-1], "He 2" ) ==0 );

	/* pointer to excited state of O+2 */
	fiddle(opac.ipo3exc[0]-1,3.855);

	/* now fix widflx array so that it is correct */
	for( i=1; i<rfield.nupper-1; ++i )
	{
		rfield.widflx[i] = ((rfield.anu[i+1] - rfield.anu[i]) + (rfield.anu[i] - rfield.anu[i-1]))/2.f;
	}

	/* these are indices for centers of B and V filters,
	 * taken from table on page 202 of Allen, AQ, 3rd ed */
	/* the B filter array offset */
	rfield.ipB_filter = ipoint( RYDLAM / WL_B_FILT );
	/* the V filter array offset */
	rfield.ipV_filter = ipoint( RYDLAM / WL_V_FILT );

	/* these are the lower and upper bounds for the G0 radiation field
	 * used by Tielens & Hollenbach in their PDR work */
	rfield.ipG0_TH85_lo =  ipoint(  6.0 / EVRYD );
	rfield.ipG0_TH85_hi =  ipoint( 13.6 / EVRYD );

	/* this is the limits for Draine & Bertoldi Habing field */
	rfield.ipG0_DB96_lo =  ipoint(  RYDLAM / 1110. );
	rfield.ipG0_DB96_hi =  ipoint( RYDLAM / 912. );

	/* this is special form of G0 that could be used in some future version, for now,
	 * use default, TH85 */
	rfield.ipG0_spec_lo = ipoint(  6.0 / EVRYD );
	rfield.ipG0_spec_hi = ipoint( RYDLAM / 912. );

	/* this index is to 1000A to obtain the extinction at 1000A */
	rfield.ip1000A = ipoint( RYDLAM / 1000. );

	/* now save current form of array */
	for( i=0; i < rfield.nupper; i++ )
	{
		rfield.AnuOrg[i] = rfield.anu[i];
		rfield.anusqr[i] = (realnum)sqrt(rfield.AnuOrg[i]);
	}

	/* following order of elements is in roughly decreasing abundance
	 * when ipShells gets a cell for the valence IP it does so through
	 * ipContEnergy, which makes sure that no two ions get the same cell
	 * earliest elements have most precise ip mapping */

	/* set up shell pointers for hydrogen */
	nelem = ipHYDROGEN;
	ion = 0;

	/* the number of shells */
	Heavy.nsShells[nelem][0] = 1;

	/*pointer to ionization threshold in energy array*/
	Heavy.ipHeavy[nelem][ion] = iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH1s];
	opac.ipElement[nelem][ion][0][0] = iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH1s];

	/* upper limit to energy integration */
	opac.ipElement[nelem][ion][0][1] = rfield.nupper;

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* set up shell pointers for helium */
		nelem = ipHELIUM;
		ion = 0;

		/* the number of shells */
		Heavy.nsShells[nelem][0] = 1;

		/*pointer to ionization threshold in energy array*/
		Heavy.ipHeavy[nelem][ion] = iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0];
		opac.ipElement[nelem][ion][0][0] = iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0];

		/* upper limit to energy integration */
		opac.ipElement[nelem][ion][0][1] = rfield.nupper;

		/* (hydrogenic) helium ion */
		ion = 1;
		/* the number of shells */
		Heavy.nsShells[nelem][1] = 1;

		/*pointer to ionization threshold in energy array*/
		Heavy.ipHeavy[nelem][ion] = iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH1s];
		opac.ipElement[nelem][ion][0][0] = iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH1s];

		/* upper limit to energy integration */
		opac.ipElement[nelem][ion][0][1] = rfield.nupper;
	}

	/* check that ionization potential of neutral carbon valence shell is
	 * positive */
	ASSERT( t_ADfA::Inst().ph1(2,5,5,0) > 0. );

	/* now fill in all sub-shell ionization array indices for elements heavier than He,
	 * this must be done after previous loop on iso.ipIsoLevNIonCon[ipH_LIKE] since hydrogenic species use
	 * iso.ipIsoLevNIonCon[ipH_LIKE] rather than ipoint in getting array index within continuum array */
	for( i=NISO; i<LIMELM; ++i )
	{
		if( dense.lgElmtOn[i])
		{
			/* i is the atomic number on the c scale, 2 for Li */
			ipShells(i);
		}
	}

	/* most of these are set in ipShells, but not h-like or he-like, so do these here */
	Heavy.Valence_IP_Ryd[ipHYDROGEN][0] = t_ADfA::Inst().ph1(0,0,ipHYDROGEN,0)/EVRYD* 0.9998787;
	Heavy.Valence_IP_Ryd[ipHELIUM][0] = t_ADfA::Inst().ph1(0,1,ipHELIUM,0)/EVRYD* 0.9998787;
	Heavy.Valence_IP_Ryd[ipHELIUM][1] = t_ADfA::Inst().ph1(0,0,ipHELIUM,0)/EVRYD* 0.9998787;
	for( nelem=2; nelem<LIMELM; ++nelem )
	{
		Heavy.Valence_IP_Ryd[nelem][nelem-1] = t_ADfA::Inst().ph1(0,1,nelem,0)/EVRYD* 0.9998787;
		Heavy.Valence_IP_Ryd[nelem][nelem] = t_ADfA::Inst().ph1(0,0,nelem,0)/EVRYD* 0.9998787;
		if( dense.lgElmtOn[nelem])
		{
			/* now confirm that all are properly set */
			for( j=0; j<=nelem; ++j )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][j] > 0.05 );
			}
			for( j=0; j<nelem; ++j )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][j] < Heavy.Valence_IP_Ryd[nelem][j+1]);
			}
		}
	}

	/* array indices to bound Compton electron recoil ionization of all elements */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem])
		{
			for( ion=0; ion<nelem+1; ++ion )
			{
				/* this is the threshold energy to Compton ionize valence shell electrons */
				energy = sqrt( Heavy.Valence_IP_Ryd[nelem][ion] * EN1RYD * ELECTRON_MASS * SPEEDLIGHT * SPEEDLIGHT ) / EN1RYD;
				/* the array index for this energy */
				ionbal.ipCompRecoil[nelem][ion] = ipoint( energy );
			}
		}
	}

	/* oxygen pointers for excited states
	 * IO3 is pointer to O++ exc state, is set above */
	oxy.i2d = ipoint(1.242);
	oxy.i2p = ipoint(1.367);
	opac.ipo1exc[0] = ipContEnergy(0.856,"O1ex");
	opac.ipo1exc[1] = ipoint(2.0);

	/* upper limit for excited state photoionization
	 * do not use ipContEnergy since only upper limit */
	opac.ipo3exc[1] = ipoint(5.0);

	/* upper level of 4363 */
	opac.ipo3exc3[0] = ipContEnergy(3.646,"O3ex");
	opac.ipo3exc3[1] = ipoint(5.0);

	/* following are various pointers for OI - Ly-beta pump problem
	 * these are delta energies for Boltzmann factors */
	/** \todo	2	this is redundant with contents of oxygen line arrays
	 * use them instead
	 * when removing this, make sure all line intensity predictions
	 * also go into oi line arrays */
	atoms.ipoiex[0] = ipoint(0.7005);
	atoms.ipoiex[1] = ipoint(0.10791);
	atoms.ipoiex[2] = ipoint(0.06925);
	atoms.ipoiex[3] = ipoint(0.01151);
	atoms.ipoiex[4] = ipoint(0.01999);

	/* >>chng 97 jan 27, move nitrogen after oxygen so that O gets the
	 * most accurate pointers
	 * Nitrogen
	 * in1(1) is thresh for photo from excited state */
	opac.in1[0] = ipContEnergy(0.893,"N1ex");

	/* upper limit */
	opac.in1[1] = ipoint(2.);

	if( (trace.lgTrace && trace.lgConBug) || (trace.lgTrace && trace.lgPointBug) )
	{
		fprintf( ioQQQ, "   ContCreatePointers:%ld energy cells used. N(1R):%4ld N(1.8):%4ld  N(4Ryd):%4ld N(O3)%4ld  N(x-ray):%5ld N(rcoil)%5ld\n", 
		  rfield.nupper, 
		  iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s], 
		  iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][ipH1s], 
		  iso.ipIsoLevNIonCon[ipH_LIKE][ipHELIUM][ipH1s], 
		  opac.ipo3exc[0], 
		  opac.ipCKshell, 
		  ionbal.ipCompRecoil[ipHYDROGEN][0] );


		fprintf( ioQQQ, "   ContCreatePointers: ipEnerGammaRay: %5ld IPPRpari produc%5ld\n", 
		  rfield.ipEnerGammaRay, opac.ippr );

		fprintf( ioQQQ, "   ContCreatePointers: H pointers;" );
		for( i=0; i <= 6; i++ )
		{
			fprintf( ioQQQ, "%4ld%4ld", i, iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][i] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "   ContCreatePointers: Oxy pnters;" );

		for( i=1; i <= 8; i++ )
		{
			fprintf( ioQQQ, "%4ld%4ld", i, Heavy.ipHeavy[ipOXYGEN][i-1] );
		}
		fprintf( ioQQQ, "\n" );
	}

	/* Magnesium
	 * following is energy for phot of MG+ from exc state producing 2798 */
	opac.ipmgex = ipoint(0.779);

	/* lower, upper edges of Ca+ excited term photoionizaiton */
	opac.ica2ex[0] = ipContEnergy(0.72,"Ca2x");
	opac.ica2ex[1] = ipoint(1.);

	/* set up factors and pointers for Fe continuum fluorescence */
	fe.ipfe10 = ipoint(2.605);

	/* following is WL(CM)**2/(8PI) * conv fac for RYD to NU *A21 */
	fe.pfe10 = (realnum)(2.00e-18/rfield.widflx[fe.ipfe10-1]);

	/* this is 353 pump, f=0.032 */
	fe.pfe11a = (realnum)(4.54e-19/rfield.widflx[fe.ipfe10-1]);

	/* this is 341.1 f=0.012 */
	fe.pfe11b = (realnum)(2.75e-19/rfield.widflx[fe.ipfe10-1]);
	fe.pfe14 = (realnum)(1.15e-18/rfield.widflx[fe.ipfe10-1]);

	/* set up energy pointers for line optical depth arrays
	 * this also increments flux, sets other parameters for lines */

	/* level 1 heavy elements line array */
	for( i=1; i <= nLevel1; i++ )
	{
		/* put null terminated line label into chLab using line array*/
		chIonLbl(chLab, &TauLines[i]);

		TauLines[i].ipCont = 
			ipLineEnergy(TauLines[i].EnergyWN * WAVNRYD, chLab ,0);
		TauLines[i].Emis->ipFine = 
			ipFineCont(TauLines[i].EnergyWN * WAVNRYD );
		/* for debugging pointer - index on f not c scale, 
		 * this will find all lines that entered a given cell 
		   if( TauLines[i].ipCont==603 )
			fprintf(ioQQQ,"( level1 %s\n", chLab);*/

		if( TauLines[i].Emis->gf > 0. )
		{
			/* the gf value was entered
			 * derive the A, call to function is gf, wl (A), g_up */
			TauLines[i].Emis->Aul = 
				(realnum)(eina(TauLines[i].Emis->gf,
			  TauLines[i].EnergyWN,TauLines[i].Hi->g));
		}
		else if( TauLines[i].Emis->Aul > 0. )
		{
			/* the Einstein A value was entered
			 * derive the gf, call to function is A, wl (A), g_up */
			TauLines[i].Emis->gf = 
				(realnum)(GetGF(TauLines[i].Emis->Aul,
			  TauLines[i].EnergyWN,TauLines[i].Hi->g));
		}
		else
		{
			fprintf( ioQQQ, " level 1 line does not have valid gf or A\n" );
			fprintf( ioQQQ, " This is ContCreatePointers\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* used to get damping constant */
		TauLines[i].Emis->dampXvel = 
			(realnum)(TauLines[i].Emis->Aul / TauLines[i].EnergyWN/PI4);

		/* derive the abs coefficient, call to function is gf, wl (A), g_low */
		TauLines[i].Emis->opacity = 
			(realnum)(abscf(TauLines[i].Emis->gf,
		  TauLines[i].EnergyWN,TauLines[i].Lo->g));
		/*fprintf(ioQQQ,"TauLinesss\t%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			i,TauLines[i].Emis->opacity , TauLines[i].Emis->gf , TauLines[i].Emis->Aul ,TauLines[i].EnergyWN,TauLines[i].Lo->g);*/

		/* excitation energy of transition in degrees kelvin */
		TauLines[i].EnergyK = 
			(realnum)(T1CM)*TauLines[i].EnergyWN;

		/* energy of photon in ergs */
		TauLines[i].EnergyErg = 
			(realnum)(ERG1CM)*TauLines[i].EnergyWN;

		/*line wavelength must be gt 0 */
		ASSERT( TauLines[i].WLAng > 0 );
	}

	/*Beginning of the dBaseLines*/
	for( i=0; i <linesAdded2; i++)
	{
		dBaseLines[i].dampXvel = (realnum)(dBaseLines[i].Aul/
					dBaseLines[i].tran->EnergyWN/PI4);
		dBaseLines[i].damp = -1000.0;
		/*Put null terminated line label into chLab*/
		strncpy(chLab,dBaseLines[i].tran->Hi->chLabel,4);
		chLab[4]='\0';

		dBaseLines[i].tran->ipCont = 
			ipLineEnergy(dBaseLines[i].tran->EnergyWN * WAVNRYD, chLab ,0);
		dBaseLines[i].ipFine = ipFineCont(dBaseLines[i].tran->EnergyWN * WAVNRYD );
		/* derive the abs coefficient, call to function is gf, wl (A), g_low */
		dBaseLines[i].opacity = 
			(realnum)(abscf(dBaseLines[i].gf,dBaseLines[i].tran->EnergyWN,
				dBaseLines[i].tran->Lo->g));
		/* excitation energy of in degrees kelvin */
		dBaseLines[i].tran->EnergyK = 
			(realnum)(T1CM)*dBaseLines[i].tran->EnergyWN;
		/* energy of photon in ergs */
		dBaseLines[i].tran->EnergyErg = 
			(realnum)(ERG1CM)*dBaseLines[i].tran->EnergyWN ;                                                           
	}
	/*end of the dBaseLines*/

	/* set the ipCont struc element for the H2 molecule */
	H2_ContPoint();

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* do remaining part of the iso sequences */
		for( nelem=2; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem])
			{
				/* generate label for this ion */
				sprintf( chLab, "%2s%2ld",elementnames.chElementSym[nelem], nelem+1-ipISO );
				/* array index for continuum edges */
				iso.ipIsoLevNIonCon[ipISO][nelem][0] = 
					ipContEnergy(iso.xIsoLevNIonRyd[ipISO][nelem][0],chLab);

				for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
				{
					/* array index for continuum edges */
					iso.ipIsoLevNIonCon[ipISO][nelem][ipHi] = ipContEnergy(iso.xIsoLevNIonRyd[ipISO][nelem][ipHi],chLab);

					/* define all line pointers */
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{

						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
							continue;

						/* some lines have negative or zero energy */
						/* >>chng 03 apr 22, this check was if less than or equal to zero,
						 * changed to lowest energy point so that ultra low energy transitions are
						 * not considered.	*/
						if( Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN * WAVNRYD < continuum.filbnd[0] )
							continue;

						Transitions[ipISO][nelem][ipHi][ipLo].ipCont = 
							ipLineEnergy(Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN * WAVNRYD , chLab ,
							iso.ipIsoLevNIonCon[ipISO][nelem][ipLo]);
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->ipFine = 
							ipFineCont(Transitions[ipISO][nelem][ipHi][ipLo].EnergyWN * WAVNRYD );
					}
				}
				iso.ipIsoLevNIonCon[ipISO][nelem][0] = ipContEnergy(iso.xIsoLevNIonRyd[ipISO][nelem][0],chLab);
			}
		}
	}
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* this will be over HI, HeII, then HeI only */
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem])
			{
				ipLo = 0;
				/* these are the extra Lyman lines */
				for( ipHi=2; ipHi < iso.nLyman_malloc[ipISO]; ipHi++ )
				{
					/* some energies are negative for inverted levels */
					/*lint -e772 chLab not initialized */
					ExtraLymanLines[ipISO][nelem][ipHi].ipCont = 
						ipLineEnergy(ExtraLymanLines[ipISO][nelem][ipHi].EnergyWN * WAVNRYD , chLab ,
						iso.ipIsoLevNIonCon[ipISO][nelem][ipLo]);

					ExtraLymanLines[ipISO][nelem][ipHi].Emis->ipFine = 
						ipFineCont(ExtraLymanLines[ipISO][nelem][ipHi].EnergyWN * WAVNRYD );
					/*lint +e772 chLab not initialized */
				}

				if( iso.lgDielRecom[ipISO] )
				{
					ASSERT( ipISO>ipH_LIKE );
					for( ipLo=0; ipLo<iso.numLevels_max[ipISO][nelem]; ipLo++ )
					{
						
						SatelliteLines[ipISO][nelem][ipLo].ipCont = ipLineEnergy(
							SatelliteLines[ipISO][nelem][ipLo].EnergyWN * WAVNRYD , chLab , 
							0);

						SatelliteLines[ipISO][nelem][ipLo].Emis->ipFine =  
							ipFineCont(SatelliteLines[ipISO][nelem][ipLo].EnergyWN * WAVNRYD );
					}
				}
			}
		}
	}

	/* some lines do not exist, the flag indicating this is ipCont == -1 */
	/* for H-like sequence, these are 2p-2s (energies degenerate) and 2s-1s, two photon */
	ipISO = ipHYDROGEN;
	/* do remaining part of the he iso sequence */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem])
		{
			/* for H-like sequence want 2p-2s (energies degenerate) and 2s-1s, two photon */
			Transitions[ipISO][nelem][ipH2p][ipH2s].ipCont = -1;
			//Transitions[ipISO][nelem][ipH2p][ipH2s].Emis->ipFine = -1;
			Transitions[ipISO][nelem][ipH2s][ipH1s].ipCont = -1;
			//Transitions[ipISO][nelem][ipH2s][ipH1s].Emis->ipFine = -1;
		}
	}

	fixit(); /* is this redundant? */
	/* for He-like sequence the majority of the transitions are bogus - A set to special value in this case */
	ipISO = ipHELIUM;
	/* do remaining part of the he iso sequence */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem])
		{
			/* this is the two photon transition in the singlets */
			Transitions[ipISO][nelem][ipHe2s1S][ipHe1s1S].ipCont = -1;
			Transitions[ipISO][nelem][ipHe2s1S][ipHe1s1S].Emis->ipFine = -1;

			for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( Transitions[ipISO][nelem][ipHi][ipLo].ipCont <= 0 )
						continue;

					if( fabs(Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul - iso.SmallA) < SMALLFLOAT )
					{
						/* iso.SmallA is value assigned to bogus transitions */
						Transitions[ipISO][nelem][ipHi][ipLo].ipCont = -1;
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->ipFine = -1;
					}
				}
			}
		}
	}

	/* inner shell transitions */
	for( i=0; i<nUTA; ++i )
	{
		if( UTALines[i].Emis->Aul > 0. )
		{
			UTALines[i].Emis->dampXvel = 
				(realnum)(UTALines[i].Emis->Aul/ UTALines[i].EnergyWN/PI4);

			/* derive the abs coefficient, call to function is gf, wl (A), g_low */
			UTALines[i].Emis->opacity = 
				(realnum)(abscf( UTALines[i].Emis->gf, UTALines[i].EnergyWN, UTALines[i].Lo->g));

			/* excitation energy of transition in degrees kelvin */
			UTALines[i].EnergyK = 
				(realnum)(T1CM*UTALines[i].EnergyWN);

			/* energy of photon in ergs */
			UTALines[i].EnergyErg = 
				(realnum)(ERG1CM*UTALines[i].EnergyWN);

			/* put null terminated line label into chLab using line array*/
			chIonLbl(chLab, &UTALines[i]);

			/* get pointer to energy in continuum mesh */
			UTALines[i].ipCont = ipLineEnergy(UTALines[i].EnergyWN * WAVNRYD , chLab,0 );
			UTALines[i].Emis->ipFine = ipFineCont(UTALines[i].EnergyWN * WAVNRYD  );

			/* find heating per absorption,
			 * first find threshold for this shell in ergs */
			/* ionization threshold in erg */
			double thresh = Heavy.Valence_IP_Ryd[UTALines[i].Hi->nelem-1][UTALines[i].Hi->IonStg-1] *EN1RYD;
			UTALines[i].Coll.heat = (UTALines[i].EnergyErg-thresh);
			ASSERT( UTALines[i].Coll.heat> 0. );
		}
	}

	/* level 2 heavy element lines */
	for( i=0; i < nWindLine; i++ )
	{

		/* derive the A, call to function is gf, wl (A), g_up */
		TauLine2[i].Emis->Aul = 
			(realnum)(eina(TauLine2[i].Emis->gf,
		  TauLine2[i].EnergyWN,TauLine2[i].Hi->g));

		/* coefficient needed for damping constant - units cm s-1 */
		TauLine2[i].Emis->dampXvel = 
			(realnum)(TauLine2[i].Emis->Aul/
		  TauLine2[i].EnergyWN/PI4);

		/* derive the abs coefficient, call to function is gf, wl (A), g_low */
		TauLine2[i].Emis->opacity = 
			(realnum)(abscf(TauLine2[i].Emis->gf,
		  TauLine2[i].EnergyWN,TauLine2[i].Lo->g));

		/* excitation energy of transition in degrees kelvin */
		TauLine2[i].EnergyK = 
			(realnum)(T1CM*TauLine2[i].EnergyWN);

		/* energy of photon in ergs */
		TauLine2[i].EnergyErg = 
			(realnum)(ERG1CM*TauLine2[i].EnergyWN);

		/* put null terminated line label into chLab using line array*/
		chIonLbl(chLab, &TauLine2[i]);

		/* get pointer to energy in continuum mesh */
		TauLine2[i].ipCont = ipLineEnergy(TauLine2[i].EnergyWN * WAVNRYD , chLab,0 );
		TauLine2[i].Emis->ipFine = ipFineCont(TauLine2[i].EnergyWN * WAVNRYD );
		/*if( TauLine2[i].ipCont==751 )
			fprintf(ioQQQ,"( atom_level2 %s\n", chLab);*/
	}

	/* hyperfine structure lines */
	for( i=0; i < nHFLines; i++ )
	{

		ASSERT( HFLines[i].Emis->Aul > 0. );
		/* coefficient needed for damping constant */
		HFLines[i].Emis->dampXvel = 
			(realnum)(HFLines[i].Emis->Aul/
			HFLines[i].EnergyWN/PI4);
		HFLines[i].Emis->damp = 1e-20f;

		/* derive the abs coefficient, call to function is gf, wl (A), g_low */
		HFLines[i].Emis->opacity = 
			(realnum)(abscf(HFLines[i].Emis->gf,
			HFLines[i].EnergyWN,
			HFLines[i].Lo->g));
		/* gf from this and 21 cm do not agree, A for HFS is 10x larger than level1 dat */
		/*fprintf(ioQQQ,"HFLinesss\t%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			i,HFLines[i].Emis->opacity , HFLines[i].Emis->gf , HFLines[i].Emis->Aul , HFLines[i].EnergyWN,HFLines[i].Lo->g);*/

		/* excitation energy of transition in degrees kelvin */
		HFLines[i].EnergyK = 
			(realnum)(T1CM*HFLines[i].EnergyWN);

		/* energy of photon in ergs */
		HFLines[i].EnergyErg = 
			(realnum)(ERG1CM*HFLines[i].EnergyWN);

		/* put null terminated line label into chLab using line array*/
		chIonLbl(chLab, &HFLines[i]);

		/* get pointer to energy in continuum mesh */
		HFLines[i].ipCont = ipLineEnergy(HFLines[i].EnergyWN * WAVNRYD , chLab,0 );
		HFLines[i].Emis->ipFine = ipFineCont(HFLines[i].EnergyWN * WAVNRYD );
	}

	/* Verner's FeII lines - do first so following labels will over write this
	 * only call if large FeII atom is turned on */
	FeIIPoint();

	/* the group of inner shell fluorescent lines */
	for( i=0; i < t_yield::Inst().nlines(); ++i )
	{
		strcpy( chLab , elementnames.chElementSym[t_yield::Inst().nelem(i)] );
		strcat( chLab , elementnames.chIonStage[t_yield::Inst().ion_emit(i)] );

		t_yield::Inst().set_ipoint( i, ipLineEnergy( t_yield::Inst().energy(i) , chLab , 0 ) );
	}

	/* ================================================================================== */
	/*        two photon two-photon 2-nu 2 nu 2 photon 2-photon                           */

	/* for induced two photon emission we need mirror image set
	 * of continuum indices for continuum between Lya and half Lya */
	/* first MALLOC LIMELM dimension of space */
	/* >>chng 02 jun 28, malloc this NISO instead of 2.	*/
	iso.ipSym2nu.reserve( NISO );

	/* now loop over the two iso-sequences with two photon continua */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		iso.ipSym2nu.reserve( ipISO, LIMELM );

		/* set up two photon emission */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				double E2nu = Transitions[ipISO][nelem][TwoS][0].EnergyWN * WAVNRYD;

				/* pointer to the Lya energy */
				iso.ipTwoPhoE[ipISO][nelem] = ipoint(E2nu);
				/* this half-energy is only used to get induced rates in two photon */
				iso.ipHalfTwoPhoE[ipISO][nelem] = ipoint(E2nu / 2.);
				while( rfield.anu[iso.ipTwoPhoE[ipISO][nelem]] > E2nu )
				{
					--iso.ipTwoPhoE[ipISO][nelem];
				}
				while( rfield.anu[iso.ipHalfTwoPhoE[ipISO][nelem]] > E2nu / 2. )
				{
					--iso.ipHalfTwoPhoE[ipISO][nelem];
				}

				iso.ipSym2nu.reserve( ipISO, nelem, iso.ipTwoPhoE[ipISO][nelem] );
			}
		}
	}

	iso.ipSym2nu.alloc();
	iso.As2nu.alloc( iso.ipSym2nu.clone() );

	/* now loop over the two iso-sequences with two photon continua */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* set up two photon emission */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				double E2nu = Transitions[ipISO][nelem][TwoS][0].EnergyWN * WAVNRYD;
				double Aul = Transitions[ipISO][nelem][TwoS][0].Emis->Aul;
				double SumShapeFunction = 0., Renorm= 0.;

				/* >>chng 02 aug 14, change upper limit to full Lya energy */
				for( i=0; i < iso.ipTwoPhoE[ipISO][nelem]; i++ )
				{
					/* energy is symmetric energy, the other side of half E2nu */
					energy = E2nu - rfield.anu[i];
					/* this is needed since mirror image of cell next to two-nu energy
					 * may be slightly negative */
					energy = MAX2( energy, rfield.anu[0] + rfield.widflx[0]/2. );
					/* find index for this symmetric energy */
					iso.ipSym2nu[ipISO][nelem][i] = ipoint(energy);
					while( rfield.anu[iso.ipSym2nu[ipISO][nelem][i]] > E2nu ||
						iso.ipSym2nu[ipISO][nelem][i] >= iso.ipTwoPhoE[ipISO][nelem])
					{
						--iso.ipSym2nu[ipISO][nelem][i];
					}
					ASSERT( iso.ipSym2nu[ipISO][nelem][i] >= 0 );
				}

				/* ipTwoPhoE is the cell holding the 2nu energy itself, and we do not want
				 * to include that in the following sum */
				for( i=0; i<iso.ipTwoPhoE[ipISO][nelem]; i++ )
				{
					double ShapeFunction;

					ASSERT( rfield.anu[i]<=E2nu );

					ShapeFunction = atmdat_2phot_shapefunction( rfield.AnuOrg[i]/E2nu, ipISO, nelem )*rfield.widflx[i]/E2nu;
					SumShapeFunction += ShapeFunction;

					/* >>refer	HI	2nu	Spitzer, L., & Greenstein, J., 1951, ApJ, 114, 407 */
					/* As2nu will add up to the A, so its units are s-1	*/ 
					iso.As2nu[ipISO][nelem][i] = (realnum)( Aul * ShapeFunction );
				}

				/* The spline function in twophoton.c causes a bit of an error in the integral of the
				 * shape function.  So we renormalize the integral to 1.	*/
				Renorm = 1./SumShapeFunction;

				for( i=0; i<iso.ipTwoPhoE[ipISO][nelem]; i++ )
				{
					iso.As2nu[ipISO][nelem][i] *= (realnum)Renorm;
				}

				/* The result should be VERY close to 1.	*/
				ASSERT( fabs( SumShapeFunction*Renorm - 1. ) < 0.00001 );
			}
		}
	}

	{
		/* this is an option to print out one of the two photon continua */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{	
#			define	NCRS	21
			double limit,ener[NCRS]={
			  0.,     0.03738,  0.07506,  0.1124,  0.1498,  0.1875,
			  0.225,  0.263,    0.300,    0.3373,  0.375,   0.4127,
			  0.4500, 0.487,    0.525,    0.5625,  0.6002,  0.6376,
			  0.6749, 0.7126,   0.75};

			nelem = ipHYDROGEN;
			ipISO = ipHYDROGEN;

			limit = iso.ipTwoPhoE[ipISO][nelem];

			for( i=0; i < NCRS; i++ )
			{
				fprintf(ioQQQ,"%.3e\t%.3e\n", ener[i] , 
					atmdat_2phot_shapefunction( ener[i]/0.75, ipISO, nelem ) );
			}

			xnew = 0.;
			/** \todo	2	what are we trying to print here?	*/
			for( i=0; i < limit; i++ )
			{
				double fach = iso.As2nu[ipISO][nelem][i]/2.*rfield.anu2[i]/rfield.widflx[i]*EN1RYD;
				fprintf(ioQQQ,"%.3e\t%.3e\t%.3e\n", 
					rfield.anu[i] , 
					iso.As2nu[ipISO][nelem][i] / rfield.widflx[i] , 
					fach );
				xnew += iso.As2nu[ipISO][nelem][i];
			}
			fprintf(ioQQQ," sum is %.3e\n", xnew );
			cdEXIT(EXIT_FAILURE);
		}
	}

	{
		/* this is an option to print out one of the two photon continua */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{	
			for( i=0; i<11; ++i )
			{
				char chLsav[10];
				TauDummy.WLAng = (realnum)(PI * pow(10.,(double)i));
				strcpy( chLsav, chLineLbl(&TauDummy) );
				fprintf(ioQQQ,"%.2f\t%s\n", TauDummy.WLAng , chLsav );
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* option to print out whole thing with "trace lines" command */
	if( trace.lgTrLine )
	{
		fprintf( ioQQQ, "       WL(Ang)   E(RYD)   IP   gl  gu      gf       A        damp     abs K\n" );
		for( i=1; i <= nLevel1; i++ )
		{
			strcpy( chLab, chLineLbl(&TauLines[i]) );
			iWL_Ang = (long)TauLines[i].WLAng;
			if( iWL_Ang > 1000000 )
			{
				iWL_Ang /= 10000;
			}
			else if( iWL_Ang > 10000 )
			{
				iWL_Ang /= 1000;
			}

			fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
			  chLab, iWL_Ang, RYDLAM/TauLines[i].WLAng, 
			  TauLines[i].ipCont, (long)(TauLines[i].Lo->g), 
			  (long)(TauLines[i].Hi->g), TauLines[i].Emis->gf, 
			  TauLines[i].Emis->Aul, TauLines[i].Emis->dampXvel, 
			  TauLines[i].Emis->opacity );
		}

		/*Atomic Or Molecular lines*/
		for( i=0; i <linesAdded2; i++)
		{
			strcpy( chLab, chLineLbl(dBaseLines[i].tran) );

			iWL_Ang = (long)dBaseLines[i].tran->WLAng;

			if( iWL_Ang > 1000000 )
			{
				iWL_Ang /= 10000;
			}
			else if( iWL_Ang > 10000 )
			{
				iWL_Ang /= 1000;
			}
			fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
			  chLab, iWL_Ang, RYDLAM/dBaseLines[i].tran->WLAng, 
			  dBaseLines[i].tran->ipCont, (long)(dBaseLines[i].tran->Lo->g), 
			  (long)(dBaseLines[i].tran->Hi->g),dBaseLines[i].gf, 
			  dBaseLines[i].Aul,dBaseLines[i].dampXvel, 
			  dBaseLines[i].opacity);
		}

		for( i=0; i < nWindLine; i++ )
		{
			strcpy( chLab, chLineLbl(&TauLine2[i]) );

			iWL_Ang = (long)TauLine2[i].WLAng;

			if( iWL_Ang > 1000000 )
			{
				iWL_Ang /= 10000;
			}
			else if( iWL_Ang > 10000 )
			{
				iWL_Ang /= 1000;
			}
			fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
			  chLab, iWL_Ang, RYDLAM/TauLine2[i].WLAng, 
			  TauLine2[i].ipCont, (long)(TauLine2[i].Lo->g), 
			  (long)(TauLine2[i].Hi->g), TauLine2[i].Emis->gf, 
			  TauLine2[i].Emis->Aul, TauLine2[i].Emis->dampXvel, 
			  TauLine2[i].Emis->opacity );
		}
		for( i=0; i < nHFLines; i++ )
		{
			strcpy( chLab, chLineLbl(&HFLines[i]) );

			iWL_Ang = (long)HFLines[i].WLAng;

			if( iWL_Ang > 1000000 )
			{
				iWL_Ang /= 10000;
			}
			else if( iWL_Ang > 10000 )
			{
				iWL_Ang /= 1000;
			}
			fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
			  chLab, iWL_Ang, RYDLAM/HFLines[i].WLAng, 
			  HFLines[i].ipCont, (long)(HFLines[i].Lo->g), 
			  (long)(HFLines[i].Hi->g), HFLines[i].Emis->gf, 
			  HFLines[i].Emis->Aul, HFLines[i].Emis->dampXvel, 
			  HFLines[i].Emis->opacity );
		}
	}

	/* this is an option to kill fine structure line optical depths */
	if( !rt.lgFstOn )
	{
		for( i=1; i <= nLevel1; i++ )
		{
			if( TauLines[i].EnergyWN < 10000. )
			{
				TauLines[i].Emis->opacity = 0.;
			}
		}

		/*Atomic or Molecular Lines-Humeshkar Nemala*/
		for( i=0; i <linesAdded2; i++)
		{
			if(dBaseLines[i].tran->EnergyWN < 10000. )
			{
				dBaseLines[i].opacity = 0.;
			}
		}
	}

	/* read in continuum bands data set */
	ContBandsCreate( "" );

	/* we're done adding lines and states to the stacks.  
	 * This flag is used to make sure we never add them again in this coreload. */
	lgLinesAdded = true;
	lgStatesAdded = true;

	return;
}

/*fiddle adjust energy bounds of cell with index ipnt so that lower energy
 * matches ionization edges exactly, called by createpoint */
/* ipnt is index on c scale */
STATIC void fiddle(long int ipnt, 
  double exact)
{
	realnum Ehi, 
	  Elo,
	  OldEner;

	DEBUG_ENTRY( "fiddle()" );

	/* make low edge of cell exact for photo integrals */
	ASSERT( ipnt >= 0 );
	ASSERT( ipnt < rfield.nupper-1 );

	/* upper edge of higher cell*/
	Ehi = rfield.anu[ipnt] + rfield.widflx[ipnt]*0.5f;
	/* lower edge of lower cell */
	Elo = rfield.anu[ipnt-1] - rfield.widflx[ipnt-1]*0.5f;

	/* >>chng 02 nov 11, do nothing if already very close,
	 * because VERY high res models had negative widths for some very close edges */
	if( fabs( Elo/exact - 1. ) < 0.001 ) 
		return;

	ASSERT( Elo <= exact );

	OldEner = rfield.anu[ipnt];

	/* centroid of desired lower energy and upper edge */
	rfield.anu[ipnt] = (realnum)((Ehi + exact)/2.);
	/* same for lower */
	rfield.anu[ipnt-1] = (realnum)((exact + Elo)/2.);

	rfield.widflx[ipnt] = (realnum)(Ehi - exact);
	rfield.widflx[ipnt-1] = (realnum)(exact - Elo);

	/* bring upper cell down a bit too, since we dont want large change in widflx */
	rfield.anu[ipnt+1] -= (OldEner - rfield.anu[ipnt])/2.f;

	/* sanity check */
	ASSERT( rfield.widflx[ipnt-1] > 0. );
	ASSERT( rfield.widflx[ipnt] > 0. );
	return;
}

/*ipShells assign continuum energy pointers to shells for all atoms,
 * called by ContCreatePointers */
STATIC void ipShells(
	/* nelem is the atomic number on the C scale, Li is 2 */
	long int nelem)
{
	char chLab[5];
	long int 
	  imax, 
	  ion, 
	  nelec, 
	  ns, 
	  nshell;
	/* following value cannot be used - will be set to proper threshold */
	double thresh=-DBL_MAX;

	DEBUG_ENTRY( "ipShells()" );

	ASSERT( nelem >= NISO);
	ASSERT( nelem < LIMELM );

	/* fills in pointers to valence shell ionization threshold
	 * PH1(a,b,c,d)
	 * a=1 => thresh, others fitting parameters
	 * b atomic number
	 * c number of electrons
	 * d shell number 7-1 */

	/* threshold in Ryd
	 * ion=0 for atom, up to nelem-1 for helium like, hydrogenic is elsewhere */
	for( ion=0; ion < nelem; ion++ )
	{
		/* generate label for ionization stage */
		/* this is short form of element name */
		strcpy( chLab, elementnames.chElementSym[nelem] );

		/* this is a number between 1 and 31 */
		strcat( chLab, elementnames.chIonStage[ion] );

		/* this is the iso sequence - must not redo sequence if done as iso */
		long int ipISO = nelem-ion;

		/* number of bound electrons */
		nelec = ipISO+1;

		/* nsShells(nelem,ion) is the number of shells for ion with nelec electrons,
		 * physical not c scale */
		imax = Heavy.nsShells[nelem][ion];

		/* loop on all inner shells, valence shell */
		for( nshell=0; nshell < imax; nshell++ )
		{
			/* ionization potential of this shell in rydbergs */
			thresh = (double)(t_ADfA::Inst().ph1(nshell,nelec-1,nelem,0)/EVRYD* 0.9998787);
			if( thresh <= 0.1 )
			{
				/* negative ip shell does not exist, set upper limit
				 * to less than lower limit so this never looped upon
				 * these are used as flags by LimitSh to check whether
				 * this is a real shell - if 1 or 2 is changed - change LimitSh!! */
				opac.ipElement[nelem][ion][nshell][0] = 2;
				opac.ipElement[nelem][ion][nshell][1] = 1;
			}
			else
			{
				/* this is lower bound to energy range for this shell */
				/* >>chng 02 may 27, change to version of ip with label, so that
				 * inner shell edges will appear */
				/*opac.ipElement[nelem][ion][nshell][0] = ipoint(thresh);*/
				opac.ipElement[nelem][ion][nshell][0] = 
					ipContEnergy( thresh , chLab );

				/* this is upper bound to energy range for this shell 
				 * LimitSh is an integer function, returns pointer
				 * to threshold of next major shell.  For k-shell it
				 * returns the values KshellLimit, default=7.35e4
				 * >>chng 96 sep 26, had been below, result zero cross sec at 
				 * many energies where opacity project did not produce state specific 
				 * cross section */
				opac.ipElement[nelem][ion][nshell][1] = 
					LimitSh(ion+1,  nshell+1,nelem+1);
				ASSERT( opac.ipElement[nelem][ion][nshell][1] > 0);
			}
		}

		ASSERT( imax > 0 && imax <= 7 );

		/* this will be index pointing to valence edge */
		/* [0] is pointer to threshold in energy array */
		opac.ipElement[nelem][ion][imax-1][0] = 
			ipContEnergy(thresh, chLab);

		/* pointer to valence electron ionization potential */
		Heavy.ipHeavy[nelem][ion] = opac.ipElement[nelem][ion][imax-1][0];

		/* ionization potential of valence shell in Ryd 
		 * thresh was evaluated above, now has last value, the valence shell */
		Heavy.Valence_IP_Ryd[nelem][ion] = thresh;

		Heavy.xLyaHeavy[nelem][ion] = 0.;
		if( ipISO >= NISO )
		{
			/* this is set of 3/4 of valence shell IP, this is important
			 * source of ots deep in cloud */
			Heavy.ipLyHeavy[nelem][ion] = 
				ipLineEnergy(thresh*0.75,chLab , 0);


			Heavy.ipBalHeavy[nelem][ion] = 
				ipLineEnergy(thresh*0.25,chLab , 0);
		}
		else
		{
			/* do not treat this simple way since done exactly with iso 
			 * sequences */
			Heavy.ipLyHeavy[nelem][ion] = -1;
			Heavy.ipBalHeavy[nelem][ion] = -1;
		}
	}

	/* above loop did up to hydrogenic, now do hydrogenic - 
	 * hydrogenic is special since arrays already set up */
	Heavy.nsShells[nelem][nelem] = 1;

	/* this is lower limit to range */
	/* hydrogenic photoionization set to special hydro array 
	 * this is pointer to threshold energy */
	/* this statement is in ContCreatePointers but has not been done when this routine called */
	/*iso.ipIsoLevNIonCon[ipH_LIKE][ipZ][ipLo] = ipContEnergy(iso.xIsoLevNIonRyd[ipH_LIKE][ipZ][ipLo],chLab);*/
	/*opac.ipElement[nelem][nelem][0][0] = iso.ipIsoLevNIonCon[ipH_LIKE][nelem][ipH1s];*/
	opac.ipElement[nelem][nelem][0][0] = ipoint( iso.xIsoLevNIonRyd[ipH_LIKE][nelem][ipH1s] );
	ASSERT( opac.ipElement[nelem][nelem][0][0] > 0 );

	/* this is the high-energy limit */
	opac.ipElement[nelem][nelem][0][1] = continuum.KshellLimit;

	Heavy.ipHeavy[nelem][nelem] = opac.ipElement[nelem][nelem][0][0];

	/* this is for backwards computability with Cambridge code */
	if( trace.lgTrace && trace.lgPointBug )
	{
		for( ion=0; ion < (nelem+1); ion++ )
		{
			fprintf( ioQQQ, "Ion:%3ld%3ld %2.2s%2.2s total shells:%3ld\n", 
			  nelem, ion+1, elementnames.chElementSym[nelem], elementnames.chIonStage[ion]
			  , Heavy.nsShells[nelem][ion] );
			for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
			{
				fprintf( ioQQQ, " shell%3ld %2.2s range eV%10.2e-%8.2e\n", 
				  ns+1, Heavy.chShell[ns], rfield.anu[opac.ipElement[nelem][ion][ns][0]-1]*
				  EVRYD, rfield.anu[opac.ipElement[nelem][ion][ns][1]-1]*EVRYD );
			}
		}
	}
	return;
}

/*LimitSh sets upper energy limit to subshell integrations */
STATIC long LimitSh(long int ion, 
  long int nshell, 
  long int nelem)
{
	long int LimitSh_v;

	DEBUG_ENTRY( "LimitSh()" );

	/* this routine returns the high-energy limit to the energy range
	 * for photoionization of a given shell
	 * */
	if( nshell == 1 )
	{
		/* this limit is high-energy limit to code unless changed with set kshell */
		LimitSh_v = continuum.KshellLimit;

	}
	else if( nshell == 2 )
	{
		/* this is 2s shell, upper limit is 1s
		 * >>chng 96 oct 08, up to high-energy limit
		 * LimitSh = ipElement(nelem,ion , 1,1)-1 */
		LimitSh_v = continuum.KshellLimit;

	}
	else if( nshell == 3 )
	{
		/* this is 2p shell, upper limit is 1s
		 * >>chng 96 oct 08, up to high-energy limit
		 * LimitSh = ipElement(nelem,ion , 1,1)-1 */
		LimitSh_v = continuum.KshellLimit;

	}
	else if( nshell == 4 )
	{
		/* this is 3s shell, upper limit is 2p
		 * >>chng 96 oct 08, up to K-shell edge
		 * LimitSh = ipElement(nelem,ion , 3,1)-1 */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][0][0] - 1;

	}
	else if( nshell == 5 )
	{
		/* this is 3p shell, upper limit is 2p
		 * >>chng 96 oct 08, up to K-shell edge
		 * LimitSh = ipElement(nelem,ion , 3,1)-1 */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][0][0] - 1;

	}
	else if( nshell == 6 )
	{
		/* this is 3d shell, upper limit is 2p
		 * >>chng 96 oct 08, up to K-shell edge
		 * LimitSh = ipElement(nelem,ion , 3,1)-1 */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][0][0] - 1;

	}
	else if( nshell == 7 )
	{
		/* this is 4s shell, upper limit is 3d */
		if( opac.ipElement[nelem-1][ion-1][5][0] < 3 )
		{
			/* this is check for empty shell 6, 3d
			 * if so then set to 3p instead */
			LimitSh_v = opac.ipElement[nelem-1][ion-1][4][0] - 
			  1;
		}
		else
		{
			LimitSh_v = opac.ipElement[nelem-1][ion-1][5][0] - 
			  1;
		}
		/* >>chng 96 sep 26, set upper limit down to 2s */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][2][0] - 1;

	}
	else
	{
		fprintf( ioQQQ, " LimitSh cannot handle nshell as large as%4ld\n", 
		  nshell );
		cdEXIT(EXIT_FAILURE);
	}
	return( LimitSh_v );
}

/*ContBandsCreate - read set of continuum bands to enter total emission into line*/
STATIC void ContBandsCreate(
	/* chFile is optional filename, if void then use default bands,
	 * if not void then use file specified,
	 * return value is 0 for success, 1 for failure */
	 const char chFile[] )
{
	char chLine[FILENAME_PATH_LENGTH_2];
	const char* chFilename;
	FILE *ioDATA;

	bool lgEOL;
	long int i,k;

	/* keep track of whether we have been called - want to be
	 * called a total of one time */
	static bool lgCalled=false;

	DEBUG_ENTRY( "ContBandsCreate()" );

	/* do nothing if second or later call*/
	if( lgCalled )
	{
		/* success */
		return;
	}
	lgCalled = true;

	/* use default filename if void string, else use file specified */
	chFilename = ( strlen(chFile) == 0 ) ? "continuum_bands.ini" : chFile;

	/* get continuum band data  */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " ContBandsCreate opening %s:", chFilename );
	}

	ioDATA = open_data( chFilename, "r" );

	/* now count how many bands are in the file */
	continuum.nContBand = 0;

	/* first line is a magic number and does not count as a band*/
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " ContBandsCreate could not read first line of %s.\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
			++continuum.nContBand;
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " ContBandsCreate could not rewind %s.\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}

	continuum.ContBandWavelength = (realnum *)MALLOC(sizeof(realnum)*(unsigned)(continuum.nContBand) );
	continuum.chContBandLabels = (char **)MALLOC(sizeof(char *)*(unsigned)(continuum.nContBand) );
	continuum.ipContBandLow = (long int *)MALLOC(sizeof(long int)*(unsigned)(continuum.nContBand) );
	continuum.ipContBandHi = (long int *)MALLOC(sizeof(long int)*(unsigned)(continuum.nContBand) );
	continuum.BandEdgeCorrLow = (realnum *)MALLOC(sizeof(realnum)*(unsigned)(continuum.nContBand) );
	continuum.BandEdgeCorrHi = (realnum *)MALLOC(sizeof(realnum)*(unsigned)(continuum.nContBand) );

	/* now make second dim, id wavelength, and lower and upper bounds */
	for( i=0; i<continuum.nContBand; ++i )
	{
		/* array of labels, each 4 long plus 0 at [4] */
		continuum.chContBandLabels[i] = (char *)MALLOC(sizeof(char)*(unsigned)(5) );
	}

	/* first line is a versioning magic number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " ContBandsCreate could not read first line of %s.\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}
	/* bands_continuum magic number here <- this string is in band_continuum.dat
	 * with comment to search for this to find magic number  */
	{
		long int m1 , m2 , m3,
			myr = 11,
			mmo = 6,
			mdy = 7;
		i = 1;
		m1 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		m2 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		m3 = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( ( m1 != myr ) ||
		    ( m2 != mmo ) ||
		    ( m3 != mdy ) )
		{
			fprintf( ioQQQ, 
				" ContBandsCreate: the version of the data file %s I found (%li %li %li)is not the current version (%li %li %li).\n", 
				chFilename ,
				m1 , m2 , m3 ,
				myr , mmo , mdy );
			fprintf( ioQQQ, 
				" ContBandsCreate: you need to update this file.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* now read in data again, but save it this time */
	k = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
		{
			/* copy 4 char label plus termination */
			strncpy( continuum.chContBandLabels[k] , chLine , 4 );
			continuum.chContBandLabels[k][4] = 0;

			/* now get central band wavelength 
			 * >>chng 06 aug 11 from 4 to 6, the first 4 char are labels and
			 * these can contain numbers, next comes a space, then the number */
			i = 6;
			continuum.ContBandWavelength[k] = (realnum)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			/* >>chng 06 feb 21, multiply by 1e4 to convert micron wavelength into Angstroms,
			 * which is assumed by the code.  before this correction the band centroid 
			 * wavelength was given in the output incorrectly listed as Angstroms.
			 * results were correct just label was wrong */
			continuum.ContBandWavelength[k] *= 1e4f;

			/* these are short and long wave limits, which are high and
			 * low energy limits - these are now wl in microns but are
			 * converted to Angstroms */
			double xHi = FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL)*1e4;
			double xLow = FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL)*1e4;
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been 3 numbers on this band line.   Sorry.\n" );
				fprintf( ioQQQ, " string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			{
				/* this is an option to print out one of the two photon continua */
				enum {DEBUG_LOC=false};
				if( DEBUG_LOC )
				{
					fprintf(ioQQQ, "READ:%s\n", chLine );
					fprintf(ioQQQ, "GOT: %s %g %g %g\n",continuum.chContBandLabels[k],
					continuum.ContBandWavelength[k] , xHi , xLow );
				}
			}

			/* make sure bands bounds are in correct order, shorter - longer wavelength*/
			if( xHi >= xLow )
			{
				fprintf( ioQQQ, " ContBandWavelength band %li "
					"edges are in improper order.\n" ,k);
				fprintf(ioQQQ,"band: %s %.3e %.3e %.3e \n",
					continuum.chContBandLabels[k],
					continuum.ContBandWavelength[k],
					xHi , 
					xLow);
				cdEXIT(EXIT_FAILURE);
			}

			// check that central wavelength is indeed between the limits
			// xHi & xLow are hi and low energy limits to band so logic reversed
			if( continuum.ContBandWavelength[k] < xHi ||
				continuum.ContBandWavelength[k] > xLow )
			{
				fprintf( ioQQQ, " ContBandWavelength band %li central "
					"wavelength not within band.\n" ,k);
				fprintf(ioQQQ,"band: %s %.3e %li %li \n",
					continuum.chContBandLabels[k],
					continuum.ContBandWavelength[k],
					continuum.ipContBandHi[k] , 
					continuum.ipContBandLow[k]);
				cdEXIT(EXIT_FAILURE);
			}

			/* get continuum index - RYDLAM is 911.6A = 1 Ryd so 1e4 converts 
			 * micron to Angstrom - xHi is high energy (not wavelength)
			 * edge of the band */
			continuum.ipContBandHi[k] = ipoint( RYDLAM / xHi );
			continuum.ipContBandLow[k] = ipoint( RYDLAM / xLow );

			/* fraction of first and last bin to include */
			continuum.BandEdgeCorrLow[k] = 
				((rfield.anu[continuum.ipContBandLow[k]-1]+
				rfield.widflx[continuum.ipContBandLow[k]-1]/2.f)-
				(realnum)(RYDLAM/xLow)) /
				rfield.widflx[continuum.ipContBandLow[k]-1];
			ASSERT( continuum.BandEdgeCorrLow[k]>=0. && continuum.BandEdgeCorrLow[k]<=1.);
			continuum.BandEdgeCorrHi[k] = ((realnum)(RYDLAM/xHi) -
				(rfield.anu[continuum.ipContBandHi[k]-1]-
				rfield.widflx[continuum.ipContBandHi[k]-1]/2.f)) /
				rfield.widflx[continuum.ipContBandHi[k]-1];
			ASSERT( continuum.BandEdgeCorrHi[k]>=0. && continuum.BandEdgeCorrHi[k]<=1.);
			/*fprintf(ioQQQ,"DEBUG bands_continuum %s %.3e %li %li \n",
				continuum.chContBandLabels[k],
				continuum.ContBandWavelength[k],
				continuum.ipContBandHi[k] , 
				continuum.ipContBandLow[k]);*/

			if( trace.lgTrace && trace.lgConBug )
			{
				if( k==0 )
					fprintf( ioQQQ, "   ContCreatePointer trace bands\n");
				fprintf( ioQQQ, 
					"     band %ld label %s low wl= %.3e low ipnt= %li "
					" hi wl= %.3e hi ipnt= %li \n", 
					k, 
					continuum.chContBandLabels[k] ,
					xLow,
					continuum.ipContBandLow[k] ,
					xHi,
					continuum.ipContBandHi[k]  );
			}
#			if 0
			// hazy table giving band properties
#			include "prt.h"
			fprintf(ioQQQ,
				"DEBUG %s & ", 
				continuum.chContBandLabels[k] );
			prt_wl( ioQQQ , continuum.ContBandWavelength[k] );
			fprintf(ioQQQ," & ");
			prt_wl( ioQQQ , xHi );
			fprintf(ioQQQ," -- ");
			prt_wl( ioQQQ , xLow );
			fprintf(ioQQQ,"\\\\ \n");
#			endif
			++k;
		}
	}
	/* now validate this incoming data */
	for( i=0; i<continuum.nContBand; ++i )
	{
		/* make sure all are positive */
		if( continuum.ContBandWavelength[i] <=0. )
		{
			fprintf( ioQQQ, " ContBandWavelength band %li has non-positive entry.\n",i );
			cdEXIT(EXIT_FAILURE);
		}
	}

	fclose(ioDATA);
	return;
}
