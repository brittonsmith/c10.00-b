/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IterStart set and save values of many variables at start of iteration after initial temperature set*/
/*IterRestart restart iteration */
#include "cddefines.h"
#include "cddrive.h"
#include "physconst.h"
#include "iso.h"
#include "taulines.h"
#include "hydrogenic.h"
#include "struc.h"
#include "dynamics.h"
#include "prt.h"
#include "hyperfine.h"
#include "dense.h"
#include "magnetic.h"
#include "continuum.h"
#include "geometry.h"
#include "h2.h"
#include "he.h"
#include "grains.h"
#include "atomfeii.h"
#include "pressure.h"
#include "stopcalc.h"
#include "conv.h"
#include "mean.h"
#include "ca.h"
#include "thermal.h"
#include "atoms.h"
#include "wind.h"
#include "opacity.h"
#include "timesc.h"
#include "trace.h"
#include "colden.h"
#include "secondaries.h"
#include "hmi.h"
#include "radius.h"
#include "phycon.h"
#include "called.h"
#include "mole.h"
#include "rfield.h"
#include "ionbal.h"
#include "atmdat.h"
#include "lines.h"
#include "molcol.h"
#include "input.h"
#include "rt.h"
#include "iterations.h"
#include "cosmology.h"
/* these are the static saved variables, set in IterStart and reset in IterRestart */

static double h2plus_heat_save, HeatH2Dish_used_save, HeatH2Dexc_used_save,
	hmihet_save, hmitot_save, H2_Solomon_dissoc_rate_used_H2g_save,deriv_HeatH2Dexc_used_save,
	H2_Solomon_dissoc_rate_used_H2s_save, H2_H2g_to_H2s_rate_used_save, H2_photodissoc_used_H2s_save,
	H2_photodissoc_used_H2g_save, 
	UV_Cont_rel2_Draine_DB96_face,
	UV_Cont_rel2_Draine_DB96_depth , UV_Cont_rel2_Habing_TH85_face ,
	**saveMoleSource , **saveMoleSink;
static realnum ***SaveMoleChTrRate;

/* arguments are atomic number, ionization stage*/
static realnum xIonFsave[LIMELM+3][LIMELM+1];
static double HeatSave[LIMELM][LIMELM];
static realnum supsav[LIMELM][LIMELM],p2nit,d5200r;
/* save values of dr and drNext */
static double drSave , drNextSave;

/* also places to save them between iterations */
static long int IonLowSave[LIMELM], 
	IonHighSave[LIMELM];

/* these are used to reset the level populations of the h and he model atoms */
/*static realnum hnsav[LIMELM][LMHLVL+1], */
static bool lgHNSAV = false;

static double  ***HOpacRatSav;
static realnum H_sum_in_CO_save;

static realnum hmsav[N_H_MOLEC];

static double ortho_save , para_save;
static double ***hnsav,
  edsav;

static realnum xMolFracsSave[LIMELM],
	gas_phase_save[LIMELM];

/*IterStart set and save values of many variables at start of iteration after initial temperature set*/
void IterStart(void)
{
	long int i, 
	  ion, 
	  ipISO,
	  n ,
	  nelem;
	double fhe1nx, 
	  ratio;

	DEBUG_ENTRY( "IterStart()" );

	/* allocate two matrices if first time in this core load 
	 * these variables are malloced here because they are file static - 
	 * used to save values at start of first zone, and reset to this value at
	 * start of new iteration */
	if( !lgHNSAV )
	{
		/* set flag so we never do this again */
		lgHNSAV = true;

		HOpacRatSav = (double ***)MALLOC(sizeof(double **)*NISO );

		hnsav = (double ***)MALLOC(sizeof(double **)*NISO );

		for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{

			HOpacRatSav[ipISO] = (double **)MALLOC(sizeof(double *)*LIMELM );

			hnsav[ipISO] = (double **)MALLOC(sizeof(double *)*LIMELM );

			/* now do the second dimension */
			for( nelem=ipISO; nelem<LIMELM; ++nelem )
			{
				HOpacRatSav[ipISO][nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(iso.numLevels_max[ipISO][nelem]+1) );

				hnsav[ipISO][nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(iso.numLevels_max[ipISO][nelem]+1) );
			}
		}
		saveMoleSource = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		saveMoleSink = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		SaveMoleChTrRate = 
			(realnum***)MALLOC(sizeof(realnum**)*(unsigned)LIMELM );
		for(nelem=0; nelem<LIMELM; ++nelem )
		{
			/* chemistry source and sink terms for ionization ladders */
			saveMoleSource[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+2) );
			saveMoleSink[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+2) );
			SaveMoleChTrRate[nelem] = 
				(realnum**)MALLOC(sizeof(realnum*)*(unsigned)(nelem+2) );
			for( ion=0; ion<nelem+2; ++ion )
			{
				SaveMoleChTrRate[nelem][ion] = 
					(realnum*)MALLOC(sizeof(realnum)*(unsigned)(nelem+2) );
			}
		}
	}

	/* IterStart is called at the start of EVERY iteration */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " IterStart called.\n" );
	}

	/* check if this is the last iteration */
	if( iteration > iterations.itermx )
		iterations.lgLastIt = true;

	/* flag set true in RT_line_one_tauinc if maser ever capped */
	rt.lgMaserCapHit = false;
	/* flag for remembering if maser ever set dr in any zone */
	rt.lgMaserSetDR = false;

	/* zero out charge transfer heating and cooling fractions */
	atmdat.HCharHeatMax = 0.;
	atmdat.HCharCoolMax = 0.;

	/* these are set false if limits of Case B tables is ever exceeded */
	for(nelem=0; nelem<HS_NZ; ++nelem )
	{
		atmdat.lgHCaseBOK[0][nelem] = true;
		atmdat.lgHCaseBOK[1][nelem] = true;
	}

	/* reset the largest number of levels in the database species */
	for( long ipSpecies=0; ipSpecies<nSpecies; ++ipSpecies )
	{
		Species[ipSpecies].numLevels_local = Species[ipSpecies].numLevels_max;
		Species[ipSpecies].CoolTotal = 0.;
	}

	/* the reason this calculation stops */
	strncpy( StopCalc.chReasonStop, "reason not specified.", 
		sizeof(StopCalc.chReasonStop) );

	/* zero fractions of He0 destruction due to 23S */
	he.nzone = 0;
	he.frac_he0dest_23S = 0.;
	he.frac_he0dest_23S_photo = 0.;

	dense.EdenMax = 0.;
	dense.EdenMin = BIGFLOAT;

	timesc.time_H2_Dest_longest = 0.;
	timesc.time_H2_Form_longest = 0.;
	/* remains neg if not evaluated */
	timesc.time_H2_Dest_here = -1.;
	timesc.time_H2_Form_here = 0.;

	for( i=0; i < mole.num_comole_calc; i++ )
	{
		timesc.AgeCOMoleDest[i] = 0.;
	}
	timesc.BigCOMoleForm = 0.;
	timesc.TimeH21cm = 0.;
	thermal.CoolHeatMax = 0.;
	hydro.HCollIonMax = 0.;

	atmdat.HIonFracMax = 0.;

	/* will save highest n=2p/1s ratio for hydrogen atom*/
	hydro.pop2mx = 0.;
	hydro.lgHiPop2 = false;

	hydro.nLyaHot = 0;
	hydro.TLyaMax = 0.;

	/* evaluate the gas and radiation pressure */
	/* this sets values of pressure.PresTotlCurr */
	pressure.PresInteg = 0.;
	pressure.pinzon = 0.;
	pressure.PresIntegElecThin = 0.;
	PresTotCurrent();

	dynamics.HeatMax = 0.;
	dynamics.CoolMax = 0.;
	if( dynamics.lgAdvection )
		DynaIterStart();

	/* reset these since we now have a valid solution */
	pressure.pbeta = 0.;
	pressure.RadBetaMax = 0.;
	pressure.lgPradCap = false;
	pressure.lgPradDen = false;
	/* this flag will say we hit the sonic point */
	pressure.lgSonicPoint = false;
	pressure.lgStrongDLimbo = false;

	/* define two timescales for equilibrium, Compton and thermal */
	timesc.tcmptn = 1.e0/(rfield.qtot*6.65e-25*dense.eden);
	timesc.time_therm_long = 1.5*dense.pden*BOLTZMANN*phycon.te/thermal.ctot;

	/* this will be total mass in computed structure */
	dense.xMassTotal = 0.;

	for( nelem=0; nelem < LIMELM; nelem++ )
	{

		if( dense.lgElmtOn[nelem] )
		{

			/* now zero out the ionic fractions */
			for( ion=0; ion < (nelem + 2); ion++ )
			{
				xIonFsave[nelem][ion] = dense.xIonDense[nelem][ion];
			}

			for( ion=0; ion < LIMELM; ion++ )
			{
				HeatSave[nelem][ion] = thermal.heating[nelem][ion];
			}
		}
	}

	/* >>chng 02 jan 28, add ipISO loop */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* >>chng 03 apr 11, from 0 to ipISO */
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( n=0; n < iso.numLevels_max[ipISO][nelem]; n++ )
				{
					HOpacRatSav[ipISO][nelem][n] = iso.ConOpacRatio[ipISO][nelem][n];
					hnsav[ipISO][nelem][n] = StatesElemNEW[nelem][nelem-ipISO][n].Pop;
				}
			}
			iso.TwoNu_induc_dn_max[ipISO][nelem] = 0.;
		}
	}

	rfield.ipEnergyBremsThin = 0;
	rfield.EnergyBremsThin = 0.;

	p2nit = atoms.p2nit;
	d5200r = atoms.d5200r;
	atoms.nNegOI = 0;
	for( i=0; i< N_OI_LEVELS; ++i )
		atoms.popoi[i] = 0.;

	/* save molecular fractions and ionization range */
	for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		/* save molecular densities */
		xMolFracsSave[nelem] = dense.xMolecules[nelem];
		gas_phase_save[nelem] = dense.gas_phase[nelem];
		IonLowSave[nelem] = dense.IonLow[nelem];
		IonHighSave[nelem] = dense.IonHigh[nelem];
	}
	/*fprintf(ioQQQ,"DEBUG IterStart save set to gas_phase hden %.4e \n",
		dense.gas_phase[0]);*/

	edsav = dense.eden;
	H_sum_in_CO_save = dense.H_sum_in_CO;

	/* save molecular densities */
	for( i=0; i<N_H_MOLEC; i++) 
	{
		hmsav[i] = hmi.Hmolec[i];
	}
	ortho_save = h2.ortho_density;
	para_save = h2.para_density;

	for( i=0; i<mole.num_comole_calc; ++i )
	{
		COmole[i]->hevmol_save = COmole[i]->hevmol;
	}

	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		/* these have one more ion than above */
		for( ion=0; ion<nelem+2; ++ion )
		{
			long int ion2;
			/* zero out the source and sink arrays */
			saveMoleSource[nelem][ion] = mole.source[nelem][ion];
			saveMoleSink[nelem][ion] = mole.sink[nelem][ion];
			/*>>chng 06 jun 27, move from here, where leak occurs, to
			* above where xMoleChTrRate first created */
			/*mole.xMoleChTrRate[nelem][ion] = (realnum*)MALLOC(sizeof(realnum)*(unsigned)(nelem+2));*/
			for( ion2=0; ion2<nelem+2; ++ion2 )
			{
				SaveMoleChTrRate[nelem][ion][ion2] = mole.xMoleChTrRate[nelem][ion][ion2];
			}
		}
	}

	hmihet_save = hmi.hmihet;
	hmitot_save = hmi.hmitot;

	h2plus_heat_save = hmi.h2plus_heat;

	HeatH2Dish_used_save = hmi.HeatH2Dish_used;
	HeatH2Dexc_used_save = hmi.HeatH2Dexc_used;

	deriv_HeatH2Dexc_used_save = hmi.deriv_HeatH2Dexc_used;
	H2_Solomon_dissoc_rate_used_H2g_save = hmi.H2_Solomon_dissoc_rate_used_H2g;
	H2_Solomon_dissoc_rate_used_H2s_save = hmi.H2_Solomon_dissoc_rate_used_H2s;

	H2_H2g_to_H2s_rate_used_save = hmi.H2_H2g_to_H2s_rate_used;

	H2_photodissoc_used_H2s_save = hmi.H2_photodissoc_used_H2s;
	H2_photodissoc_used_H2g_save = hmi.H2_photodissoc_used_H2g;

	UV_Cont_rel2_Draine_DB96_face = hmi.UV_Cont_rel2_Draine_DB96_face;
	UV_Cont_rel2_Draine_DB96_depth = hmi.UV_Cont_rel2_Draine_DB96_depth;
	UV_Cont_rel2_Habing_TH85_face = hmi.UV_Cont_rel2_Habing_TH85_face;

	/* save zone thickness, and next zone thickness */
	drSave = radius.drad;
	drNextSave = radius.drNext;
	/* remember largest and smallest dr over iteration */
	radius.dr_min_last_iter = BIGFLOAT;
	radius.dr_max_last_iter = 0.;

	geometry.nprint = iterations.IterPrnt[iteration-1];

	colden.TotMassColl = 0.;
	colden.tmas = 0.;
	colden.wmas = 0.;
	colden.rjnmin = FLT_MAX;
	colden.ajmmin = FLT_MAX;

	thermal.nUnstable = 0;
	thermal.lgUnstable = false;

	rfield.extin_mag_B_point = 0.;
	rfield.extin_mag_V_point = 0.;
	rfield.extin_mag_B_extended = 0.;
	rfield.extin_mag_V_extended = 0.;

	/* plus 1 is to zero out point where unit integration occurs */
	for( i=0; i < rfield.nupper+1; i++ )
	{
		/* diffuse fields continuum */
		/* >>chng 03 nov 08, recover SummedDif */
		rfield.SummedDifSave[i] = rfield.SummedDif[i];
		rfield.ConEmitReflec[0][i] = 0.;
		rfield.ConEmitOut[0][i] = 0.;
		rfield.ConInterOut[i] = 0.;
		/* save OTS continuum for next iteration */
		rfield.otssav[i][0] = rfield.otscon[i];
		rfield.otssav[i][1] = rfield.otslin[i];
		rfield.outlin[0][i] = 0.;
		rfield.outlin_noplot[i] = 0.;
		rfield.ConOTS_local_OTS_rate[i] = 0.;
		rfield.ConOTS_local_photons[i] = 0.;
		rfield.DiffuseEscape[i] = 0.;

		/* save initial absorption, scattering opacities for next iteration */
		opac.opacity_abs_savzon1[i] = opac.opacity_abs[i];
		opac.opacity_sct_savzon1[i] = opac.opacity_sct[i];

		/* will accumulate optical depths through cloud */
		opac.TauAbsFace[i] = opac.taumin;
		opac.TauScatFace[i] = opac.taumin;
		/* >>chng 99 dec 04, having exactly 1 zone thickness for first zone caused discontinuity
		 * for heating in very high T models in func_map.in.  zone 1 and 2 were 20% different,
		 * since tau in is 1e-20, e2 is 0.9999, and so some H ots
		 * these were not here at all - changed to dr/2 */
		/* attenuation of flux by optical depths IN THIS ZONE 
		 * DirectionalCosin is 1/COS(theta), is usually 1, reset with illuminate command,
		 * option for illumination of slab at an angle */
		/* >>chng 04 oct 09, from drad to radius.drad_x_fillfac - include fill fac, PvH */
		opac.ExpZone[i] = sexp(opac.opacity_abs[i]*radius.drad_x_fillfac/2.*geometry.DirectionalCosin);

		/* e(-tau) in inward direction, up to illuminated face */
		opac.ExpmTau[i] = (realnum)opac.ExpZone[i];

		/* e2(tau) in inward direction, up to illuminated face, this is used to get the
		 * recombination escape probability */
		opac.E2TauAbsFace[i] = (realnum)e2(opac.TauAbsFace[i] );
	}

	t_fe2ovr_la::Inst().zero_opacity();

	/* this zeros out arrays to hold mean ionization fractions
	 * later entered by mean
	 * read out by setlim */
	mean.MeanZero();

	/* zero out the column densities */
	for( i=0; i < NCOLD; i++ )
	{
		colden.colden[i] = 0.;
	}
	colden.He123S = 0.;
	colden.H0_ov_Tspin = 0.;
	colden.OH_ov_Tspin = 0.;
	colden.coldenH2_ov_vel = 0.;

	/* upper and lower levels of H0 1s */
	colden.H0_21cm_upper =0;
	colden.H0_21cm_lower =0;

	h2.ortho_colden = 0.;
	h2.para_colden = 0.;

	for( i=0; i < mole.num_comole_calc; i++ )
	{
		COmole[i]->hevcol = 0.;
	}
	for( i=0; i < 5; i++ )
	{
		colden.C2Pops[i] = 0.;
		colden.C2Colden[i] = 0.;
		/* pops and column density for SiII atom */
		colden.Si2Pops[i] = 0.;
		colden.Si2Colden[i] = 0.;
	}
	for( i=0; i < 3; i++ )
	{
		/* pops and column density for CI atom */
		colden.C1Pops[i] = 0.;
		colden.C1Colden[i] = 0.;
		/* pops and column density for OI atom */
		colden.O1Pops[i] = 0.;
		colden.O1Colden[i] = 0.;
		/* pops and column density for CIII atom */
		colden.C3Pops[i] = 0.;
	}
	for( i=0; i < 4; i++ )
	{
		/* pops and column density for CIII atom */
		colden.C3Colden[i] = 0.;
	}

	/* these are some line of sight emission measures */
	colden.dlnenp = 0.;
	colden.dlnenHep = 0.;
	colden.dlnenHepp = 0.;
	colden.dlnenCp = 0.;
	colden.H0_ov_Tspin = 0.;
	colden.OH_ov_Tspin = 0.;

	/* now zero heavy element molecules */
	molcol("ZERO",ioQQQ);

	/* this will be sum of all free free heating over model */
	thermal.FreeFreeTotHeat = 0.;

	thermal.thist = 0.;
	thermal.tlowst = 1e20f;

	wind.AccelAver = 0.;
	wind.acldr = 0.;
	ionbal.ilt = 0;
	ionbal.iltln = 0;
	ionbal.ilthn = 0;
	ionbal.ihthn = 0;
	ionbal.ifail = 0;

	secondaries.SecHIonMax = 0.;
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		for( ion=0; ion<nelem+1; ++ion )
		{
			supsav[nelem][ion] = secondaries.csupra[nelem][ion];
		}
	}
	secondaries.x12sav = secondaries.x12tot;
	secondaries.hetsav = secondaries.HeatEfficPrimary;
	secondaries.savefi = secondaries.SecIon2PrimaryErg;
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		for( ion=0; ion<nelem+1; ++ion )
		{
			ionbal.CompRecoilHeatRateSave[nelem][ion] = ionbal.CompRecoilHeatRate[nelem][ion];
			ionbal.CompRecoilIonRateSave[nelem][ion] = ionbal.CompRecoilIonRate[nelem][ion];
		}
	}

	/* these will keep track of the number of convergence failures that occur */
	conv.nTeFail = 0;
	conv.nTotalFailures = 0;
	conv.nPreFail = 0;
	conv.failmx = 0.;
	conv.nIonFail = 0;
	conv.nPopFail = 0;
	conv.nNeFail = 0;
	conv.nGrainFail = 0;
	conv.dCmHdT = 0.;

	/* abort flag */
	lgAbort = false;

	GrainStartIter();

	rfield.comtot = 0.;

	co.codfrc = 0.;
	co.codtot = 0.;

	hmi.HeatH2DexcMax = 0.;
	hmi.CoolH2DexcMax = 0.;
	hmi.h2dfrc = 0.;
	hmi.h2line_cool_frac = 0.;
	hmi.h2dtot = 0.;
	thermal.HeatLineMax = 0.;
	thermal.GBarMax = 0.;
	hyperfine.cooling_max = 0.;
	hydro.cintot = 0.;
	geometry.lgZoneTrp = false;

	hmi.h2pmax = 0.;

	/************************************************************************
	 *
	 * allocate space for lines arrays 
	 *
	 ************************************************************************/

	if( iteration == 1 )
	{
		/* there are three types of calls to lines() 
		 * ipass = -1, first call, only count number of lines
		 * ipass =  0, second pass, save labels and wavelengths
		 * ipass =  1, integrate intensity*/
		LineSave.ipass = -1;
		lines();
		ASSERT( LineSave.nsum > 0 );

		/* in a grid or MPI run this may not be the first time here,
		 * return old memory and grab new appropriate for this size,
		 * since number of lines to be stored can change */
		if( LineSv != NULL )
			free( LineSv );
		if( LineSvSortWL != NULL )
			free( LineSvSortWL );

		/* this is the large main line array */
		LineSv = (LinSv*)MALLOC((unsigned)LineSave.nsum*sizeof( LinSv ) );
		LineSvSortWL = (LinSv*)MALLOC((unsigned)LineSave.nsum*sizeof( LinSv ));

		/* this will be used as sanity check in future iterations and grid points */
		LineSave.nsumAllocated = LineSave.nsum;

		/* this is only done on first iteration since will integrate over time */
		for( i=0; i < LineSave.nsum; i++ )
		{
			/* this is array of wavelengths, only defined for some lines,
			 * negative is sentinel saying whether line should be added into total continuum */
			LineSv[i].SumLine[2] = 0;
			LineSv[i].SumLine[3] = 0;
		}

		/* there are three calls to lines() 
		 * first call with ipass = -1 was above, only counted number
		 * of lines and malloced space.  This is second call and will do
		 * one-time creation of line labels */
		LineSave.ipass = 0;
		lines();
		/* has to be positive */
		ASSERT( LineSave.nsum > 0);
		/* in the future calls to lines will result in integrations */
		LineSave.ipass = 1;

		if( trace.lgTrace )
			fprintf( ioQQQ, "%7ld lines printed in main line array\n", 
			  LineSave.nsum );
	}

	/* this was set in call to lines above */
	ASSERT( LineSave.nsum > 0);
	ASSERT( LineSave.nsumAllocated == LineSave.nsum );

	/* zero emission line arrays - this has to be done on every iteration */
	for( i=0; i < LineSave.nsum; i++ )
	{
		LineSv[i].SumLine[0] = 0;
		LineSv[i].SumLine[1] = 0;
	}

	/* now zero out some variables set by last call to LINES */
	hydro.cintot = 0.;
	thermal.totcol = 0.;
	rfield.comtot = 0.;
	thermal.FreeFreeTotHeat = 0.;
	thermal.power = 0.;
	thermal.HeatLineMax = 0.;
	thermal.GBarMax = 0.;
	hyperfine.cooling_max = 0.;

	hmi.h2pmax = 0.;

	co.codfrc = 0.;
	co.codtot = 0.;

	hmi.h2dfrc = 0.;
	hmi.h2line_cool_frac = 0.;
	hmi.h2dtot = 0.;
	timesc.sound = 0.;

	if( LineSave.lgNormSet )
	{
		/* set normalization line index to zero from negative initial value so that
		 * we pass the assert in cdLine, in order to get the norm line index */
		LineSave.ipNormWavL = 0;
		LineSave.ipNormWavL = cdLine( LineSave.chNormLab , LineSave.WavLNorm , &fhe1nx , &ratio );
		if( LineSave.ipNormWavL < 0 )
		{
			/* did not find the line if return is negative */
			fprintf( ioQQQ, "PROBLEM could not find the normalisation line.\n");
			fprintf( ioQQQ, 
			"IterStart could not find a line with a wavelength of %f and a label of %s\n", 
			  LineSave.WavLNorm, LineSave.chNormLab );
			fprintf( ioQQQ, "Please check the emission line output to find the correct line identification.\n");
			fprintf( ioQQQ, "Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		/* set normalization line index to zero from negative initial value so that
		 * we pass the assert in cdLine, in order to get the norm line index */
		LineSave.ipNormWavL = 0;
		LineSave.ipNormWavL = cdLine( "TOTL" , 4861.36f , &fhe1nx , &ratio );
		if( LineSave.ipNormWavL < 0 )
		{
			/* did not find the line if return is negative */
			fprintf( ioQQQ, "PROBLEM could not find the default normalisation line.\n");
			fprintf( ioQQQ, 
			"IterStart could not find a line with a wavelength of 4861 and a label of TOTL\n" );
			fprintf( ioQQQ, "Please check the emission line output to find the correct line identification.\n");
			fprintf( ioQQQ, "Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* set up stop line command on first iteration 
	 * find index for lines and save for future iterations 
	 * StopCalc.nstpl is zero (false) if no stop line commands entered */
	if( iteration == 1 && StopCalc.nstpl )
	{
		/* nstpl is number of stop line commands, 0 if none entered */
		for( long int nStopLine=0; nStopLine < StopCalc.nstpl; nStopLine++ )
		{
			double relint, absint ;
			/* returns array index for line in array stack if we found the line, 
			 * return negative of total number of lines as debugging aid if line not found */
			StopCalc.ipStopLin1[nStopLine] = cdLine( StopCalc.chStopLabel1[nStopLine], 
				/* wavelength of line in angstroms, not format printed by code */
				StopCalc.StopLineWl1[nStopLine], &relint, &absint );

			if( StopCalc.ipStopLin1[nStopLine]<0 )
			{
				fprintf( ioQQQ, 
				  " IterStart could not find first line in STOP LINE command, number %ld with label *%s* and wl %.1f\n", 
				  StopCalc.ipStopLin1[nStopLine] ,
				  StopCalc.chStopLabel1[nStopLine],
				  StopCalc.StopLineWl1[nStopLine]);
				cdEXIT(EXIT_FAILURE);
			}

			StopCalc.ipStopLin2[nStopLine] = cdLine( StopCalc.chStopLabel2[nStopLine], 
				/* wavelength of line in angstroms, not format printed by code */
				StopCalc.StopLineWl2[nStopLine], &relint, &absint );

			if( StopCalc.ipStopLin2[nStopLine] < 0 )
			{
				fprintf( ioQQQ, 
				" IterStart could not find second line in STOP LINE command, line number %ld with label *%s* and wl %.1f\n", 
				StopCalc.ipStopLin1[nStopLine] ,
				StopCalc.chStopLabel1[nStopLine],
				StopCalc.StopLineWl1[nStopLine]);
				cdEXIT(EXIT_FAILURE);
			}

			if( trace.lgTrace )
			{
				fprintf( ioQQQ, 
					" stop line 1 is number %5ld wavelength is %f label is %4.4s\n", 
					StopCalc.ipStopLin1[nStopLine], 
					LineSv[StopCalc.ipStopLin1[nStopLine]].wavelength, 
					LineSv[StopCalc.ipStopLin1[nStopLine]].chALab );
				fprintf( ioQQQ, 
					" stop line 2 is number %5ld wavelength is %f label is %4.4s\n", 
				  StopCalc.ipStopLin2[nStopLine], 
				  LineSv[StopCalc.ipStopLin2[nStopLine]].wavelength, 
				  LineSv[StopCalc.ipStopLin2[nStopLine]].chALab  );
			}
		}
	}

	/* option to only print last iteration */
	if( prt.lgPrtLastIt )
	{
		if( iteration == 1 )
		{
			called.lgTalk = false;
		}

		/* initial condition of TALK may be off if optimization used or not master rank
		 * sec part is for print last command 
		 * lgTalkForcedOff is set true when cdTalk is called with false
		 * to turn off printing */
		if( iterations.lgLastIt && !prt.lgPrtStart && !called.lgTalkForcedOff )
		{
			called.lgTalk = called.lgTalkSave;
		}
	}

	if( trace.lgTrace && trace.lgOptcBug )
	{
		if( iteration > 1 )
		{
			fprintf( ioQQQ, " Outer optical depths defined.\n" );
		}
		else
		{
			fprintf( ioQQQ, " Outer optical depths NOT defined.\n" );
		}

		/* heavy element lines */
		for( i=0; i < nLevel1; i++ )
		{
			fprintf( ioQQQ, "%6f:%8.1e%8.1e%8.1e", 
				TauLines[i].WLAng, 
			  TauLines[i].Emis->TauIn, 
			  TauLines[i].Emis->TauTot, 
			  TauLines[i].Emis->Pesc );
		}

		/*Atomic and molecular lines-Humeshkar Nemala*/
		for( i=0; i <linesAdded2; i++)
		{
			fprintf( ioQQQ, "%6f:%8.1e%8.1e%8.1e", 
				 dBaseLines[i].tran->WLAng, 
			  dBaseLines[i].TauIn, 
			  dBaseLines[i].TauTot, 
			  dBaseLines[i].Pesc);
		}

		fprintf( ioQQQ, "\n" );
	}

	if( opac.lgCaseB )
	{
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " IterStart does not change mid-zone optical depth since CASE B\n" );
		}
	}

	/* check if induced recombination can be ignored */
	hydro.FracInd = 0.;
	hydro.fbul = 0.;

	/* remember some things about effects of induced rec on H only
	 * don't do ground state since SPHERE turns it off */
	for( i=ipH2p; i < (iso.numLevels_max[ipH_LIKE][ipHYDROGEN] - 1); i++ )
	{
		if( Transitions[ipH_LIKE][ipHYDROGEN][i+1][i].Emis->Aul <= iso.SmallA )
			continue;

		ratio = iso.RecomInducRate[ipH_LIKE][ipHYDROGEN][i]*iso.PopLTE[ipH_LIKE][ipHYDROGEN][i]/
			SDIV(iso.RecomInducRate[ipH_LIKE][ipHYDROGEN][i]*iso.PopLTE[ipH_LIKE][ipHYDROGEN][i] + 
			iso.RadRecomb[ipH_LIKE][ipHYDROGEN][i][ipRecRad]*
			iso.RadRecomb[ipH_LIKE][ipHYDROGEN][i][ipRecNetEsc]);
		if( ratio > hydro.FracInd )
		{
			hydro.FracInd = (realnum)ratio;
			hydro.ndclev = i;
		}

		ratio = Transitions[ipH_LIKE][ipHYDROGEN][i+1][i].Emis->pump/
			(Transitions[ipH_LIKE][ipHYDROGEN][i+1][i].Emis->pump + 
			Transitions[ipH_LIKE][ipHYDROGEN][i+1][i].Emis->Aul);

		if( ratio > hydro.fbul )
		{
			hydro.fbul = (realnum)ratio;
			hydro.nbul = i;
		}
	}

	if( trace.lgTrace )
		fprintf( ioQQQ, " IterStart returns.\n" );
	return;
}

/*IterRestart reset many variables at the start of a new iteration 
 * called by cloudy after calculation is completed, when more iterations 
 * are needed - the iteration has been incremented when this routine is
 * called so iteration == 2 after first iteration, we are starting
 * the second iteration */
void IterRestart(void)
{
	long int i, 
	  ion, 
	  ipISO ,
	  n, 
	  nelem;
	double SumOTS;

	DEBUG_ENTRY( "IterRestart()" );

	/* this is case where temperature floor has been set, if it was hit
	 * then we did a constant temperature calculation, and must go back to 
	 * a thermal solution 
	 * test on thermal.lgTemperatureConstantCommandParsed distinguishes
	 * from temperature floor option, so not reset if constant temperature
	 * was actually set */
	if( StopCalc.TeFloor > 0. && !thermal.lgTemperatureConstantCommandParsed )
	{
		thermal.lgTemperatureConstant = false;
		thermal.ConstTemp = 0.;
	}

	lgAbort = false;

	/* reset some parameters needed for magnetic field */
	Magnetic_reinit();

	opac.stimax[0] = 0.;
	opac.stimax[1] = 0.;

	H2_Reset();

	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		/* these have one more ion than above */
		for( ion=0; ion<nelem+2; ++ion )
		{
			long int ion2;
			/* zero out the source and sink arrays */
			mole.source[nelem][ion] = saveMoleSource[nelem][ion];
			mole.sink[nelem][ion] = saveMoleSink[nelem][ion];
			/*>>chng 06 jun 27, move from here, where leak occurs, to
			* above where xMoleChTrRate first created */
			/*mole.xMoleChTrRate[nelem][ion] = (realnum*)MALLOC(sizeof(realnum)*(unsigned)(nelem+2));*/
			for( ion2=0; ion2<nelem+2; ++ion2 )
			{
				mole.xMoleChTrRate[nelem][ion][ion2] = SaveMoleChTrRate[nelem][ion][ion2];
			}
		}
	}

	/* reset molecular abundances */
	for( i=0; i<N_H_MOLEC; i++) 
	{
		hmi.Hmolec[i] = hmsav[i];
	}
	hmi.H2_total = hmi.Hmolec[ipMH2g] + hmi.Hmolec[ipMH2s];
	/*fprintf(ioQQQ," IterRestar sets H2 total to %.2e\n",hmi.H2_total );*/
	h2.ortho_density = ortho_save;
	h2.para_density = para_save;

	/* zero out the column densities */
	for( i=0; i < NCOLD; i++ )
	{
		colden.colden[i] = 0.;
	}
	colden.He123S = 0.;
	colden.coldenH2_ov_vel = 0.;

	for( i=0; i < mole.num_comole_calc; i++ )
	{
		/* column densities */
		COmole[i]->hevcol = 0.;
		/* largest fraction of atoms in molecules */
		COmole[i]->xMoleFracMax = 0.;
	}

	for( i=0; i< mole.num_comole_calc; ++i )
	{
		COmole[i]->hevmol = COmole[i]->hevmol_save;

		/* check for NaN */
		ASSERT( !isnan( COmole[i]->hevmol ) );
	}

	/* close out this iteration if dynamics or time dependence is enabled */
	if( dynamics.lgAdvection )
		DynaIterEnd();

	rfield.extin_mag_B_point = 0.;
	rfield.extin_mag_V_point = 0.;
	rfield.extin_mag_B_extended = 0.;
	rfield.extin_mag_V_extended = 0.;

	hmi.hmihet = hmihet_save;
	hmi.hmitot = hmitot_save;

	hmi.BiggestH2 = -1.f;
	hmi.h2plus_heat = h2plus_heat_save;
	hmi.HeatH2Dish_used = HeatH2Dish_used_save;
	hmi.HeatH2Dexc_used = HeatH2Dexc_used_save;

	hmi.deriv_HeatH2Dexc_used = (realnum)deriv_HeatH2Dexc_used_save;
	hmi.H2_Solomon_dissoc_rate_used_H2g = H2_Solomon_dissoc_rate_used_H2g_save;
	hmi.H2_Solomon_dissoc_rate_used_H2s = H2_Solomon_dissoc_rate_used_H2s_save;
	hmi.H2_H2g_to_H2s_rate_used = H2_H2g_to_H2s_rate_used_save;
	hmi.H2_photodissoc_used_H2s = H2_photodissoc_used_H2s_save;
	hmi.H2_photodissoc_used_H2g = H2_photodissoc_used_H2g_save;

	hmi.UV_Cont_rel2_Draine_DB96_face = (realnum)UV_Cont_rel2_Draine_DB96_face;
	hmi.UV_Cont_rel2_Draine_DB96_depth = (realnum)UV_Cont_rel2_Draine_DB96_depth;
	hmi.UV_Cont_rel2_Habing_TH85_face = (realnum)UV_Cont_rel2_Habing_TH85_face;

	rfield.ipEnergyBremsThin = 0;
	rfield.EnergyBremsThin = 0.;
	rfield.lgUSphON = false;

	radius.glbdst = 0.;
	/* zone thickness, and next zone thickness */
	radius.drad = drSave;
	radius.drad_mid_zone = drSave/2.;
	radius.drNext = drNextSave;

	/* was min dr ever used? */
	radius.lgDrMinUsed = false;

	continuum.cn4861 = continuum.sv4861;
	continuum.cn1216 = continuum.sv1216;

	/* total luminosity */
	continuum.totlsv = continuum.TotalLumin;

	/* set debugger on now if NZONE desired is 0 */
	if( (trace.nznbug == 0 && iteration >= trace.npsbug) && trace.lgTrOvrd )
	{
		if( trace.nTrConvg==0 )
		{
			trace.lgTrace = true;
		}
		else
			/* trace convergence entered = but with negative flag = make positive,
			 * abs and not mult by -1 since may trigger more than one time */
			trace.nTrConvg = abs( trace.nTrConvg );

		fprintf( ioQQQ, " IterRestart called.\n" );
	}
	else
	{
		trace.lgTrace = false;
	}

	/* zero secondary suprathermals variables for first ionization attem */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		for( ion=0; ion<nelem+1; ++ion )
		{
			secondaries.csupra[nelem][ion] = supsav[nelem][ion];
		}
	}
	secondaries.x12tot = secondaries.x12sav;
	secondaries.HeatEfficPrimary = secondaries.hetsav;
	secondaries.SecIon2PrimaryErg = secondaries.savefi;
	for( nelem=0; nelem<LIMELM; ++nelem)
	{
		for( ion=0; ion<nelem+1; ++ion )
		{
			ionbal.CompRecoilHeatRate[nelem][ion] = ionbal.CompRecoilHeatRateSave[nelem][ion];
			ionbal.CompRecoilIonRate[nelem][ion] = ionbal.CompRecoilIonRateSave[nelem][ion];
		}
	}

	wind.lgVelPos = true;
	wind.AccelMax = 0.;
	wind.AccelAver = 0.;
	wind.acldr = 0.;
	wind.windv = wind.windv0;

	thermal.nUnstable = 0;
	thermal.lgUnstable = false;

	pressure.pbeta = 0.;
	pressure.RadBetaMax = 0.;
	pressure.lgPradCap = false;
	pressure.lgPradDen = false;
	/* this flag will say we hit the sonic point */
	pressure.lgSonicPoint = false;
	pressure.lgStrongDLimbo = false;

	pressure.PresInteg = 0.;
	pressure.pinzon = 0.;
	pressure.PresIntegElecThin = 0.;

	pressure.RhoGravity_dark = 0.;
	pressure.RhoGravity_self = 0.;
	pressure.RhoGravity_external = 0.;
	pressure.RhoGravity = 0.;
	pressure.IntegRhoGravity = 0.;

	dense.eden = edsav;
	dense.EdenHCorr = edsav;
	dense.EdenTrue = edsav;
	dense.H_sum_in_CO = H_sum_in_CO_save;

	for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		/* reset molecular densities */
		dense.gas_phase[nelem] = gas_phase_save[nelem];
		dense.xMolecules[nelem] = xMolFracsSave[nelem];
		dense.IonLow[nelem] = IonLowSave[nelem];
		dense.IonHigh[nelem] = IonHighSave[nelem];
	}
	/*fprintf(ioQQQ,"DEBUG IterRestart gas_phase set to save  hden %.4e\n",
		dense.gas_phase[0]);*/

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( n=ipH1s; n < iso.numLevels_max[ipISO][nelem]; n++ )
				{
					iso.ConOpacRatio[ipISO][nelem][n] = HOpacRatSav[ipISO][nelem][n];
					StatesElemNEW[nelem][nelem-ipISO][n].Pop = hnsav[ipISO][nelem][n];
				}
			}
		}
	}

	if( trace.lgTrace && trace.lgNeBug )
	{
		fprintf( ioQQQ, "     EDEN set to%12.4e by IterRestart.\n", 
		  dense.eden );
	}

#	if 0
	for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* reset total gas phase abundance to the original set */
			dense.gas_phase[nelem] = abund.solar[nelem]*dense.gas_phase[ipHYDROGEN];
		}
		else
		{
			/* >>chng 04 apr 20, set to zero if element is off */
			dense.gas_phase[nelem] = 0.;
		}
	}
#	endif
	for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		for( ion=0; ion < (nelem + 2); ion++ )
		{
			dense.xIonDense[nelem][ion] = xIonFsave[nelem][ion];
		}
		for( i=0; i < LIMELM; i++ )
		{
			thermal.heating[nelem][i] = HeatSave[nelem][i];
		}
	}

	GrainRestartIter();

	/* continuum was saved in flux_total_incident */
	for( i=0; i < rfield.nupper; i++ )
	{
		/* time-constant part of beamed continuum */
		rfield.flux_beam_const[i] = rfield.flux_beam_const_save[i];
		/* continuum flux_time_beam_save has the initial value of the
		 * time-dependent beamed continuum */
		rfield.flux_beam_time[i] = rfield.flux_time_beam_save[i]*
			rfield.time_continuum_scale;

		if( cosmology.lgDo )
		{
			double slope_ratio;
			double fac_old = TE1RYD*rfield.anu[i]/(CMB_TEMP*(1. + cosmology.redshift_start  ));
			double fac_new = TE1RYD*rfield.anu[i]/(CMB_TEMP*(1. + cosmology.redshift_current));

			if( fac_old > log10(DBL_MAX) )
			{
				slope_ratio = 0.;
			}
			else if( fac_new > log10(DBL_MAX) )
			{
				slope_ratio = 0.;
			}
			else if( fac_old > 1.e-5 )
			{
				slope_ratio = (exp(fac_old) - 1.)/(exp(fac_new) - 1.);
			}
			else
			{
				slope_ratio = (fac_old*(1. + fac_old/2.))/(fac_new*(1. + fac_new/2.));
			}

			rfield.flux_isotropic[i] = rfield.flux_isotropic_save[i] * 
				pow( ((realnum)1.f + cosmology.redshift_current)/((realnum)1.f + cosmology.redshift_start), (realnum)4.f ) * (realnum)slope_ratio;
		}
		else
		{
			/* the isotropic continuum source is time steady */
			rfield.flux_isotropic[i] = rfield.flux_isotropic_save[i];
		}

		rfield.flux[0][i] = rfield.flux_isotropic[i] + rfield.flux_beam_time[i] +
			rfield.flux_beam_const[i];

		rfield.SummedDif[i] = rfield.SummedDifSave[i];
		rfield.SummedCon[i] = rfield.flux[0][i] + rfield.SummedDif[i];
		rfield.trans_coef_zone[i] = 1.f;
		rfield.trans_coef_total[i] = 1.f;

		rfield.OccNumbIncidCont[i] = rfield.flux[0][i]*rfield.convoc[i];
		rfield.otscon[i] = rfield.otssav[i][0];
		rfield.otslin[i] = rfield.otssav[i][1];
		rfield.ConInterOut[i] = 0.;
		rfield.OccNumbBremsCont[i] = 0.;
		rfield.OccNumbDiffCont[i] = 0.;
		rfield.OccNumbContEmitOut[i] = 0.;
		rfield.outlin[0][i] = 0.;
		rfield.outlin_noplot[i] = 0.;
		rfield.ConOTS_local_OTS_rate[i] = 0.;
		rfield.ConOTS_local_photons[i] = 0.;
		rfield.DiffuseEscape[i] = 0.;

		opac.opacity_abs[i] = opac.opacity_abs_savzon1[i];
		opac.OldOpacSave[i] = opac.opacity_abs_savzon1[i];
		opac.opacity_sct[i] = opac.opacity_sct_savzon1[i];
		opac.albedo[i] = 
			opac.opacity_sct[i]/MAX2(1e-30,opac.opacity_sct[i] + opac.opacity_abs[i]);
		opac.tmn[i] = 1.;
		/* >>chng 99 dec 04, having exactly 1 for first zone caused discontinuity
		 * for heating in very high T models in func_map.in.  zone 1 and 2 were 20% different,
		 * since tau in is 1e-20, e2 is 0.9999, and so some H ots
		opac.ExpmTau[i] = 1.;
		opac.E2TauAbsFace[i] = 1.;*/
		/* >>chng 99 dec 04, having exactly 1 for first zone caused discontinuity
		 * for heating in very high T models in func_map.in.  zone 1 and 2 were 20% different,
		 * since tau in is 1e-20, e2 is 0.9999, and so some H ots
		 * these were not here at all*/
		/* attenuation of flux by optical depths IN THIS ZONE 
		 * DirectionalCosin is 1/COS(theta), is usually 1, reset with illuminate command,
		 * option for illumination of slab at an angle */
		opac.ExpZone[i] = sexp(opac.opacity_abs[i]*radius.drad/2.*geometry.DirectionalCosin);

		/* e(-tau) in inward direction, up to illuminated face */
		opac.ExpmTau[i] = (realnum)opac.ExpZone[i];

		/* e2(tau) in inward direction, up to illuminated face */
		opac.E2TauAbsFace[i] = (realnum)e2(opac.TauAbsFace[i]);
		rfield.reflin[0][i] = 0.;
		rfield.ConEmitReflec[0][i] = 0.;
		rfield.ConEmitOut[0][i] = 0.;
		rfield.ConRefIncid[0][i] = 0.;

		/* escape in the outward direction
		 * on second and later iterations define outward E2 */
		if( iteration > 1 )
		{
			/* e2 from current position to outer edge of shell */
			realnum tau = MAX2(SMALLFLOAT , opac.TauAbsTotal[i] - opac.TauAbsFace[i] );
			opac.E2TauAbsOut[i] = (realnum)e2( tau );
			ASSERT( opac.E2TauAbsOut[i]>=0. && opac.E2TauAbsOut[i]<=1. );
		}
		else
			opac.E2TauAbsOut[i] = 1.;

	}

	/* update continuum */
	RT_OTS_Update(&SumOTS);

	thermal.FreeFreeTotHeat = 0.;
	atoms.p2nit = p2nit;
	atoms.d5200r = d5200r;

	if( called.lgTalk )
	{
		fprintf( ioQQQ, "\f\n          Start Iteration Number %ld   %75.75s\n\n\n", 
		  iteration, input.chTitle );
	}

	/* reset variable to do with FeII atom, in FeIILevelPops */
	FeIIReset();
	return;
}

/* do some work with ending the iteration */
void IterEnd(void)
{

	DEBUG_ENTRY( "IterEnd()" );

	/* give indication of geometry */
	double fac = radius.depth/radius.rinner;
	if( fac < 0.1 )
	{
		geometry.lgGeoPP = true;
	}
	else
	{
		geometry.lgGeoPP = false;
	}

	if( iteration > dynamics.n_initial_relax )
	{
	    // save cumulative lines
	    for( long n=0; n<LineSave.nsum; ++n )
	    {
			 for( long nEmType=0; nEmType<2; ++nEmType )
			 {
				  LineSv[n].SumLine[nEmType+2] +=  LineSv[n].SumLine[nEmType]*dynamics.timestep;
			 }
	    }
		 // save cumulative continua
		  for( long n=0; n<rfield.nflux; ++n)
		  {

				rfield.flux[1][n] += rfield.flux[0][n]*(realnum)dynamics.timestep;
				rfield.ConEmitReflec[1][n] += rfield.ConEmitReflec[0][n]*(realnum)dynamics.timestep;
				rfield.ConEmitOut[1][n] += rfield.ConEmitOut[0][n]*(realnum)dynamics.timestep;
				rfield.ConRefIncid[1][n] += rfield.ConRefIncid[0][n]*(realnum)dynamics.timestep;
				rfield.flux_total_incident[1][n] += rfield.flux_total_incident[0][n]*(realnum)dynamics.timestep;
				rfield.reflin[1][n] += rfield.reflin[0][n]*(realnum)dynamics.timestep;
				rfield.outlin[1][n] += rfield.outlin[0][n]*(realnum)dynamics.timestep;
		  }
	}

	
	/* save previous iteration's results */
	struc.nzonePreviousIteration = nzone;
	for( long i=0; i<struc.nzonePreviousIteration; ++i )
	{
		struc.depth_last[i] = struc.depth[i];
		struc.drad_last[i] = struc.drad[i];
	}

	/* all continue were attenuated in last zone in radius_increment to represent the intensity
	 * in the middle of the next zone - this was too much since we now want
	 * intensity emergent from outer edge of last zone */
	for( long i=0; i < rfield.nflux; i++ )
	{
		{
			enum{DEBUG_LOC=false};
			if( DEBUG_LOC)
			{
				fprintf(ioQQQ,"i=%li opac %.2e \n", i, 
					(double)opac.opacity_abs[i]*radius.drad_x_fillfac/2.*(double)geometry.DirectionalCosin );
			}
		}
		double tau = opac.opacity_abs[i]*radius.drad_x_fillfac/2.*geometry.DirectionalCosin;
		ASSERT( tau > 0. );
		fac = sexp( tau );

		/* >>chng 02 dec 14, add first test to see whether product of ratio is within
		 * range of a float - ConInterOut can be finite and fac just above a small float,
		 * so ratio exceeds largest size of a float */
		/*if( fac > SMALLFLOAT )*/
		if( (realnum)(fac/SDIV(rfield.ConInterOut[i]))>SMALLFLOAT && fac > SMALLFLOAT )
		{
			rfield.ConInterOut[i] /= (realnum)fac;
			rfield.outlin[0][i] /= (realnum)fac;
			rfield.outlin_noplot[i] /= (realnum)fac;
		}
	}

	/* remember thickness of previous iteration */
	radius.StopThickness[iteration-1] = radius.depth;
	return;
}
