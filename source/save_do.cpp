/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveDo produce save output during calculation,
 * chTime is 'MIDL' during calculation, 'LAST' at the end */
/*SaveNewContinuum produce the 'save new continuum' output */
/*SaveLineStuff save optical depths or source functions for all transferred lines */
/*Save1Line called by SaveLineStuff to produce output for one line */
/*SaveNewContinuum produce the 'save new continuum' output */
/*SaveLineIntensity produce the 'save lines intensity' output */
/* save h emission, for AGN3 chapter 4, routine is below */
/*SaveResults1Line do single line of output for the save results and save line intensity commands */
/* the number of emission lines across one line of printout */
/*SaveSpecial generate output for the save special command */
/*SaveResults save results from save results command */
/*SaveResults1Line do single line of output for the save results and save line intensity commands */
/*SaveGaunts called by save gaunts command to output gaunt factors */
/*SaveFeII_cont  return appropriate FeII banded continuum, inward, outward, or total */
/*FindStrongestLineLabels find strongest lines contributing to point in continuum energy mesh, output in some save commands */
#include "cddefines.h"
#include "cddrive.h"
#include "physconst.h"
#include "mean.h"
#include "taulines.h"
#include "struc.h"
#include "iso.h"
#include "mole.h"
#include "hyperfine.h"
#include "rt.h"
#include "lines_service.h"
#include "doppvel.h"
#include "dense.h"
#include "h2.h"
#include "magnetic.h"
#include "hydrogenic.h"
#include "secondaries.h"
#include "grainvar.h"
#include "lines.h"
#include "dynamics.h"
#include "colden.h"
#include "continuum.h"
#include "ionbal.h"
#include "yield.h"
#include "prt.h"
#include "iterations.h"
#include "heavy.h"
#include "conv.h"
#include "geometry.h"
#include "called.h"
#include "helike.h"
#include "opacity.h"
#include "rfield.h"
#include "phycon.h"
#include "timesc.h"
#include "radius.h"
#include "atomfeii.h"
#include "monitor_results.h"
#include "thermal.h"
#include "wind.h"
#include "hmi.h"
#include "pressure.h"
#include "elementnames.h"
#include "ipoint.h"
#include "gammas.h"
#include "atmdat.h"
#include "hcmap.h"
#include "input.h"
#include "save.h"
#include "optimize.h"
#include "warnings.h"
#include "grid.h"

// this needs C linkage since it is passed as a pointer to a C library routine
extern "C" int wavelength_compare (const void * a, const void * b)
{
	LinSv a1 = *(LinSv*)a;
	LinSv b1 = *(LinSv*)b;
	/* comparison is b-a so we get inverse wavelength order (increasing energy order) */
	if( b1.wavelength > a1.wavelength )
		return 1;
	else if( b1.wavelength < a1.wavelength )
		return -1;
	else
		return 0;
}

// find strongest lines contributing to point in continuum energy mesh, output in some save commands
STATIC void FindStrongestLineLabels( void )
{
	long low_index=0;
	long high_index=0;
	long j_min = 0;
	double MaxFlux = 0.;
	long ipMaxFlux = 0;
	long j = 0;

	ASSERT( LineSave.ipass==1 );
	// copy and sort
	memcpy(LineSvSortWL, LineSv, (unsigned)LineSave.nsum*sizeof( LinSv ));
	qsort((void *)LineSvSortWL,(size_t)LineSave.nsum,sizeof( LinSv ),wavelength_compare);

	while( rfield.anu[j_min]+0.5*rfield.widflx[j_min] < RYDLAM/LineSvSortWL[0].wavelength)
		j_min++;

	for( j=0; j<rfield.nflux; j++ )
	{
		if( j < j_min )
		{
			strcpy( rfield.chLineLabel[j], "    " );
			continue;
		}

		ASSERT( LineSvSortWL[low_index].wavelength != 0. );

		while( RYDLAM/LineSvSortWL[low_index].wavelength < rfield.anu[j]-0.5*rfield.widflx[j] && low_index < LineSave.nsum-1 )
		{
			low_index++;
			if( LineSvSortWL[low_index].wavelength == 0. )
			{
				// hit the end of real wavelengths.  Pad rest of labels with spaces
				for( long j1=j; j1<rfield.nflux; j1++ )
					strcpy( rfield.chLineLabel[j1], "    " );
				return;
			}
		}

		high_index = low_index;
		ASSERT( LineSvSortWL[high_index].wavelength != 0. );

		while( RYDLAM/LineSvSortWL[high_index].wavelength < rfield.anu[j]+0.5*rfield.widflx[j] && high_index < LineSave.nsum-1 )
		{
			high_index++;
			if( LineSvSortWL[high_index].wavelength == 0. )
			{
				high_index--;
				break;
			}
		}
		// while loop found first one greater than j bin, decrement again to get back into j bin
		high_index--;

		ASSERT( LineSvSortWL[low_index].wavelength > 0. );
		ASSERT( LineSvSortWL[high_index].wavelength > 0. );
		ASSERT( RYDLAM/LineSvSortWL[low_index].wavelength >= rfield.anu[j]-0.5*rfield.widflx[j] );
		ASSERT( RYDLAM/LineSvSortWL[high_index].wavelength <= rfield.anu[j]+0.5*rfield.widflx[j] );

		MaxFlux = 0.;
		ipMaxFlux = 0;

		for( long k = low_index; k <= high_index; k++ )
		{
			if( strcmp( LineSvSortWL[k].chALab, "Coll" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Heat" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Pump" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "pump" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "nInu" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "nFnu" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "InwT" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "InwC" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Inwd" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Ca A" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Ca B" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "+Col" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Pho+" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Pcon" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Q(H)" )==0 ||
				strcmp( LineSvSortWL[k].chALab, "Unit" )==0 )
				continue;

			if( LineSvSortWL[k].SumLine[0] > MaxFlux )
			{
				MaxFlux = LineSvSortWL[k].SumLine[0];
				ipMaxFlux = k;
			}
		}

		/* line	label */
		if( ipMaxFlux > 0 )
			strcpy( rfield.chLineLabel[j], LineSvSortWL[ipMaxFlux].chALab );
	}

	return;
}

// implements the absorption option on the
// set save line width command
inline realnum PrettyTranmission(long j, realnum transmission)
{
	if( save.ResolutionAbs < realnum(0.) )
		// option to conserve energy
		return transmission;
	else
	{
		realnum corr = save.ResolutionAbs*rfield.widflx[j]/rfield.anu[j];
		return realnum(max(0., 1. - (1.-transmission)*corr ));
	}
}

/*SaveResults1Line do single line of output for the save results and save line intensity commands */
/* the number of emission lines across one line of printout */
STATIC void SaveResults1Line(
  FILE * ioPUN, 
  /* 4 char + null string */
  const char *chLab, 
  realnum wl, 
  double xInten, 
  const char *chFunction);/* null terminated string saying what to do */

/*SaveGaunts called by save gaunts command to output gaunt factors */
STATIC void SaveGaunts(FILE* ioPUN);

/*SaveResults save results from save results command */
/*SaveResults1Line do single line of output for the save results and save line intensity commands */
STATIC void SaveResults(FILE* ioPUN);

STATIC void SaveLineStuff(
  FILE * ioPUN,
  const char *chJob , 
  realnum xLimit);

/* save h emission, for chapter 4, routine is below */
STATIC void AGN_Hemis(FILE *ioPUN );

/*SaveNewContinuum produce the 'save new continuum' output */
STATIC void SaveNewContinuum(FILE * ioPUN );

/*SaveLineIntensity produce the 'save lines intensity' output */
STATIC void SaveLineIntensity(FILE * ioPUN , long int ipPun);

char *chDummy;

/*SaveFeII_cont  return appropriate FeII banded continuum, inward, outward, or total */
STATIC realnum SaveFeII_cont( long int ipCont , long ipFeII_Cont_type )
{

	DEBUG_ENTRY( "SaveFeII_cont()" );

	// [0]=wavelength, [1] inward; [2] outward, 3=> total
	if( ipFeII_Cont_type == 3 )
		return( FeII_Cont[ipCont][1]+FeII_Cont[ipCont][2] );
	else
		return(FeII_Cont[ipCont][ipFeII_Cont_type]);
}

void SaveDo(
	/* chTime is null terminated 4 char string, either "MIDL" or "LAST" */
	const char *chTime) 
{
	bool lgTlkSav;
	long int
	  ipPun, 
	  i,
	  i1, 
	  ion, 
	  ipConMax, 
	  ipHi,
	  ipLinMax, 
	  ipLo,
	  ips, 
	  j, 
	  jj, 
	  limit, 
	  nelem;
	double ConMax, 
	  RateInter, 
	  av, 
	  conem, 
	  eps, 
	  flxatt, 
	  flxcor, 
	  flxin, 
	  flxref, 
	  flxtrn, 
	  fout, 
	  fref, 
	  fsum, 
	  opConSum, 
	  opLinSum, 
	  stage, 
	  sum, 
	  texc, 
	  xLinMax;

	DEBUG_ENTRY( "SaveDo()" );

	/* 
	 * the "last" option on save command, to save on last iteration,
	 * is parsed at the top of the loop in only one place.  
	 * no further action is needed at all for save last to work
	 * ok throughout this routine 
	 */

	/* 
	 * each branch can have a test whether chTime is or is not "LAST"
	 *
	 * if( strcmp(chTime,"LAST") == 0 )  <== print after iteration is complete 
	 *
	 * if "LAST" then this is last call to routine after iteration complete
	 * save only if "LAST" when results at end of iteration are needed
	 *
	 * if( strcmp(chTime,"LAST") != 0 )  <== print for every zone 
	 *
	 * test for .not."LAST" is for every zone result, where you do not
	 * want to save last zone twice
	 */

	/* return if no save to do */
	if( save.nsave < 1 )
	{ 
		return;
	}

	/* during a grid calculation this routine saves grid points after
	 * cloudy is called.  we may output it below */
	if( grid.lgGrid )
	{
		if( chTime[0]=='L' )
			GridGatherInCloudy();
	}

	// sort line labels if this is last call, this avoids multiple calls if several
	// output options need sorted labels and is safer since labels will be sorted in
	// case new code is added that reports the strong lines.  The disadvantage is that
	// we sort even if the labels are not used
	if( strcmp(chTime,"LAST") == 0 )
	{
		// sort emission line intensities so strongest lines are reported
		FindStrongestLineLabels();
	}

	for( ipPun=0; ipPun < save.nsave; ipPun++ )
	{
		if( save.lgPunHeader[ipPun] && !nMatch(save.chHeader[ipPun],save.chNONSENSE) )
		{
			fprintf( save.ipPnunit[ipPun], "%s", save.chHeader[ipPun] );
			save.lgPunHeader[ipPun] = false;
		}

		/* this global variable to remember where in the save stack we are */
		save.ipConPun = ipPun;

		/* iterations.lgLastIt is true if this is last iteration
		 * lgPunLstIter set true if 'last' key occurred on save command
		 * normally is false.  This will skip saving if last set and
		 * this is not last iteration */
		if( iterations.lgLastIt || (!save.lgPunLstIter[ipPun]) )
		{

			if( strcmp(save.chSave[ipPun],"ABUN") == 0 )
			{
				/* save abundances vs depth */
				if( strcmp(chTime,"LAST") != 0 )
				{
					fprintf( save.ipPnunit[ipPun], "%.2f", 
						log10(MAX2(SMALLFLOAT,dense.gas_phase[ipHYDROGEN])) );
					for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
					{
						/* >>chng 05 feb 03, protect against non-positive abundances,
						 * bug caught by Marcelo Castellanos */
						fprintf( save.ipPnunit[ipPun], "\t%.2f", 
						  log10(MAX2(SMALLFLOAT,dense.gas_phase[nelem])) );
					}
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"21CM") == 0 )
			{
				/* save information about 21 cm line */
				if( strcmp(chTime,"LAST") != 0 )
				{
					fprintf( save.ipPnunit[ipPun], 
					  "%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
					  /* depth, cm */
					  radius.depth_mid_zone,
					  hyperfine.Tspin21cm ,
					  phycon.te ,
					  /* temperature from Lya - 21 cm optical depth ratio */
					  3.84e-7* Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauCon /
					  SDIV( HFLines[0].Emis->TauCon ),
					  /*TexcLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] ),*/
					  HFLines[0].Lo->Pop ,
					  HFLines[0].Hi->Pop ,
					  OccupationNumberLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] ),
					  HFLines[0].Emis->TauCon , 
					  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauCon,
					  HFLines[0].Emis->PopOpc,
					  /* term in () is density (cm-3) of 1s, so this is n(1s) / Ts */
					  (StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop)/SDIV( hyperfine.Tspin21cm),
					  /* why was above multiplied by this following term? */
					  /* *HFLines[0].EnergyErg/BOLTZMANN/4.,*/
					  HFLines[0].Emis->TauIn,
					  TexcLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] ) ,
					  colden.H0_ov_Tspin,
					  /*>>chng 27 mar, GS, integrated 21cm spin temperature*/
					  colden.H0_21cm_lower,
					  colden.H0_21cm_upper,
					  -0.068/log((colden.H0_21cm_upper/3.)/colden.H0_21cm_lower)
					  );
				}
			}

			else if( strcmp(save.chSave[ipPun],"AGES") == 0 )
			{
				/* save timescales vs depth */
				if( strcmp(chTime,"LAST") != 0 )
				{
					fprintf( save.ipPnunit[ipPun], "%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
					  /* depth, cm */
					  radius.depth_mid_zone,
					  /* cooling timescale */
					  dense.pden*BOLTZMANN*1.5*phycon.te/ thermal.htot, 
					  /* H2 destruction timescale */
					  timesc.time_H2_Dest_here, 
					  /* CO destruction timescale */
					  timesc.AgeCOMoleDest[findspecies("CO")->index], 
					  /* OH destruction timescale */
					  timesc.AgeCOMoleDest[findspecies("OH")->index], 
					  /* H recombination timescale */
					  1./(dense.eden*2.90e-10/(phycon.te70*phycon.te10/phycon.te03)) );
				}
			}

			else if( strcmp(save.chSave[ipPun]," AGN") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					if( strcmp( save.chSaveArgs[ipPun], "HECS" ) == 0 )
					{
						/* this routine is in helike.c */
						AGN_He1_CS(save.ipPnunit[ipPun]);
					}
					if( strcmp( save.chSaveArgs[ipPun], "HEMI" ) == 0 )
					{
						/* save h emiss, for chapter 4, routine is below */
						AGN_Hemis(save.ipPnunit[ipPun]);
					}
					else
					{
						fprintf( ioQQQ, " SaveDo does not recognize flag %4.4s set for AGN save.  This is impossible.\n", 
						  save.chSave[ipPun] );
						ShowMe();
						cdEXIT(EXIT_FAILURE);
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"MONI") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save the monitor output */
					lgCheckMonitors(save.ipPnunit[ipPun]);
				}
			}

			else if( strcmp(save.chSave[ipPun],"AVER") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save the averages output */
					save_average( ipPun );
				}
			}

			else if( strncmp(save.chSave[ipPun],"CHA",3) == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* one of the charge transfer options, all in chargtran.c */
					ChargTranPun( save.ipPnunit[ipPun] , save.chSave[ipPun] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"COOL") == 0 )
			{
				/* save cooling, routine in file of same name */
				if( strcmp(chTime,"LAST") != 0 )
					CoolSave(save.ipPnunit[ipPun]);
			}

			else if( strncmp(save.chSave[ipPun],"DYN" , 3) == 0 )
			{
				/* save dynamics xxx, information about dynamical solutions */
				if( strcmp(chTime,"LAST") != 0 )
					DynaSave(save.ipPnunit[ipPun] ,save.chSave[ipPun][3] );
			}

			else if( strcmp(save.chSave[ipPun],"ENTH") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
					fprintf( save.ipPnunit[ipPun],
						"%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
						radius.depth_mid_zone,
						phycon.EnthalpyDensity,
						phycon.EnergyExcitation,
						phycon.EnergyIonization,
						phycon.EnergyBinding ,
						0.5*POW2(wind.windv)*dense.xMassDensity ,	/* KE */
						5./2.*pressure.PresGasCurr ,				/* thermal plus PdV work */
						magnetic.EnthalpyDensity);						/* magnetic terms */
			}

			else if( strcmp(save.chSave[ipPun],"COLU") == 0 )
			{
				/* save column densities */
				if( strcmp(chTime,"LAST") == 0 )
				{
					PrtColumns(save.ipPnunit[ipPun]);
				}
			}

			else if( strcmp(save.chSave[ipPun],"COLS") == 0 )
			{
				/* save some column densities */
				if( strcmp(chTime,"LAST") == 0 )
				{
					save_colden( save.ipPnunit[ipPun] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"COMP") == 0 )
			{
				/* Compton energy exchange coefficients */
				if( strcmp(chTime,"LAST") != 0 )
				{
					for( jj=0; jj<rfield.nflux; jj = jj + save.ncSaveSkip)
					{
						fprintf( save.ipPnunit[ipPun], "%10.2e%10.2e%10.2e\n", 
						  rfield.anu[jj], rfield.comup[jj]/rfield.widflx[jj], 
						  rfield.comdn[jj]/rfield.widflx[jj] );
					}
				}
			}

			/* save continuum commands */
			else if( strcmp(save.chSave[ipPun],"CON ") == 0 )
			{
				/* this is the usual "save continuum" case */
				/* >>chng 06 apr 03, add every option to do every zone */
				/* if save every is set then nSaveEveryZone must be positive 
				 * was init to -1 */
				bool lgPrintThis =false;
				if( save.lgSaveEveryZone[ipPun] )
				{
					/* this branch, every option is on line so want to print every n zone */
					if( strcmp(chTime,"LAST") != 0 )
					{
						/* not last zone - print first and intermediate cases */
						if( nzone==1 )
						{
							lgPrintThis = true;
						}
						else if( nzone%save.nSaveEveryZone[ipPun]==0 )
						{
							lgPrintThis = true;
						}
					}
					else
					{
						/* this is last zone, print only if did not trip on above */
						if( nzone!=1 && nzone%save.nSaveEveryZone[ipPun]!=0 )
						{
							lgPrintThis = true;
						}
					}
				}
				else
				{
					/* this branch, not "every", so only print the last zone */
					if( strcmp(chTime,"LAST") == 0 )
						lgPrintThis = true;
				}
				ASSERT( !save.lgSaveEveryZone[ipPun] || save.nSaveEveryZone[ipPun]>0 );
				if( lgPrintThis )
				{
					if( save.lgSaveEveryZone[ipPun] && nzone!=1)
						fprintf( save.ipPnunit[ipPun], "%s\n",
							save.chHashString );

					/* option to also print same arrays but for cumulative arrays */
					int nEmType = (int)save.punarg[ipPun][0];
					ASSERT( nEmType==0 || nEmType==1 );

					for( j=0; j<rfield.nflux; j = j+save.ncSaveSkip)
					{
						/* four continua predicted here;
						 * incident, attenuated incident, emitted,
						 * then attenuated incident + emitted, last reflected
						 * reflected continuum is stored relative to inner radius
						 * others are stored for this radius */

						/* NB this code also used in save emitted,
						 * transmitted continuum commands */

						/* the incident continuum */
						flxin = rfield.flux_total_incident[nEmType][j]*rfield.anu2[j]*
						  EN1RYD/rfield.widflx[j];

						// a value < 0. indicates that energy should be conserved
						realnum resolution = ( save.Resolution < realnum(0.) ) ?
							rfield.anu[j]/rfield.widflx[j] : save.Resolution;

						/* the reflected continuum */
						flxref = (rfield.anu2[j]*(rfield.ConRefIncid[nEmType][j]+rfield.ConEmitReflec[nEmType][j])/rfield.widflx[j] +
							rfield.anu[j]*resolution*rfield.reflin[nEmType][j])*EN1RYD;

						/* the attenuated incident continuum */
						flxatt = rfield.flux[nEmType][j]*rfield.anu2[j]*EN1RYD/
						  rfield.widflx[j]*radius.r1r0sq * 
							PrettyTranmission( j, rfield.trans_coef_total[j] );

						/* the outward emitted continuum */
						conem = (rfield.ConEmitOut[nEmType][j]/
						  rfield.widflx[j]*rfield.anu2[j] + resolution*
						  rfield.outlin[nEmType][j]*rfield.anu[j])*radius.r1r0sq*
						  EN1RYD*geometry.covgeo;

						/* sum of emitted and transmitted continua */
						flxtrn = conem + flxatt;

						/* photon energy in appropriate energy or wavelength units */
						fprintf( save.ipPnunit[ipPun],"%.5e\t", AnuUnit(rfield.AnuOrg[j]) );
						/* incident continuum */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", flxin ); 
						/* trans cont */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", flxatt ); 
						/* DiffOut cont */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", conem ); 
						/* net trans cont */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", flxtrn ); 
						/* reflected cont */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", flxref ); 
						/* total cont */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", flxref + flxtrn );

						/* reflected lines */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", rfield.anu[j]*resolution*rfield.reflin[nEmType][j]*EN1RYD );

						/* outward lines */
						fprintf( save.ipPnunit[ipPun],"%.3e\t", rfield.anu[j]*resolution*rfield.outlin[nEmType][j]*EN1RYD*
							radius.r1r0sq*geometry.covgeo );

						fprintf( save.ipPnunit[ipPun], "%4.4s\t%4.4s\t", 
						/* line	label */
						  rfield.chLineLabel[j] ,
						/* cont label*/
						  rfield.chContLabel[j] );
						/* number of lines within that cell over cell width
						 * save raw continuum has number by itself */
						fprintf( save.ipPnunit[ipPun], "%.2f\n", rfield.line_count[j]/rfield.widflx[j]*rfield.anu[j] );
					}
				}
			}

			/* this is the save spectrum command, 
			 * the new continuum command that will replace the previous one */
			else if( strcmp(save.chSave[ipPun],"CONN") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
					SaveNewContinuum( save.ipPnunit[ipPun]);
			}

			else if( strcmp(save.chSave[ipPun],"CONC") == 0 )
			{
				/* save incident continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* incident continuum */
					for( j=0; j<rfield.nflux; j = j + save.ncSaveSkip)
					{
						flxin = rfield.flux_total_incident[0][j]*rfield.anu2[j]*
						  EN1RYD/rfield.widflx[j];
						/* >>chng 96 oct 22, format of anu to .5 to resolve energy mesh near 1 Ryd */
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.3e\n", 
						  AnuUnit(rfield.AnuOrg[j]), flxin );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONG") == 0 )
			{
				/* save emitted grain continuum in optically thin limit */
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* GrainMakeDiffuse broke out emission into types 
					 * according to matType */
					for( j=0; j<rfield.nflux; j = j + save.ncSaveSkip)
					{
						fsum = gv.GraphiteEmission[j]*rfield.anu2[j]*
						  EN1RYD/rfield.widflx[j] *radius.r1r0sq*geometry.covgeo;

						fout = gv.SilicateEmission[j]*rfield.anu2[j]*
						  EN1RYD/rfield.widflx[j] *radius.r1r0sq*geometry.covgeo;

						/* anu is .5e format to resolve energy mesh near 1 Ryd 
						 * AnuUnit gives anu in whatever units were set with units option */
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.3e\t%.3e\t%.3e\n", 
						  AnuUnit(rfield.AnuOrg[j]) , fsum , fout , 
						  /* total emission per unit volume defined in GrainMakeDiffuse
						   * used in RT_diffuse to form total grain emission */
						  gv.GrainEmission[j]*rfield.anu2[j]*
						  EN1RYD/rfield.widflx[j] *radius.r1r0sq*geometry.covgeo );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONR") == 0 )
			{
				/* save reflected continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */
				if( strcmp(chTime,"LAST") == 0 )
				{
					if( geometry.lgSphere )
					{
						fprintf( save.ipPnunit[ipPun], " Reflected continuum not predicted when SPHERE is set.\n" );
						fprintf( ioQQQ , 
							"\n\n>>>>>>>>>>>>>\n Reflected continuum not predicted when SPHERE is set.\n" );
						cdEXIT(EXIT_FAILURE);
					}

					for( j=0; j<rfield.nflux; j = j + save.ncSaveSkip)
					{
						// a value < 0. indicates that energy should be conserved
						realnum resolution = ( save.Resolution < realnum(0.) ) ?
							rfield.anu[j]/rfield.widflx[j] : save.Resolution;

						/* reflected continuum */
						flxref = rfield.anu2[j]*(rfield.ConRefIncid[0][j]+rfield.ConEmitReflec[0][j])/
						  rfield.widflx[j]*EN1RYD;
						/* reflected lines */
						fref = rfield.anu[j]*resolution*rfield.reflin[0][j]*EN1RYD;
						/* ratio of reflected to incident continuum, the albedo */
						if( rfield.flux_total_incident[0][j] > 1e-25 )
						{
							av = rfield.ConRefIncid[0][j]/rfield.flux_total_incident[0][j];
						}
						else
						{
							av = 0.;
						}
						/* >>chng 96 oct 22, format of anu to .5 to resolve energy mesh near 1 Ryd */
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4s\n", 
						  AnuUnit(rfield.AnuOrg[j]), flxref, fref, flxref + fref, 
						  av, rfield.chContLabel[j] );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CNVE") == 0 )
			{
				/* the save convergence error command */
				if( strcmp(chTime,"LAST") != 0 )
				{
					fprintf( save.ipPnunit[ipPun], 
						"%.5e\t%li\t%.4e\t%.4e\t%.4f\t%.4e\t%.4e\t%.3f\t%.4e\t%.4e\t%.4f\n", 
						radius.depth_mid_zone, 
						conv.nPres2Ioniz,
						pressure.PresTotlCorrect, 
						pressure.PresTotlCurr, 
						(pressure.PresTotlCorrect - pressure.PresTotlCurr)*100./pressure.PresTotlCorrect, 
						dense.EdenTrue,
						dense.eden,
						(dense.EdenTrue - dense.eden)*100./dense.EdenTrue,
						thermal.htot,
						thermal.ctot,
						(thermal.htot - thermal.ctot)*100./thermal.htot );
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONB") == 0 )
			{
				/* save continuum bins binning */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */
				if( strcmp(chTime,"LAST") != 0 )
				{
					for( j=0; j<rfield.nupper; j = j + save.ncSaveSkip)
					{
						fprintf( save.ipPnunit[ipPun], "%14.5e %14.5e %14.5e\n", 
						  AnuUnit(rfield.AnuOrg[j]), rfield.anu[j], rfield.widflx[j] );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"COND") == 0 )
			{
				ASSERT( fp_equal( save.punarg[ipPun][0] , (realnum)2. ) ||
					fp_equal( save.punarg[ipPun][0] , (realnum)1.) );

				/* save diffuse continuum the local line and continuous emission */
				if( strcmp(chTime,"LAST") == 0 && 
					fp_equal(save.punarg[ipPun][0] , (realnum)1.) )
				{
					/* this option to save diffuse emission for last zone
					 * as series of columns */
					for( j=0; j<rfield.nflux; j = j+save.ncSaveSkip)
					{
						// a value < 0. indicates that energy should be conserved
						realnum resolution = ( save.Resolution < realnum(0.) ) ?
							rfield.anu[j]/rfield.widflx[j] : save.Resolution;

						double EmisLin = resolution*EN1RYD*
							rfield.DiffuseLineEmission[j]*rfield.anu[j];
						double EmisCon = rfield.ConEmitLocal[nzone][j]*
							rfield.anu2[j]*EN1RYD/rfield.widflx[j]; 
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.5e\t%.5e\t%.5e\n", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  EmisCon ,
						  EmisLin , 
						  EmisCon+EmisLin);
					}
				}
				else if( strcmp(chTime,"LAST") != 0 && 
					fp_equal(save.punarg[ipPun][0] , (realnum)2.) )
				{
					/* diffuse emission per zone, with each a long row */
					static bool lgMustPrintHeader=true;
					if( lgMustPrintHeader )
					{
						lgMustPrintHeader = false;
						fprintf( save.ipPnunit[ipPun], "%.5e", 
							AnuUnit(rfield.AnuOrg[0]) );
						for( j=1; j<rfield.nflux; j = j+save.ncSaveSkip)
						{
							fprintf( save.ipPnunit[ipPun], "\t%.5e", 
								AnuUnit(rfield.AnuOrg[j]) );
						}
						fprintf( save.ipPnunit[ipPun], "\n" );
					}
					// a value < 0. indicates that energy should be conserved
					realnum resolution = ( save.Resolution < realnum(0.) ) ?
						rfield.anu[0]/rfield.widflx[0] : save.Resolution;
					double EmisLin = resolution*EN1RYD*
						rfield.DiffuseLineEmission[0]*rfield.anu[0];
					double EmisCon = rfield.ConEmitLocal[nzone][0]*
						rfield.anu2[0]*EN1RYD/rfield.widflx[0]; 
					fprintf( save.ipPnunit[ipPun], "%.5e", 
						EmisCon+EmisLin);
					for( j=1; j<rfield.nflux; j = j+save.ncSaveSkip)
					{
						// a value < 0. indicates that energy should be conserved
						resolution = ( save.Resolution < realnum(0.) ) ?
							rfield.anu[j]/rfield.widflx[j] : save.Resolution;
						double EmisLin = resolution*EN1RYD*
							rfield.DiffuseLineEmission[j]*rfield.anu[j];
						double EmisCon = rfield.ConEmitLocal[nzone][j]*
							rfield.anu2[j]*EN1RYD/rfield.widflx[j]; 
						fprintf( save.ipPnunit[ipPun], "\t%.5e", 
							EmisCon+EmisLin);
					}
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONE") == 0 )
			{
				/* save emitted continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save emitted continuum */
					for( j=0; j<rfield.nflux;j = j +save.ncSaveSkip)
					{
						// a value < 0. indicates that energy should be conserved
						realnum resolution = ( save.Resolution < realnum(0.) ) ?
							rfield.anu[j]/rfield.widflx[j] : save.Resolution;

						/* this is the reflected component */
						flxref = (rfield.anu2[j]*(rfield.ConRefIncid[0][j]+rfield.ConEmitReflec[0][j])/
						  rfield.widflx[j] + rfield.anu[j]*resolution*
						  rfield.reflin[0][j])*EN1RYD;

						/* this is the total emission in the outward direction */
						conem = (rfield.ConEmitOut[0][j])/rfield.widflx[j]*rfield.anu2[j] + 
							resolution*rfield.outlin[0][j]*rfield.anu[j];

						conem *= radius.r1r0sq*EN1RYD*geometry.covgeo;

						/* output: photon energy, reflected, outward, total emission
						 *  >>chng 96 oct 22, format of anu to .5e to resolve energy mesh near 1 Ryd */
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.3e\t%.3e\t%.3e\t%4.4s\t%4.4s\n", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  flxref, 
						  conem, 
						  flxref + conem, 
						  rfield.chLineLabel[j], 
						  rfield.chContLabel[j]
						   );
					}
				}
			}

			/* save fine continuum command */
			else if( strcmp(save.chSave[ipPun],"CONf") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					long nu_hi , nskip;
					if( save.punarg[ipPun][0] > 0. )
						/* possible lower bounds to energy range - 
						 * 0 if not set with range option*/
						j = ipFineCont( save.punarg[ipPun][0] );
					else
						j = 0;

					/* upper limit set with range option */
					if( save.punarg[ipPun][1]> 0. )
						nu_hi = ipFineCont( save.punarg[ipPun][1]);
					else
						nu_hi = rfield.nfine;

					/* number of cells to bring together, default is 10 */
					nskip = (long)save.punarg[ipPun][2];

					do
					{
						realnum sum1 = rfield.fine_opt_depth[j];
						realnum xnu = rfield.fine_anu[j];
						for( jj=1; jj<nskip; ++jj )
						{
							xnu += rfield.fine_anu[j+jj];
							sum1 += rfield.fine_opt_depth[j+jj];
						}
						fprintf( save.ipPnunit[ipPun], 
							"%.6e\t%.3e\n", 
							AnuUnit(xnu/nskip), 
							sexp(sum1/nskip) );
						j += nskip;
					} while( j < nu_hi );
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONi") == 0 )
			{
				/* save continuum interactions */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */

				/* continuum interactions */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* this is option to set lowest energy */
					if( save.punarg[ipPun][0] <= 0. )
					{
						i1 = 1;
					}
					else if( save.punarg[ipPun][0] < 100. )
					{
						i1 = ipoint(save.punarg[ipPun][0]);
					}
					else
					{
						i1 = (long int)save.punarg[ipPun][0];
					}

					fref = 0.;
					fout = 0.;
					fsum = 0.;
					sum = 0.;
					flxin = 0.;

					for( j=i1-1; j < rfield.nflux; j++ )
					{
						fref += rfield.flux[0][j]*opac.opacity_abs[j];
						fout += rfield.otslin[j]*opac.opacity_abs[j];
						fsum += rfield.otscon[j]*opac.opacity_abs[j];
						sum += rfield.ConInterOut[j]*opac.opacity_abs[j];
						flxin += (rfield.outlin[0][j] + rfield.outlin_noplot[j])*opac.opacity_abs[j];
					}
					fprintf( save.ipPnunit[ipPun], "%10.2e%10.2e%10.2e%10.2e%10.2e\n", 
					  fref, fout, fsum, sum, flxin );
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONI") == 0 )
			{
				/* save ionizing continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */

				if( save.lgSaveEveryZone[ipPun] || (strcmp(chTime,"LAST") == 0) )
				{
					/* this flag will remember whether we have ever printed anything */
					bool lgPrt=false;
					if( save.lgSaveEveryZone[ipPun] )
						fprintf(save.ipPnunit[ipPun],"#save every zone %li\n", nzone);

					/* save ionizing continuum command
					 * this is option to set lowest energy,
					 * if no number was entered then this was zero */
					if( save.punarg[ipPun][0] <= 0. )
						i1 = 1;
					else if( save.punarg[ipPun][0] < 100. )
						i1 = ipoint(save.punarg[ipPun][0]);
					else
						i1 = (long int)save.punarg[ipPun][0];

					sum = 0.;
					for( j=i1-1; j < rfield.nflux; j++ )
					{
						flxcor = rfield.flux[0][j] + 
						  rfield.otslin[j] + 
						  rfield.otscon[j] + 
						  rfield.ConInterOut[j] +
						  rfield.outlin[0][j] + rfield.outlin_noplot[j];

						sum += flxcor*opac.opacity_abs[j];
					}

					if( sum > 0. )
						sum = 1./sum;
					else
						sum = 1.;

					fsum = 0.;

					for( j=i1-1; j<rfield.nflux; ++j)
					{
						flxcor = rfield.flux[0][j] + 
						  rfield.otslin[j] + 
						  rfield.otscon[j] + 
						  rfield.ConInterOut[j]+
						  rfield.outlin[0][j] + rfield.outlin_noplot[j];

						fsum += flxcor*opac.opacity_abs[j];

						/* punched quantities are freq, flux, flux*cross sec,
						 * fraction of total, integral fraction of total */
						RateInter = flxcor*opac.opacity_abs[j]*sum;

						/* punage(ipPun,2) is lowest interaction rate to consider, def=0.01 (1 percent) */
						/* >>chng 01 nov 22, format to c-friendly */
						if( (RateInter >= save.punarg[ipPun][1]) && (flxcor > SMALLFLOAT) )
						{
							lgPrt = true;
							/* >>chng 96 oct 22, format of anu to 11.5 to resolve energy mesh near 1 Ryd */
							fprintf( save.ipPnunit[ipPun], 
								"%li\t%.5e\t%.2e\t%.2e\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.4s\t%.4s\n", 
							  j,
							  AnuUnit(rfield.AnuOrg[j]), 
							  flxcor, 
							  flxcor*opac.opacity_abs[j], 
							  rfield.flux[0][j]/flxcor, 
							  rfield.otslin[j]/flxcor, 
							  rfield.otscon[j]/flxcor, 
							  (rfield.outlin[0][j] + rfield.outlin_noplot[j])/flxcor, 
							  rfield.ConInterOut[j]/flxcor, 
							  RateInter, 
							  fsum*sum, 
							  rfield.chLineLabel[j], 
							  rfield.chContLabel[j] );
						}
					}
					if( !lgPrt )
					{
						/* entered logical block but did not print anything */
						fprintf(save.ipPnunit[ipPun],
							" punchdo, the PUNCH IONIZING CONTINUUM command "
							"did not find a strong point, sum and fsum were %.2e %.2e\n",
							sum,fsum);
						fprintf(save.ipPnunit[ipPun],
							" punchdo, the low-frequency energy was %.5e Ryd\n",
							rfield.anu[i1-1]);
						fprintf(save.ipPnunit[ipPun],
							" You can reset the threshold for the lowest fractional "
							"interaction to print with the second number of the save command\n"
							" The fraction was %.3f and this was too large.\n",
							save.punarg[ipPun][1]);
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONS") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					// continuum volume emissivity and opacity as a function of radius
					// command was "save continuum emissivity" */
					if( save.ipEmisFreq < 0 )
						save.ipEmisFreq = ipoint(save.emisfreq.Ryd());
					j = save.ipEmisFreq-1;

					fprintf( save.ipPnunit[ipPun], 
						 "%.14e\t%.14e\t%.5e\t%.5e\t%.5e\n", 
						 radius.Radius_mid_zone,
						 radius.depth_mid_zone,
						 rfield.anu2[j]*rfield.ConEmitLocal[nzone][j]/rfield.widflx[j]*EN1RYD, 
						 opac.opacity_abs[j], 
						 opac.opacity_sct[j] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"CORA") == 0 )
			{
				/* save raw continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */

				if( strcmp(chTime,"LAST") == 0 )
				{
					/* this option to save all raw ionizing continuum */
					for( j=0;j<rfield.nflux;j = j + save.ncSaveSkip)
					{
						fprintf( save.ipPnunit[ipPun], 
							"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%4.4s\t%4.4s\t", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  rfield.flux[0][j], 
						  rfield.otslin[j], 
						  rfield.otscon[j], 
						  rfield.ConRefIncid[0][j],
						  rfield.ConEmitReflec[0][j], 
						  rfield.ConInterOut[j],
						  rfield.outlin[0][j]+rfield.outlin_noplot[j], 
						  rfield.ConEmitOut[0][j] ,
						  rfield.chLineLabel[j], 
						  rfield.chContLabel[j]
						  );
						/* number of lines within that cell */
						fprintf( save.ipPnunit[ipPun], "%li\n", rfield.line_count[j] );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONo") == 0 )
			{
				/* save outward local continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */

				if( strcmp(chTime,"LAST") == 0 )
				{
					/* option to save out outward continuum here */
					for( j=0;j<rfield.nflux; j = j + save.ncSaveSkip)
					{
						fprintf( save.ipPnunit[ipPun], "%11.5e%10.2e%10.2e\n", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  rfield.ConEmitOut[0][j]+ rfield.outlin[0][j] + rfield.outlin_noplot[j], 
						  (rfield.flux[0][j] + rfield.otscon[j] + rfield.otslin[j] + 
						  rfield.ConInterOut[j])*opac.opacity_abs[j]*
						  rfield.anu[j] );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONO") == 0 )
			{
				/* save outward continuum */
				/* set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */

				if( strcmp(chTime,"LAST") == 0 )
				{
					/* option to save out outward continuum */
					for( j=0; j<rfield.nflux; j = j + save.ncSaveSkip)
					{
						// a value < 0. indicates that energy should be conserved
						realnum resolution = ( save.Resolution < realnum(0.) ) ?
							rfield.anu[j]/rfield.widflx[j] : save.Resolution;

						fprintf( save.ipPnunit[ipPun], "%11.5e%10.2e%10.2e%10.2e%10.2e\n", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  rfield.flux[0][j]*rfield.anu2[j]* EN1RYD/rfield.widflx[j]*radius.r1r0sq, 
						  rfield.ConInterOut[j]/rfield.widflx[j]*rfield.anu2[j]*radius.r1r0sq*EN1RYD, 
						  resolution*(rfield.outlin[0][j]+rfield.outlin_noplot[j])*rfield.anu[j]*radius.r1r0sq*EN1RYD, 
						  (rfield.ConInterOut[j]/rfield.widflx[j]*
						  rfield.anu2[j] + resolution*(rfield.outlin[0][j]+rfield.outlin_noplot[j])*
						  rfield.anu[j])*radius.r1r0sq*EN1RYD );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CONT") == 0 )
			{
				/* save transmitted continuum - this is not the main "save continuum"
				 * command - search on "CON " above 
				 * set pointer for possible change in units of energy in continuum
				 * AnuUnit will give anu in whatever units were set with save units */

				if( strcmp(chTime,"LAST") == 0 )
				{
					fprintf( save.ipPnunit[ipPun], "#\n" );
					fprintf( save.ipPnunit[ipPun], "%32ld # file format version number\n",
						 VERSION_TRNCON );
					fprintf( save.ipPnunit[ipPun], "%s # check 1\n",
						 continuum.mesh_md5sum.c_str() );
					union {
						double x;
						uint32 i[2];
					} u;
					u.x = double(rfield.emm);
					if( cpu.big_endian() )
						fprintf( save.ipPnunit[ipPun], "%23.8x %8.8x # check 2\n",
							 u.i[0], u.i[1] );
					else
						fprintf( save.ipPnunit[ipPun], "%23.8x %8.8x # check 2\n",
							 u.i[1], u.i[0] );
					u.x = double(rfield.egamry);
					if( cpu.big_endian() )
						fprintf( save.ipPnunit[ipPun], "%23.8x %8.8x # check 3\n",
							 u.i[0], u.i[1] );
					else
						fprintf( save.ipPnunit[ipPun], "%23.8x %8.8x # check 3\n",
							 u.i[1], u.i[0] );
					u.x = continuum.ResolutionScaleFactor;
					if( cpu.big_endian() )
						fprintf( save.ipPnunit[ipPun], "%23.8x %8.8x # check 4\n",
							 u.i[0], u.i[1] );
					else
						fprintf( save.ipPnunit[ipPun], "%23.8x %8.8x # check 4\n",
							 u.i[1], u.i[0] );
					fprintf( save.ipPnunit[ipPun], "%32ld # nflux\n",
						 (rfield.nflux+save.ncSaveSkip-1)/save.ncSaveSkip );
					fprintf( save.ipPnunit[ipPun], "#\n" );

					/* this option to save transmitted continuum */
					for( j=0; j < rfield.nflux; j += save.ncSaveSkip )
					{
						/* attenuated incident continuum
						 * >>chng 97 jul 10, remove SaveLWidth from this one only since
						 * we must conserve energy even in lines 
						 * >>chng 07 apr 26 include transmission coefficient */
						flxatt = rfield.flux[0][j]*rfield.anu2[j]*EN1RYD/
						  rfield.widflx[j]*radius.r1r0sq*rfield.trans_coef_total[j];

						/*conem = (rfield.ConOutNoInter[j] + rfield.ConInterOut[j]+rfield.outlin[0][j])*
						  rfield.anu2[j];
						conem *= radius.r1r0sq*EN1RYD*geometry.covgeo;*/
						/* >>chng 00 jan 03, above did not include all contributors.  
						 * Pasted in below from usual
						 * save continuum command */
						/* >>chng 04 jul 15, removed factor of save.SaveLWidth -
						 * this should not be there to conserve energy, as explained in hazy
						 * where command was documented, and in comment above.  caught by PvH */
						/* >>chng 04 jul 23, incorrect use of outlin - before multiplied by an2,
						 * quantity should be photons per Ryd, since init quantity is
						 * photons per cell.  Must div by widflx.  caught by PvH  */
						/*conem = (rfield.ConEmitOut[0][j]/rfield.widflx[j]*rfield.anu2[j] +
						  rfield.outlin[0][j]*rfield.anu[j])*radius.r1r0sq*
						  EN1RYD*geometry.covgeo;*/
						conem = (rfield.ConEmitOut[0][j] + rfield.outlin[0][j]) / rfield.widflx[j]*
							rfield.anu2[j]*radius.r1r0sq*EN1RYD*geometry.covgeo;

						flxtrn = conem + flxatt;

						/* use AnuOrg here instead of anu since probably
						 * going to be used by table read
						 * and we want good anu array for sanity check*/
						fprintf( save.ipPnunit[ipPun],"%.5e\t%.3e\t%.3e\n",
							 AnuUnit(rfield.AnuOrg[j]), flxtrn,
							 rfield.trans_coef_total[j] );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"CON2") == 0 )
			{
				/* save total two-pohoton continuum  */
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* this option to save diffuse continuum */
					for( j=0; j<rfield.nflux; j = j+save.ncSaveSkip)
					{
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.5e\t%.5e\n", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  rfield.TotDiff2Pht[j]/rfield.widflx[j] , 
						  rfield.TotDiff2Pht[j]*rfield.anu2[j]*EN1RYD/rfield.widflx[j]);
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSE") == 0 )
			{
				/* save grain extinction - includes only grain opacity, not total */
				if( strcmp(chTime,"LAST") != 0 )
				{
					fprintf( save.ipPnunit[ipPun], " %.5e\t", 
						radius.depth_mid_zone );

					/* visual extinction of an extended source (like a PDR)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\t" , rfield.extin_mag_V_extended);

					/* visual extinction of point source (star)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\n" , rfield.extin_mag_V_point);
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSO") == 0 )
			{
				/* save grain opacity */
				if( strcmp(chTime,"LAST") == 0 )
				{
					for( j=0; j < rfield.nflux; j++ )
					{
						double scat;
						fprintf( save.ipPnunit[ipPun], 
						  "%.5e\t%.2e\t%.2e\t%.2e\t", 
						  /* photon energy or wavelength */
						  AnuUnit(rfield.AnuOrg[j]), 
						  /* total opacity, discount forward scattering */
						  gv.dstab[j] + gv.dstsc[j], 
						  /* absorption opacity */
						  gv.dstab[j], 
						  /* scatter, with forward discounted */
						  gv.dstsc[j] );
						/* add together total scattering, discounting 1-g */
						scat = 0.;
						/* sum over all grain species */
						for( size_t nd=0; nd < gv.bin.size(); nd++ )
						{
							scat += gv.bin[nd]->pure_sc1[j]*gv.bin[nd]->dstAbund;
						}
						/* finally, scattering including effects of forward scattering */
						fprintf( save.ipPnunit[ipPun], 
						  "%.2e\t", scat );
						fprintf( save.ipPnunit[ipPun], 
						  "%.2e\n", gv.dstsc[j]/(gv.dstab[j] + gv.dstsc[j]) );
					}
				}
			}

			/* save grain abundance and save grain D/G ratio commands */
			else if( strcmp(save.chSave[ipPun],"DUSA") == 0 ||
				 strcmp(save.chSave[ipPun],"DUSD") == 0 )
			{
				bool lgDGRatio = ( strcmp(save.chSave[ipPun],"DUSD") == 0 );

				/* grain abundance */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeaderDRRatio = true,
						lgMustPrtHeaderAbundance=true;
					/* print grain header first if this has not yet been done */
					if( ( lgMustPrtHeaderDRRatio && lgDGRatio ) || 
					    ( lgMustPrtHeaderAbundance && !lgDGRatio ) )
					{
						/* only print one header for each case, but must
						 * track separately if both used in same sim */
						if( lgMustPrtHeaderDRRatio && lgDGRatio )
							lgMustPrtHeaderDRRatio = false;
						else if( lgMustPrtHeaderAbundance &&!lgDGRatio )
							lgMustPrtHeaderAbundance = false;

						fprintf( save.ipPnunit[ipPun], "#Depth" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\ttotal\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\txxxx\n" );
					}
					fprintf( save.ipPnunit[ipPun], " %.5e", 
						radius.depth_mid_zone );
					/* grain abundance per bin in g/cm^3 */
					double total = 0.;
					for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
					{
						double abund = gv.bin[nd]->IntVol*gv.bin[nd]->dustp[0]*
							gv.bin[nd]->cnv_H_pCM3;
						if( lgDGRatio )
							abund /= dense.xMassDensity;
						fprintf( save.ipPnunit[ipPun], "\t%.3e", abund );
						total += abund;
					}
					fprintf( save.ipPnunit[ipPun], "\t%.3e\n", total );
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSP") == 0 )
			{
				/* grain potential */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeader = true;
					/* do labels first if this is first zone */
					if( lgMustPrtHeader )
					{
						/* first print string giving grain id */
						fprintf( save.ipPnunit[ipPun], "#Depth" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* don't need to do this, ever again */
						lgMustPrtHeader = false;
					}
					fprintf( save.ipPnunit[ipPun], " %.5e", 
						radius.depth_mid_zone );
					/* grain potential in eV */
					for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
						fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->dstpot*EVRYD );
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSR") == 0 )
			{
				/* grain H2 formation rates */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeader = true;
					/* do labels first if this is first zone */
					if( lgMustPrtHeader )
					{
						/* first print string giving grain id */
						fprintf( save.ipPnunit[ipPun], "#Depth" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* don't need to do this, ever again */
						lgMustPrtHeader = false;
					}
					fprintf( save.ipPnunit[ipPun], " %.5e", 
						radius.depth_mid_zone );
					/* grain formation rate for H2 */
					for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
						fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->rate_h2_form_grains_used );
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUST") == 0 )
			{
				/* grain temperatures - K*/
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeader = true;
					/* do labels first if this is first zone */
					if( lgMustPrtHeader )
					{
						/* first print string giving grain id */
						fprintf( save.ipPnunit[ipPun], "#Depth" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* don't need to do this, ever again */
						lgMustPrtHeader = false;
					}
					fprintf( save.ipPnunit[ipPun], " %.5e", 
						radius.depth_mid_zone );
					for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
						fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->tedust );
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSC") == 0 )
			{
				/* save grain charge - eden from grains and 
				 * charge per grain in electrons / grain */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeader = true;
					/* do labels first if this is first zone */
					if( lgMustPrtHeader )
					{
						/* first print string giving grain id */
						fprintf( save.ipPnunit[ipPun], "#Depth\tne(grn)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)\tne\t" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* don't need to do this, ever again */
						lgMustPrtHeader = false;
					}

					fprintf( save.ipPnunit[ipPun], " %.5e\t%.4e", 
						radius.depth_mid_zone ,
						/* electron density contributed by grains, in e/cm^3, 
						 * positive number means grain supplied free electrons */
						gv.TotalEden );

					/* average charge per grain in electrons */
					for( size_t nd=0; nd < gv.bin.size(); ++nd )
					{
						fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AveDustZ );
					}
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSH") == 0 )
			{
				/* grain heating */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeader = true;
					/* save grain charge, but do labels first if this is first zone */
					if( lgMustPrtHeader )
					{
						/* first print string giving grain id */
						fprintf( save.ipPnunit[ipPun], "#Depth" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* don't need to do this, ever again */
						lgMustPrtHeader = false;
					}
					fprintf( save.ipPnunit[ipPun], " %.5e", 
						radius.depth_mid_zone );
					/* grain heating */
					for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
						fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->GasHeatPhotoEl );
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp(save.chSave[ipPun],"DUSV") == 0 )
			{
				/* grain drift velocities */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* used to print header exactly one time */
					static bool lgMustPrtHeader = true;
					/* save grain velocity, but do labels first if this is first zone */
					if( lgMustPrtHeader )
					{
						/* first print string giving grain id */
						fprintf( save.ipPnunit[ipPun], "#Depth" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%s", gv.bin[nd]->chDstLab );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* now print grain radius converting from cm to microns */
						fprintf( save.ipPnunit[ipPun], "#grn rad (mic)" );
						for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
							fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->AvRadius*1e4 );
						fprintf( save.ipPnunit[ipPun], "\n" );

						/* don't need to do this, ever again */
						lgMustPrtHeader = false;
					}
					fprintf( save.ipPnunit[ipPun], " %.5e", 
						radius.depth_mid_zone );
					/* grain drift velocity in km/s */
					for( size_t nd=0; nd < gv.bin.size(); ++nd ) 
						fprintf( save.ipPnunit[ipPun], "\t%.3e", gv.bin[nd]->DustDftVel*1e-5 );
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			/* >>chng 02 dec 30, separated scattering cross section and asymmetry factor, PvH */
			else if( strcmp(save.chSave[ipPun],"DUSQ") == 0 )
			{
				/* save grain Qs */
				if( strcmp(chTime,"LAST") == 0 )
				{
					for( j=0; j < rfield.nflux; j++ )
					{
						fprintf( save.ipPnunit[ipPun], "%11.4e", 
						  rfield.anu[j] );
						for( size_t nd=0; nd < gv.bin.size(); nd++ )
						{
							fprintf( save.ipPnunit[ipPun], "%9.1e%9.1e", 
							   gv.bin[nd]->dstab1[j]*4./gv.bin[nd]->IntArea,
							   gv.bin[nd]->pure_sc1[j]*gv.bin[nd]->asym[j]*4./gv.bin[nd]->IntArea );
						}
						fprintf( save.ipPnunit[ipPun], "\n" );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"ELEM") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					realnum renorm = 1.f;

					/* this is the index for the atomic number on the physical scale */
					/* >>chng 04 nov 23, use c scale throughout */
					nelem = (long int)save.punarg[ipPun][0];
					ASSERT( nelem >= ipHYDROGEN );

					/* don't do this if element is not turned on */
					if( dense.lgElmtOn[nelem] )
					{
						/* >>chng 04 nov 23, add density option, leave as cm-3 
						* default is still norm to total of that element */
						if( save.punarg[ipPun][1] == 0 )
							renorm = dense.gas_phase[nelem];

						fprintf( save.ipPnunit[ipPun], " %.5e", radius.depth_mid_zone );

						for( j=0; j <= (nelem + 1); ++j)
						{
							fprintf( save.ipPnunit[ipPun], "\t%.2e", 
							dense.xIonDense[nelem][j]/renorm );
						}
						if( nelem==ipHYDROGEN )
						{
							/* H2 */
							fprintf( save.ipPnunit[ipPun], "\t%.2e", 
								hmi.H2_total/renorm );
						}
						/* >>chng 04 nov 23 add C and O fine structure pops */
						else if( nelem==ipCARBON )
						{
							fprintf( save.ipPnunit[ipPun], "\t%.2e\t%.2e\t%.2e", 
								colden.C1Pops[0]/renorm, colden.C1Pops[1]/renorm, colden.C1Pops[2]/renorm);
							fprintf( save.ipPnunit[ipPun], "\t%.2e\t%.2e", 
								colden.C2Pops[0]/renorm, colden.C2Pops[1]/renorm);
							fprintf( save.ipPnunit[ipPun], "\t%.2e", 
								findspecies("CO")->hevmol/renorm );
						}
						else if( nelem==ipOXYGEN )
						{
							fprintf( save.ipPnunit[ipPun], "\t%.2e\t%.2e\t%.2e", 
								colden.O1Pops[0]/renorm, colden.O1Pops[1]/renorm, colden.O1Pops[2]/renorm);
						}
						fprintf( save.ipPnunit[ipPun], "\n" );
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"RECA") == 0 )
			{
				/* this will create table for AGN3 then exit, 
				 * routine is in makerecom.c */
				ion_recombAGN( save.ipPnunit[ipPun] );
				cdEXIT(EXIT_FAILURE);
			}

			else if( strcmp(save.chSave[ipPun],"RECE") == 0 )
			{
				/* save recombination efficiencies, 
				 * option turned on with the  "save recombination efficiencies" command
				 * output for the save recombination coefficients command is actually
				 * produced by a series of routines, as they generate the recombination
				 * coefficients.  these include 
				 * dielsupres, helium, hydrorecom, iibod, and makerecomb*/
				fprintf( save.ipPnunit[ipPun], 
					"%12.4e %12.4e %12.4e %12.4e\n", 
				  iso.RadRecomb[ipH_LIKE][ipHYDROGEN][0][ipRecRad], 
				  iso.RadRecomb[ipH_LIKE][ipHYDROGEN][0][ipRecNetEsc] ,
				  iso.RadRecomb[ipH_LIKE][ipHYDROGEN][2][ipRecRad],
				  iso.RadRecomb[ipH_LIKE][ipHYDROGEN][2][ipRecNetEsc]);
			}

			else if( strncmp(save.chSave[ipPun],"FEc" , 3 ) == 0 )
			{
				/* FeII continuum */
				if( strcmp(chTime,"LAST") == 0 )
				{

					if( save.punarg[ipPun][0] == 1 )
					{
						// row format used for grids - one time print of wavelength grid first
						char chTempHeader[100];
						long ipFeII_Cont_type;
						if( strcmp(save.chSave[ipPun],"FEcI") == 0 )
						{
							strcpy(chTempHeader , "#FeII inward: Wl(A)\tInt[erg cm-2 s-1]\n");
							ipFeII_Cont_type = 1;
						}
						else if( strcmp(save.chSave[ipPun],"FEcO") == 0 )
						{
							strcpy(chTempHeader , "#FeII outward: Wl(A)\tInt[erg cm-2 s-1]\n");
							ipFeII_Cont_type = 2;
						}
						else if( strcmp(save.chSave[ipPun],"FEcT") == 0 )
						{
							strcpy(chTempHeader , "#FeII total: Wl(A)\tInt[erg cm-2 s-1]\n");
							ipFeII_Cont_type = 3;
						}
						else
							TotalInsanity();
						
						if( save.lgPunHeader[ipPun] && 
							(!grid.lgGrid || optimize.nOptimiz==0) )
						{
							// one time print of FeII continuum header
							j = 0;
							fprintf( save.ipPnunit[ipPun], "%s%.5f", 
								chTempHeader ,
								AnuUnit((realnum)RYDLAM/((FeII_Cont[j][0]+FeII_Cont[j+1][0])/2.f) ));
							for( j=1; j < nFeIIConBins; j++ )
							{
								fprintf( save.ipPnunit[ipPun], "\t%.5e", 
									AnuUnit((realnum)RYDLAM/((FeII_Cont[j][0]+FeII_Cont[j+1][0])/2.f) ));
							}
							fprintf( save.ipPnunit[ipPun], "\n");
							save.lgPunHeader[ipPun] = false;
						}
						// give intensities
						fprintf( save.ipPnunit[ipPun], "%.2f", 
							SaveFeII_cont( 0 , ipFeII_Cont_type ));
						for( j=1; j < nFeIIConBins; j++ )
						{
								fprintf( save.ipPnunit[ipPun], "\t%e",
									SaveFeII_cont( j , ipFeII_Cont_type ) );
						}
						fprintf( save.ipPnunit[ipPun], "\n");
					}
					else if( save.punarg[ipPun][0] == 2 )
					{
						// the default, four columns, wl, total, inward, outward
						// one time print of header
						if( save.lgPunHeader[ipPun] && 
							(!grid.lgGrid || optimize.nOptimiz==0) )
						{
							fprintf( save.ipPnunit[ipPun], "FeiI wl, total, inward, outward\n");
							save.lgPunHeader[ipPun] = false;
						}
						for( j=0; j < nFeIIConBins; j++ )
						{
							fprintf( save.ipPnunit[ipPun], "%.5f",
								AnuUnit((realnum)RYDLAM/((FeII_Cont[j][0]+FeII_Cont[j+1][0])/2.f) ));
							fprintf( save.ipPnunit[ipPun], "\t%e",
								SaveFeII_cont( j , 3 ) );// total
							fprintf( save.ipPnunit[ipPun], "\t%e",
								SaveFeII_cont( j , 1 ) );// inward
							fprintf( save.ipPnunit[ipPun], "\t%e",
								SaveFeII_cont( j , 2 ) );// outward
							fprintf( save.ipPnunit[ipPun], "\n");
						}
					}
				}
			}

			/* save column densities */
			else if( strcmp(save.chSave[ipPun],"FENl") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save FeII level energies and stat weights, followed by column density */
					FeIIPunchColden( save.ipPnunit[ipPun] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"FE2l") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save FeII level energies and stat weights */
					FeIIPunchLevels( save.ipPnunit[ipPun] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"FEli") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save line intensities, routine is in atom_feii.c */
					FeIISaveLines( save.ipPnunit[ipPun] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"FE2o") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save line optical depths, routine is in atom_feii.c */
					FeIIPunchOpticalDepth( save.ipPnunit[ipPun] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"FRED") == 0 )
			{
				/* set with save Fred command, this punches some stuff from
				 * Fred Hamann's dynamics project */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* Fred's list */
					fprintf( save.ipPnunit[ipPun], "%.5e\t%.5e\t%.3e\t%.3e\t%.3e"
						"\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e"
						"\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e"
						"\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e"
						"\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e"
						"\t%.3e\t%.3e\n",
						radius.Radius, radius.depth ,wind.windv/1e5,
						wind.dvdr,
						dense.gas_phase[ipHYDROGEN], dense.eden , phycon.te,
						wind.AccelLine , wind.AccelCont ,
						wind.fmul , 
						// acceleration in this zone due to electron scattering,
						// if incident SED was not attenuated
						pressure.pinzon_PresIntegElecThin/dense.xMassDensity/radius.drad_x_fillfac ,
						mean.xIonMean[0][ipHYDROGEN][0][0] , mean.xIonMean[0][ipHYDROGEN][1][0] ,
						mean.xIonMean[0][ipHELIUM][0][0] , mean.xIonMean[0][ipHELIUM][1][0] ,
						mean.xIonMean[0][ipHELIUM][2][0] ,
						mean.xIonMean[0][ipCARBON][1][0] , mean.xIonMean[0][ipCARBON][2][0] ,
						mean.xIonMean[0][ipCARBON][3][0] ,
						mean.xIonMean[0][ipOXYGEN][0][0] , mean.xIonMean[0][ipOXYGEN][1][0] ,
						mean.xIonMean[0][ipOXYGEN][2][0] , mean.xIonMean[0][ipOXYGEN][3][0] ,
						mean.xIonMean[0][ipOXYGEN][4][0] , mean.xIonMean[0][ipOXYGEN][5][0] ,
						mean.xIonMean[0][ipOXYGEN][6][0] , mean.xIonMean[0][ipOXYGEN][7][0] ,
						dense.xIonDense[ipHYDROGEN][0] , dense.xIonDense[ipHYDROGEN][1] ,
						dense.xIonDense[ipHELIUM][0] , dense.xIonDense[ipHELIUM][1] ,
						dense.xIonDense[ipHELIUM][2] ,
						dense.xIonDense[ipCARBON][1] , dense.xIonDense[ipCARBON][2] ,
						dense.xIonDense[ipCARBON][3] ,
						dense.xIonDense[ipOXYGEN][0] , dense.xIonDense[ipOXYGEN][1] ,
						dense.xIonDense[ipOXYGEN][2] , dense.xIonDense[ipOXYGEN][3] ,
						dense.xIonDense[ipOXYGEN][4] , dense.xIonDense[ipOXYGEN][5] ,
						dense.xIonDense[ipOXYGEN][6] , dense.xIonDense[ipOXYGEN][7] ,
						mean.xIonMean[0][ipMAGNESIUM][1][0] , dense.xIonDense[ipMAGNESIUM][1],
						TauLines[ipT1032].Emis->TauIn ,
						TauLines[ipT1032].Emis->TauCon );
				}
			}

			else if( strcmp(save.chSave[ipPun],"FE2d") == 0 )
			{
				/* save some departure coefficients for large FeII atom */
				if( strcmp(chTime,"LAST") != 0 )
					FeIIPunDepart(save.ipPnunit[ipPun],false);
			}

			else if( strcmp(save.chSave[ipPun],"FE2D") == 0 )
			{
				/* save all departure coefficients for large FeII atom */
				if( strcmp(chTime,"LAST") != 0 )
					FeIIPunDepart(save.ipPnunit[ipPun],true);
			}

			else if( strcmp(save.chSave[ipPun],"FE2p") == 0 )
			{
				bool lgFlag = false;
				if( save.punarg[ipPun][2] )
					lgFlag = true;
				/* save small subset of level populations for large FeII atom */
				if( strcmp(chTime,"LAST") != 0 )
					FeIIPunPop(save.ipPnunit[ipPun],false,0,0,
					lgFlag);
			}

			else if( strcmp(save.chSave[ipPun],"FE2P") == 0 )
			{
				bool lgFlag = false;
				if( save.punarg[ipPun][2] )
					lgFlag = true;
				/* save range of level populations for large FeII atom */
				if( strcmp(chTime,"LAST") != 0 )
					FeIIPunPop(save.ipPnunit[ipPun],
					true,
					(long int)save.punarg[ipPun][0],
					(long int)save.punarg[ipPun][1],
					lgFlag);
			}

			/* save spectra in fits format */
			else if( strcmp(save.chSave[ipPun],"FITS") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					saveFITSfile( save.ipPnunit[ipPun], NUM_OUTPUT_TYPES );
				}
			}
			/* save gammas (but without element) */
			else if( strcmp(save.chSave[ipPun],"GAMt") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					long ns;
					/* save photoionization rates, with the PUNCH GAMMAS command */
					for( nelem=0; nelem < LIMELM; nelem++ )
					{
						for( ion=0; ion <= nelem; ion++ )
						{
							for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
							{
								fprintf( save.ipPnunit[ipPun], "%3ld%3ld%3ld%10.2e%10.2e%10.2e", 
									nelem+1, ion+1, ns+1, 
									ionbal.PhotoRate_Shell[nelem][ion][ns][0], 
									ionbal.PhotoRate_Shell[nelem][ion][ns][1] ,
									ionbal.PhotoRate_Shell[nelem][ion][ns][2] );

								for( j=0; j < t_yield::Inst().nelec_eject(nelem,ion,ns); j++ )
								{
									fprintf( save.ipPnunit[ipPun], "%5.2f",
										 t_yield::Inst().elec_eject_frac(nelem,ion,ns,j) );
								}
								fprintf( save.ipPnunit[ipPun], "\n" );
							}
						}
					}
				}
			}

			/* save gammas element, ion */
			else if( strcmp(save.chSave[ipPun],"GAMe") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					int ns;
					nelem = (long)save.punarg[ipPun][0];
					ion = (long)save.punarg[ipPun][1];
					/* valence shell */
					ns = Heavy.nsShells[nelem][ion]-1;
					/* show what some of the ionization sources are */
					GammaPrt( 
						opac.ipElement[nelem][ion][ns][0] , 
						opac.ipElement[nelem][ion][ns][1] , 
						opac.ipElement[nelem][ion][ns][2] , 
						save.ipPnunit[ipPun], 
						ionbal.PhotoRate_Shell[nelem][ion][ns][0] , 
						ionbal.PhotoRate_Shell[nelem][ion][ns][0]*0.1 );
				}
			}

			else if( strcmp(save.chSave[ipPun],"GAUN") == 0 )
			{
				/* save gaunt factors */
				if( strcmp(chTime,"LAST") != 0 )
					SaveGaunts(save.ipPnunit[ipPun]);
			}

			/* save parameters of grid calculation */
			else if( strcmp(save.chSave[ipPun],"GRID") == 0 )
			{
				/* save of all grid parameters */
				if( grid.lgGrid && strcmp(chTime,"LAST") == 0 )
				{
					if( optimize.nOptimiz == 0 )
					{
						/* start of line gives abort and warning summary */	
						fprintf(save.ipPnunit[ipPun], "#Index\tAbort?\tWarnings?" );
						/* print start of each variable command line */
						for( i=0; i< grid.nintparm; i++ )
						{
							char chStr[10];
							strncpy( chStr , optimize.chVarFmt[i] , 9 );
							/* make sure this small bit of string is terminated */
							chStr[9] = '\0';
							fprintf(save.ipPnunit[ipPun], "\t%s", chStr );
						}
						fprintf(save.ipPnunit[ipPun], "\tgrid parameter string\n" );
					}
					/* abort / warning summary for this sim */
					fprintf(save.ipPnunit[ipPun], "%9.9ld\t%c\t%c",
						optimize.nOptimiz,
						TorF(lgAbort) , 
						TorF(warnings.lgWarngs) );
					/* the grid parameters */
					char chGridParam[INPUT_LINE_LENGTH];
					char chStringHold[100];
					sprintf( chStringHold ,"%f",grid.interpParameters[optimize.nOptimiz][0] );
					strcpy( chGridParam , chStringHold );
					for( j=0; j< grid.nintparm; j++ )
					{
						if( j>0 )
						{
							sprintf( chStringHold ,", %f",grid.interpParameters[optimize.nOptimiz][j] );
							strcat( chGridParam , chStringHold );
						}
						fprintf(save.ipPnunit[ipPun], "\t%f",
							grid.interpParameters[optimize.nOptimiz][j] );
					}
					fprintf(save.ipPnunit[ipPun], "\t%s\n", chGridParam );
				}
			}

			else if( strcmp(save.chSave[ipPun],"HISp") == 0 )
			{
				/* save pressure history of current zone */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* note if pressure convergence failure occurred in history that follows */
					if( !conv.lgConvPres )
					{
						fprintf( save.ipPnunit[ipPun], 
							"#PROBLEM  Pressure not converged iter %li zone %li density-pressure follows:\n",
							iteration , nzone );
					}
					/* note if temperature convergence failure occurred in history that follows */
					if( !conv.lgConvTemp )
					{
						fprintf( save.ipPnunit[ipPun], 
							"#PROBLEM  Temperature not converged iter %li zone %li density-pressure follows:\n",
							iteration , nzone );
					}
					for( unsigned long k=0; k < conv.hist_pres_density.size(); ++k )
					{
						/* save history of density - pressure, with correct pressure */
						fprintf( save.ipPnunit[ipPun] , "%2li %4li\t%.5e\t%.5e\t%.5e\n",
							iteration,
							nzone,
							conv.hist_pres_density[k],
							conv.hist_pres_current[k],
							conv.hist_pres_correct[k]);
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"HISt") == 0 )
			{
				/* save temperature history of current zone */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* note if pressure convergence failure occurred in history that follows */
					if( !conv.lgConvPres )
					{
						fprintf( save.ipPnunit[ipPun], 
							"#PROBLEM  Pressure not converged iter %li zone %li temp heat cool follows:\n",
							iteration , nzone );
					}
					/* note if temperature convergence failure occurred in history that follows */
					if( !conv.lgConvTemp )
					{
						fprintf( save.ipPnunit[ipPun], 
							"#PROBLEM  Temperature not converged iter %li zone %li temp heat cool follows:\n",
							iteration , nzone );
					}
					for( unsigned long k=0; k < conv.hist_temp_temp.size(); ++k )
					{
						/* save history of density - pressure, with correct pressure */
						fprintf( save.ipPnunit[ipPun] , "%2li %4li\t%.5e\t%.5e\t%.5e\n",
							iteration,
							nzone,
							conv.hist_temp_temp[k],
							conv.hist_temp_heat[k],
							conv.hist_temp_cool[k]);
					}
				}
			}

			else if( strncmp(save.chSave[ipPun],"H2",2) == 0 )
			{
				/* all save info on large H2 molecule include H2 PDR pdr */
				H2_PunchDo( save.ipPnunit[ipPun] , save.chSave[ipPun] , chTime, ipPun );
			}

			else if( strcmp(save.chSave[ipPun],"HEAT") == 0 )
			{
				/* save heating, routine in file of same name */
				if( strcmp(chTime,"LAST") != 0 )
					SaveHeat(save.ipPnunit[ipPun]);
			}

			else if( strncmp(save.chSave[ipPun],"HE",2) == 0 )
			{
				/* various save helium commands */
				/* save helium line wavelengths */
				if( strcmp(save.chSave[ipPun] , "HELW") == 0 )
				{
					if( strcmp(chTime,"LAST") == 0 )
					{
						/* save helium & he-like wavelengths, first header */
						fprintf( save.ipPnunit[ipPun], 
							"Z\tElem\t2 1P->1 1S\t2 3P1->1 1S\t2 3P2->1 1S"
							"\t2 3S->1 1S\t2 3P2->2 3S\t2 3P1->2 3S\t2 3P0->2 3S" );
						fprintf( save.ipPnunit[ipPun], "\n" );
						for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
						{
							/* print element name, nuclear charge */
							fprintf( save.ipPnunit[ipPun], "%li\t%s", 
								nelem+1 , elementnames.chElementSym[nelem] );
							/*prt_wl print floating wavelength in Angstroms, in output format */
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2p1P][ipHe1s1S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2p3P1][ipHe1s1S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2p3P2][ipHe1s1S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2s3S][ipHe1s1S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2p3P2][ipHe2s3S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2p3P1][ipHe2s3S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\t" );
							prt_wl( save.ipPnunit[ipPun] , 
								Transitions[ipHE_LIKE][nelem][ipHe2p3P0][ipHe2s3S].WLAng );
							fprintf( save.ipPnunit[ipPun], "\n"); 
						}
					}
				}
				else
					TotalInsanity();
			}

			/* save hummer, results needed for Lya transport, to feed into David's routine */
			else if( strcmp(save.chSave[ipPun],"HUMM") == 0 )
			{
				eps = Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Aul/
				  (Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Coll.ColUL *dense.eden);
				fprintf( save.ipPnunit[ipPun], 
					" %.5e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
				  radius.depth_mid_zone, 
				  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauIn, 
				  StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop, 
				  StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop, 
				  phycon.te, Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->damp, eps );
			}

			else if( strncmp( save.chSave[ipPun] , "HYD", 3 ) == 0 )
			{
				/* various save hydrogen commands */
				if( strcmp(save.chSave[ipPun],"HYDc") == 0 )
				{
					if( strcmp(chTime,"LAST") != 0 )
					{
						/* save hydrogen physical conditions */
						fprintf( save.ipPnunit[ipPun], 
						  " %.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
						  radius.depth_mid_zone, phycon.te, dense.gas_phase[ipHYDROGEN], dense.eden, 
						  dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN], 
						  dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN], 
						  hmi.H2_total/dense.gas_phase[ipHYDROGEN], 
						  hmi.Hmolec[ipMH2p]/dense.gas_phase[ipHYDROGEN], 
						  hmi.Hmolec[ipMH3p]/dense.gas_phase[ipHYDROGEN], 
						  hmi.Hmolec[ipMHm]/dense.gas_phase[ipHYDROGEN] );
					}
				}

				else if( strcmp(save.chSave[ipPun],"HYDi") == 0 )
				{
					if( strcmp(chTime,"LAST") != 0 )
					{
						/* save hydrogen ionization
						 * this will be total decays to ground */
						RateInter = 0.;
						stage = iso.ColIoniz[ipH_LIKE][ipHYDROGEN][0]*dense.eden*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;
						fref = iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;
						fout = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop;
						/* 06 aug 28, from numLevels_max to _local. */
						for( ion=ipH2s; ion < iso.numLevels_local[ipH_LIKE][ipHYDROGEN]; ion++ )
						{
							/* this is total decays to ground */
							RateInter += 
								Transitions[ipH_LIKE][ipHYDROGEN][ion][ipH1s].Emis->Aul*
							  (Transitions[ipH_LIKE][ipHYDROGEN][ion][ipH1s].Emis->Pesc + 
							  Transitions[ipH_LIKE][ipHYDROGEN][ion][ipH1s].Emis->Pelec_esc + 
							  Transitions[ipH_LIKE][ipHYDROGEN][ion][ipH1s].Emis->Pdest);
							/* total photo from all levels */
							fref += iso.gamnc[ipH_LIKE][ipHYDROGEN][ion]*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ion].Pop;
							/* total col ion from all levels */
							stage += iso.ColIoniz[ipH_LIKE][ipHYDROGEN][ion]*dense.eden*
							  StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ion].Pop;
							fout += StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ion].Pop;
						}
						
						/* make these relative to parent ion */
						stage /= dense.xIonDense[ipHYDROGEN][1];
						fref /= dense.xIonDense[ipHYDROGEN][1];
						fout /= dense.xIonDense[ipHYDROGEN][1];

						fprintf( save.ipPnunit[ipPun], "hion\t%4ld\t%.2e\t%.2e\t%.2e", 
						  nzone, 
						  /* photo and collision ion rates have units s-1 */
						  iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s], 
						  iso.ColIoniz[ipH_LIKE][ipHYDROGEN][ipH1s]* dense.EdenHCorr,
						  ionbal.RateRecomTot[ipHYDROGEN][0] );

						fprintf( save.ipPnunit[ipPun], "\t%.2e", 
							iso.RadRec_caseB[ipH_LIKE][ipHYDROGEN] );

						fprintf( save.ipPnunit[ipPun], 
							"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
						  dense.xIonDense[ipHYDROGEN][1]/dense.xIonDense[ipHYDROGEN][0], 
						  iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]/(ionbal.RateRecomTot[ipHYDROGEN][0]), 
						  iso.RadRecomb[ipH_LIKE][ipHYDROGEN][1][ipRecEsc], 
						  RateInter, 
						  fref/MAX2(1e-37,fout), 
						  stage/MAX2(1e-37,fout), 
						  /* simple H+ */
						  safe_div( iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]*dense.xIonDense[ipHYDROGEN][0], dense.eden*dense.xIonDense[ipHYDROGEN][1] ),
						  secondaries.csupra[ipHYDROGEN][0]);

						GammaPrt(iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][0],rfield.nflux,iso.ipOpac[ipH_LIKE][ipHYDROGEN][ipH1s],
						  save.ipPnunit[ipPun],iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s],iso.gamnc[ipH_LIKE][ipHYDROGEN][ipH1s]*
						  0.05);
					}
				}

				else if( strcmp(save.chSave[ipPun],"HYDp") == 0 )
				{
					if( strcmp(chTime,"LAST") != 0 )
					{
						/* save hydrogen populations 
						 * first give total atom and ion density [cm-3]*/
						fprintf( save.ipPnunit[ipPun], "%.5e\t%.2e\t%.2e", 
						  radius.depth_mid_zone, 
						  dense.xIonDense[ipHYDROGEN][0], 
						  dense.xIonDense[ipHYDROGEN][1] );

						/* next give state-specific densities [cm-3] */
						for( j=ipH1s; j < iso.numLevels_local[ipH_LIKE][ipHYDROGEN]-1; j++ )
						{
							fprintf( save.ipPnunit[ipPun], "\t%.2e", 
								StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][j].Pop );
						}
						fprintf( save.ipPnunit[ipPun], "\n" );
					}
				}

				else if( strcmp(save.chSave[ipPun],"HYDl") == 0 )
				{
					if( strcmp(chTime,"LAST") == 0 )
					{
						/* save hydrogen line 
						 * gives intensities and optical depths */
						/** \todo	1	these hardwired numbers for indices
						 * can't be right */
						for( ipHi=4; ipHi<iso.numLevels_local[ipH_LIKE][ipHYDROGEN]-1; ++ipHi )
						{
							for( ipLo=ipHi-5; ipLo<ipHi; ++ipLo )
							{
								if( ipLo<0 )
									continue;
								fprintf(save.ipPnunit[ipPun], "%li\t%li\t%.4e\t%.2e\n",
									ipHi, ipLo,
									/* wl in cm */
									1./Transitions[ipH_LIKE][ipHYDROGEN][ipHi][ipLo].EnergyWN,
									Transitions[ipH_LIKE][ipHYDROGEN][ipHi][ipLo].Emis->TauIn );
							}
						}
					}
				}

				/* save hydrogen Lya - some details about Lya */
				else if( strcmp(save.chSave[ipPun],"HYDL") == 0 )
				{
					if( strcmp(chTime,"LAST") != 0 )
					{
						/* the population ratio for Lya */
						double popul = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop/SDIV(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop);
						/* the excitation temperature of Lya */
						texc = TexcLine( &Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s] );
						fprintf( save.ipPnunit[ipPun], 
						  "%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
						  radius.depth_mid_zone,
						  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauIn, 
						  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauTot, 
						  popul, 
						  texc, 
						  phycon.te, 
						  texc/phycon.te ,
						  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Pesc, 
						  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Pdest, 
						  Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->pump, 
						  opac.opacity_abs[Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].ipCont-1],
						  opac.albedo[Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].ipCont-1] );
					}
				}

				else if( strcmp(save.chSave[ipPun],"HYDr") == 0 )
				{
					/* save hydrogen recc - recombination cooling for AGN3 */
					TempChange(2500.f, false);
					while( phycon.te <= 20000. )
					{
						double r1;
						double ThinCoolingCaseB; 

						r1 = HydroRecCool(1,0);
						ThinCoolingCaseB = pow(10.,((-25.859117 + 
						  0.16229407*phycon.telogn[0] + 
						  0.34912863*phycon.telogn[1] - 
						  0.10615964*phycon.telogn[2])/(1. + 
						  0.050866793*phycon.telogn[0] - 
						  0.014118924*phycon.telogn[1] + 
						  0.0044980897*phycon.telogn[2] + 
						  6.0969594e-5*phycon.telogn[3])))/phycon.te;

						fprintf( save.ipPnunit[ipPun], " %10.2e\t", 
							phycon.te);
						fprintf( save.ipPnunit[ipPun], " %10.2e\t", 
							(r1+ThinCoolingCaseB)/(BOLTZMANN*phycon.te) );

						fprintf( save.ipPnunit[ipPun], " %10.2e\t", 
							r1/(BOLTZMANN*phycon.te));

						fprintf( save.ipPnunit[ipPun], " %10.2e\n", 
							ThinCoolingCaseB/(BOLTZMANN*phycon.te));

						TempChange(phycon.te *2.f , false);
					}
					/* must exit since we have disturbed the solution */
					fprintf(ioQQQ , "save agn now exits since solution is disturbed.\n");
					cdEXIT( EXIT_SUCCESS );
				}
				else
					TotalInsanity();
			}

			else if( strcmp(save.chSave[ipPun],"IONI") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save mean ionization distribution */
					PrtMeanIon( 'i', false , save.ipPnunit[ipPun] );
				}
			}

			/* save ionization rates */
			else if( strcmp(save.chSave[ipPun],"IONR") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* this is element number */
					nelem = (long)save.punarg[ipPun][0];
					fprintf( save.ipPnunit[ipPun], 
						"%.5e\t%.4e\t%.4e", 
						radius.depth_mid_zone,
						dense.eden ,
						dynamics.Rate);
					/* >>chng 04 oct 15, from nelem+2 to nelem+1 - array over read -
					 * caught by PnH */
					for( ion=0; ion<nelem+1; ++ion )
					{
						fprintf( save.ipPnunit[ipPun], 
							"\t%.4e\t%.4e\t%.4e\t%.4e", 
							dense.xIonDense[nelem][ion] ,
							ionbal.RateIonizTot(nelem,ion) ,
							ionbal.RateRecomTot[nelem][ion] ,
							dynamics.Source[nelem][ion] );
					}
					fprintf( save.ipPnunit[ipPun], "\n");
				}
			}

			else if( strcmp(save.chSave[ipPun]," IP ") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save valence shell ip's */
					for( nelem=0; nelem < LIMELM; nelem++ )
					{
						int ion_big;
						double energy;

						/* this is the largest number of ion stages per line */
						const int NELEM_LINE = 10;
						/* this loop in case all ions do not fit across page */
						for( ion_big=0; ion_big<=nelem; ion_big += NELEM_LINE )
						{
							int ion_limit = MIN2(ion_big+NELEM_LINE-1,nelem);

							/* new line then element name */
							fprintf( save.ipPnunit[ipPun], 
								"\n%2.2s", elementnames.chElementSym[nelem]);

							/* print ion stages across line */
							for( ion=ion_big; ion <= ion_limit; ++ion )
							{
								fprintf( save.ipPnunit[ipPun], "\t%4ld", ion+1 );
							}
							fprintf( save.ipPnunit[ipPun], "\n" );

							/* this loop is over all shells */
							ASSERT( ion_limit < LIMELM );
							/* upper limit is number of shells in atom */
							for( ips=0; ips < Heavy.nsShells[nelem][ion_big]; ips++ )
							{

								/* print shell label */
								fprintf( save.ipPnunit[ipPun], "%2.2s", Heavy.chShell[ips]);

								/* loop over possible ions */
								for( ion=ion_big; ion<=ion_limit; ++ion )
								{

									/* does this subshell exist for this ion? break if it does not*/
									/*if( Heavy.nsShells[nelem][ion]<Heavy.nsShells[nelem][0] )*/
									if( ips >= Heavy.nsShells[nelem][ion] )
										break;

									/* array elements are shell, numb of electrons, element, 0 */
									energy = t_ADfA::Inst().ph1(ips,nelem-ion,nelem,0);

									/* now print threshold with correct format */
									if( energy < 10. )
									{
										fprintf( save.ipPnunit[ipPun], "\t%6.3f", energy );
									}
									else if( energy < 100. )
									{
										fprintf( save.ipPnunit[ipPun], "\t%6.2f", energy );
									}
									else if( energy < 1000. )
									{
										fprintf( save.ipPnunit[ipPun], "\t%6.1f", energy );
									}
									else
									{
										fprintf( save.ipPnunit[ipPun], "\t%6ld",  (long)(energy) );
									}
								}

								/* put cs at end of long line */
								fprintf( save.ipPnunit[ipPun], "\n" );
							}
						}
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"LINC") == 0 )
			{
				/* save line cumulative */
				if( strcmp(chTime,"LAST") != 0 )
				{
					save_line(save.ipPnunit[ipPun],"PUNC",
						save.lgEmergent[ipPun]); 
				}
			}

			else if( strcmp(save.chSave[ipPun],"LIND") == 0 )
			{
				/* save line data, then stop */
				SaveLineData(save.ipPnunit[ipPun]);
			}

			else if( strcmp(save.chSave[ipPun],"LINL") == 0 )
			{
				/* save line labels */
				bool lgPrintAll=false;
				/* LONG keyword on save line labels command sets this to 1 */
				if( save.punarg[ipPun][0]>0. )
					lgPrintAll = true;
				prt_LineLabels(save.ipPnunit[ipPun] , lgPrintAll );
			}

			else if( strcmp(save.chSave[ipPun],"LINO") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save line optical depths, routine is below, file static */
					SaveLineStuff(save.ipPnunit[ipPun],"optical" , save.punarg[ipPun][0]);
				}
			}

			else if( strcmp(save.chSave[ipPun],"LINP") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					static bool lgFirst=true;
					/* save line populations, need to do this twice if very first
					 * call since first call to SaveLineStuff generates atomic parameters
					 * rather than level pops, routine is below, file static */
					SaveLineStuff(save.ipPnunit[ipPun],"populat" , save.punarg[ipPun][0]);
					if( lgFirst )
					{
						lgFirst = false;
						SaveLineStuff(save.ipPnunit[ipPun],"populat" , save.punarg[ipPun][0]);
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"LINS") == 0 )
			{
				/* save line emissivity */
				if( strcmp(chTime,"LAST") != 0 )
				{
					save_line(save.ipPnunit[ipPun],"PUNS",
						save.lgEmergent[ipPun]);
				}
			}

			else if( strcmp(save.chSave[ipPun],"LINR") == 0 )
			{
				/* save line RT */
				if( strcmp(chTime,"LAST") != 0 )
					Save_Line_RT( save.ipPnunit[ipPun]);
			}

			else if( strcmp(save.chSave[ipPun],"LINA") == 0 )
			{
				/* save line array */
				if( strcmp(chTime,"LAST") == 0 )
				{
					/* save out all lines with energies */
					for( j=0; j < LineSave.nsum; j++ )
					{
						if( LineSv[j].wavelength > 0. && 
							LineSv[j].SumLine[0] > 0. )
						{
							/* line energy, in units set with units option */
							fprintf( save.ipPnunit[ipPun], "%12.5e", 
							  AnuUnit((realnum)RYDLAM/LineSv[j].wavelength) );
							/* line label */
							fprintf( save.ipPnunit[ipPun], "\t%4.4s ",
								LineSv[j].chALab );
							/* wavelength */
							prt_wl( save.ipPnunit[ipPun], LineSv[j].wavelength );
							/* intrinsic intensity */
							fprintf( save.ipPnunit[ipPun], "\t%8.3f", 
								log10(SDIV(LineSv[j].SumLine[0]) ) + radius.Conv2PrtInten );
							/* emergent line intensity, r recombination  */
							fprintf( save.ipPnunit[ipPun], "\t%8.3f", 
								log10(SDIV(LineSv[j].SumLine[1]) ) + radius.Conv2PrtInten );
							/* type of line, i for info, etc */
							fprintf( save.ipPnunit[ipPun], " \t%c\n", 
							  LineSv[j].chSumTyp);
						}
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"LINI") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 && 
					(nzone/save.LinEvery)*save.LinEvery != nzone )
				{
					/* this is the last zone
					 * save line intensities - but do not do last zone twice */
					SaveLineIntensity(save.ipPnunit[ipPun] , ipPun );
				}
				else if( strcmp(chTime,"LAST") != 0 )
				{
					/* following so we only save first zone if LinEvery reset */
					if( (save.lgLinEvery && nzone == 1) || 
					  (nzone/save.LinEvery)*save.LinEvery == 
					  nzone )
					{
						/* this is middle of calculation
						 * save line intensities */
						SaveLineIntensity(save.ipPnunit[ipPun] , ipPun );
					}
				}
			}

			else if( strcmp( save.chSave[ipPun],"LEIL") == 0)
			{
				/* some line intensities for the Leiden PDR,
				 * but only do this when calculation is complete */
				if( strcmp(chTime,"LAST") == 0 )
				{
					double absval , rel;
					long int n;
					/* the lines we will find,
					 * for a sample list of PDR lines look at LineList_PDR_H2.dat
					 * in the cloudy data dir */
					/* the number of H2 lines */
					const int NLINE_H2 = 31; 
					/* the number of lines which are not H2 */
					const int NLINE_NOTH_H2 = 5; 
					/* the labels and wavelengths for the lines that are not H2 */
					char chLabel[NLINE_NOTH_H2][5]=
					{"C  2", "O  1", "O  1","C  1", "C  1" };
					double Wl[NLINE_NOTH_H2]=
					{157.6 , 63.17 , 145.5 ,609.2 , 369.7 ,  };
					/* these are wavelengths in microns, conv to Angstroms before call */
					/* >>chng 05 sep 06, many of following wavelengths updated to agree
					 * with output - apparently not updated when energies changed */
					double Wl_H2[NLINE_H2]=
					{2.121 ,
					 28.21 , 17.03 , 12.28 , 9.662 , 8.024 , 6.907 , 6.107 , 5.510 , 5.051 , 4.693 ,	
					 4.408 , 4.180 , 3.996 , 3.845 , 3.724 , 3.626 , 3.547 , 3.485 , 3.437 , 3.403 ,
					 3.381 , 3.368 , 3.365 , 3.371 , 3.387 , 3.410 , 3.441 , 3.485 , 3.542 , 3.603};
					/* print a header for the lines */
					for( n=0; n<NLINE_NOTH_H2; ++n )
					{
						fprintf(save.ipPnunit[ipPun], "%s\t%.2f",chLabel[n] , Wl[n]);
						/* get the line, non positive return says didn't find it */
						/* arguments are 4-char label, wavelength, return log total intensity, linear rel inten */
						if( cdLine( chLabel[n] , (realnum)(Wl[n]*1e4) , &absval , &rel ) <= 0 )
						{
							fprintf(save.ipPnunit[ipPun], " did not find\n");
						}
						else
						{
							fprintf(save.ipPnunit[ipPun], "\t%.3e\t%.3e\n",pow(10.,rel),absval);
						}
					}
					fprintf(save.ipPnunit[ipPun], "\n\n\n");

					/* only print the H2 lines if the big molecule is turned on */
					if( h2.lgH2ON )
					{
						fprintf(save.ipPnunit[ipPun], 
							"Here are some of the H2 Intensities, The first one is the\n"
							"1-0 S(0) line and the following ones are the 0-0 S(X)\n"
							"lines where X goes from 0 to 29\n\n");
						for( n=0; n<NLINE_H2; ++n )
						{
							fprintf(save.ipPnunit[ipPun], "%s\t%.3f","H2  " , Wl_H2[n]);
							/* get the line, non positive return says didn't find it */
							if( cdLine( "H2  " , (realnum)(Wl_H2[n]*1e4) , &absval , &rel ) <= 0 )
							{
								fprintf(save.ipPnunit[ipPun], " did not find\n");
							}
							else
							{
								fprintf(save.ipPnunit[ipPun], "\t%.3e\t%.3e\n",pow(10.,rel),absval);
							}
						}
					}
				}
			}

			else if( strcmp( save.chSave[ipPun],"LEIS") == 0)
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* get some column densities we shall need */
					double col_ci , col_oi , col_cii, col_heii;
					if( cdColm("carb" , 1 , &col_ci ) )
						TotalInsanity();
					if( cdColm("carb" , 2 , &col_cii ) )
						TotalInsanity();
					if( cdColm("oxyg" , 1 , &col_oi ) )
						TotalInsanity();
					if( cdColm("heli" , 2 , &col_heii ) )
						TotalInsanity();
					/* save Leiden structure - some numbers for the Leiden PDR model comparisons */
					fprintf( save.ipPnunit[ipPun], 
					"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t"
					"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t"
					"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t" 
					"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t"
					"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
					/* depth of this point */
					radius.depth_mid_zone,
					/* A_V for an extended source */
					0.00,	
					/* A_V for a point source */
					rfield.extin_mag_V_point,
					/* temperature */
					phycon.te ,
					dense.xIonDense[ipHYDROGEN][0],
					hmi.H2_total,
					dense.xIonDense[ipCARBON][0],
					dense.xIonDense[ipCARBON][1],
					dense.xIonDense[ipOXYGEN][0],
					findspecies("CO")->hevmol,
					findspecies("O2")->hevmol,
					findspecies("CH")->hevmol,
					findspecies("OH")->hevmol,
					dense.eden,
					dense.xIonDense[ipHELIUM][1],
					dense.xIonDense[ipHYDROGEN][1],
					hmi.Hmolec[ipMH3p],
					colden.colden[ipCOL_H0],
					colden.colden[ipCOL_H2g]+colden.colden[ipCOL_H2s],
					col_ci,
					col_cii,
					col_oi,
					findspecies("CO")->hevcol,
					findspecies("O2")->hevcol,
					findspecies("CH")->hevcol,
					findspecies("OH")->hevcol,
					colden.colden[ipCOL_elec],
					col_heii,
					colden.colden[ipCOL_Hp],
					colden.colden[ipCOL_H3p],
					hmi.H2_Solomon_dissoc_rate_used_H2g ,
					gv.rate_h2_form_grains_used_total,
					hmi.UV_Cont_rel2_Draine_DB96_depth,
					/* CO and C dissociation rate */
					CO_findrk("PHOTON,CO=>C,O"),
					/* total CI ionization rate */
					ionbal.PhotoRate_Shell[ipCARBON][0][2][0],
					/* total heating, erg cm-3 s-1 */
					thermal.htot,
					/* total cooling, erg cm-3 s-1 */
					thermal.ctot,
					/* GrnP grain photo heating */
					thermal.heating[0][13],
					/* grain collisional cooling */
					MAX2(0.,gv.GasCoolColl),	
					/* grain collisional heating */
					-1.*MIN2(0.,gv.GasCoolColl),	
					/* COds - CO dissociation heating */
					thermal.heating[0][9],
					/* H2dH-Heating due to H2 dissociation */
					hmi.HeatH2Dish_used,
					/* H2vH-Heating due to collisions with H2 */
					hmi.HeatH2Dexc_used ,
					/* ChaT - charge transfer heating */
					thermal.heating[0][24] ,
					/* cosmic ray heating */
					thermal.heating[1][6] ,
					/* heating due to atoms of various heavy elements */
					thermal.heating[ipMAGNESIUM][0],
					thermal.heating[ipSULPHUR][0],
					thermal.heating[ipSILICON][0],
					thermal.heating[ipIRON][0],
					thermal.heating[ipSODIUM][0],
					thermal.heating[ipALUMINIUM][0],
					thermal.heating[ipCARBON][0],
					TauLines[ipT610].Coll.cool,
					TauLines[ipT370].Coll.cool,
					TauLines[ipT157].Coll.cool,
					TauLines[ipT63].Coll.cool,
					TauLines[ipT146].Coll.cool );
				}
			}

			else if( strcmp( save.chSave[ipPun],"LLST") == 0)
			{
				/* save linelist command - do on last iteration */
				if( strcmp(chTime,"LAST") == 0 )
				{
					fprintf( save.ipPnunit[ipPun], 
						"%li" , iteration );

					/* -1 is flag saying that this save command was not set */
					if( save.nLineList[ipPun] < 0 )
						TotalInsanity();

					/* loop over all lines in the file we read */
					for( j=0; j<save.nLineList[ipPun]; ++j )
					{
						double relative , absolute, PrtQuantity;
						if( (cdLine( save.chLineListLabel[ipPun][j] , 
							save.wlLineList[ipPun][j]  , 
							&relative , &absolute , save.lgEmergent[ipPun]) ) <=0 )
						{
							if( !h2.lgH2ON && strncmp( save.chLineListLabel[ipPun][j] , "H2  " , 4 )==0 )
							{
								static bool lgMustPrintFirstTime = true;
								if( lgMustPrintFirstTime )
								{
									/* it's an H2 line and H2 is not being done - ignore it */
									fprintf( ioQQQ,"Did not find an H2 line, the large model is not "
										"included, so I will ignore it.  Log intensity set to -30.\n" );
									fprintf( ioQQQ,"I will totally ignore any future missed H2 lines\n");
									lgMustPrintFirstTime = false;
								}
								relative = -30.f;
								absolute = -30.f;
							}
							else if( lgAbort )
							{
								/* we are in abort mode */
								relative = -30.f;
								absolute = -30.f;
							}
							else
							{
								fprintf(ioQQQ,"DISASTER - did not find a line in the Line List table\n");
								cdEXIT(EXIT_FAILURE);
							}
						}

						/* options to do either relative or absolute intensity
						 * default is relative, is absolute keyword on line then
						 * punarg set to 1 */
						/* straight line intensities */
						if( save.punarg[ipPun][0] > 0 )
							PrtQuantity = absolute;
						else
							PrtQuantity = relative;
						/* if taking ratio print every other line as ratio
						 * with previous line */
						if( save.lgLineListRatio[ipPun] )
						{
							/* do line pair ratios */
							static double SaveQuantity = 0;
							if( is_odd(j) )
								fprintf( save.ipPnunit[ipPun], "\t%.4e" , 
									SaveQuantity / SDIV( PrtQuantity ) );
							else
								SaveQuantity = PrtQuantity;
						}
						else
						{
							fprintf( save.ipPnunit[ipPun], "\t%.4e" , PrtQuantity );
						}
					}
					fprintf( save.ipPnunit[ipPun], "\n" );
				}
			}

			else if( strcmp( save.chSave[ipPun],"CHRT") == 0)
			{
				if ('\0' == save.chHeader[ipPun][1]) 
				{
					CO_punch_mol(NULL,save.optname[ipPun].c_str(),save.chHeader[ipPun],0.0);
					fprintf( save.ipPnunit[ipPun], "%s", save.chHeader[ipPun]+1 ); // Leave out '#' already printed
				}
				/* save chemistry CO rates command */
				if( strcmp(chTime,"LAST") != 0 )
				{
					CO_punch_mol(save.ipPnunit[ipPun],save.optname[ipPun].c_str(),NULL,radius.depth_mid_zone);
				}
			}

			else if( strcmp( save.chSave[ipPun],"RCO ") == 0)
			{
				/* save chemistry CO rates command */
				if( strcmp(chTime,"LAST") != 0 )
				{
					CO_punch_mol(save.ipPnunit[ipPun],"CO",NULL,radius.depth_mid_zone);
				}
			}

			/*>>chng 06 jul 13, Nick Abel add OH rates */
			else if( strcmp( save.chSave[ipPun],"ROH ") == 0)
			{
				/* save chemistry OH rates command */
				if( strcmp(chTime,"LAST") != 0 )
				{
					CO_punch_mol(save.ipPnunit[ipPun],"OH",NULL,radius.depth_mid_zone);
				}
			}
			else if( strcmp(save.chSave[ipPun],"MAP ") == 0 )
			{
				/* do the map now if we are at the zone, or if this
				 * is the LAST call to this routine and map not done yet */
				if(  !hcmap.lgMapDone &&
					(nzone == hcmap.MapZone  ||  strcmp(chTime,"LAST") == 0 ) )
				{
					lgTlkSav = called.lgTalk;
					called.lgTalk = cpu.lgMPI_talk();
					hcmap.lgMapBeingDone = true;
					map_do(save.ipPnunit[ipPun] , " map");
					called.lgTalk = lgTlkSav;
				}
			}

			else if( strcmp(save.chSave[ipPun],"MOLE") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					static bool lgMustPrintHeader = true;
					if( lgMustPrintHeader )
					{
						lgMustPrintHeader = false;
						fprintf( save.ipPnunit[ipPun], 
							"#depth\tAV(point)\tAV(extend)\tCO diss rate\tC recom rate");

						for(i=0; i<N_H_MOLEC; ++i )
						{
							fprintf( save.ipPnunit[ipPun], "\t%s", hmi.chLab[i] );
						}
						for(i=0; i<mole.num_comole_calc; ++i )
						{
							fprintf( save.ipPnunit[ipPun], "\t%s",COmole[i]->label );
						}
						fprintf( save.ipPnunit[ipPun], "\n");
					}
					/* molecules, especially for PDR, first give radius */
					fprintf( save.ipPnunit[ipPun], "%.5e\t" , radius.depth_mid_zone );

					/* visual extinction of point source (star)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\t" , rfield.extin_mag_V_point);

					/* visual extinction of an extended source (like a PDR)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\t" , rfield.extin_mag_V_extended);

					/* carbon monoxide photodissociation rate */
					fprintf( save.ipPnunit[ipPun], "%.5e\t" , CO_findrk("PHOTON,CO=>C,O") );

					/* carbon recombination rate */
					fprintf( save.ipPnunit[ipPun], "%.5e" , ionbal.RateRecomTot[ipCARBON][0] );

					/* now do all the H molecules */
					for(j=0; j<N_H_MOLEC; ++j )
					{
						fprintf(save.ipPnunit[ipPun],"\t%.2e",hmi.Hmolec[j] );
					}

					/* now all the C molecules */
					for(j=0; j<mole.num_comole_calc; ++j )
					{
						fprintf(save.ipPnunit[ipPun],"\t%.2e",COmole[j]->hevmol );
					}
					fprintf(save.ipPnunit[ipPun],"\n");
				}
			}

			else if( strcmp(save.chSave[ipPun],"OPAC") == 0 )
			{
				/* save opacity- routine will parse which type of opacity save to do */
				if(  save.lgSaveEveryZone[ipPun] || strcmp(chTime,"LAST") == 0 )
					save_opacity(save.ipPnunit[ipPun],ipPun);
			}

			/* save coarse optical depths command */
			else if( strcmp(save.chSave[ipPun],"OPTc") == 0 )
			{
				if( save.lgSaveEveryZone[ipPun] || strcmp(chTime,"LAST") == 0 )
				{
					for( j=0; j < rfield.nflux; j++ )
					{
						fprintf( save.ipPnunit[ipPun], 
							"%13.5e\t%.3e\t%12.4e\t%.3e\n", 
						  AnuUnit(rfield.AnuOrg[j]), 
						  opac.TauAbsFace[j]+opac.TauScatFace[j], 
						  opac.TauAbsFace[j], 
						  opac.TauScatFace[j] );
					}
				}
			}

			/* save fine optical depths command */
			else if( strcmp(save.chSave[ipPun],"OPTf") == 0 )
			{
				if( save.lgSaveEveryZone[ipPun] || strcmp(chTime,"LAST") == 0 )
				{
					long nu_hi , nskip;
					if( save.punarg[ipPun][0] > 0. )
						/* possible lower bounds to energy range - will be zero if not set */
						j = ipFineCont( save.punarg[ipPun][0] );
					else
						j = 0;

					/* upper limit */
					if( save.punarg[ipPun][1]> 0. )
						nu_hi = ipFineCont( save.punarg[ipPun][1]);
					else
						nu_hi = rfield.nfine;

					/* we will bring nskip cells together into one printed
					 * number to make output smaller - default is 10 */
					nskip = (long)save.punarg[ipPun][2];
					do
					{
						realnum sum1 = rfield.fine_opt_depth[j];
						realnum sum2 = rfield.fine_opac_zone[j];
						/* want to report the central wavelength of the cell */
						realnum xnu = rfield.fine_anu[j];
						for( jj=1; jj<nskip; ++jj )
						{
							sum1 += rfield.fine_opt_depth[j+jj];
							sum2 += rfield.fine_opac_zone[j+jj];
							xnu += rfield.fine_anu[j+jj];
						}
						if( sum2 > 0. )
							fprintf( save.ipPnunit[ipPun], 
							  "%12.6e\t%.3e\t%.3e\n", 
							  AnuUnit(xnu/nskip), 
							  sum1/nskip , 
							  sum2/nskip);
						j += nskip;
					}while( j < nu_hi );
				}
			}

			else if( strcmp(save.chSave[ipPun]," OTS") == 0 )
			{
				ConMax = 0.;
				xLinMax = 0.;
				opConSum = 0.;
				opLinSum = 0.;
				ipLinMax = 1;
				ipConMax = 1;

				for( j=0; j < rfield.nflux; j++ )
				{
					opConSum += rfield.otscon[j]*opac.opacity_abs[j];
					opLinSum += rfield.otslin[j]*opac.opacity_abs[j];
					if( rfield.otslin[j]*opac.opacity_abs[j] > xLinMax )
					{
						xLinMax = rfield.otslin[j]*opac.opacity_abs[j];
						ipLinMax = j+1;
					}
					if( rfield.otscon[j]*opac.opacity_abs[j] > ConMax )
					{
						ConMax = rfield.otscon[j]*opac.opacity_abs[j];
						ipConMax = j+1;
					}
				}
				fprintf( save.ipPnunit[ipPun], 
				  "tot con lin=%.2e%.2e lin=%.4s%.4e%.2e con=%.4s%.4e%.2e\n", 
				  opConSum, opLinSum, rfield.chLineLabel[ipLinMax-1]
				  , rfield.anu[ipLinMax-1], xLinMax, rfield.chContLabel[ipConMax-1]
				  , rfield.anu[ipConMax-1], ConMax );
			}

			else if( strcmp(save.chSave[ipPun],"OVER") == 0 )
			{
				/* save overview
				 * this is the floor for the smallest ionization fractions printed */
				double toosmall = SMALLFLOAT ,
					hold;

				/* overview of model results,
				 * depth, te, hden, eden, ion fractions H, He, c, O */
				if( strcmp(chTime,"LAST") != 0 )
				{

					/* print the depth */
					fprintf( save.ipPnunit[ipPun], "%.5e\t", radius.depth_mid_zone );

					/* temperature, heating */
					if(dynamics.Cool() > dynamics.Heat()) 
					{
						fprintf( save.ipPnunit[ipPun], "%.4f\t%.3f", 
						  log10(phycon.te), log10(SDIV(thermal.htot-dynamics.Heat())) ); 
					}
					else
					{
						double diff = fabs(thermal.htot-dynamics.Cool());
						fprintf( save.ipPnunit[ipPun], "%.4f\t%.3f", 
						  log10(phycon.te), log10(SDIV(diff)) ); 
					}

					/* hydrogen and electron densities */
					fprintf( save.ipPnunit[ipPun], "\t%.4f\t%.4f", 
					  log10(dense.gas_phase[ipHYDROGEN]), log10(dense.eden) );

					/* molecular fraction of hydrogen */
					fprintf( save.ipPnunit[ipPun], "\t%.4f", 
					  /*log10(MAX2(toosmall,2.*hmi.Hmolec[ipMH2g]/dense.gas_phase[ipHYDROGEN])) );*/
					  log10(MAX2(toosmall,2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN])) );

					/* ionization fractions of hydrogen */
					fprintf( save.ipPnunit[ipPun], "\t%.4f\t%.4f", 
					  log10(MAX2(toosmall,dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN])), 
					  log10(MAX2(toosmall,dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN])) );

					/* ionization fractions of helium */
					for( j=1; j <= 3; j++ )
					{
						double arg1 = SDIV(dense.gas_phase[ipHELIUM]);
						arg1 = MAX2(toosmall,dense.xIonDense[ipHELIUM][j-1]/arg1 );
						fprintf( save.ipPnunit[ipPun], "\t%.4f", 
						  log10(arg1) );
					}

					/* carbon monoxide molecular fraction of CO */
					hold = SDIV(dense.gas_phase[ipCARBON]);
					hold = findspecies("CO")->hevmol/hold;
					hold = MAX2(toosmall, hold );
					fprintf( save.ipPnunit[ipPun], "\t%.4f", log10(hold) );

					/* ionization fractions of carbon */
					for( j=1; j <= 4; j++ )
					{
						hold = SDIV(dense.gas_phase[ipCARBON]);
						hold = MAX2(toosmall,dense.xIonDense[ipCARBON][j-1]/hold);
						fprintf( save.ipPnunit[ipPun], "\t%.4f", 
						  log10(hold) );
					}

					/* ionization fractions of oxygen */
					for( j=1; j <= 6; j++ )
					{
						hold = SDIV(dense.gas_phase[ipOXYGEN]);
						hold = MAX2(toosmall,dense.xIonDense[ipOXYGEN][j-1]/hold);
						fprintf( save.ipPnunit[ipPun], "\t%.4f", 
						  log10(hold) );
					}

					/* visual extinction of point source (star)*/
					fprintf( save.ipPnunit[ipPun], "\t%.2e" , rfield.extin_mag_V_point);

					/* visual extinction of an extended source (like a PDR)*/
					fprintf( save.ipPnunit[ipPun], "\t%.2e\n" , rfield.extin_mag_V_extended);
				}
			}

			else if( strcmp(save.chSave[ipPun]," PDR") == 0 )
			{
				/* this is the save PDR command */
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* convert optical depth at wavelength of V filter
					 * into magnitudes of extinction */
					/* >>chyng 03 feb 25, report extinction to illuminated face,
					 * rather than total extinction which included far side when
					 * sphere was set */
					/*av = opac.TauTotalGeo[0][rfield.ipV_filter-1]*1.08574;*/

					fprintf( save.ipPnunit[ipPun], 
						"%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t", 
					  radius.depth_mid_zone, 
					  /* total hydrogen column density, all forms */
					  colden.colden[ipCOL_HTOT], 
					  phycon.te, 
					  /* fraction of H that is atomic */
					  dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN], 
					  /* ratio of n(H2) to total H, == 0.5 when fully molecular */
					  2.*hmi.Hmolec[ipMH2g]/dense.gas_phase[ipHYDROGEN], 
					  2.*hmi.Hmolec[ipMH2s]/dense.gas_phase[ipHYDROGEN], 
					  /* atomic to total carbon */
					  dense.xIonDense[ipCARBON][0]/SDIV(dense.gas_phase[ipCARBON]), 
					  findspecies("CO")->hevmol/SDIV(dense.gas_phase[ipCARBON]), 
					  findspecies("H2O")->hevmol/SDIV(dense.gas_phase[ipOXYGEN]),
					  /* hmi.UV_Cont_rel2_Habing_TH85 is field relative to Habing background, dimensionless */
					  hmi.UV_Cont_rel2_Habing_TH85_depth);

					/* visual extinction due to dust alone, of point source (star)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\t" , rfield.extin_mag_V_point);

					/* visual extinction due to dust alone,  of an extended source (like a PDR)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\t" , rfield.extin_mag_V_extended);

					/* visual extinction (all sources) of a point source (like a PDR)*/
					fprintf( save.ipPnunit[ipPun], "%.2e\n" , opac.TauAbsGeo[0][rfield.ipV_filter] );
				}
			}

			else if( strcmp(save.chSave[ipPun],"PHYS") == 0 )
			{
				if( strcmp(chTime,"LAST") != 0 )
				{
					/* save physical conditions */
					fprintf( save.ipPnunit[ipPun], "%.5e\t%.4e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
					  radius.depth_mid_zone, phycon.te, dense.gas_phase[ipHYDROGEN], 
					  dense.eden, thermal.htot, wind.AccelTotalOutward, geometry.FillFac );
				}
			}

			else if( strcmp(save.chSave[ipPun],"PRES") == 0 )
			{
				/* the save pressure command */
				if( strcmp(chTime,"LAST") != 0 )
				{
					fprintf( save.ipPnunit[ipPun], 
					  "%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%c\n", 
					  /*A 1 #P depth */
					  radius.depth_mid_zone, 
					  /*B 2 Pcorrect */
					  pressure.PresTotlCorrect, 
					  /*C 3 Pcurrent */
					  pressure.PresTotlCurr, 
					  /*D 4 Pln + pintg 
					   * >>chng 06 apr 19, subtract pinzon the acceleration added in this zone
					   * since is not total at outer edge of zone, above is at inner edge */
					  pressure.PresTotlInit + pressure.PresInteg - pressure.pinzon, 
					  /*E 5 pgas (0) */
					  pressure.PresTotlInit, 
					  /*F 6 Pgas */
					  pressure.PresGasCurr, 
					  /*G 7 Pram */
					  pressure.PresRamCurr, 
					  /*H 8 P rad in lines */
					  pressure.pres_radiation_lines_curr, 
					  /*I 9 Pinteg subtract continuum rad pres which has already been added on */
					  pressure.PresInteg - pressure.pinzon, 
					  /*J 10 V(wind km/s) wind speed in km/s */
					  wind.windv/1e5,
					  /*K cad(km/s) sound speed in km/s */
					  timesc.sound_speed_adiabatic/1e5,
					  /* the magnetic pressure */
					  magnetic.pressure ,
					  /* the local turbulent velocity in km/s */
					  DoppVel.TurbVel/1e5 ,
					  /* turbulent pressure */
					  pressure.PresTurbCurr*DoppVel.lgTurb_pressure ,
					  /* gravitational pressure */
					  pressure.IntegRhoGravity,
					  // the integral of electron scattering acceleration in
					  // the absence of any absorptio, minus acceleration in current
					  // zone, which has been added in - done this way, result is
					  // zero in first zonen
					  pressure.PresIntegElecThin-pressure.pinzon_PresIntegElecThin,
					  // is this converged?
					  TorF(conv.lgConvPres) );
				}
			}
			else if( strcmp(save.chSave[ipPun],"PREL") == 0 )
			{
				/* line pressure contributors */
				fprintf( save.ipPnunit[ipPun], 
					"%.5e\t%.3e\t%.3e\t", 
					/*A 1 #P depth */
					radius.depth_mid_zone ,
					pressure.PresTotlCurr,
					pressure.pres_radiation_lines_curr/SDIV(pressure.PresTotlCurr) );
				PrtLinePres(save.ipPnunit[ipPun]);

			}

			else if( save.chSave[ipPun][0]=='R' )
			{
				/* work around internal limits to Microsoft vs compiler */
				if( strcmp(save.chSave[ipPun],"RADI") == 0 )
				{
					/* save radius information for all zones */
					if( strcmp(chTime,"LAST") != 0 )
					{
						fprintf( save.ipPnunit[ipPun], "%ld\t%.5e\t%.4e\t%.4e\n", 
						  nzone, radius.Radius_mid_zone, radius.depth_mid_zone, 
						  radius.drad );
					}
				}

				else if( strcmp(save.chSave[ipPun],"RADO") == 0 )
				{
					/* save radius information for only the last zone */
					if( strcmp(chTime,"LAST") == 0 )
					{
						fprintf( save.ipPnunit[ipPun], "%ld\t%.5e\t%.4e\t%.4e\n", 
							nzone, radius.Radius_mid_zone, radius.depth_mid_zone, 
							radius.drad );
					}
				}

				else if( strcmp(save.chSave[ipPun],"RESU") == 0 )
				{
					/*  save results of the calculation */
					if( strcmp(chTime,"LAST") == 0 )
						SaveResults(save.ipPnunit[ipPun]);
				}
				else
				{
					/* this can't happen */
					TotalInsanity();
				}
			}

			else if( strcmp(save.chSave[ipPun],"SECO") == 0 )
			{
				/*  save secondary ionization */
				if( strcmp(chTime,"LAST") != 0 )
					fprintf(save.ipPnunit[ipPun],
					"%.5e\t%.3e\t%.3e\t%.3e\n",
					radius.depth ,
					secondaries.csupra[ipHYDROGEN][0],
					secondaries.csupra[ipHYDROGEN][0]*2.02,
					secondaries.x12tot );
			}

			else if( strcmp(save.chSave[ipPun],"SOUS") == 0 )
			{
				/* full spectrum of continuum source function at 1 depth
				 *  command was "save source spectrum" */
				if( strcmp(chTime,"LAST") != 0 )
				{
					limit = MIN2(rfield.ipMaxBolt,rfield.nflux);
					for( j=0; j < limit; j++ )
					{
						fprintf( save.ipPnunit[ipPun], 
							"%.5e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", 
						  AnuUnit(rfield.AnuOrg[j]),
						  rfield.ConEmitLocal[nzone][j]/rfield.widflx[j], 
						  opac.opacity_abs[j], 
						  rfield.ConSourceFcnLocal[nzone][j], 
						  rfield.ConSourceFcnLocal[nzone][j]/plankf(j) ,
						  safe_div(rfield.ConSourceFcnLocal[nzone][j],rfield.flux[0][j]));
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"SOUD") == 0 )
			{
				/* parts of continuum source function vs depth
				 * command was save source function depth */
				j = iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s] + 2;
				fprintf( save.ipPnunit[ipPun], 
					"%.4e\t%.4e\t%.4e\t%.4e\n", 
				  opac.TauAbsFace[j-1], 
				  rfield.ConEmitLocal[nzone][j-1]/rfield.widflx[j-1]/MAX2(1e-35,opac.opacity_abs[j-1]), 
				  rfield.otscon[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1], 
				  rfield.otscon[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][0]-1]/opac.opacity_abs[iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1] );
			}

			/* this is save special option */
			else if( strcmp(save.chSave[ipPun],"SPEC") == 0 )
			{
				SaveSpecial(save.ipPnunit[ipPun],chTime);
			}

			/* save species information - number of levels */
			else if( strcmp(save.chSave[ipPun],"SPCn") == 0 )
			{
				static bool lgMustPrintHeader=true;
				if( lgMustPrintHeader )
				{
					lgMustPrintHeader = false;
					fprintf( save.ipPnunit[ipPun], "Levels\t%s",
						Species[0].chLabel );
					for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
					{
						fprintf( save.ipPnunit[ipPun], "\t%s",
								Species[ipSpecies].chLabel );
					}
					fprintf( save.ipPnunit[ipPun], "\n" );
					fprintf( save.ipPnunit[ipPun], "Max\t%li" ,
						Species[0].numLevels_max );
					for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
					{
						fprintf( save.ipPnunit[ipPun], "\t%li" ,
							Species[ipSpecies].numLevels_max );
					}
					fprintf( save.ipPnunit[ipPun], "\n");
				}
				fprintf( save.ipPnunit[ipPun], "Local\t%li" ,
					Species[0].numLevels_local );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					fprintf( save.ipPnunit[ipPun], "\t%li" ,
						Species[ipSpecies].numLevels_local );
				}
				fprintf( save.ipPnunit[ipPun], "\n");
			}
			else if( strcmp(save.chSave[ipPun],"SPCp") == 0 )
			{

				// print level populations for all species, for debugging
				// future - should do only for one species, identified as in "H2O"
				long mostlevels = -1;
				for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
				{
					if( Species[ipSpecies].lgActive )
						mostlevels = MAX2(mostlevels,Species[ipSpecies].numLevels_local);
				}

				static bool lgMustPrintHeader=true;
				if( lgMustPrintHeader )
				{
					lgMustPrintHeader = false;
					fprintf( save.ipPnunit[ipPun], "Species");
					for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
					{
						fprintf( save.ipPnunit[ipPun], "\t%s", Species[ipSpecies].chLabel );
					}
					fprintf( save.ipPnunit[ipPun], "\n" );
				}

				if( mostlevels < 0 )
				{
					fprintf(save.ipPnunit[ipPun],"No species were active\n");
				}
				else
				{
					for (long j=0;j<mostlevels;j++)
					{
						fprintf( save.ipPnunit[ipPun], "Level %li",j+1);

						for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
						{
							if ( Species[ipSpecies].numLevels_local > j && dBaseStates[ipSpecies][j].Pop != 0.)
							{
								fprintf( save.ipPnunit[ipPun], "\t%1.2E" ,dBaseStates[ipSpecies][j].Pop);
							}
							else
								fprintf( save.ipPnunit[ipPun], "\t X");
						}
						fprintf( save.ipPnunit[ipPun], "\n");
					}
				}
			}

			else if( strcmp(save.chSave[ipPun],"TEMP") == 0 )
			{
				static double deriv_old=-1;
				double deriv=-1. , deriv_sec;
				/* temperature and its derivatives */
				fprintf( save.ipPnunit[ipPun], "%.5e\t%.4e\t%.2e", 
					radius.depth_mid_zone, 
					phycon.te, 
					thermal.dCooldT );
				/* if second zone then have one deriv */
				if( nzone >1 )
				{
					deriv = (phycon.te - struc.testr[nzone-2])/ radius.drad;
					fprintf( save.ipPnunit[ipPun], "\t%.2e", deriv );
					/* if third zone then have second deriv */
					if( nzone > 2 )
					{
						deriv_sec = (deriv-deriv_old)/ radius.drad;
						fprintf( save.ipPnunit[ipPun], "\t%.2e", 
						  deriv_sec );
					}
					deriv_old = deriv;
				}
				fprintf( save.ipPnunit[ipPun], "\n");
			}

			/* time dependent model */
			else if( strcmp(save.chSave[ipPun],"TIMD") == 0 )
			{
				if( strcmp(chTime,"LAST") == 0 )
					DynaPunchTimeDep( save.ipPnunit[ipPun] , "END" );
			}

			/* execution time per zone */
			else if( strcmp(save.chSave[ipPun],"XTIM") == 0 )
			{
				static double ElapsedTime , ZoneTime;
				if( nzone<=1 )
				{
					ElapsedTime = cdExecTime();
					ZoneTime = 0.;
				}
				else
				{
					double t = cdExecTime();
					ZoneTime = t - ElapsedTime;
					ElapsedTime = t;
				}

				/* zone, time for this zone, elapsed time */
				fprintf( save.ipPnunit[ipPun], " %ld\t%.3f\t%.2f\n", 
				  nzone,  ZoneTime , ElapsedTime );
			}

			else if( strcmp(save.chSave[ipPun],"TPRE") == 0 )
			{
				/* temperature and its predictors, turned on with save tprid */
				fprintf( save.ipPnunit[ipPun], "%5ld %11.4e %11.4e %11.4e %g\n", 
				  nzone, phycon.TeInit, phycon.TeProp, phycon.te, 
				  (phycon.TeProp- phycon.te)/phycon.te );
			}

			else if( strcmp(save.chSave[ipPun],"WIND") == 0 )
			{
				/* wind velocity, radiative acceleration, and ratio total
				 * to electron scattering acceleration */
				/* first test only save last zone */
				if( (save.punarg[ipPun][0] == 0 && strcmp(chTime,"LAST") == 0)
					||
					/* this test save all zones */
					(save.punarg[ipPun][0] == 1  && strcmp(chTime,"LAST") != 0 ) )
				{
					fprintf( save.ipPnunit[ipPun], 
						"%.5e\t%.5e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n", 
						radius.Radius_mid_zone, 
						radius.depth_mid_zone, 
						wind.windv, 
						wind.AccelTotalOutward, 
						wind.AccelLine,
						wind.AccelCont ,
						wind.fmul ,
						wind.AccelGravity );
				}
			}

			else if( strcmp(save.chSave[ipPun],"XATT") == 0 )
			{
				/* attenuated incident continuum */
				ASSERT( grid.lgOutputTypeOn[2] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 2 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XRFI") == 0 )
			{
				/* reflected incident continuum */
				ASSERT( grid.lgOutputTypeOn[3] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 3 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XINC") == 0 )
			{
				/* incident continuum */
				ASSERT( grid.lgOutputTypeOn[1] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 1 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XDFR") == 0 )
			{
				/* reflected diffuse continuous emission */
				ASSERT( grid.lgOutputTypeOn[5] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 5 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XDFO") == 0 )
			{
				/* diffuse continuous emission outward */
				ASSERT( grid.lgOutputTypeOn[4] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 4 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XLNR") == 0 )
			{
				/* reflected lines */
				ASSERT( grid.lgOutputTypeOn[7] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 7 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XLNO") == 0 )
			{
				/* outward lines */
				ASSERT( grid.lgOutputTypeOn[6] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 6 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XREF") == 0 )
			{
				/* total reflected, lines and continuum */
				ASSERT( grid.lgOutputTypeOn[9] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 9 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XTOT") == 0 )
			{
				/* total spectrum, reflected plus transmitted */
				ASSERT( grid.lgOutputTypeOn[0] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 0 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XTRN") == 0 )
			{
				/* total outward, lines and continuum */
				ASSERT( grid.lgOutputTypeOn[8] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 8 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			else if( strcmp(save.chSave[ipPun],"XSPM") == 0 )
			{
				/* exp(-tau) to the illuminated face */
				ASSERT( grid.lgOutputTypeOn[10] );

				if( strcmp(chTime,"LAST") == 0 )
				{
					if( grid.lgGrid )
						saveFITSfile( save.ipPnunit[ipPun], 10 );
					else
					{
						fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
			/* there are a few "save" commands that are handled elsewhere
			 * save dr is an example.  These will have lgRealSave set false */
			else if( save.lgRealSave[ipPun] )
			{
				/* this is insanity, internal flag set in ParseSave not seen here */
				fprintf( ioQQQ, " PROBLEM DISASTER SaveDo does not recognize flag %4.4s set by ParseSave.  This is impossible.\n", 
				  save.chSave[ipPun] );
				TotalInsanity();
			}

			/* print special hash string to separate out various iterations 
			 * chTime is LAST on last iteration
			 * save.lgHashEndIter flag is true by default, set false
			 * with "no hash" keyword on save command
			 * save.lg_separate_iterations is true by default, set false
			 * when save time dependent calcs since do not want special
			 * character between time steps 
			 * grid.lgGrid is only true when doing a grid of calculations */
			if( strcmp(chTime,"LAST") == 0 && 
				!(iterations.lgLastIt && !grid.lgGrid ) && 
				save.lgHashEndIter[ipPun] &&
				save.lg_separate_iterations[ipPun] &&
				!save.lgFITS[ipPun] )
			{
				if( dynamics.lgTimeDependentStatic && strcmp( save.chHashString , "TIME_DEP" )==0 )
				{
					fprintf( save.ipPnunit[ipPun], "\"time=%f\n",
						dynamics.time_elapsed );
				}
				else
				{
					fprintf( save.ipPnunit[ipPun], "%s\n",
						save.chHashString );
				}
			}
			if( save.lgFLUSH )
				fflush( save.ipPnunit[ipPun] );
		}
	}
	return;
}

/*SaveLineIntensity produce the 'save lines intensity' output */
STATIC void SaveLineIntensity(FILE * ioPUN, long int ipPun )
{
	long int i;

	DEBUG_ENTRY( "SaveLineIntensity()" );

	/* used to save out all the emission line intensities
	 * first initialize the line image reader */

	fprintf( ioPUN, "**********************************************************************************************************************************\n" );
	input.echo(ioPUN);

	/* now print any cautions or warnings */
	cdWarnings( ioPUN);
	cdCautions( ioPUN);
	fprintf( ioPUN, "zone=%5ld\n", nzone );
	fprintf( ioPUN, "**********************************************************************************************************************************\n" );
	fprintf( ioPUN, "begin emission lines\n" );

	/* only save non-zero intensities */
	SaveResults1Line(ioPUN,"    ",0,0.,"Start");

	// check whether intrinsic or emergent line emissivity
	bool lgEmergent = false;
	if( save.punarg[ipPun][0] > 0 )
		lgEmergent = true;

	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSv[i].SumLine[lgEmergent] > 0. )
		{
			SaveResults1Line(ioPUN,(char*)LineSv[i].chALab,LineSv[i].wavelength,
			  LineSv[i].SumLine[save.lgEmergent[ipPun]], "Line ");
		}
	}

	SaveResults1Line(ioPUN,"    ",0,0.,"Flush");

	fprintf( ioPUN, "     \n" );
	fprintf( ioPUN, "**********************************************************************************************************************************\n" );

	return;
}

/* lgSaveOpticalDepths true says save optical depths */
static bool lgPopsFirstCall , lgSaveOpticalDepths;

/*SaveLineStuff save optical depths or source functions for all transferred lines */
STATIC void SaveLineStuff(
  FILE * ioPUN,
  const char *chJob , 
  realnum xLimit )
{

	long int nelem , ipLo , ipHi , i , ipISO;
	long index = 0;
	static bool lgFirst=true;

	DEBUG_ENTRY( "SaveLineStuff()" );

	/*find out which job this is and set a flag to use later */
	if( strcmp( &*chJob , "optical" ) == 0 )
	{
		/* save line optical depths */
		lgSaveOpticalDepths = true;
		lgPopsFirstCall = false;
	}
	else if( strcmp( &*chJob , "populat" ) == 0 )
	{
		lgSaveOpticalDepths = false;
		/* level population information */
		if( lgFirst )
		{
			lgPopsFirstCall = true;
			fprintf(ioPUN,"index\tAn.ion\tgLo\tgUp\tE(wn)\tgf\n");
			lgFirst = false;
		}
		else
		{
			lgPopsFirstCall = false;
		}
	}
	else
	{
		fprintf( ioQQQ, " insane job in SaveLineStuff =%s\n", 
		  &*chJob );
		cdEXIT(EXIT_FAILURE);
	}

	/* loop over all lines, calling put1Line to create info (routine located below) */
	/* hydrogen like lines */
	/* >>chng 02 may 16, had been explicit H and He-like loops */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem]  )
			{
				/* 06 aug 28, from numLevels_max to _local. */
				for( ipHi=1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
				{
					for( ipLo=0; ipLo <ipHi; ipLo++ )
					{
						if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
							continue;

						++index;
						Save1Line( &Transitions[ipISO][nelem][ipHi][ipLo], ioPUN, xLimit, index, GetDopplerWidth(dense.AtomicWeight[nelem]) );
					}
				}
				/* also do extra Lyman lines if optical depths are to be done,
				 * these are line that are included only for absorption, not in the
				 * model atoms */
				if( lgSaveOpticalDepths )
				{
					/* >>chng 02 aug 23, for he-like, had starting on much too high a level since
					 * index was number of levels - caught by Adrian Turner */
					/* now output extra line lines, starting one above those already done above */
					/*for( ipHi=iso.numLevels_max[ipISO][nelem]; ipHi < iso.nLyman[ipISO]; ipHi++ )*/
					/* 06 aug 28, from numLevels_max to _local. */
					for( ipHi=StatesElemNEW[nelem][nelem-ipISO][iso.numLevels_local[ipISO][nelem]-1].n+1; ipHi < iso.nLyman[ipISO]; ipHi++ )
					{
						++index;
						Save1Line( &ExtraLymanLines[ipISO][nelem][ipHi], ioPUN, xLimit, index, GetDopplerWidth(dense.AtomicWeight[nelem]) );
					}
				}
			}
		}
	}

	/* index from 1 to skip over dummy line */
	for( i=1; i < nLevel1; i++ )
	{
		++index;
		Save1Line( &TauLines[i], ioPUN, xLimit, index, GetDopplerWidth(dense.AtomicWeight[TauLines[i].Hi->nelem-1]) );
	}

	for( i=0; i < nWindLine; i++ )
	{
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			++index;
			Save1Line( &TauLine2[i], ioPUN, xLimit, index, GetDopplerWidth(dense.AtomicWeight[TauLine2[i].Hi->nelem-1]) );
		}
	}

	for( i=0; i < nUTA; i++ )
	{
		++index;
		Save1Line( &UTALines[i], ioPUN, xLimit, index, GetDopplerWidth(dense.AtomicWeight[UTALines[i].Hi->nelem-1]) );
	}


	/* do optical depths of FeII lines, if large atom is turned on */
	FeIIPunchLineStuff( ioPUN , xLimit  , index);

	/* save optical depths of H2 lines */
	H2_PunchLineStuff( ioPUN , xLimit  , index);

	/*fprintf(ioPUN, "##################################\n"); */
	fprintf( ioPUN , "%s\n",save.chHashString );
	return;
}

/*Save1Line called by SaveLineStuff to produce output for one line */
void Save1Line( transition * t , FILE * ioPUN , realnum xLimit  , long index, realnum DopplerWidth )
{

	if( lgSaveOpticalDepths )
	{
		/* optical depths, no special first time, only print them */
		if( t->Emis->TauIn >= xLimit )
		{
			/* label like "C  4" or "H  1" */
			fprintf( ioPUN , "%2.2s%2.2s\t", 
			  elementnames.chElementSym[t->Hi->nelem-1] ,
			  elementnames.chIonStage[t->Hi->IonStg-1]  );

			/* print wavelengths, either line in main printout labels, 
			 * or in various units in exponential notation - prt_wl is in prt.c */
			if( strcmp( save.chConPunEnr[save.ipConPun], "labl" )== 0 )
			{
				prt_wl( ioPUN , t->WLAng );
			}
			else
			{
				/* this converts energy in Rydbergs into any of the other units */
				fprintf( ioPUN , "%.7e", AnuUnit((realnum)(t->EnergyWN * WAVNRYD)) );
			}
			/* print the optical depth */
			fprintf( ioPUN , "\t%.3f", t->Emis->TauIn  );
			/* damping constant */
			fprintf(ioPUN, "\t%.3e", 
				t->Emis->dampXvel / DopplerWidth );
			fprintf(ioPUN, "\n");
		}
	}
	else if( lgPopsFirstCall )
	{
		char chLbl[11];
		/* first call to line populations, print atomic parameters and indices */
		strcpy( chLbl, chLineLbl(t) );
		fprintf(ioPUN, "%li\t%s" , index , chLbl );
		/* stat weights */
		fprintf(ioPUN, "\t%.0f\t%.0f", 
			t->Lo->g ,t->Hi->g);
		/* energy difference, gf */
		fprintf(ioPUN, "\t%.2f\t%.3e", 
			t->EnergyWN ,t->Emis->gf);
		fprintf(ioPUN, "\n");
	}
	else
	{
		/* not first call, so do level populations and indices defined above */
		if( t->Hi->Pop > xLimit )
		{
			/* >>chng 05 may 08, add abundances, which for iso-seq species is
			 * the density of the parent ion, for other lines, is unity.
			 * had not been included so pops for iso seq were rel to parent ion.
			 * caught by John Everett */
			/* multiplication by abundance no longer necessary since iso pops now denormalized */
			fprintf(ioPUN,"%li\t%.2e\t%.2e\n", index, t->Lo->Pop, t->Hi->Pop );
		}
	}
}

/*SaveNewContinuum produce the 'save new continuum' output */
STATIC void SaveNewContinuum(FILE * ioPUN )
{
	long int ipLo, ipHi,
		j ,
		ncells;

	double wllo, wlhi;

	double *energy, 
		*cont_incid,
		*cont_atten,
		*diffuse_in,
		*diffuse_out,
		*emis_lines_out,
		*emis_lines_in;

	/* get the low limit */
	wllo = rfield.anu[0];
	/* get high-energy limit */
	wlhi = rfield.anu[rfield.nflux-1];
	/* use native continuum mesh */
	ipLo = ipoint(wllo)-1;
	ipHi = ipoint(wlhi)-1;
	ncells = ipHi - ipLo + 1;

	/* now allocate the space */
	energy = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	cont_incid = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	cont_atten = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	diffuse_in = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	diffuse_out = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	emis_lines_out = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	emis_lines_in = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double) );
	/*emis_lines_pump_out = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double));
	  emis_lines_pump_in = (double *)MALLOC( (size_t)(ncells+1)*sizeof(double));*/

	/* now convert to rydbergs */
	for( j=0; j<ncells; ++j )
	{
		energy[j] = rfield.AnuOrg[j+ipLo] - rfield.widflx[j+ipLo]/2.;
	}

	fixit();	
	/* all of these should be rerouted to cdSPEC2, but units not right
	 * need ergs multiplied for one, and continuum and lines may not be added correctly.
	 * Goal is to abandon current cdSPEC and replace it with current cdSPEC2 */

#if	1
	/* for cdSPEC the energy vector is the lower edge of the energy cell */
	/* get incident continuum */
	cdSPEC( 1 , energy , ncells , cont_incid );
	/* get attenuated incident continuum */
	cdSPEC( 2 , energy , ncells , cont_atten );
	/* get diffuse continuous emission, reflected */
	cdSPEC( 5 , energy , ncells , diffuse_in );
	/* get continuous emission outward direction */
	cdSPEC( 4 , energy , ncells , diffuse_out );
	/** \todo	2	- NB - if continuum resolution changed the lines WILL NOT WORK */
	/* get all outward lines */
	cdSPEC( 6 , energy , ncells , emis_lines_out );
	/* get all reflected lines */
	cdSPEC( 7 , energy , ncells , emis_lines_in );
#else
	cdSPEC2( 1, rfield.nflux, 0, rfield.nflux - 1, cont_incid );
	/* get attenuated incident continuum */
	cdSPEC2( 2, rfield.nflux, 0, rfield.nflux - 1, cont_atten );
	/* get diffuse continuous emission, reflected */
	cdSPEC2( 5, rfield.nflux, 0, rfield.nflux - 1, diffuse_in );
	/* get continuous emission outward direction */
	cdSPEC2( 4, rfield.nflux, 0, rfield.nflux - 1, diffuse_out );
	/* get all outward lines */
	cdSPEC2( 6, rfield.nflux, 0, rfield.nflux - 1, emis_lines_out );
	/* get all reflected lines */
	cdSPEC2( 7, rfield.nflux, 0, rfield.nflux - 1, emis_lines_in );
#endif

	/* for this example we will do a wavelength range */
	for( j=0; j<ncells-1; ++j )
	{
		/* photon energy in appropriate energy or wavelength units */
		fprintf( ioPUN,"%.5e\t", AnuUnit((realnum)(energy[j]+rfield.widflx[j+ipLo]/2.) ) );
		fprintf( ioPUN,"%.3e\t", cont_incid[j] );
		fprintf( ioPUN,"%.3e\t", cont_atten[j] );
		fprintf( ioPUN,"%.3e\t", diffuse_in[j]+diffuse_out[j] );
		fprintf( ioPUN,"%.3e", 
			emis_lines_out[j]+emis_lines_in[j]/*+emis_lines_pump_out[j]+emis_lines_pump_in[j]*/ );
		fprintf( ioPUN,"\n" );
	}

	free(energy);
	free(cont_incid);
	free(diffuse_in);
	free(diffuse_out);
	free(cont_atten);
	free(emis_lines_out);
	free(emis_lines_in);
	/*free(emis_lines_pump_out);
	free(emis_lines_pump_in );*/
}

/* save AGN3 hemiss, for Chapter 4, routine is below */
STATIC void AGN_Hemis(FILE *ioPUN )
{
	const int NTE = 4;
	realnum te[NTE] = {5000., 10000., 15000., 20000.};
	realnum *agn_continuum[NTE];
	double TempSave = phycon.te;
	long i , j;

	DEBUG_ENTRY( "AGN_Hemis()" );

	/* make table of continuous emission at various temperatuers */
	/* first allocate space */
	for( i=0;i<NTE; ++i)
	{
		agn_continuum[i] = (realnum *)MALLOC((unsigned)rfield.nflux*sizeof(realnum) );

		/* set the next temperature */
		/* recompute everything at this new temp */
		TempChange(te[i] , true);
		/* converge the pressure-temperature-ionization solution for this zone */
		ConvPresTempEdenIoniz();

		/* now get the thermal emission */
		RT_diffuse();
		for(j=0;j<rfield.nflux; ++j )
		{
			agn_continuum[i][j] = rfield.ConEmitLocal[nzone][j]/(realnum)dense.eden/
				(dense.xIonDense[ipHYDROGEN][1] + dense.xIonDense[ipHELIUM][1] + dense.xIonDense[ipHELIUM][2] );
		}
	}

	/* print title for line */
	fprintf(ioPUN,"wl");
	for( i=0;i<NTE; ++i)
	{
		fprintf(ioPUN,"\tT=%.0f",te[i]);
	}
	fprintf( ioPUN , "\tcont\n"); 

	/* not print all n temperatures across a line */
	for(j=0;j<rfield.nflux; ++j )
	{
		fprintf( ioPUN , "%12.5e", 
		  AnuUnit(rfield.AnuOrg[j]) );
		/* loop over the temperatures, and for each, calculate a continuum */
		for( i=0;i<NTE; ++i)
		{
			fprintf(ioPUN,"\t%.3e",agn_continuum[i][j]*rfield.anu2[j]*EN1RYD/rfield.widflx[j]);
		}
		/* cont label and end of line*/
		fprintf( ioPUN , "\t%s\n" , rfield.chContLabel[j]); 
	}

	/* now free the continua */
	for( i=0;i<NTE; ++i)
	{
		free( agn_continuum[i] );
	}

	/* Restore temperature stored before this routine was called	*/
	/* and force update */
	TempChange(TempSave , true);

	fprintf( ioQQQ, "AGN_Hemis - result of save AGN3 hemis - I have left the code in a disturbed state, and will now exit.\n");
	cdEXIT(EXIT_FAILURE);
}

/*SaveResults save results from save results command */
/*SaveResults1Line do single line of output for the save results and save line intensity commands */
STATIC void SaveResults(FILE* ioPUN)
{
	long int i , nelem , ion;

	DEBUG_ENTRY( "SaveResults()" );

	/* used to save out line intensities, optical depths,
	 * and column densities */

	fprintf( ioPUN, "**********************************************************************************************************************************\n" );
	input.echo(ioPUN);

	/* first print any cautions or warnings */
	cdWarnings(ioPUN);
	cdCautions(ioPUN);
	fprintf( ioPUN, "**********************************************************************************************************************************\n" );

	fprintf( ioPUN, "C*OPTICAL DEPTHS ELECTRON=%10.3e\n", opac.telec );

	fprintf( ioPUN, "BEGIN EMISSION LINES\n" );
	SaveResults1Line(ioPUN,"    ",0,0.,"Start");

	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSv[i].SumLine[0] > 0. )
		{
			SaveResults1Line(ioPUN,(char*)LineSv[i].chALab,LineSv[i].wavelength,
			  LineSv[i].SumLine[0],"Line ");
		}
	}

	SaveResults1Line(ioPUN,"    ",0,0.,"Flush");

	fprintf( ioPUN, "     \n" );

	fprintf( ioPUN, "BEGIN COLUMN DENSITIES\n" );

	/* this dumps out the whole array,*/
	/* following loop relies on LIMELM being 30, assert it here in case
	 * this is ever changed */
	ASSERT( LIMELM == 30 );
	/* this order of indices is to keep 30 as the fastest variable,
	 * and the 32 (LIMELM+1) as the slower one */
	for( nelem=0; nelem<LIMELM; nelem++ )
	{
		for(ion=0; ion < nelem+1; ion++)
		{
			fprintf( ioPUN, " %10.3e", mean.xIonMean[0][nelem][ion][0] );
			/* throw line feed every 10 numbers */
			if( nelem==9|| nelem==19 || nelem==29 )
			{
				fprintf( ioPUN, "\n" );
			}
		}
	}

	fprintf( ioPUN, "END OF RESULTS\n" );
	fprintf( ioPUN, "**********************************************************************************************************************************\n" );
	return;
}

static const int LINEWIDTH = 6;

/*SaveResults1Line do single line of output for the save results and save line intensity commands */
/* the number of emission lines across one line of printout */
STATIC void SaveResults1Line(
  FILE * ioPUN, 
  /* 4 char + null string */
  const char *chLab, 
  realnum wl, 
  double xInten, 
  /* null terminated string saying what to do */
  const char *chFunction)
{

	long int i;
	static realnum wavelength[LINEWIDTH];
	static long ipLine;
	static double xIntenSave[LINEWIDTH];
	static char chLabSave[LINEWIDTH][5];

	DEBUG_ENTRY( "SaveResults1Line()" );

	/* if LineWidth is changed then change format in write too */

	if( strcmp(chFunction,"Start") == 0 )
	{
		ipLine = 0;
	}
	else if( strcmp(chFunction,"Line ") == 0 )
	{
		/* save results in array so that they can be printed when done */
		wavelength[ipLine] = wl;
		strcpy( chLabSave[ipLine], chLab );
		xIntenSave[ipLine] = xInten;

		/* now increment the counter and then check if we have filled the line, 
		 * and so should print it */
		++ipLine;
		/* do print now if we are in column mode (one line per line) or if we have filled up
		 * the line */
		if( ( strcmp(save.chPunRltType,"column") == 0 ) || ipLine == LINEWIDTH )
		{
			/* "array " is usual array 6 wide, "column" is one line per line */
			for( i=0; i < ipLine; i++ )
			{
				fprintf( ioPUN, " %4.4s ", chLabSave[i] );
				prt_wl( ioPUN, wavelength[i] );
				fprintf( ioPUN,"\t%.3e", xIntenSave[i]);
				/* >>chng 02 apr 24, do not print type */
				/* single column for input into data base */
				if( strcmp(save.chPunRltType,"column") == 0 )				
					fprintf( ioPUN, "\n" );
			}
			/* only put cr if we did not just put one */
			if( strcmp(save.chPunRltType,"array ") == 0 )				
				fprintf( ioPUN, " \n" );
			ipLine = 0;
		}
	}
	else if( strcmp(chFunction,"Flush") == 0 )
	{
		if( ipLine > 0 )
		{
			/* this is an option to print many emission lines across an output line,
			 * the array option, or a single column of numbers, the column option
			 * that appears on the "save results" and "save intensity" commands
			 */
			/* usual array 6 wide */
			for( i=0; i < ipLine; i++ )
			{
				fprintf( ioPUN, " %4.4s", chLabSave[i] );
				prt_wl( ioPUN, wavelength[i] );
				fprintf( ioPUN,"\t%.3e", xIntenSave[i]);
				/* >>chng 02 apr 24, do not print type */
				/* single column for input into data base */
				if( strcmp(save.chPunRltType,"column") == 0 )				
					fprintf( ioPUN, "\n" );
			}
			if( strcmp(save.chPunRltType,"array ") == 0 )
				fprintf( ioPUN, " \n" );
		}
	}
	else
	{
		fprintf( ioQQQ, " SaveResults1Line called with insane request=%5.5s\n",
		  chFunction );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

static const int NENR_GAUNT = 37;
static const int NTE_GAUNT = 21;

/*SaveGaunts called by save gaunts command to output gaunt factors */
STATIC void SaveGaunts(FILE* ioPUN)
{
	long int i, 
	  charge,
	  ite, 
	  j;

	realnum ener[NENR_GAUNT], 
	  ste[NTE_GAUNT],
	  z;
	double tempsave;
    double g[NENR_GAUNT][NTE_GAUNT], gfac;


	DEBUG_ENTRY( "SaveGaunts()" );

	/* this routine is called from the PUNCH GAUNTS command
	 * it drives the gaunt factor routine to save gaunts over full range */
	tempsave = phycon.te;

	for( i=0; i < NTE_GAUNT; i++ )
	{
		ste[i] = 0.5f*i;
	}

	for( i=0; i < NENR_GAUNT; i++ )
	{
		ener[i] = 0.5f*i - 8.f;
	}

	for( charge = 1; charge<LIMELM; charge++ )
	{
		/* nuclear charge */
		z = (realnum)log10((double)charge);

		/* energy is log of energy */
		for( ite=0; ite < NTE_GAUNT; ite++ )
		{
			phycon.alogte = ste[ite];
			phycon.te = pow(10.,phycon.alogte);

			for( j=0; j < NENR_GAUNT; j++ )
			{
				gfac = cont_gaunt_calc( phycon.te, (double)charge, pow(10.,(double)ener[j]) );
				/*fprintf(ioQQQ, "z %.2e ener %.2e temp %.2e gfac %.3e \n",
					pow(10.,z), pow(10.,ener[j]), (double)phycon.te, gfac );*/
				g[j][ite] = gfac;
			}
		}

		/* now save out the results */
		fprintf( ioPUN, "      " );
		for( i=1; i <= NTE_GAUNT; i++ )
		{
			fprintf( ioPUN, "\t%6.1f", ste[i-1] );
		}
		fprintf( ioPUN, "\n" );

		for( j=0; j < NENR_GAUNT; j++ )
		{
			fprintf( ioPUN, "\t%6.2f", ener[j] );
			for( ite=0; ite < NTE_GAUNT; ite++ )
			{
				fprintf( ioPUN, "\t%6.2f", log10(g[j][ite]) );
			}
			fprintf( ioPUN, "\n" );
		}

		fprintf( ioPUN, "      " );
		for( i=0; i < NTE_GAUNT; i++ )
		{
			fprintf( ioPUN, "\t%6.1f", ste[i] );
		}
		fprintf( ioPUN, "\n" );

		/* now do the same thing as triplets */
		fprintf( ioPUN, "      " );
		for( i=0; i < NTE_GAUNT; i++ )
		{
			fprintf( ioPUN, "\t%6.1f", ste[i] );
		}
		fprintf( ioPUN, "\n" );

		for( i=0; i < NTE_GAUNT; i++ )
		{
			for( j=0; j < NENR_GAUNT; j++ )
			{
				fprintf( ioPUN, "\t%6.2f\t%6.2f\t%6.2e\n", ste[i], ener[j], 
				  log10(g[j][i]) );
			}
		}

		fprintf( ioPUN, "Below is log(u), log(gamma**2), gff\n" );
		/* do the same thing as above, but this time print log(u) and log(gamma2) instead of temp and energy.	*/
		for( i=0; i < NTE_GAUNT; i++ )
		{
			for( j=0; j < NENR_GAUNT; j++ )
			{
				fprintf( ioPUN, "\t%6.2f\t%6.2f\t%6.2e\n", 2.*z + log10(TE1RYD) - ste[i] , log10(TE1RYD)+ener[j]-ste[i], 
				log10(g[j][i]) );
			}
		}
		fprintf( ioPUN, "end of charge = %li\n", charge);
		fprintf( ioPUN, "****************************\n");
	}

	phycon.te = tempsave;
	phycon.alogte = log10( phycon.te );
	return;
}
