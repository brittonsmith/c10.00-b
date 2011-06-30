/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_hydro put H-like iso sequence into line intensity stack */
#include "cddefines.h"
#include "atmdat.h"
#include "dense.h"
#include "prt.h"
#include "hydrogenic.h"
#include "iso.h"
#include "rfield.h"
#include "geometry.h"
#include "lines.h"
#include "lines_service.h"
#include "phycon.h"
#include "radius.h"
#include "secondaries.h"
#include "taulines.h"
#include "trace.h"

void lines_hydro(void)
{
	long ipISO = ipH_LIKE;
	long int i, nelem, ipHi, ipLo;
	char chLabel[5]="    ";

	double hbetab, 
		em , 
		pump ,
		caseb;

	DEBUG_ENTRY( "lines_hydro()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   lines_hydro called\n" );

	// this can be changed with the atom levels command but must be at least 3
	ASSERT( iso.n_HighestResolved_max[ipH_LIKE][ipHYDROGEN] >= 3 );
	ASSERT( iso.n_HighestResolved_max[ipH_LIKE][ipHELIUM] >= 3 );

	i = StuffComment( "H-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start H -like iso sequence ");

	linadd(MAX2(0.,iso.xLineTotCool[ipH_LIKE][ipHYDROGEN]),912,"Clin",'c',
		"  total collisional cooling due to all hydrogen lines ");

	linadd(MAX2(0.,-iso.xLineTotCool[ipH_LIKE][ipHYDROGEN]),912,"Hlin",'h'	,
		"  total collisional heating due to all hydrogen lines ");

	linadd(MAX2(0.,iso.xLineTotCool[ipH_LIKE][ipHELIUM]),228,"Clin",'c',
		"  total collisional cooling due to all HeII lines ");

	linadd(MAX2(0.,-iso.xLineTotCool[ipH_LIKE][ipHELIUM]),228,"Hlin",'h'	,
		"  total collisional heating due to all HeII lines ");

	/*fprintf(ioQQQ," debugg\t%.2e\t%.2e\t%.2e\n", 
		radius.drad,
		iso.xLineTotCool[ipH_LIKE][ipHYDROGEN] , 
		iso.cLya_cool[ipH_LIKE][ipHYDROGEN]);*/

	/* >>chng 95 jun 25 changed from info to cooling to pick this up in primal.in   */
	linadd(MAX2(0.,iso.cLya_cool[ipH_LIKE][ipHYDROGEN]),1216,"Cool",'i',
		"collisionally excited La cooling ");

	linadd(MAX2(0.,-iso.cLya_cool[ipH_LIKE][ipHYDROGEN]),1216,"Heat",'i',
		"  collisionally de-excited La heating ");

	linadd(MAX2(0.,iso.cLyrest_cool[ipH_LIKE][ipHYDROGEN]),960,"Crst",'i',
		"  cooling due to n>2 Lyman lines ");

	linadd(MAX2(0.,-iso.cLyrest_cool[ipH_LIKE][ipHYDROGEN]),960,"Hrst",'i',
		"  heating due to n>2 Lyman lines ");

	linadd(MAX2(0.,iso.cBal_cool[ipH_LIKE][ipHYDROGEN]),4861,"Crst",'i',
		"  cooling due to n>3 Balmer lines ");

	linadd(MAX2(0.,-iso.cBal_cool[ipH_LIKE][ipHYDROGEN]),4861,"Hrst",'i',
		"  heating due to n>3 Balmer lines ");

	linadd(MAX2(0.,iso.cRest_cool[ipH_LIKE][ipHYDROGEN]),0,"Crst",'i',
		"  cooling due to higher Paschen lines ");

	linadd(MAX2(0.,-iso.cRest_cool[ipH_LIKE][ipHYDROGEN]),0,"Hrst",'i',
		"  heating due to higher Paschen lines ");

	/* remember largest fractional ionization of H due to secondaries */
	secondaries.SecHIonMax = MAX2( secondaries.SecHIonMax , secondaries.sec2total );

	/* remember fraction of H ionizations due to ct */
	atmdat.HIonFracMax = MAX2( atmdat.HIonFracMax, atmdat.HIonFrac);

	/* remember largest fraction of thermal collisional ionization of H ground state */
	hydro.HCollIonMax = 
		(realnum)MAX2( hydro.HCollIonMax , hydro.H_ion_frac_collis );

	linadd(secondaries.x12tot*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop*1.634e-11,1216,"LA X" ,'i',
		"Lyaa contribution from suprathermal secondaries from ground ");

	/* factor of 0.4836 is ratio of A(4-2)/(A(4-3)+A(4-2))
	 * the IPLNPUMP is the actual pumping rate per atom */
	/* H-beta produced by continuum pumping in optically thin ld limit */
	pump = (double)(Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH1s].Emis->pump*StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop*4.09e-12*0.4836);
	linadd(pump,4861,"Pump",'r',
		   "part of Hbeta formed by continuum pumping");

	linadd(MAX2(0.,iso.coll_ion[ipH_LIKE][ipHYDROGEN]),0,"CION",'c',
		"collision ionization cooling of hydrogen ");

	linadd(MAX2(-iso.coll_ion[ipH_LIKE][ipHYDROGEN],0.),0,"3bHt",'h',
		"  this is the heating due to 3-body recombination ");

	linadd(MAX2(0.,iso.coll_ion[ipH_LIKE][ipHELIUM]),0,"He2C",'c',
		"collision ionization cooling of He+ ");

	linadd(MAX2(-iso.coll_ion[ipH_LIKE][ipHELIUM],0.),0,"He2H",'h',
		"  this is the heating due to 3-body recombination onto He+");

	fixit();  //why is there a zero here?
	linadd(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop*0.*iso.pestrk[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s]*1.634e-11,1216,"Strk",'i',
	  "  Stark broadening contribution to line ");

	linadd(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH3s].Pop*iso.pestrk[ipH_LIKE][ipHYDROGEN][ipH3s][ipH2p]*3.025e-12,
	  6563,"Strk",'i',
	  "  Stark broadening contribution to line ");

	linadd(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH4s].Pop*iso.pestrk[ipH_LIKE][ipHYDROGEN][ipH4s][ipH2p]*4.084e-12,
	  4861,"Strk",'i',
	  "Stark broadening contribution to line ");

	linadd(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH4p].Pop*iso.pestrk[ipH_LIKE][ipHYDROGEN][ipH4p][ipH3s]*1.059e-12,
	  18751,"Strk",'i',
		   " Stark broadening contribution to line ");

	/* pestrk[5,4] is A[4,5]*pest[4,5] 
	 * Stark broadening contribution to line */
	if( iso.n_HighestResolved_max[ipH_LIKE][ipHYDROGEN] >= 5 )
	{
		long ip5p = iso.QuantumNumbers2Index[ipH_LIKE][ipHYDROGEN][5][1][2];
		linadd(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ip5p].Pop*iso.pestrk[ipH_LIKE][ipHYDROGEN][ip5p][ipH4s]*4.900e-13,40512,"Strk",'i',
			"Stark broadening part of line");
	}
	/* this can fail if RT_line_all never updates the ots rates, a logic error,
	 * but only assert this during actual calculation (ipass>0), */
	ASSERT( LineSave.ipass  <1 ||
		Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->ots>= 0.);

	linadd(Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->ots*Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].EnergyErg, 1216,"Dest",'i',
		"  portion of line lost due to absorp by background opacity ");

	/* portion of line lost due to absorb by background opacity */
	linadd(Transitions[ipH_LIKE][ipHYDROGEN][ipH3p][ipH2s].Emis->ots*Transitions[ipH_LIKE][ipHYDROGEN][ipH3p][ipH2s].EnergyErg, 6563,"Dest",'i',
		"Ha destroyed by background opacity");

	/* portion of line lost due to absorp by background opacity */
	if( iso.n_HighestResolved_max[ipH_LIKE][ipHYDROGEN] >= 5 )
	{
		long ip5p = iso.QuantumNumbers2Index[ipH_LIKE][ipHYDROGEN][5][1][2];
		linadd(Transitions[ipH_LIKE][ipHYDROGEN][ip5p][ipH4s].Emis->ots*Transitions[ipH_LIKE][ipHYDROGEN][ip5p][ipH4s].EnergyErg,40516, "Dest",'i',
			"portion of line lost due to absorb by background opacity");
	}

	/* portion of line lost due to absorb by background opacity */
	if( iso.numLevels_max[ipH_LIKE][ipHYDROGEN] > ipH4p )
		linadd(Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH2s].Emis->ots*Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH2s].EnergyErg, 4861,"Dest",'i',
			"portion of line lost due to absorb by background opacity");

	/* portion of line lost due to absorb by background opacity */
	if( iso.numLevels_max[ipH_LIKE][ipHYDROGEN] > ipH4p )
		linadd(Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH3s].Emis->ots*Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH3s].EnergyErg ,18751, "Dest",'i',
			"portion of line lost due to absorb by background opacity");

	linadd(StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop*Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Aul*
		hydro.dstfe2lya*Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].EnergyErg , 1216 , "Fe 2" , 'i',
		"Ly-alpha destroyed by overlap with FeII " );

	linadd(iso.RadRec_caseB[ipH_LIKE][ipHYDROGEN]*dense.xIonDense[ipHYDROGEN][1]*dense.eden * 1.64e-11,1216,"Ca B",'i',
		" simple high-density case b intensity of Ly-alpha, no two photon ");

	/* these entries only work correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		/* H-beta computed from Q(H) and specified covering factor */
		if( nzone == 1 )
		{
			/* evaluate the case b emissivity by interpolating on the hummer & storey tables */
			caseb = rfield.qhtot*
				atmdat_HS_caseB( 4 , 2 , 1 , phycon.te , dense.eden, 'b' ) / iso.RadRec_caseB[ipH_LIKE][ipHYDROGEN];
			/* the atmdat_HS_caseB returned -1 if the physical conditions were outside range of validity.  
			 * In this case use simple approximation with no temperature or density dependence */
			if( caseb < 0 )
			{
				caseb = rfield.qhtot*4.75e-13;
			}
			LineSv[LineSave.nsum].SumLine[0] = 0.;
			LineSv[LineSave.nsum].SumLine[1] = 0.;
		}
		else
		{
			caseb = 0.;
		}
		/* H-beta computed from Q(H) and specified covering factor */
		linadd( caseb/radius.dVeffAper*geometry.covgeo , 4861 , "Q(H)" , 'i' ,
			"Case B H-beta computed from Q(H) and specified covering factor");

		if( nzone == 1 )
		{
			caseb = rfield.qhtot*Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].EnergyErg;
			LineSv[LineSave.nsum].SumLine[0] = 0.;
			LineSv[LineSave.nsum].SumLine[1] = 0.;
		}
		else
		{
			caseb = 0.;
		}
		/* >>chng 02 nov 05, better approximation for Lya for temperature of first zone */
		linadd( caseb/radius.dVeffAper*geometry.covgeo , 1216 , "Q(H)" , 'i',
			"Ly-alpha from Q(H), high-dens lim, specified covering factor" );
	}

	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			ASSERT( iso.n_HighestResolved_max[ipH_LIKE][nelem] >= 3 );

			for( ipHi=1; ipHi < iso.numLevels_max[ipISO][nelem]; ipHi++ )
			{
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
						continue;

					/* this is in real units not emissivity*/
					Transitions[ipISO][nelem][ipHi][ipLo].Emis->phots = 
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->Aul*
						StatesElemNEW[nelem][nelem-ipISO][ipHi].Pop*
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->Pesc;

					/* now find line intensity  */
					Transitions[ipISO][nelem][ipHi][ipLo].Emis->xIntensity = 
						Transitions[ipISO][nelem][ipHi][ipLo].Emis->phots*
						Transitions[ipISO][nelem][ipHi][ipLo].EnergyErg;
				}
			}
		}
	}

	/* create emissivity or intensity for hydrogenic species,
	 * first combine/bring balmer series together */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.IonHigh[nelem] == nelem + 1 )
		{
			/* bring nL - n'L' emission together as n-n' emission. */
			for( ipHi=3; ipHi < iso.numLevels_max[ipH_LIKE][nelem]; ipHi++ )
			{
				long index_of_nHi_P;

				/* is ipHi is collapsed level, index_of_nHi_P is ipHi */
				if( N_(ipHi) > iso.n_HighestResolved_max[ipH_LIKE][nelem] )
					index_of_nHi_P = ipHi;
				else
					index_of_nHi_P = iso.QuantumNumbers2Index[ipH_LIKE][nelem][ N_(ipHi) ][1][2];

				/* only need to consider resolved lower level here */
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{	
					long index_of_nLo_S = iso.QuantumNumbers2Index[ipH_LIKE][nelem][ N_(ipLo) ][0][2];

					/* jump out if ipLo is collapsed 
					 * NB this must be up to n_HighestResolved_local and not n_HighestResolved_max */
					if( N_(ipLo) > iso.n_HighestResolved_local[ipH_LIKE][nelem] || N_(ipLo) == N_(ipHi) )
						break;

					if( Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Aul <= iso.SmallA )
						continue;

					/* add everything into nP - n'S, skip if current indices are those levels. */
					if( ipHi == index_of_nHi_P && ipLo == index_of_nLo_S )
						continue;
					else
					{
						/* add resolved line to nP - n'S */
						Transitions[ipH_LIKE][nelem][index_of_nHi_P][index_of_nLo_S].Emis->xIntensity +=
							Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->xIntensity;
						/* zero out the resolved line */
						Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->xIntensity = 0;
						//ASSERT( Transitions[ipH_LIKE][nelem][index_of_nHi_P][index_of_nLo_S].Emis->xIntensity > 0. );
					}
				}
			}
		}
	}

	/* H beta recombination, assuming old case B */
	hbetab = (double)((pow(10.,-20.89 - 0.10612*POW2(phycon.alogte - 4.4)))/
	  phycon.te);
	/* need to pass this assert if CaBo is to have valid array indices for ipCont */
	/* 06 aug 28, from numLevels_max to _local. */
	/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
	ASSERT( iso.numLevels_max[ipH_LIKE][ipHYDROGEN] > 4 );
	hbetab *= dense.xIonDense[ipHYDROGEN][1]*dense.eden;

	lindst(hbetab, -4861 ,"CaBo",
		1 ,'i',false,
		" this is old case b based on Ferland (1980) PASP ");

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* need to pass this assert if CaBo is to have valid array indices for ipCont */
		/* 06 aug 28, from numLevels_max to _local. */
		/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
		ASSERT( iso.numLevels_max[ipH_LIKE][ipHELIUM] > 4 );
		/* 1640 1640 1640 */
		em = 2.03e-20/(phycon.te70*phycon.te10*phycon.te03);
		em *= dense.xIonDense[ipHELIUM][2]*dense.eden;

		lindst(em,-1640,"CaBo",
			1,'i',false,
			" old prediction of He II 1640, Case B at low densities");

		/* hydrogenic helium */
		/* old prediction of He II 4686, case B */
		em = 2.52e-20/(pow(phycon.te,1.05881));
		em *= dense.xIonDense[ipHELIUM][2]*dense.eden;

		lindst(em,-4686,"CaBo",	1,'i',false,
			   " old prediction of He II 4686, Case B at low densities");
	}

	/* predict case b intensities of hydrogen lines */
	if( LineSave.ipass <= 0 )
	{
		for(nelem=0; nelem<HS_NZ; ++nelem )
		{
			atmdat.lgHCaseBOK[0][nelem] = true;
			atmdat.lgHCaseBOK[1][nelem] = true;
		}
	}
	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* HS_NZ is limit to charge of elements in HS predictions, now 8 == oxygen */
			/* but don't do the minor elements - these were not read in and so should not be
			 * printed - remove equivalent if statement in createdata to read them in */
			if( nelem < HS_NZ && (nelem<2 || nelem>4) )
			{
				int iCase;
				for( iCase=0; iCase<2; ++iCase )
				{
					char chAB[2]={'A','B'};
					char chLab[5]="Ca  ";

					/* adding iCase means start from n=1 for case A, n=2 for Case B,
					 * note that principal quantum number is on physics scale, not C */
					/* 06 aug 28, both of these from numLevels_max to _local. */
					/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
					for( ipLo=1+iCase; ipLo<MIN2(10,iso.n_HighestResolved_max[ipH_LIKE][nelem] + iso.nCollapsed_max[ipH_LIKE][nelem]); ++ipLo )
					{
						for( ipHi=ipLo+1; ipHi< MIN2(25,iso.n_HighestResolved_max[ipH_LIKE][nelem] + iso.nCollapsed_max[ipH_LIKE][nelem]+1); ++ipHi )
						{
							realnum wl;
							double case_b_Intensity;
							long int ipCHi , ipCLo;
							/* Put case b predictions into line stack
							 * NB NB NB each Hummer & Storey case b line must be 
							 * explicitly clobbered by hand in routine final if 
							 * atmdat.lgHCaseBOK[iCase][nelem] flag is set false
							 * since this indicates that we exceeded bounds of table,
							 * DO NOT want to print lines in that case */

							/* first do case b emissivity of balmer lines */

							/* get HS predictions */
							case_b_Intensity = atmdat_HS_caseB( ipHi,ipLo , nelem+1, phycon.te , dense.eden, chAB[iCase] );
							if( case_b_Intensity<=0. )
							{
								atmdat.lgHCaseBOK[iCase][nelem] = false;
								case_b_Intensity = 0.;
							}

							case_b_Intensity *= dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;

							if( iCase==0 && ipLo==1 )
							{
								/* get physical scal prin quant numbers onto cloudy c scale */
								ipCHi = ipHi;
								ipCLo = 0;
							}
							else
							{
								/* get physical scal prin quant numbers onto cloudy c scale */
								ipCHi = ipHi;
								ipCLo = ipLo;
							}

							/* make label either Ca A or Ca B */
							chLab[3] = chAB[iCase];

							/* new treatment is different from old for indices greater than 2. */
							if( ipCHi > 2 )
							{
								if( ipCLo >  2 )
								{
									/* if both indices above two, just treat as nP to n'S transition. */
									ipCHi = iso.QuantumNumbers2Index[ipH_LIKE][nelem][ipCHi][1][2];
									ipCLo = iso.QuantumNumbers2Index[ipH_LIKE][nelem][ipCLo][0][2];
								}
								else if(  ipCLo ==  2 )
								{
									/* treat as nS to 2P transition. */
									ipCHi = iso.QuantumNumbers2Index[ipH_LIKE][nelem][ipCHi][0][2];
								}
								else if( ipCLo ==  1 || ipCLo == 0 )
								{
									/* treat as nP to n'S transition. */
									ipCHi = iso.QuantumNumbers2Index[ipH_LIKE][nelem][ipCHi][1][2];
								}
							}

							/* this is wavelength of interpolated case b from HS tables */
							wl = Transitions[ipH_LIKE][nelem][ipCHi][ipCLo].WLAng;
							atmdat.WaveLengthCaseB[nelem][ipHi][ipLo] = wl;

							lindst(case_b_Intensity,wl,chLab,Transitions[ipH_LIKE][nelem][ipCHi][ipCLo].ipCont,'i',false,
								" case a or case b from Hummer & Storey tables" );
						}
					}
				}
			}

			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */

			/* >>chng 01 sep 10, do two-photon here */
			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				/* two photon emission */
				chIonLbl(chLabel, &Transitions[ipH_LIKE][nelem][ipH2s][ipH1s]);
			}
			linadd( Transitions[ipH_LIKE][nelem][ipH2s][ipH1s].Emis->xIntensity , 0,chLabel,'r',
				"two-photon emission");

			linadd(
				StatesElemNEW[nelem][nelem-ipH_LIKE][ipH2s].Pop*
				iso.TwoNu_induc_dn[ipH_LIKE][nelem]*
				Transitions[ipH_LIKE][nelem][ipH2s][ipH1s].EnergyErg,
				22, chLabel ,'i',
				"induced two-photon emission ");

			/* 06 aug 28, from numLevels_max to _local. */
			/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
			for( ipLo=ipH1s; ipLo < iso.numLevels_max[ipH_LIKE][nelem]-1; ipLo++ )
			{
				/* don't bother with decays to 2p since we set them to zero above */
				if( ipLo==ipH2p )
					continue;

				/* set number of levels we want to print, first is default,
				 * only print real levels, second is set with "print line
				 * iso collapsed" command */
				long int nLoop  = iso.numLevels_max[ipH_LIKE][nelem] - iso.nCollapsed_max[ipH_LIKE][nelem];
				if( prt.lgPrnIsoCollapsed )
					nLoop  = iso.numLevels_max[ipH_LIKE][nelem];

				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* skip the 2s-2p and two photon transitions */
					/* >>chng 02 feb 10, had been explicit levels, change to general ipCont */
					if( Transitions[ipH_LIKE][nelem][ipHi][ipLo].ipCont < 1 )
						continue;

					char chComment[23];
					GenerateTransitionConfiguration( &Transitions[ipH_LIKE][nelem][ipHi][ipLo], chComment);
					fixit(); // put chComment instead of the below, cant just punt chComment there.
					PutLine(&Transitions[ipH_LIKE][nelem][ipHi][ipLo],
						"predicted line, all processes included");
				}
			}
		}
	}
	
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_hydro returns\n" );
	}
	return;
}
