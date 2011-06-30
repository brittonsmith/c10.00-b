/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_general put general information and energetics into line intensity stack */
/*GetMaxhLine find the strongest heating line */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "hydrogenic.h"
#include "dense.h"
#include "thermal.h"
#include "continuum.h"
#include "geometry.h"
#include "dynamics.h"
#include "rt.h"
#include "iso.h"
#include "rfield.h"
#include "trace.h"
#include "ionbal.h"
#include "lines_service.h"
#include "radius.h"
#include "lines.h"
/*GetMaxhLine find the strongest heating line */
STATIC void GetMaxhLine(void);

void lines_general(void)
{
	long int i, 
	  ipHi, 
	  ipLo, 
	  nelem, 
	  ipnt;

	double 
	  hbetac, 
	  HeatMetal ,
	  ee511, 
	  hlalph;

	DEBUG_ENTRY( "lines_general()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_general called\n" );
	}

	i = StuffComment( "general properties" );
	linadd( 0., (realnum)i , "####", 'i',
		" start of general properties");

	/* total H-beta from multi-level atom */
	nelem = ipHYDROGEN;
	ipLo = ipH2p;

	// this can be changed with the atom levels command but must be at
	// least 3
	ASSERT( iso.n_HighestResolved_max[ipH_LIKE][nelem] >=3 );

	if( iso.n_HighestResolved_max[ipH_LIKE][nelem] >= 4 )
	{
		ipHi = ipH4p;
		hbetac = 
		(Transitions[ipH_LIKE][nelem][ipH4p][ipH2s].Emis->Aul * 
		Transitions[ipH_LIKE][nelem][ipH4p][ipH2s].Emis->Pesc *
		StatesElemNEW[nelem][nelem-ipH_LIKE][ipH4p].Pop +
		Transitions[ipH_LIKE][nelem][ipH4s][ipH2p].Emis->Aul *
		Transitions[ipH_LIKE][nelem][ipH4s][ipH2p].Emis->Pesc *
		StatesElemNEW[nelem][nelem-ipH_LIKE][ipH4s].Pop +
		Transitions[ipH_LIKE][nelem][ipH4d][ipH2p].Emis->Aul *
		Transitions[ipH_LIKE][nelem][ipH4d][ipH2p].Emis->Pesc *
		StatesElemNEW[nelem][nelem-ipH_LIKE][ipH4d].Pop) *
		Transitions[ipH_LIKE][nelem][ipHi][ipLo].EnergyErg;
	}
	else
	{
		// atom levels command does not allow < 3
		ASSERT( iso.n_HighestResolved_max[ipH_LIKE][nelem] == 3 );
		ipHi = 6;
		hbetac = 
			(Transitions[ipH_LIKE][nelem][ipHi][ipH2s].Emis->Aul * 
			Transitions[ipH_LIKE][nelem][ipHi][ipH2s].Emis->Pesc +
			Transitions[ipH_LIKE][nelem][ipHi][ipH2p].Emis->Aul *
			Transitions[ipH_LIKE][nelem][ipHi][ipH2p].Emis->Pesc ) *
			StatesElemNEW[nelem][nelem-ipH_LIKE][ipHi].Pop *
			Transitions[ipH_LIKE][nelem][ipHi][ipLo].EnergyErg;
	}

	/* these lines added to outlin in metdif - following must be false 
	 * this passes array index for line energy in continuum mesh - in rest
	 * of code this is set by a previous call to PntForLine, this index
	 * is on the f not c scale */
	rt.fracin = Transitions[ipH_LIKE][nelem][ipHi][ipH2s].Emis->FracInwd;
	lindst(hbetac,Transitions[ipH_LIKE][nelem][ipHi][ipH2s].WLAng,"TOTL",
		Transitions[ipH_LIKE][nelem][ipHi][ipH2s].ipCont,'i',false,
		   " H I Balmer beta predicted by model atom " );
	rt.fracin = 0.5;

	if( iso.n_HighestResolved_max[ipH_LIKE][nelem] < 4 )
	{
		// we need to have something for Hb "H  1" and "Inwd"
		lindst(hbetac,Transitions[ipH_LIKE][nelem][ipHi][ipH2s].WLAng,"H  1",
			Transitions[ipH_LIKE][nelem][ipHi][ipH2s].ipCont,'i',false,
			" H I Balmer beta predicted by model atom " );
		// we need to have something for Hb "H  1" and "Inwd"
		lindst(hbetac/2.,Transitions[ipH_LIKE][nelem][ipHi][ipH2s].WLAng,"Inwd",
			Transitions[ipH_LIKE][nelem][ipHi][ipH2s].ipCont,'i',false,
			" H I Balmer beta predicted by model atom " );
	}

	/* total Ly-a from multi-level atom */
	ipHi = ipH2p;
	ipLo = ipH1s;
	hlalph = 
		Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Aul* 
		StatesElemNEW[nelem][nelem-ipH_LIKE][ipHi].Pop*
		Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->Pesc*
		Transitions[ipH_LIKE][nelem][ipHi][ipLo].EnergyErg;

	rt.fracin = Transitions[ipH_LIKE][nelem][ipHi][ipLo].Emis->FracInwd;
	lindst(hlalph,Transitions[ipH_LIKE][nelem][ipHi][ipLo].WLAng,"TOTL",
		Transitions[ipH_LIKE][nelem][ipHi][ipLo].ipCont,'i',false ,
		" H I Lya predicted from model atom ");
	rt.fracin = 0.5;

	/* this entry only works correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		linadd(continuum.totlsv/radius.dVeffAper,0,"Inci",'i',
		       "total luminosity in incident continuum");
		/* ipass is flag to indicate whether to only set up line array
		 * (ipass=0) or actually evaluate lines intensities (ipass=1) */
		if( LineSave.ipass > 0 )
		{
			continuum.totlsv = 0.;
		}
	}

	linadd(thermal.htot,0,"TotH",'i',
		"  total heating, all forms, information since individuals added later ");

	linadd(thermal.ctot,0,"TotC",'i',
		"  total cooling, all forms, information since individuals added later ");

	linadd(thermal.heating[0][0],0,"BFH1",'h',
		"  hydrogen photoionization heating, ground state only ");

	linadd(thermal.heating[0][1],0,"BFHx",'h',
		"  net hydrogen photoionization heating less rec cooling, all excited states normally zero, positive if excited states are net heating ");

	linadd(thermal.heating[0][22],0,"Line",'h',
		"  heating due to induced lines absorption of continuum ");
	if( thermal.htot > 0. )
	{
		if( thermal.heating[0][22]/thermal.htot > thermal.HeatLineMax )
		{
			thermal.HeatLineMax = (realnum)(thermal.heating[0][22]/thermal.htot);
			/* finds the strongest heating line */
			GetMaxhLine();
		}
	}

	linadd(thermal.heating[1][0]+thermal.heating[1][1]+thermal.heating[1][2],0,"BFHe",'h',
	  "  total helium photoionization heating, all stages ");

	HeatMetal = 0.;
	/* some sums that will be printed in the stack */
	for( nelem=2; nelem<LIMELM; ++nelem)
	{
		/* we now have final solution for this element */
		for( i=dense.IonLow[nelem]; i < dense.IonHigh[nelem]; i++ )
		{
			ASSERT( i < LIMELM );
			/* total metal photo heating for LINES */
			HeatMetal += thermal.heating[nelem][i];
		}
	}

	linadd(HeatMetal,0,"TotM",'h',
		"  total heavy element photoionization heating, all stages ");

	linadd(thermal.heating[0][21],0,"pair",'h',
		"  heating due to pair production ");

	/* ipass is flag to indicate whether to only set up line array
	 * (ipass=0) or actually evaluate lines intensities (ipass=1) */
	if( LineSave.ipass > 0 )
	{
		/* this will be max local heating due to bound compton */
		ionbal.CompHeating_Max = MAX2( ionbal.CompHeating_Max , ionbal.CompRecoilHeatLocal/thermal.htot);
	}
	else
	{
		ionbal.CompHeating_Max = 0.;
	}

	linadd(ionbal.CompRecoilHeatLocal,0,"Cbnd",'h',
		"  heating due to bound compton scattering ");

	linadd(rfield.cmheat,0,"ComH",'h',
		"  Compton heating ");

	linadd(CoolHeavy.tccool,0,"ComC",'c',
		"  total Compton cooling ");

	/* record max local heating due to advection */
	dynamics.HeatMax = MAX2( dynamics.HeatMax , dynamics.Heat() /thermal.htot );
	/* record max local cooling due to advection */
	dynamics.CoolMax = MAX2( dynamics.CoolMax , dynamics.Cool() /thermal.htot );

	linadd(dynamics.Cool()  , 0 , "advC" , 'i',
		"  cooling due to advection " );

	linadd(dynamics.Heat() , 0 , "advH" , 'i' ,
		"  heating due to advection ");

	linadd( thermal.char_tran_heat ,0,"CT H",'h',
		" heating due to charge transfer ");

	linadd( thermal.char_tran_cool ,0,"CT C",'c',
		" cooling due to charge transfer ");

	linadd(thermal.heating[1][6],0,"CR H",'h',
		" cosmic ray heating ");

	linadd(thermal.heating[0][20],0,"extH",'h',
		" extra heat added to this zone, from HEXTRA command ");

	linadd(CoolHeavy.cextxx,0,"extC",'c',
		" extra cooling added to this zone, from CEXTRA command ");

	// 511 keV annihilation line, counts as recombination line since
	// neither heating nor cooling, but does remove energy
	ee511 = (dense.gas_phase[ipHYDROGEN] + 4.*dense.gas_phase[ipHELIUM])*ionbal.PairProducPhotoRate[0]*2.*8.20e-7;
	PntForLine(2.427e-2,"e-e+",&ipnt);
	lindst(ee511,(realnum)2.427e-2,"e-e+",ipnt,'r',true,
		" 511keV annihilation line " );

	linadd(CoolHeavy.expans,0,"Expn",'c',
		"  expansion cooling, only non-zero for wind ");

	linadd(iso.RadRecCool[ipH_LIKE][ipHYDROGEN],0,"H FB",'i',
		"  H radiative recombination cooling ");

	linadd(MAX2(0.,iso.FreeBnd_net_Cool_Rate[ipH_LIKE][ipHYDROGEN]),0,"HFBc",'c',
		"  net free-bound cooling ");

	linadd(MAX2(0.,-iso.FreeBnd_net_Cool_Rate[ipH_LIKE][ipHYDROGEN]),0,"HFBh",'h',
		"  net free-bound heating ");

	linadd(iso.RecomInducCool_Rate[ipH_LIKE][ipHYDROGEN],0,"Hind",'c',
		"  cooling due to induced rec of hydrogen ");

	linadd(CoolHeavy.cyntrn,0,"Cycn",'c',
		"  cyclotron cooling ");

	// cooling due to database species
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		// Species label may be too long for linadd
		char chLabel[5];
		strncpy( chLabel , Species[ipSpecies].chLabel , 4 );
		chLabel[4] = '\0';
		// this is information, 'i', since individual lines
		// have been added as cooling or heating
		linadd(Species[ipSpecies].CoolTotal,0, chLabel,'i',
			" net cooling due to database species");
	}

	return;
}

/*GetMaxhLine find the strongest heating line */
STATIC void GetMaxhLine(void)
{
	long int i;
	double strong;

	DEBUG_ENTRY( "GetMaxhLine()" );

	/* routine called to find which is the strongest heating line */
	strong = 0.;

	/* possible for levlmax to remain 0 if induced processes turned off */
	thermal.levlmax = 0;

	for( i=1; i <= nLevel1; i++ )
	{
		if( TauLines[i].Coll.heat > strong )
		{
			strong = TauLines[i].Coll.heat;
			thermal.levlmax = 1;
			thermal.ipHeatlmax = i;
		}
	}

	for( i=0; i < nWindLine; i++ )
	{
		if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
		{
			if( TauLine2[i].Coll.heat > strong )
			{
				strong = TauLine2[i].Coll.heat;
				thermal.levlmax = 2;
				thermal.ipHeatlmax = i;
			}
		}
	}

	for( i=0; i < nHFLines; i++ )
	{
		if( HFLines[i].Coll.heat > strong )
		{
			strong = HFLines[i].Coll.heat;
			thermal.levlmax = 3;
			thermal.ipHeatlmax = i;
		}
	}

	/* external database lines */
	for( i=0; i < linesAdded2; i++ )
	{
		if(dBaseLines[i].tran->Coll.heat > strong )
		{
			strong = dBaseLines[i].tran->Coll.heat;
			thermal.levlmax = 4;
			thermal.ipHeatlmax = i;
		}
	}

	fixit();  // all other line stacks need to be included here.
	// can we just sweep over line stack?  Is that ready yet?

	return;
}
