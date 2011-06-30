/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_lv1_k_zn place lines of elements potassium and heavier into lines storage stack */
#include "cddefines.h"
#include "cddrive.h"
#include "coolheavy.h"
#include "ca.h"
#include "fe.h"
#include "rfield.h"
#include "dense.h"
#include "phycon.h"
#include "radius.h"
#include "taulines.h"
#include "trace.h"
#include "lines_service.h"
#include "rt.h"
#include "atomfeii.h"
#include "lines.h"

void lines_lv1_k_zn(void)
{
	long int i, 
	  ipnt,
	  ilo,
	  ihi;

	double c10, 
	  c14, 
	  eff, 
	  fela, 
	  r14;

	DEBUG_ENTRY( "lines_lv1_k_zn()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_k_zn called\n" );
	}

	PutLine(&TauLines[ipKI7745],
		"  potassium K I 7745 ");

	PutLine(&TauLines[ipxK03462],
		"  [K III] 4.62 microns ");

	PutLine(&TauLines[ipxK04598],
		"  [KIV] 5.983 min  ");

	PutLine(&TauLines[ipxK04154],
		"  [KIV] 15.39 micron ");

	PutLine(&TauLines[ipxK06882],
		"  [KVI] 8.823 micron ");

	PutLine(&TauLines[ipxK06557],
		"  [KVI]  5.575 micron ");

	PutLine(&TauLines[ipxK07319],
		"  [K VII] 3.189 microns ");

	PutLine(&TauLines[ipxK11425],
		" K 11 4249.99A ");

	PutLine(&TauLines[ipCaI4228],
		"  calcium Ca I 4228 ");

	linadd(ca.Cakh,3933,"Ca 2",'i',
		" coll excited calcium k+h " );

	linadd(ca.Cair,8579,"Ca 2",'i' ,
		" infrared triplet ");

	linadd(ca.c7306,7306,"Ca 2",'i',
		" forbidden lines, 7291+7324 together " );

	linadd(ca.dCakh,3933,"Phot",'i' ,
		" fraction H Ly-alpha destruction of excited levels ");

	linadd(ca.dCaf12,7306,"Phot",'i' ,
		" fraction H Ly-alpha destruction of excited levels ");

	PntForLine(3934.,"Ca2K",&ipnt);
	lindst(ca.Cak,3934,"Ca2K",ipnt,'c',true,
		" individual lines from five level atom");


	PntForLine(3969.,"Ca2H",&ipnt);
	lindst(ca.Cah,3969,"Ca2H",ipnt,'c',true,
		" individual lines from five level atom" );


	PntForLine(8498.,"Ca2X",&ipnt);
	lindst(ca.Cax,8498,"Ca2X",ipnt,'c',true,
		" individual lines from five level atom " );


	PntForLine(8542.,"Ca2Y",&ipnt);
	lindst(ca.Cay,8542,"Ca2Y",ipnt,'c',true,
		"  individual lines from five level atom" );


	PntForLine(8662.,"Ca2Z",&ipnt);
	lindst(ca.Caz,8662,"Ca2Z",ipnt,'c',true,
		" individual lines from five level atom" );


	PntForLine(7291.,"CaF1",&ipnt);
	lindst(ca.Caf1,7291,"CaF1",ipnt,'c',true,
		" individual lines from five level atom" );


	PntForLine(7324.,"CaF2",&ipnt);
	lindst(ca.Caf2,7324,"CaF2",ipnt,'c',true,
		" individual lines from five level atom" );

	eff = dense.eden*dense.xIonDense[ipCALCIUM][2]*5.4e-21/(phycon.te/
	  phycon.te10/phycon.te10);
	linadd(eff,3933,"Rec ",'i',
		" recombination contribution to CaII emission" );

	PutLine(&TauLines[ipTCa3],
		"  Ca IV 3.2 micron ");

	PutLine(&TauLines[ipTCa4],
		"  Ca V 4.16, 11.47 micron");

	PutLine(&TauLines[ipTCa12],
		"  Ca V 4.16, 11.47 micron ");

	PntForLine(6087.,"Ca 5",&ipnt);
	lindst(ca.Ca6087,6087,"Ca 5",ipnt,'c',true ,
		"  Ca V optical and UV lines, collisional excitation, 3-level atom");

	PntForLine(5311.,"Ca 5",&ipnt);
	lindst(ca.c5311,5311,"Ca 5",ipnt,'c',true ,
		" Ca V optical and UV lines, collisional excitation, 3-level atom");

	PntForLine(2414.,"Ca 5",&ipnt);
	lindst(ca.c2414,2414,"Ca 5",ipnt,'c',true ,
		"  Ca V optical and UV lines, collisional excitation, 3-level atom");


	PntForLine(3997.,"Ca 5",&ipnt);
	lindst(ca.c3997,3997,"Ca 5",ipnt,'c',true,
		" Ca V optical and UV lines, collisional excitation, 3-level atom" );

	PutLine(&TauLines[ipCa0741],
		"  [Ca VII] 4.09 microns" );

	PutLine(&TauLines[ipCa0761],
		"  [Ca VII] 6.15 microns " );


	PntForLine(5620.,"Ca 7",&ipnt);
	lindst(ca.Ca5620,5620,"Ca 7",ipnt,'c',true,
		" Ca VII optical and UV lines, collisional excitation, 3-level atom" );


	PntForLine(4941.,"Ca 7",&ipnt);
	lindst(ca.Ca4941,4941,"Ca 7",ipnt,'c',true,
		" Ca VII optical and UV lines, collisional excitation, 3-level atom" );


	PntForLine(2112.,"Ca 7",&ipnt);
	lindst(ca.Ca2112,2112,"Ca 7",ipnt,'c',true ,
		" Ca VII optical and UV lines, collisional excitation, 3-level atom");


	PntForLine(3688.,"Ca 7",&ipnt);
	lindst(ca.Ca3688,3688,"Ca 7",ipnt,'c',true,
		" Ca VII optical and UV lines, collisional excitation, 3-level atom" );

	PutLine(&TauLines[ipCa08232],
		"  [Ca VIII]  2.32 microns, A Saraph and Strey ");

	PutLine(&TauLines[ipCa12333],
		"  [Ca 12] 3328.78A ");

	PutLine(&TauLines[ipTCa302],
		"  Ca 18 Li seq 2s2p, 302, 345 separate ");

	PutLine(&TauLines[ipTCa345],
		"  Ca 18 Li seq 2s2p, 302, 345 separate ");

	PutLine(&TauLines[ipTCa19],
		"  Ca 18 Li seq 2s3p,  ");


	PntForLine(22.08e4,"Sc 2",&ipnt);
	lindst(CoolHeavy.Sc22p08m,22.08e4,"Sc 2",ipnt,'c',true,
		" Sc II 2.08 (1-3) " );


	PntForLine(24.1e4,"Sc 2",&ipnt);
	lindst(CoolHeavy.Sc24p1m,24.1e4,"Sc 2",ipnt,'c',true,
		" Sc II 4.1 micron (1-2)" );


	PntForLine(24.2e4,"Sc 2",&ipnt);
	lindst(CoolHeavy.Sc24p2m,24.2e4,"Sc 2",ipnt,'c',true,
		"  Sc II 4.22 (2-3)" );


	PntForLine(3933.,"Sc 3",&ipnt);
	lindst(CoolHeavy.Sc33936,3933,"Sc 3",ipnt,'c',true,
		" Sc III 3936" );

	PutLine(&TauLines[ipSc05231],
		"  [Sc V] 1.46 microns ");


	PntForLine(5054.,"Sc 6",&ipnt);
	lindst(CoolHeavy.Sc45058,5054,"Sc 6",ipnt,'c',true ,
		" Sc VI 5054 (1-2)");


	PntForLine(3592.,"Sc 6",&ipnt);
	lindst(CoolHeavy.Sc43595,3592,"Sc 6",ipnt,'c',true,
		" Sc VI 3595 (2-3)" );


	PntForLine(2100.,"Sc 6",&ipnt);
	lindst(CoolHeavy.Sc42100,2100,"Sc 6",ipnt,'c',true,
		"  Sc VI 2100 (1-3)" );

	PutLine(&TauLines[ipSc13264],
		"  [Sc 13] 2637.97A");


	PntForLine(1.21e4,"Ti 3",&ipnt);
	lindst(CoolHeavy.Ti3l21,1.21e4,"Ti 3",ipnt,'c',true,
		" Ti III 1.21 micron, (actually multiplet) 2-1 transition from model atom " );


	PntForLine(9594.,"Ti 3",&ipnt);
	lindst(CoolHeavy.Ti3l31,9594,"Ti 3",ipnt,'c',true,
		" Ti III 9594, 3-1 transition, (actually multiplet) from model atom" );

	PntForLine(4.57e4,"Ti 3",&ipnt);
	lindst(CoolHeavy.Ti3l32,4.57e4,"Ti 3",ipnt,'c',true,
		" Ti III 4.57 micron, 3-2 transition, (actually multiplet) from model atom" );

	PutLine(&TauLines[ipTi06172],
		"  [Ti VI] 1.72 microns ");

	PutLine(&TauLines[ipTi14212],
		"  [Ti XIV] 2117.79 ");


	PntForLine(8823.,"V  3",&ipnt);
	lindst(CoolHeavy.V38830,8823,"V  3",ipnt,'c',true ,
		"  V III 8823 ");


	PntForLine(8507.,"V  3",&ipnt);
	lindst(CoolHeavy.V38507,8507,"V  3",ipnt,'c',true,
		"  V III 8507" );


	PntForLine(7735.,"V  4",&ipnt);
	lindst(CoolHeavy.V47741,7735,"V  4",ipnt,'c',true,
		"  V IV 7741 1-3" );


	PntForLine(9489.,"V  4",&ipnt);
	lindst(CoolHeavy.V49496,9489,"V  4",ipnt,'c',true,
		" V IV 9496 2-1 " );


	PntForLine(4.19e4,"V  4",&ipnt);
	lindst(CoolHeavy.V44p2m,4.19e4,"V  4",ipnt,'c',true,
		"  V IV 4.19 micron 3-2" );

	PutLine(&TauLines[ipVa07130],
		"  [V VII] 1.304 microns ");

	PutLine(&TauLines[ipVa15172],
		" [V 15] 1721.38 ");

	PntForLine(5828.,"Cr 3",&ipnt);
	lindst(CoolHeavy.Cr3l21,5828,"Cr 3",ipnt,'c',true,
		" [CrIII] multiplet blend at 5828A" );

	PntForLine(7267.,"Cr 4",&ipnt);
	lindst(CoolHeavy.Cr4l21,7267,"Cr 4",ipnt,'c',true,
		" [CrIV] 2 - 1 multiplet blend at 7272" );


	PntForLine(6801.,"Cr 4",&ipnt);
	lindst(CoolHeavy.Cr4l31,6801,"Cr 4",ipnt,'c',true,
		" [CrIV] 3 - 1 multiplet blend at 6806" );


	PntForLine(7979.,"Cr 5",&ipnt);
	lindst(CoolHeavy.Cr5l21,7979,"Cr 5",ipnt,'c',true,
		" [CrV] 2 - 1 multiplet blend at 7985" );

	PntForLine(6577.,"Cr 5",&ipnt);
	lindst(CoolHeavy.Cr5l31,6577,"Cr 5",ipnt,'c',true,
		"  [CrV] 3 - 1 multiplet blend at 6582" );


	PntForLine(3.75e4,"Cr 5",&ipnt);
	lindst(CoolHeavy.Cr5l32,3.75e4,"Cr 5",ipnt,'c',true,
		" [CrV] 3 - 2 multiplet blend at 3.75 microns " );

	PutLine(&TauLines[ipCr08101],
		"  [Cr VIII] 1.01 microns ");

	PutLine(&TauLines[ipCr16141],
		"  [Cr 16] 1410.60 ");

	PutLine(&TauLines[ipxMn0979],
		" [Mn IX] 7968.5 A ");

	PutLine(&TauLines[ipxMn1712],
		" [Mn 17] 1169.59 ");

	/* bob Rubin's UV line
	 * f2 = dense.xIonDense(26,4)*sexp(50 764./te)*0.45*cdsqte/6.*7.01e-12
	 * call linadd( f2 , 2837 , 'BobR' , 'i')
	 * f2 = dense.xIonDense(26,4)*sexp(55 989./te)*0.384*cdsqte/6.*7.74e-12
	 * call linadd( f2 , 2568 , 'BobR' , 'i') */

	/* iron */

	PutLine(&TauLines[ipFe1_24m],
		"  Fe 1 24m ");

	PutLine(&TauLines[ipFe1_35m],
		"  Fe 1 35m ");

	PutLine(&TauLines[ipFe1_54m],
		"  Fe 1 54m ");

	PutLine(&TauLines[ipFe1_111m],
		"  Fe 1 111m ");

	PutLine(&TauLines[ipFeI3884],
		"  Fe 1 3884 ");

	PutLine(&TauLines[ipFeI3729],
		"  Fe 1 3729 ");

	PutLine(&TauLines[ipFeI3457],
		"  Fe 1 3457 ");

	PutLine(&TauLines[ipFeI3021],
		"  Fe 1 3021 ");

	PutLine(&TauLines[ipFeI2966],
		"  Fe 1 2966 ");

	linadd(MAX2(0.,FeII.Fe2_large_cool+FeII.Fe2_UVsimp_cool),0,"Fe2c",'c' ,
		"total of all Fe 2 cooling, both simple UV and large atom together ");

	linadd(MAX2(0.,-FeII.Fe2_large_cool-FeII.Fe2_UVsimp_cool),0,"Fe2h",'h' ,
		"total of all Fe 2 heating, both simple UV and large atom together ");

	linadd(FeII.for7,4300,"Fe 2",'i' ,
		" Fe 2 forbidden 2-1 transition from Netzer's atom ");

	PutLine(&TauLines[ipTuv3],
		" 2400 in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTr48],
		" 6200 in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTFe16],
		" 1080 in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTFe26],
		" 1500  in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTFe34],
		" 11500 in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTFe35],
		" 2500 in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTFe46],
		" 2300 in simple Wills, Netzer, Wills FeII");
	PutLine(&TauLines[ipTFe56],
		" 8900 in simple Wills, Netzer, Wills FeII");

	/* option to save all intensities predicted by large FeII atom,
	 * code is in FeIILevelPops */
	FeIIAddLines();
	/* we were called by lines, and we want to zero out Fe2SavN */
	for( long ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( long ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* only evaluate real transitions */
			if( Fe2LevN[ipHi][ipLo].ipCont > 0 ) 
				PutLine( &Fe2LevN[ipHi][ipLo] ," Fe II emission" );
		}
	}

	/* emission bands from the model Fe II atom */
	for( i=0; i < nFeIIBands; i++ )
	{
		double SumBandInward;
		/* [i][0] is center wavelength, [i][1] and [i][2] are upper and
		 * lower bounds in Angstroms.  These are set in FeIIZero 
		 * units are erg s-1 cm-3 */
		eff = FeIISumBand(FeII_Bands[i][1],FeII_Bands[i][2],
			&SumBandInward);

		linadd(eff,FeII_Bands[i][0],"Fe2b",'i' ,
			" total Fe II emission in Fe II bands, as defined in FeII_bands.ini ");
		linadd(SumBandInward,FeII_Bands[i][0],"Inwd",'i' ,
			" inward Fe II emission in Fe II bands, as defined in FeII_bands.ini ");
	}

	// integrate the pseudo continuum of FeII emission
	if( LineSave.ipass > 0 )
	{
		// initialize
		if( nzone == 1 )
		{
			for( i=0; i < nFeIIConBins; i++ )
			{
				/* initialize arrays */
				FeII_Cont[i][1] = 0.;
				FeII_Cont[i][2] = 0.;
			}
		}

		// integrate
		for( i=0; i < nFeIIConBins; i++ )
		{
			double SumBandInward;
			/* [i][0] is total intensity in cell, [i][1] and [i][2] are lower and
			 * upper bounds in Angstroms.  these are set in FeIIZero *
			 * find total emission from large FeII atom, integrated over band */
			double TotalFeII = FeIISumBand(FeII_Cont[i][0],FeII_Cont[i+1][0],
				&SumBandInward);
			FeII_Cont[i][1] += (realnum)(SumBandInward*radius.dVeffAper);
			FeII_Cont[i][2] += (realnum)(MAX2(0.,TotalFeII-SumBandInward)*radius.dVeffAper);
			/*fprintf(ioQQQ,"DEBUG feii\t%li\t%.2e\n", i, FeII_Cont[i][0]);*/
		}
	}

	PutLine(&TauLines[ipT191],
		"  anomalous Fe 2 transition at 1787, RMT 191");

	linadd(fe.Fe3CoolTot,0,"Fe3c",'c' ,
		" chng 05 dec 16, FeIII code created by Kevin Blagrave  Fe3c 0 - total cooling due to 14-level Fe 3 atom ");

	/* Fe 3 14-level atom 
	 * following from print statements within loop */
	/* Fe 3 22.92m from Blagrave 14-level atom */
 	/* Fe 3 13.53m from Blagrave 14-level atom */
 	/* Fe 3 33.03m from Blagrave 14-level atom */
 	/* Fe 3 10.72m from Blagrave 14-level atom */
 	/* Fe 3 20.15m from Blagrave 14-level atom */
 	/* Fe 3 51.67m from Blagrave 14-level atom */
 	/* Fe 3 9.732m from Blagrave 14-level atom */
 	/* Fe 3 16.91m from Blagrave 14-level atom */
 	/* Fe 3 34.66m from Blagrave 14-level atom */
 	/* Fe 3 105.3m from Blagrave 14-level atom */
 	/* Fe 3  5152A from Blagrave 14-level atom */
 	/* Fe 3  5271A from Blagrave 14-level atom */
 	/* Fe 3  5356A from Blagrave 14-level atom */
 	/* Fe 3  5412A from Blagrave 14-level atom */
 	/* Fe 3  5440A from Blagrave 14-level atom */
 	/* Fe 3  4986A from Blagrave 14-level atom */
 	/* Fe 3  5097A from Blagrave 14-level atom */
 	/* Fe 3  5177A from Blagrave 14-level atom */
 	/* Fe 3  5230A from Blagrave 14-level atom */
 	/* Fe 3  5256A from Blagrave 14-level atom */
 	/* Fe 3 15.47m from Blagrave 14-level atom */
 	/* Fe 3  4925A from Blagrave 14-level atom */
 	/* Fe 3  5033A from Blagrave 14-level atom */
 	/* Fe 3  5111A from Blagrave 14-level atom */
 	/* Fe 3  5162A from Blagrave 14-level atom */
 	/* Fe 3  5188A from Blagrave 14-level atom */
 	/* Fe 3 11.16m from Blagrave 14-level atom */
 	/* Fe 3 40.04m from Blagrave 14-level atom */
 	/* Fe 3  4881A from Blagrave 14-level atom */
 	/* Fe 3  4988A from Blagrave 14-level atom */
 	/* Fe 3  5064A from Blagrave 14-level atom */
 	/* Fe 3  5114A from Blagrave 14-level atom */
 	/* Fe 3  5139A from Blagrave 14-level atom */
 	/* Fe 3 9.282m from Blagrave 14-level atom */
 	/* Fe 3 23.21m from Blagrave 14-level atom */
 	/* Fe 3 55.20m from Blagrave 14-level atom */
 	/* Fe 3  4833A from Blagrave 14-level atom */
 	/* Fe 3  4937A from Blagrave 14-level atom */
 	/* Fe 3  5012A from Blagrave 14-level atom */
 	/* Fe 3  5061A from Blagrave 14-level atom */
 	/* Fe 3  5085A from Blagrave 14-level atom */
 	/* Fe 3 7.789m from Blagrave 14-level atom */
 	/* Fe 3 15.69m from Blagrave 14-level atom */
 	/* Fe 3 25.79m from Blagrave 14-level atom */
 	/* Fe 3 48.41m from Blagrave 14-level atom */
 	/* Fe 3  4714A from Blagrave 14-level atom */
 	/* Fe 3  4813A from Blagrave 14-level atom */
 	/* Fe 3  4884A from Blagrave 14-level atom */
 	/* Fe 3  4931A from Blagrave 14-level atom */
 	/* Fe 3  4954A from Blagrave 14-level atom */
 	/* Fe 3 5.543m from Blagrave 14-level atom */
 	/* Fe 3 8.638m from Blagrave 14-level atom */
 	/* Fe 3 11.01m from Blagrave 14-level atom */
 	/* Fe 3 13.76m from Blagrave 14-level atom */
 	/* Fe 3 19.22m from Blagrave 14-level atom */
 	/* Fe 3  4659A from Blagrave 14-level atom */
 	/* Fe 3  4755A from Blagrave 14-level atom */
 	/* Fe 3  4825A from Blagrave 14-level atom */
 	/* Fe 3  4870A from Blagrave 14-level atom */
 	/* Fe 3  4893A from Blagrave 14-level atom */
 	/* Fe 3 4.859m from Blagrave 14-level atom */
 	/* Fe 3 7.085m from Blagrave 14-level atom */
 	/* Fe 3 8.608m from Blagrave 14-level atom */
 	/* Fe 3 10.20m from Blagrave 14-level atom */
 	/* Fe 3 12.92m from Blagrave 14-level atom */
 	/* Fe 3 39.41m from Blagrave 14-level atom */
 	/* Fe 3  4608A from Blagrave 14-level atom */
 	/* Fe 3  4702A from Blagrave 14-level atom */
 	/* Fe 3  4770A from Blagrave 14-level atom */
 	/* Fe 3  4814A from Blagrave 14-level atom */
 	/* Fe 3  4836A from Blagrave 14-level atom */
 	/* Fe 3 4.356m from Blagrave 14-level atom */
 	/* Fe 3 6.063m from Blagrave 14-level atom */
 	/* Fe 3 7.146m from Blagrave 14-level atom */
 	/* Fe 3 8.208m from Blagrave 14-level atom */
 	/* Fe 3 9.884m from Blagrave 14-level atom */
 	/* Fe 3 20.34m from Blagrave 14-level atom */
 	/* Fe 3 42.06m from Blagrave 14-level atom */
 	/* Fe 3  4574A from Blagrave 14-level atom */
 	/* Fe 3  4668A from Blagrave 14-level atom */
 	/* Fe 3  4734A from Blagrave 14-level atom */
 	/* Fe 3  4778A from Blagrave 14-level atom */
 	/* Fe 3  4800A from Blagrave 14-level atom */
 	/* Fe 3 4.077m from Blagrave 14-level atom */
 	/* Fe 3 5.535m from Blagrave 14-level atom */
 	/* Fe 3 6.423m from Blagrave 14-level atom */
 	/* Fe 3 7.269m from Blagrave 14-level atom */
 	/* Fe 3 8.554m from Blagrave 14-level atom */
 	/* Fe 3 15.41m from Blagrave 14-level atom */
 	/* Fe 3 25.31m from Blagrave 14-level atom */
 	/* Fe 3 63.56m from Blagrave 14-level atom */
	for( ihi=1; ihi<NLFE3; ++ihi )
	{
		for( ilo=0; ilo<ihi; ++ilo )
		{
			/* emission in these lines */
			PntForLine(fe.Fe3_wl[ihi][ilo],"Fe 3",&ipnt);
#			if 0
			fprintf( ioQQQ,"\t/* FeIII ");
			prt_wl( ioQQQ , (realnum)(fe.Fe3_wl[ihi][ilo]+0.5) );
			fprintf( ioQQQ," from Blagrave 14-level atom */\n" );
#			endif
			lindst( fe.Fe3_emiss[ihi][ilo] , (realnum)(fe.Fe3_wl[ihi][ilo]+0.5) , "Fe 3",ipnt,'c',true,
				" " );
		}
	}

	/*>>chng 05 dec 18, following are now in the above */
	/* sum of 3p and 3g states together */
	/*	linadd(CoolHeavy.c5270,0,"Fe 3",'c' ); */

	/* Fe 3 5270, predictions from Garstang et al 78
	PntForLine(5270.,"Fe 3",&ipnt);
	lindst(CoolHeavy.c5270*0.2090,5270,"Fe 3",ipnt,'c',true );*/

	/* Fe 3 5270, predictions from Garstang et al 78 
	PntForLine(4658.,"Fe 3",&ipnt);
	lindst(CoolHeavy.c5270*0.3667,4658,"Fe 3",ipnt,'c',true ); */

	PutLine(&TauLines[ipT1122]," Fe 3 1122 entire multiplet");

	linadd(fe.Fe4CoolTot,0,"Fe4c",'i',
		" Fe4c 0 - total cooling due to 12-level Fe 4 atom " );


	PntForLine(3096.,"Fe 4",&ipnt);
	lindst(fe.fe40401,3096,"Fe 4",ipnt,'c',true,
		" Fe 4 3096.A, 4-1 and 5-1 transitions together"  );


	PntForLine(2836.,"Fe 4",&ipnt);
	lindst(fe.fe42836,2836,"Fe 4",ipnt,'c',true,
		" Fe 4 2835.7A, 6-1 transition, 4P5/2 - 6S5/2 "  );


	PntForLine(2829.,"Fe 4",&ipnt);
	lindst(fe.fe42829,2829,"Fe 4",ipnt,'c',true,
		"   Fe 4 2829.4A, 7-1 transition, 4P3/2 - 6S5/2"  );


	PntForLine(2567.,"Fe 4",&ipnt);
	lindst(fe.fe42567,2567,"Fe 4",ipnt,'c',true,
		"  Fe 4 2567.6+ 2567.4. 11-1 and 12-1 transitions"  );


	PntForLine(2.774e4,"Fe 4",&ipnt);
	lindst(fe.fe41207,2.774e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.774 microns 12-7 transition "  );


	PntForLine(2.714e4,"Fe 4",&ipnt);
	lindst(fe.fe41206,2.714e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.714 microns 12-6 transition "  );


	PntForLine(2.716e4,"Fe 4",&ipnt);
	lindst(fe.fe41106,2.716e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.716 microns 11-6 transition"  );


	PntForLine(2.806e4,"Fe 4",&ipnt);
	lindst(fe.fe41007,2.806e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.806 microns 10-7 transition "  );


	PntForLine(2.865e4,"Fe 4",&ipnt);
	lindst(fe.fe41008,2.865e4,"Fe 4",ipnt,'c',true ,
		"  Fe 4 2.865 microns 10-8 transition");


	PntForLine(2.836e4,"Fe 4",&ipnt);
	lindst(fe.fe40906,2.836e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.836 microns 9-6 transition" );


	PntForLine(3892.,"Fe 5",&ipnt);
	lindst(CoolHeavy.c3892,3892,"Fe 5",ipnt,'c',true,
		" Fe 5  3892+3839" );

	linadd(CoolHeavy.c5177,0,"Fe 6",'c' ,
		" all of 2G lines together first ");


	PntForLine(5177.,"Fe 6",&ipnt);
	lindst(CoolHeavy.c5177*0.354,5177,"Fe 6",ipnt,'c',true,
		" Fe 6 5177, approximate correct " );

	linadd(fe.Fe7CoolTot,0,"Fe7c",'c' ,
		" Fe7c 0 - total cooling due to n-level Fe 7 atom ");

	/* >>chng 04 nov 04, move to multi-level system */
	for( ilo=0; ilo<NLFE7-1; ++ilo )
	{
		/* must not do 1-0 or 2-1, which are transferred lines */
		for( ihi=MAX2(3,ilo+1); ihi<NLFE7; ++ihi )
		{

			PntForLine(fe.Fe7_wl[ihi][ilo],"Fe 7",&ipnt);
			lindst( fe.Fe7_emiss[ihi][ilo] , (realnum)(fe.Fe7_wl[ihi][ilo]+0.5) , "Fe 7",ipnt,'c',true,
				" emission in these lines" );
		}
	}
#	if 0
	PntForLine(5721.,"Fe 7",&ipnt);
	lindst( fe.Fe7_5721 , 5721 , "Fe 7",ipnt,'c',true,
		" " );

	PntForLine(6601.,"Fe 7",&ipnt);
	lindst( fe.Fe7_6601 , 6601 , "Fe 7",ipnt,'c',true,
		" " );

	PntForLine(3760.,"Fe 7",&ipnt);
	lindst( fe.Fe7_3760 , 3760 , "Fe 7",ipnt,'c',true,
		" " );

	PntForLine(3588.,"Fe 7",&ipnt);
	lindst( fe.Fe7_3588 , 3588 , "Fe 7",ipnt,'c',true,
		" " );
#	endif

	PutLine(&TauLines[ipFe0795],
		"   [Fe 7] 9.51 micron ");

	PutLine(&TauLines[ipFe0778],
		"  [Fe 7] 7.81 micron ");

	/* [Fe 7] 6087 
	PntForLine(6087.,"Fe 7",&ipnt);
	lindst(CoolHeavy.c6087,6087,"Fe 7",ipnt,'c',true );*/

	/* [Fe 7] 5722 
	PntForLine(5722.,"Fe 7",&ipnt);
	lindst(CoolHeavy.Fe5722,5722,"Fe 7",ipnt,'c',true );*/

	PutLine(&TauLines[ipT245],
		"   Be-seq lines ");

	PntForLine(242.,"Fe 7",&ipnt);
	lindst(CoolHeavy.c242,242,"Fe 7",ipnt,'c',true,
		"  Fe 9 242 j=1 slower decay");

	PutLine(&TauLines[ipT352],
		" the E1 transition that can pump [Fe X] ");

	/* optically thin Fe X pumping */
	eff = 69.4/(69.4 + 0.27*dense.cdsqte);

	/* coll excitation of 352 which decays to excited state of 6374
	 * assumes 17/56 (ratio of A's) go to excited state */
	c10 = TauLines[ipT352].Emis->phots*eff*.01676*TauLines[ipT352].EnergyErg*352/6374.;

	/* Fe 10 and Fe 14 from Mason 75
	 * total (coll, pumped)
	 * call linadd( C6374+C10 , 6374 , 'Fe10','i') */
	PutExtra( c10 );


	PutLine(&TauLines[ipFe106375],
		" [Fe 10] 6375, collisions with pumping too ");

	/* collisional contribution
	 * call linadd( C6374 , 6374 , 'Coll','c')
	 * collisions of E1 line, plus pumped by continuum fluorescence
	 * call linadd( C10 , 6374 , ' 352','c')
	 * Fe XI 7892, 6.08 micron 
	PutLine(&TauLines[ipTFe07]);*/

	/** \todo	2	put this line back in!
	 * EFF = 43.6 / (43.6 + 0.27*COLFAC)
	 * R11 = FE(11)*FLUX(IPFE10)*3.122E-12*EFF *
	 *  1 ( PFE11A*ESCINC(T353(1),1E-4)/(1.+17.0/5.3*T353(3) )  +
	 *  2 PFE11B*ESCINC(T353(1)/3.,1E-4)/(1.+11.0/12.3*T353(3) ) )
	 * contribution to Fe11 from continuum fluorescence
	 * call linadd( R11    , 7892 , 'Pump','i')
	 *
	 * [Fe 11] 6.08 microns
	PutLine(&TauLines[ipTFe61]); */

	/* Fe 11 2649 collisional excitation
	PntForLine(2649.,"Fe11",&ipnt);
	lindst(CoolHeavy.c2649,2649,"Fe11",ipnt,'c',true ); */

	/*  Fe 11 1467 collisional excitation 
	PntForLine(1467.,"Fe11",&ipnt);
	lindst(CoolHeavy.c1467,1467,"Fe11",ipnt,'c',true );*/

	linadd(fe.Fe11CoolTot,0,"Fe11",'c' ,
		" >>chng 05 dec 18, add Fe 11  Fe11 0 - total cooling due to 5-level Fe 11 atom ");

	/* Fe 11 5-level atom */
	for( ihi=1; ihi<NLFE11; ++ihi )
	{
		for( ilo=0; ilo<ihi; ++ilo )
		{
			PntForLine(fe.Fe11_wl[ihi][ilo],"Fe11",&ipnt);
			lindst( fe.Fe11_emiss[ihi][ilo] , (realnum)(fe.Fe11_wl[ihi][ilo]+0.5) , "Fe11",ipnt,'c',true,
				" emission in these lines" );
		}
	}

	PntForLine(1242.,"Fe12",&ipnt);
	lindst(CoolHeavy.c1242,1242,"Fe12",ipnt,'c',true ,
		" Fe 12, 1242, 1349 together, collisional excitation");

	PntForLine(2170.,"Fe12",&ipnt);
	lindst(CoolHeavy.c2170,2170,"Fe12",ipnt,'c',true ,
		" Fe 12, 2170, 2406 together, collisional excitation");


	PntForLine(2568.,"Fe12",&ipnt);
	lindst(CoolHeavy.c2568,2568,"Fe12",ipnt,'c',true,
		"  Fe12 2904, 2567, 3567, 3073 together, collisional excitation" );

	/* >>chng 05 dec 18, add Fe 13  */
	/* Fe13 0 - total cooling due to 5-level Fe 13 atom */
	linadd(fe.Fe13CoolTot,0,"Fe13",'c' ,
		   "total cooling due to Fe 13 model atom ");

	/* Fe 13 5-level atom */
	for( ihi=1; ihi<NLFE13; ++ihi )
	{
		for( ilo=0; ilo<ihi; ++ilo )
		{

			PntForLine(fe.Fe13_wl[ihi][ilo],"Fe13",&ipnt);
			lindst( fe.Fe13_emiss[ihi][ilo] , (realnum)(fe.Fe13_wl[ihi][ilo]+0.5) , "Fe13",ipnt,'c',true ,
				" Fe 13 emission");
		}
	}

	/* Fe 14 optically thin in line 344 */
	eff = 60.3/(60.3 + 0.23*dense.cdsqte/4.);
	r14 = dense.xIonDense[ipIRON][14-1]*fe.pfe14*rfield.flux[0][fe.ipfe10-1]*
	  3.75e-12*eff/(1. + 24./.63*
	  TauLines[ipT347].Emis->Pesc)*esc_PRD_1side(TauLines[ipT347].Emis->TauIn,1e-4);

	linadd(CoolHeavy.c5303+r14,5303,"Fe14",'i',
		" total emission in Fe 14 5304");

	PntForLine(5303.,"Fe14",&ipnt);
	lindst(CoolHeavy.c5303,5303,"Coll",ipnt,'c' ,true ,
		" Fe 14 5304  contribution from collisional excitation ");

	lindst(r14,5303,"Pump",ipnt,'r' ,true ,
		" Fe 14 5304  continuum fluorescense ");

	/** \todo	2	put this in */

	c14 = 0.;
	linadd(c14,5303," 347",'c' ,
		" collisional excitation of E1 line ");

	PutLine(&TauLines[ipFe17_17],
		" Fe 17 17.1A M2 line");

	PutLine(&TauLines[ipFe18975],
		" Fe 18 974.86A ");

	PntForLine(7047.,"Fe19",&ipnt);
	lindst(CoolHeavy.c7082,7047,"Fe19",ipnt,'c',true,
		"  O-like Fe19, 3P ground term, 7046.72A vacuum wl, 1328.90A  >>chng 01 aug 10, updated wavelengths   Fe 19 7047 '85 " );


	PntForLine(1328.,"Fe19",&ipnt);
	lindst(CoolHeavy.c1328,1329,"Fe19",ipnt,'c',true,
		" Fe 19 1329" );


	PntForLine(592.,"Fe19",&ipnt);
	lindst(CoolHeavy.c592,592,"Fe19",ipnt,'c',true,
		"  Fe 19 from loulergue et al '85" );


	PntForLine(1118.,"Fe19",&ipnt);
	lindst(CoolHeavy.c1118,1118,"Fe19",ipnt,'c',true ,
		" Fe 19 from loulergue et al '85");

	PutLine(&TauLines[ipTFe13],
		"   next two 3p ground state lines, collisional excitation ");

	PutLine(&TauLines[ipTFe23],
		"  collisional excitation ");

	PutLine(&TauLines[ipTFe20_578],
		"  Fe20 721.40A, 578");
	PutLine(&TauLines[ipTFe20_721],
		" ");

	linadd(
		TauLines[ipFe22_247].Emis->xIntensity+
		TauLines[ipFe22_217].Emis->xIntensity+
		TauLines[ipFe22_348].Emis->xIntensity+
		TauLines[ipFe22_292].Emis->xIntensity+
		TauLines[ipFe22_253].Emis->xIntensity,
		260,"TOTL",'i',
		" Fe 22 845.6A  total intensity of Fe22, all lines in the multiplet " );
	PutLine(&TauLines[ipFe22_247],
		"Fe 22 247");
	PutLine(&TauLines[ipFe22_217],
		"Fe 22 217");
	PutLine(&TauLines[ipFe22_348],
		"Fe 22 348");
	PutLine(&TauLines[ipFe22_292],
		"Fe 22 292");
	PutLine(&TauLines[ipFe22_253],
		"Fe 22 253");

	/*  Fe 23 1909-like 262.6 */
	PntForLine(263.,"Fe23",&ipnt);
	lindst(CoolHeavy.c263,263,"Fe23",ipnt,'c',true,
		"Fe 23 1909-like 262.6" );


	PutLine(&TauLines[ipT192],
		" Fe 24 only 192 of 255, 192 Li seq doublet, collisional excitation ");

	PutLine(&TauLines[ipT255],
		"  Fe 24 255 of 255, 192 Li seq doublet, collisional excitation ");

	PutLine(&TauLines[ipT11],
		"  Fe 24 Li seq 2s3p collisional excitation ");

	/* recombination Ka */
	if( dense.lgElmtOn[ipIRON] )
	{
		/* these lines added to outlin in metdif - following must be false
		 * fela = xLyaHeavy(nelem,nelem)*dense.xIonDense(nelem,nelem+1) */
		fela = Transitions[ipH_LIKE][ipIRON][ipH2p][ipH1s].Emis->xIntensity;
	}
	else
	{
		fela = 0.;
	}

	/* >>chng 02 jan 14, add grain fe to this sum */
	/* total intensity of K-alpha line */
	/*linadd((fe.fekcld+fe.fegrain)*1.03e-8+(fe.fekhot+fela)*1.11e-8,2,"FeKa",'i' );*/
	if( dense.lgElmtOn[ipIRON] )
	{
		lindst((fe.fekcld+fe.fegrain)*1.03e-8+(fe.fekhot+fela)*1.11e-8,1.78f,"FeKa",
			Transitions[ipH_LIKE][ipIRON][ipH2p][ipH1s].ipCont,'i',false,
			   "total intensity of K-alpha line" );
	}

	linadd(fela*1.11e-8,2,"FeLr",'i' ,
		" recombination from fully stripped ion ");

	/* >>chng 03 aug 14, label changed from TotH to AugH to be like rest total hot iron Ka; */
	linadd((fe.fekhot+fela)*1.11e-8,2,"AugH",'i' ,
		"  Auger hot iron, assumes case b for H and He-like ");

	linadd(fe.fekcld*1.03e-8,2,"AugC",'i',
		" Auger production of cold iron, less than or 17 times ionized " );

	linadd(fe.fegrain*1.03e-8,2,"AugG",'i' ,
		" grain production of cold iron ");

	PutLine(&TauLines[ipCo11527],
		"  [Co XI] 5168. A ");

	PutLine(&TauLines[ipNi1_7m],
		"  nickel  [Ni I] 7m ");

	/* nickel*/


	PutLine(&TauLines[ipNi1_11m],
		"  [Ni I] 11m ");

	PutLine(&TauLines[ipxNi1242],
		" [Ni XII] 4230.8 A ");

	/* copper */

	/* zinc */
	 PutLine(&TauLines[ipZn04363],
		 "zinc iv 3.625 microns, cs and A just made up ");

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_k_zn returns\n" );
	}
	return;
}
