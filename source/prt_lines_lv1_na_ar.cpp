/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_lv1_na_ar place lines of elements sodium through argon into lines storage stack */
#include "cddefines.h"
#include "coolheavy.h"
#include "sil.h"
#include "phycon.h"
#include "embesq.h"
#include "taulines.h"
#include "dense.h"
#include "ionbal.h"
#include "trace.h"
#include "lines_service.h"
#include "lines.h"

void lines_lv1_na_ar(void)
{
	long int ipnt;
	double drec, 
	  fac, 
	  rec, 
	  sum, 
	  t4;

	DEBUG_ENTRY( "lines_lv1_na_ar()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_na_ar called\n" );
	}

	t4 = phycon.te/1e4;

	PutLine(&TauLines[ipT5895],
		" sodium sum of Na D lines");

	PutLine(&TauLines[ipfsNa373],
		" [NaIII] 7.319 micron ");

	PutLine(&TauLines[ipfsNa490],
		" [NaIV] 9.048 micron ");

	PutLine(&TauLines[ipfsNa421],
		" [NaIV] 21.29 micron ");

	PntForLine(1365.,"Na 5",&ipnt);
	lindst(CoolHeavy.c1365,1365,"Na 5",ipnt,'c',true ,
		" [NaV] 1365, sum of 1365.1+1365.8; cs only guess ");

	PntForLine(2067.,"Na 5",&ipnt);
	lindst(CoolHeavy.c2067,2067,"Na 5",ipnt,'c',true ,
		" [NaV] 2067, sum of 2066.9+2068.4; cs only guess ");

	PntForLine(4017.,"Na 5",&ipnt);
	lindst(CoolHeavy.c4017,4017,"Na 5",ipnt,'c',true,
		" [NaV] 4017, sum of 4010.9+4016.7+4022.7; cs only guess " );

	PntForLine(2569.,"Na 6",&ipnt);
	lindst(CoolHeavy.c2569,2569,"Na 6",ipnt,'c',true,
		" [Na VI] 2568.9 " );

	PntForLine(1357.,"Na 6",&ipnt);
	lindst(CoolHeavy.c1357,1357,"Na 6",ipnt,'c',true ,
		" [Na VI] 1356.6 ");

	PntForLine(2972.,"Na 6",&ipnt);
	lindst(CoolHeavy.c2972/(1.+1./3.02),2972,"Na 6",ipnt,'c',true ,
		" [Na VI] 2971.9 ");

	PntForLine(2872.,"Na 6",&ipnt);
	lindst(CoolHeavy.c2972/(1.+3.02),2872,"Na 6",ipnt,'c',true,
		" [Na VI] 2872.7 ");

	PutLine(&TauLines[ipxNa6143],
		" [NaVI] 14.32 micron ");

	PutLine(&TauLines[ipxNa6862],
		" [NaVI] 8.62 micron ");

	PutLine(&TauLines[ipxNa0746],
		" [NaVII] 4.675 micron ");

	PutLine(&TauLines[ipT4561],
		"Magnesium I 4571, O I data for coll strength and trans prob ");

	PutLine(&TauLines[ipMgI2853],
		" Mg I 2853 ");

	PutLine(&TauLines[ipMgI2026],
		" Mg I 2026 ");

	linadd(TauLines[ipT2796].Emis->xIntensity+TauLines[ipT2804].Emis->xIntensity,2798,
		"TOTL",'i',"Mg II 2798 add both lines of multiplet together " );

	/* sum of inward fracs of lines */
	sum =  TauLines[ipT2796].Emis->xIntensity*TauLines[ipT2796].Emis->FracInwd + 
		TauLines[ipT2804].Emis->xIntensity*TauLines[ipT2804].Emis->FracInwd;

	linadd(sum,2798,"Inwd",'i',
		"inward part of Mg II 2798" );

	PutLine(&TauLines[ipT2796],
		"one member of Mg II multiplet");

	PutLine(&TauLines[ipT2804],
		"one member of Mg II multiplet");

	PutLine(&TauLines[ipTMg4],
		" Mg IV 4.5 micron ");

	PutLine(&TauLines[ipTMg14],
		" Mg V  13.5 micron emission");

	PutLine(&TauLines[ipTMg6],
		" Mg V 5.6 micron emission ");

	PutLine(&TauLines[ipxMg52855],
		" [Mg 5] 2571, 2893 ");

	PutLine(&TauLines[ipxMg52417],
		" [Mg 5] 2417.5A, 3-2 in model atom" );

	PutLine(&TauLines[ipxMg51325],
		" [Mg 5] 1324.58A, 3-1 in model atom ");

	PntForLine(1806.,"Mg 6",&ipnt);
	lindst(CoolHeavy.c1806,1806,"Mg 6",ipnt,'c',true ,
		" MG VI");

	PutLine(&TauLines[ipxMg72569],
		" [Mg 7] 2510, 2629, 2-1 transitions, together");

	PutLine(&TauLines[ipxMg72261],
		" [Mg 7] 3-2 transition, 2261.5 ");

	PutLine(&TauLines[ipxMg71190],
		" [Mg 7] 3-1 transition, 1189.82A ");

	PutLine(&TauLines[ipfsMg755] ,
		" Mg 7 IR line 5.50 microns ");

	PutLine(&TauLines[ipfsMg790],
		"  Mg 7 IR line 9.03 microns ");

	PutLine(&TauLines[ipxMg08303],
		" [Mg 8] 3.03 micron ");

	PutLine(&TauLines[ipT705]," Mg 9 704.5 ");

	linadd(TauLines[ipTMg610].Emis->xIntensity+TauLines[ipTMg625].Emis->xIntensity,615,"TOTL",'i',
		" Mg 10 614.9 both of doublet, li seq 2s 2p" );

	PutLine(&TauLines[ipTMg610],
		"");
	PutLine(&TauLines[ipTMg625],
		"");

	PutLine(&TauLines[ipT58],
		" part of Mg 10 destroyed by background opacity Mg 10 58 li seq 2s 3p ");
	linadd(TauLines[ipTMg610].Emis->ots*TauLines[ipTMg610].EnergyErg+
		TauLines[ipTMg625].Emis->ots*TauLines[ipTMg625].EnergyErg,615,"dest",'i',
		"" );

	PutLine(&TauLines[ipAlI3957],
		" Aluminum Al I 3957 ");

	PutLine(&TauLines[ipAlI3090],
		" Al I 3090 ");

	linadd(embesq.em2669+TauLines[ipT2670].Emis->xIntensity,2665,"totl",'i',
		"Al II 1671 total emission in Al II] 2669.7, 2660 doublet" );
	PutLine(&TauLines[ipT2670],
		" ");

	linadd(embesq.em2669,2660,"Al 2",'i',
		"emission in Al II] 2669 alone" );

	linadd(TauLines[ipT1855].Emis->xIntensity+TauLines[ipT1863].Emis->xIntensity,1860,"TOTL",'i',
		" Al III" );
	sum = TauLines[ipT1855].Emis->xIntensity*TauLines[ipT1855].Emis->FracInwd + 
		TauLines[ipT1863].Emis->xIntensity* TauLines[ipT1863].Emis->FracInwd;

	linadd(sum,1860,"Inwd",'i',
		" inward part of AlIII line" );
	PutLine(&TauLines[ipT1855],
		"");
	PutLine(&TauLines[ipT1863],
		"");

	PutLine(&TauLines[ipAl529],
		" [Al V] 2.905 micron ");

	PutLine(&TauLines[ipAl6366],
		" [Al VI] 3.66 micron ");

	PutLine(&TauLines[ipAl6912],
		" [Al VI] 9.12 micron");

	PntForLine(2428.,"Al 6",&ipnt);
	lindst(CoolHeavy.c2428/(1.+1./3.73),2428,"Al 6",ipnt,'c',true,
		" [Al VI] 2428.4 " );

	PntForLine(2601.,"Al 6",&ipnt);
	lindst(CoolHeavy.c2428/(1.+3.73),2601,"Al 6",ipnt,'c',true,
		"  [Al VI] 2601.0");

	PntForLine(1170.,"Al 6",&ipnt);
	lindst(CoolHeavy.c1170,1170,"Al 6",ipnt,'c',true,
		" [Al VI] 1169.86 ");

	PntForLine(2125.,"Al 6",&ipnt);
	lindst(CoolHeavy.c2125,2125,"Al 6",ipnt,'c',true ,
		"  [Al VI] 2124.95" );

	PutLine(&TauLines[ipAl8575],
		" [Al VIII] 5.75 micron " );

	PutLine(&TauLines[ipAl8370],
		" [Al VIII] 3.70 micron");

	PutLine(&TauLines[ipAl09204],
		"  [Al IX] 2.04 micron, no collision strength, A NIST ");

	PutLine(&TauLines[ipT639]," Al X ");

	linadd(TauLines[ipTAl550].Emis->xIntensity+TauLines[ipTAl568].Emis->xIntensity,556,"TOTL",'i',
		" Al 11, Li seq 2s2p" );

	PutLine(&TauLines[ipTAl550],
		"");

	PutLine(&TauLines[ipTAl568],
		"");

	PutLine(&TauLines[ipTAl48],
		"  Al 11, Li seq 2s3p ");

	PutLine(&TauLines[ipSi1_130m],
		" silicon Silicon  Si I 130m ");

	 PutLine(&TauLines[ipSi1_68m],
		 " Si I 68m");

	PutLine(&TauLines[ipSii2518],
		"  Si I 2518A ");

	 PutLine(&TauLines[ipSii2215],
		 " Si I 2215A ");

	PutLine(&TauLines[ipTSi35],
		" Silicon II 35 micron ");

	linadd(
		TauLines[ipSi2_2334].Emis->xIntensity+
		TauLines[ipSi2_2329].Emis->xIntensity+
		TauLines[ipSi2_2350].Emis->xIntensity+
		TauLines[ipSi2_2344].Emis->xIntensity+
		TauLines[ipSi2_2336].Emis->xIntensity,
		2335,"TOTL",'i',
		"total intensity of S IV] 1406, all lines in the multiplet" );
	PutLine(&TauLines[ipSi2_2334],
		" ");
	PutLine(&TauLines[ipSi2_2329],
		" ");
	PutLine(&TauLines[ipSi2_2350],
		" ");
	PutLine(&TauLines[ipSi2_2344],
		" ");
	PutLine(&TauLines[ipSi2_2336],
		" ");

	PutLine(&TauLines[ipT1808],
		" SI II 1808, permitted resonance line, collisionally excited ");

	PutLine(&TauLines[ipT1527],
		"  SI II 1527, permitted resonance line, collisionally excited ");

 	PutLine(&TauLines[ipT1305],
		" SI II 1305, permitted resonance line, collisionally excited ");

	PutLine(&TauLines[ipT1260],
		"  SI II 1260, permitted resonance line, collisionally excited ");

	/* SI II 1260, rough guess of dielec contribution */
	drec = dense.xIonDense[13][2]*dense.eden*7.6e-7/phycon.te32*
	  1.57e-11*ionbal.DielSupprs[1][2];

	fac = emit_frac(&TauLines[ipT1260]);
	PntForLine(1260.,"Si 2",&ipnt);
	lindst(drec*fac,1260,"diel",ipnt,'r',true,
		" fac = (1.-TauLines[ipT1260].ColOvTot);" );


	PntForLine(1909.,"Si 2",&ipnt);
	lindst(drec*1260./1909.,1909,"diel",ipnt,'r',true,
		" dielectronic recombination SiII 1909" );

	rec = 1e-12*(0.1152/t4 - 0.3082 + 4.4734*t4 + 0.0207*t4*t4)/pow(t4,1.5)*
	  sexp(0.2981/t4);

	/* >>chng 96 jul 8, added dielectronic recombination contribution to 1207 */
	/*rec *= dense.xIonDense[13][3]*dense.eden*(1.-TauLines[ipT1207].ColOvTot)**/
	rec *= dense.xIonDense[13][3]*dense.eden*emit_frac(&TauLines[ipT1207])*
	  1.65e-11*ionbal.DielSupprs[1][2];

	PutExtra(MAX2(0.,rec));

	PutLine(&TauLines[ipT1207],
		" SI III 1207, collisional excitation and dielectronic recombination ");

	linadd(MAX2(0.,rec),1207,"rec ",'i',
		" Si III 1207, dielectronic recombination only" );

	linadd(embesq.em1895+TauLines[ipT1895].Emis->xIntensity,1888,"TOTL",'i',
		" Si III] 1892+1883, total intensity of both lines" );
	PutLine(&TauLines[ipT1895],
		" ");

	PntForLine(1883.,"Si 3",&ipnt);
	lindst(embesq.em1895,1883,"Si 3",ipnt,'c',true ,
		" Si III] 1883 by itself");

	/*fac = (1.-TauLines[ipT1895].ColOvTot);*/
	fac = emit_frac(&TauLines[ipT1895]);

	double p1895 = ionbal.PhotoRate_Shell[ipSILICON][1][2][0]*
		dense.xIonDense[ipSILICON][1]*0.85;
	linadd( p1895*1.05e-11*fac,1895,"PHOT",'i',
		" photoproduction by inner shell removal" );

	linadd(TauLines[ipT1403].Emis->xIntensity+TauLines[ipT1394].Emis->xIntensity,1397,"TOTL",'i',
		" Si IV 1397, collisionally excited " );

	sum = TauLines[ipT1403].Emis->xIntensity*TauLines[ipT1403].Emis->FracInwd + 
		TauLines[ipT1394].Emis->xIntensity* TauLines[ipT1394].Emis->FracInwd;

	linadd(sum,1397,"Inwd",'i',
		" inward part of SiIV 1397" );
	PutLine(&TauLines[ipT1403],
		" ");
	PutLine(&TauLines[ipT1394],
		" ");

	PutLine(&TauLines[ipSi619],
		"  SI VI 1.9641 micron ");

	PntForLine(2148.,"Si 7",&ipnt);
	lindst(sil.c2148,2148,"Si 7",ipnt,'c',true,
		" SI VII, 2148, O III like, collisionally excited" );

	PutLine(&TauLines[ipTSi25],
		"  Si VII 2.48, 6.49 micron, collisionally excited ");

	PutLine(&TauLines[ipTSi65],
		"  Si VII 2.48, 6.49 micron, collisionally excited ");

	PntForLine(1446.,"Si 8",&ipnt);
	lindst(sil.c1446,1446,"Si 8",ipnt,'c',true,
		" SI VIII 1446, OIII like, collisionally excited" );

	PntForLine(1985.,"Si 9",&ipnt);
	lindst(sil.c1985,1985,"Si 9",ipnt,'c',true,
		" SI IX 1985, 2150, collisionally excited" );


	PntForLine(949.,"Si 9",&ipnt);
	lindst(sil.c949,949,"Si 9",ipnt,'c',true,
		" collisionally excited" );

	PntForLine(1815.,"Si 9",&ipnt);
	lindst(sil.c1815,1815,"Si 9",ipnt,'c',true,
		" collisionally excited " );

	PutLine(&TauLines[ipTSi4],
		"  SI 9, 3.86, 2.84 3P fine structure lines ");

	PutLine(&TauLines[ipTSi3],
		"  SI 9, 3.86, 2.84 3P fine structure lines ");

	PntForLine(691.,"Si 9",&ipnt);
	lindst(sil.c691,691,"Si 9",ipnt,'c',true ,
		"  both components of 5S-3P doublet");

	PutLine(&TauLines[ipSi10_606],
		" SI 10 606A, actually group of 4 intercombination lines ");

	 PutLine(&TauLines[ipSi10143],
		 " [Si 10] 1.43 micron, collisionally excited ");

	PntForLine(583.,"Si11",&ipnt);
	lindst(sil.c583,581,"Si11",ipnt,'c',true,
		"  Si 11 582.9, collisionally excited  >>chng 01 may 23, wavelength from 583 to 581" );

	linadd(TauLines[ipTSi499].Emis->xIntensity+TauLines[ipTSi521].Emis->xIntensity,506,"TOTL",'i' ,
		"emission total Si 12 506 + 499 ");

	PutLine(&TauLines[ipTSi499],
		"  Si 12 506 li seq 2s 2p ");
	PutLine(&TauLines[ipTSi521],
		" ");

	PutLine(&TauLines[ipTSi41],
		"  Si 12 40.9A, li seq 2s 3p ");

	PutLine(&TauLines[ipP0260],
		" phosphorus  [P II] 60.64 micron ");

	PutLine(&TauLines[ipP0233],
		"  [P II] 32.87 micron ");

	PntForLine(16400.,"P  2",&ipnt);
	lindst(CoolHeavy.p2_32,16400,"P  2",ipnt,'c',true ,
		" 3-2 1.64 micron ");

	PntForLine(4669.,"P  2",&ipnt);
	lindst(CoolHeavy.p2_31*0.75,4669,"P  2",ipnt,'c',true ,
		"  >>chng 01 may 15, add these lines  P 2 3-1 4670, 4738 vac wl, 4669, 4737 air");

	PntForLine(4737.,"P  2",&ipnt);
	lindst(CoolHeavy.p2_31*0.25,4737,"P  2",ipnt,'c',true ,
		" P 2 3-1 4670, 4738 vac wl, 4737, 4737 air");

	PntForLine(11890.,"P  2",&ipnt);
	lindst(CoolHeavy.p2_21*0.75,11890,"P  2",ipnt,'c',true,
		" 2-1 1.147, 1.189 micron" );

	PntForLine(11470.,"P  2",&ipnt);
	lindst(CoolHeavy.p2_21*0.25,11470,"P  2",ipnt,'c',true ,
		"  [P II] 1.14 micron");

	PutLine(&TauLines[ipP0318],
		"  [P III] 17.885 micron ");

	PutLine(&TauLines[ipP713],
		"  [P VII] 1.3745 micron ");

	PutLine(&TauLines[ipP848],
		" [P VIII] 4.85 micron ");

	PutLine(&TauLines[ipP817],
		"  [P VIII] 1.735 microns ");

	PutLine(&TauLines[ipP1027],
		"  [P X] 2.708 micron ");

	PutLine(&TauLines[ipP1018],
		"  [P X] 1.868 microns ");

	PutLine(&TauLines[ipS1_25m],
		"  sulphur  S I 25m ");

	PutLine(&TauLines[ipS1_56m],
		"  S I 56m ");

	linadd(1.75e-22*dense.eden*dense.xIonDense[ipSULPHUR][1]/phycon.sqrte/
	  phycon.te10/phycon.te03,1807,"S 1R",'i',
	  " guesstimate of Sulphur I triplet excited state recombination rate.this is to check whether photoexcit of S II is ever important S I 1807 recombination " );

	linadd(CoolHeavy.c6731,6720,"S  2",'i',
		" S II 6731 + 6716 together " );

	linadd(CoolHeavy.S4070+CoolHeavy.S4078,4074,
		"S  2",'i',"S II 4070 +4078 together" );

	PntForLine(10330.,"S  2",&ipnt);
	lindst(CoolHeavy.c10330,10330,"S  2",ipnt,'c',true,
		" S II N=3 lines, all four lines together " );

	PntForLine(6731.,"S II",&ipnt);
	lindst(CoolHeavy.S6733,6731,"S II",ipnt,'c',true,
		" individual line from five level atom" );

	PntForLine(6716.,"S II",&ipnt);
	lindst(CoolHeavy.S6718,6716,"S II",ipnt,'c',true,
		" individual line from five level atom" );

	PntForLine(4070.,"S II",&ipnt);
	lindst(CoolHeavy.S4070,4070,"S II",ipnt,'c',true,
		" individual line from five level atom " );

	PntForLine(4078.,"S II",&ipnt);
	lindst(CoolHeavy.S4078,4078,"S II",ipnt,'c',true,
		" individual line from five level atom" );

	linadd(CoolHeavy.S10323,10323,"S II",'c',
		" individual line from five level atom " );

	linadd(CoolHeavy.S10289,10289,"S II",'c',
		" individual line from five level atom" );

	linadd(CoolHeavy.S10373,10373,"S II",'c' ,
		" individual line from five level atom ");

	linadd(CoolHeavy.S10339,10339,"S II",'c',
		" individual line from five level atom " );

	 PutLine(&TauLines[ipT1256],
		 " resonance line near NV, collisionally excited ");

	PutLine(&TauLines[ipTS19],
		" S III fine structure 18.7 ");

	PutLine(&TauLines[ipTS34],
		"  S III fine structure 34 ");

	PutLine(&TauLines[ipTS1720],
		"  S III] 1713.12, 1728.94 ");

	PntForLine(9532.,"S  3",&ipnt);
	lindst(CoolHeavy.c9532/(1.+1./2.48),9532,"S  3",ipnt,'c',true ,
		" [S III] 9532 alone ");

	PntForLine(9069.,"S  3",&ipnt);
	lindst(CoolHeavy.c9532/(1.+2.48),9069,"S  3",ipnt,'c',true,
		" [S III] 9069 alone" );

	PntForLine(6312.,"S  3",&ipnt);
	lindst(CoolHeavy.c6312,6312,"S  3",ipnt,'c',true ,
		" [S III] 6312, trans-auroral temperature sensitive ");

	PntForLine(3722.,"S  3",&ipnt);
	lindst(CoolHeavy.c6312*0.59,3722,"S  3",ipnt,'c',true,
		" [S III] 3722, same upper level as 6312" );

	PutLine(&TauLines[ipT1194],
		"  WL, other data, from Ho + Henry Ap.J. 1984 ");

	PutLine(&TauLines[ipTS11],
		"  S IV 10.5 micron, collisionally excited (label is 105) ");

	linadd(
		TauLines[ipS4_1405].Emis->xIntensity+
		TauLines[ipS4_1398].Emis->xIntensity+
		TauLines[ipS4_1424].Emis->xIntensity+
		TauLines[ipS4_1417].Emis->xIntensity+
		TauLines[ipS4_1407].Emis->xIntensity,
		1406,"TOTL",'i',
		" total intensity of S IV] 1406, all lines in the multiplet " );
	PutLine(&TauLines[ipS4_1405],
		" ");
	PutLine(&TauLines[ipS4_1398],
		" ");
	PutLine(&TauLines[ipS4_1424],
		" ");
	PutLine(&TauLines[ipS4_1417],
		" ");
	PutLine(&TauLines[ipS4_1407],
		" ");

		linadd(embesq.em1198+TauLines[ipT1198].Emis->xIntensity,1198,"TOTL",'i',
			" S V 1198] both lines together " );

	PutLine(&TauLines[ipT1198],
		"  S V 1198] the stronger transition ");

	linadd(embesq.em1198,1188,"S  5",'i',
		" Be seq, weaker of the two transitions" );

	PutLine(&TauLines[ipT786],
		"  S V 786.5, collisionally excited ");

	linadd(TauLines[ipT944].Emis->xIntensity+TauLines[ipT933].Emis->xIntensity,933,"TOTL",'i',
		" S VI 933+944, collisionally excited total S VI 933+944 " );
	PutLine(&TauLines[ipT944],
		" ");
	PutLine(&TauLines[ipT933],
		" ");

	PutLine(&TauLines[ipfsS810],
		" S VIII 9913 ");

	PutLine(&TauLines[ipfsS912],
		"  S IX 1.25 microns ");

	PutLine(&TauLines[ipfsS938],
		"  S IX  3.78 microns ");

	PutLine(&TauLines[ipfsS1119],
		"  S XI 1.91 microns ");

	PutLine(&TauLines[ipfsS1114],
		"  S XI 1.40 microns ");

	PntForLine(1715.,"S  9",&ipnt);
	lindst(CoolHeavy.c1715,1715,"S  9",ipnt,'c',true ,
		" S IX 1715, 1987, collisionally excited");

	PntForLine(1213.,"S 10",&ipnt);
	lindst(CoolHeavy.c1213,1213,"S 10",ipnt,'c',true,
		" S X 1213, 1197, collisionally excited " );

	PntForLine(1826.,"S 11",&ipnt);
	lindst(CoolHeavy.c1826,1826,"S 11",ipnt,'c',true,
		"  S XI 1615, 1826, collisionally excited" );

	PutLine(&TauLines[ipfsS1207],
		"  S XII 7611.2 ");

	/* group of four intercombination lines
	PntForLine(520.,"S 12",&ipnt);
	lindst(CoolHeavy.c520,520,"S 12",ipnt,'c',true ); */
	PutLine(&TauLines[ipS12_520],
		"  >>chng 03 sep 27, move to atom_level2 from expion");

	PntForLine(488.,"S 13",&ipnt);
	lindst(CoolHeavy.c488,488,"S 13",ipnt,'c',true,
		" S XIII 488.4, 1909 like, collisionally excited" );

	linadd(TauLines[ipTSu418].Emis->xIntensity+TauLines[ipTSu446].Emis->xIntensity,427,"TOTL",'i',
		" S 14 506 li seq 2s2p " );
	PutLine(&TauLines[ipTSu418],
		" ");
	PutLine(&TauLines[ipTSu446],
		" ");

	PutLine(&TauLines[ipTSu30],
		"  S 14 30 li seq 2s3p ");

	PutLine(&TauLines[ipCl1_11m],
		"  chlorine lines  [Cl I] 11 micron ");

	PutLine(&TauLines[ipfsCl233],
		"  [Cl II] 33.281 micron ");

	PutLine(&TauLines[ipfsCl214],
		"  [Cl II] 14.3678 micron ");

	PntForLine(8578.7,"Cl 2",&ipnt);
	lindst(CoolHeavy.c8579*0.791,8579,"Cl 2",ipnt,'c',true,
		" Chlorine II 8578.7, 9123.6 doublet");

	PntForLine(9123.6,"Cl 2",&ipnt);
	lindst(CoolHeavy.c8579*0.209,9124,"Cl 2",ipnt,'c',true,
		"  Chlorine II 8578.7, 9123.6 doublet" );

	PntForLine(6161.8,"Cl 2",&ipnt);
	lindst(CoolHeavy.c6164,6162,"Cl 2",ipnt,'c',true,
		"  Chlorine II 6161.8 auroral line  >>chng 03 feb 24, change wavelength from 6164 to correct 6161.8 " );

	PntForLine(3677.9,"Cl 2",&ipnt);
	lindst(CoolHeavy.c3679,3678,"Cl 2",ipnt,'c',true,
		"  Chlorine II 3677.9 auroral line >>chng 03 feb 24, to correct wavelength " );

	linadd(CoolHeavy.c5525,5525,"TOTL",'i',
		" Cl III 5519, 5539 doublet, both together " );

	linadd(CoolHeavy.c3350,3350,"TOTL",'i',
		" Cl III 3354, 3344 doublet, both together " );

	linadd(CoolHeavy.c8494,8494,"TOTL",'i',
		" Cl III 8504, 8436, 8552, 8483 multiplet, all together " );

	PntForLine(5538.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl5539,5538,"Cl 3",ipnt,'c',true,
		" Cl III 5538  " );

	PntForLine(5518.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl5519,5518,"Cl 3",ipnt,'c',true,
		" Cl III 5518" );

	PntForLine(3354.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl3354,3354,"Cl 3",ipnt,'c',true,
		" Cl III 3354  " );

	PntForLine(3344.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl3344,3344,"Cl 3",ipnt,'c',true,
		" Cl III 3344 " );

	PntForLine(8504.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl8504,8504,"Cl 3",ipnt,'c',true ,
		" Cl III 8504  ");

	PntForLine(8436.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl8436,8436,"Cl 3",ipnt,'c',true,
		" Cl III 8436" );

	PntForLine(8552.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl8552,8552,"Cl 3",ipnt,'c',true,
		" Cl III 8552 " );

	PntForLine(8483.,"Cl 3",&ipnt);
	lindst(CoolHeavy.Cl8483,8483,"Cl 3",ipnt,'c',true,
		" Cl III 8483" );

	PutLine(&TauLines[ipCl04203],
		" [Cl IV] fine structure line 20.354 microns");

	PutLine(&TauLines[ipCl04117],
		" [Cl IV] fine structure line 11.741 microns ");

	PntForLine(8047.,"Cl 4",&ipnt);
	lindst(CoolHeavy.c8047*0.667,8047,"Cl 4",ipnt,'c',true,
		"  ClIV 8047" );

	PntForLine(7532.,"Cl 4",&ipnt);
	lindst(CoolHeavy.c8047*0.333,7532,"Cl 4",ipnt,'c',true,
		"  ClIV 7532" );

	PntForLine(3119.,"Cl 4",&ipnt);
	lindst(CoolHeavy.c3119,3119,"Cl 4",ipnt,'c',true,
		" ClIV 3119" );

	PntForLine(5324.,"Cl 4",&ipnt);
	lindst(CoolHeavy.c5324,5324,"Cl 4",ipnt,'c',true,
		" ClIV 5324" );

	PutLine(&TauLines[ipCl973],
		"  Cl IX 7334A ");

	PutLine(&TauLines[ipCl1030],
		"  Cl X 3.05 microns ");

	PutLine(&TauLines[ipCl1092],
		"  Cl X 9223A ");

	PutLine(&TauLines[ipTAr7],
		" Argon II 7 micron ");

	PntForLine(7135.,"Ar 3",&ipnt);
	lindst(CoolHeavy.c7136/(1.+1./4.144),7135,"Ar 3",ipnt,'c',true ,
		" Argon III 7135");

	PntForLine(7751.,"Ar 3",&ipnt);
	lindst(CoolHeavy.c7136/(1.+4.144),7751,"Ar 3",ipnt,'c',true,
		" Argon III 7751" );

	PntForLine(5192.,"Ar 3",&ipnt);
	lindst(CoolHeavy.c5192,5192,"Ar 3",ipnt,'c',true,
		" Argon III 5192" );

	PntForLine(3109.,"Ar 3",&ipnt);
	lindst(CoolHeavy.c3109*0.9894,3109,"Ar 3",ipnt,'c',true,
		"  Argon III 3109" );

	PntForLine(3005.,"Ar 3",&ipnt);
	lindst(CoolHeavy.c3109*(1.-0.9894),3005,"Ar 3",ipnt,'c',true ,
		" Argon III 3005 ");

	PutLine(&TauLines[ipTAr22],
		"  Argon III 21.8, 9 micron lines");

	PutLine(&TauLines[ipTAr9],
		" Argon III 21.8, 9 micron lines ");

	linadd(CoolHeavy.Ar4740+CoolHeavy.Ar4711,4725,"TOTL",'i',
		" Argon IV 4711 + 4740 together, 4740=90%" );

	linadd(CoolHeavy.Ar2868+CoolHeavy.Ar2854,2860,"TOTL",'i',
		" [AvIV] 2868, 2854 together " );

	linadd(CoolHeavy.Ar7237+CoolHeavy.Ar7331+CoolHeavy.Ar7171+CoolHeavy.Ar7263,7250,"TOTL",'i',
	  " [AvIV] auroral lines, 7237, 7331, 7171, 7263 " );

	PntForLine(4740.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar4740,4740,"Ar 4",ipnt,'c',true,
		" [Ar IV] 4740" );

	PntForLine(4711.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar4711,4711,"Ar 4",ipnt,'c',true,
		"  [Ar IV] 4711" );

	PntForLine(2868.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar2868,2868,"Ar 4",ipnt,'c',true,
		" [Ar IV] 2868" );

	PntForLine(2854.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar2854,2854,"Ar 4",ipnt,'c',true,
		"  [Ar IV] 2854" );


	PntForLine(7263.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar7263,7263,"Ar 4",ipnt,'c',true,
		"  [Ar IV] 7263" );

	PntForLine(7171.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar7171,7171,"Ar 4",ipnt,'c',true,
		"  [Ar IV] 7171" );

	PntForLine(7331.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar7331,7331,"Ar 4",ipnt,'c',true,
		"  [Ar IV] 7331" );

	PntForLine(7237.,"Ar 4",&ipnt);
	lindst(CoolHeavy.Ar7237,7237,"Ar 4",ipnt,'c',true,
		" [Ar IV] 7237" );

	PntForLine(7005.,"Ar 5",&ipnt);
	lindst(CoolHeavy.c7007/(1.+1./2.143),7005,"Ar 5",ipnt,'c',true,
		" Argon V, 3P lines, 7005, collisionally excited" );

	PntForLine(6435.,"Ar 5",&ipnt);
	lindst(CoolHeavy.c7007/(1.+2.143),6435,"Ar 5",ipnt,'c',true ,
		" Argon V, 3P lines, 6435, collisionally excited");

	PntForLine(4626.,"Ar 5",&ipnt);
	lindst( CoolHeavy.c4626 , 4626 ,"Ar 5",ipnt,'c',true,
		" >>chng 01 mar 10, add following two lines  Argon V, 4626" );

	PntForLine(2691.,"Ar 5",&ipnt);
	lindst( CoolHeavy.c2691 , 2691 ,"Ar 5",ipnt,'c',true ,
		" Argon V, 2691");

	PutLine(&TauLines[ipTAr13],
		"  Argon V fine structure lines, 13.09, 7.903 micron line ");

	PutLine(&TauLines[ipTAr8],
		"  Argon V fine structure lines, 13.09, 7.903 micron line ");

	PutLine(&TauLines[ipAr06453],
		"  [Ar VI] 4.53 micron ");

	PutLine(&TauLines[ipAr1055],
		"  Ar X 5533.4A ");

	PutLine(&TauLines[ipAr1126],
		" Ar XI 2.60 microns ");

	PutLine(&TauLines[ipAr1178],
		" Ar XI 6917A ");

	PntForLine(4413.,"Ar14",&ipnt);
	lindst(CoolHeavy.fs4413,4413,"Ar14",ipnt,'c',true,
		" Ar XIV 4413, predicted lambda, not observed(?)" );

	PntForLine(409.,"Ar15",&ipnt);
	lindst(CoolHeavy.c409,409,"Ar15",ipnt,'c',true,
		" collisionally excited " );

	PutLine(&TauLines[ipT354],
		" collisionally excited Ar16 354, Li seq 2s-2p");

	PutLine(&TauLines[ipT389],
		"  collisionally excited Ar16 389, Li seq 2s-2p ");

	PutLine(&TauLines[ipT25],
		" collisionally excited Ar16 365, Li seq 2s-3p ");

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_na_ar returns\n" );
	}
	return;
}
