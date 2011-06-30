/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_setup convert level 1 and level 2 line parameters and pointers into internal form 
 * used by code, line data were read in by atmdat_readin */
#include "cddefines.h"
#include "physconst.h"
#include "lines_service.h"
#include "prt.h"
#include "taulines.h"
#include "opacity.h"
#include "lines.h"

/* following used for generating array indices to level1 lines,
 * start process, use it, then end */
STATIC void initFindLevLine( void );
STATIC long ipFindLevLine( realnum , long , long );
STATIC void endFindLevLine( void );

void lines_setup(void)
{
	long int i;

	static bool lgFirst = true;
	bool lgSane;

	DEBUG_ENTRY( "lines_setup()" );

	/* this routine takes the line parameters in the wind block data and sorts
	 * them into what is needed for the actual line optical depth arrays */

	/** \todo	2	initialization already done at this point */

	/* this is the dummy line */
	TauLines[0].WLAng = 0.;
	TauLines[0].Lo->g = 0;
	TauLines[0].Hi->g = 0.;
	TauLines[0].Emis->gf = 0.;
	TauLines[0].EnergyWN = 0.;
	TauLines[0].Hi->IonStg = 0;
	TauLines[0].Hi->nelem = 0;
	/* this is an impossible value of iRedisFun */
	TauLines[0].Emis->iRedisFun = 0;
	TauLines[0].Emis->Aul = 0.;

	/* the first valid line is [0] since zero is the dummy */
	if( TauLines[1].EnergyWN <= 0. )
	{
		fprintf( ioQQQ, " PROBLEM Insane value for TauLines array.\n" );
		fprintf( ioQQQ, " Was block data LineData linked in??\n" );
		fprintf( ioQQQ, " Check that it compiled OK (it probably did not).\n" );
		TotalInsanity();
	}

	/* check that all lines have valid data */
	lgSane = true;
	for( i=1; i <= nLevel1; i++ )
	{

		if( TauLines[i].Lo->g <= 0. )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane lower stat wght\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		if( TauLines[i].Hi->g <= 0. )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane upper stat wght\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		if( TauLines[i].EnergyWN <= 0. )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane energy WN\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		if( TauLines[i].Hi->IonStg <= 0 )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane ioniz stage\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		if( TauLines[i].Hi->nelem <= 0 || TauLines[i].Hi->nelem > (int)LIMELM )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane Nelem\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		if( TauLines[i].Hi->IonStg > TauLines[i].Hi->nelem )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane IonStg>Nelem\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		if( TauLines[i].Emis->iRedisFun == 0 )
		{
			fprintf( ioQQQ, "  routine lines_setup, insane line redis fcn\n" );
			fprintf( ioQQQ, " line index is %5ld\n", i );
			lgSane = false;
		}

		/* use energies for wavelengths in air if wl not forced with wl number on line */
		/* >>chng 03 oct 07, only make correction for index ref if
		 * if wl was not already set - this is an option to allow
		 * the printed wl to be specified in the level1.dat file */
		if( TauLines[i].WLAng <= 0. )
		{
			/* make following an air wavelength */
			TauLines[i].WLAng = 
				(realnum)(1.0e8/
			TauLines[i].EnergyWN/
			RefIndex( TauLines[i].EnergyWN));
		}
		{
			/*@-redef@*/
			enum{DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC  )
			{
				char chString[10];
				chIonLbl(chString,&TauLines[i]);
				fprintf( ioQQQ,"%s ", chString );
				prt_wl( ioQQQ , TauLines[i].WLAng );
				fprintf(ioQQQ,"\n");
			}
		}
	}

	if( !lgSane )
	{
		fprintf( ioQQQ, " Insane value for line arrays encountered.\n" );
		fprintf( ioQQQ, " Was block data lines linked in??\n" );
		fprintf( ioQQQ, " Were errors intreoducted into the line array?\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* set up array to store hits for each line */
	initFindLevLine( );

	/* in following calls to ipFindLevLine the numbers are the integer wavelength
	 * used in the printout, the ion stage, and the atomic number */

	/* carbon line optical depth data */
	ipT1656 = ipFindLevLine( 1656 , 1 , 6 );

	ipT9830 = ipFindLevLine( 9830 , 1 , 6 );

	ipT8727 = ipFindLevLine( 8727 , 1 , 6 );

	ipC2_2325 = ipFindLevLine( 2325 , 2 , 6 );
	ipC2_2324 = ipFindLevLine( 2324 , 2 , 6 );
	ipC2_2329 = ipFindLevLine( 2329 , 2 , 6 );
	ipC2_2328 = ipFindLevLine( 2328 , 2 , 6 );
	ipC2_2327 = ipFindLevLine( 2327 , 2 , 6 );

	ipT1335 = ipFindLevLine( 1335 , 2 , 6 );

	ipT1909 = ipFindLevLine( 1910 , 3 , 6 );

	ipT977 = ipFindLevLine( 977 , 3 , 6 );

	ipT1550 = ipFindLevLine( 1551 , 4 , 6 );

	ipT1548 = ipFindLevLine( 1548 , 4 , 6 );

	ipT386 = ipFindLevLine( 386 , 3 , 6 );

	ipT310 = ipFindLevLine( 310 , 3 , 6 );

	/*CIII* 1175, lower level = upper level of 1909*/
	ipc31175 = ipFindLevLine( 1176 , 3 , 6 );

	ipT291 = ipFindLevLine( 291 , 3 , 6 );

	ipT280 = ipFindLevLine( 280 , 3 , 6 );

	ipT274 = ipFindLevLine( 274 , 3 , 6 );

	ipT270 = ipFindLevLine( 270 , 3 , 6 );

	ipT312 = ipFindLevLine( 312 , 4 , 6 );

	/*carbon fine structure lines added by Jim Kingdon*/
	ipT610 = ipFindLevLine( 6092000 , 1 , 6 );

	ipT370 = ipFindLevLine( 3697000 , 1 , 6 );

	ipT157 = ipFindLevLine( 1576000 , 2 , 6 );

	/*nitrogen line optical depth data*/
	ipNI_pumpIndirect  = ipFindLevLine( 953.9f , 1 , 7 );
	ipNI_pumpDirect[0] = ipFindLevLine( 954.1f , 1 , 7 );
	ipNI_pumpDirect[1] = ipFindLevLine( 951.0f , 1 , 7 );
	ipNI_pumpDirect[2] = ipFindLevLine( 955.8f , 1 , 7 );
	ipNI_pumpDirect[3] = ipFindLevLine( 959.4f , 1 , 7 );
	ipNI_pumpDirect[4] = ipFindLevLine( 955.5f , 1 , 7 );
	ipNI_pumpDirect[5] = ipFindLevLine( 951.2f , 1 , 7 );
	ipNI_pumpDirect[6] = ipFindLevLine( 960.2f , 1 , 7 );
	ipNI_pumpDirect[7] = ipFindLevLine(1159.8f , 1 , 7 );
	ipNI_pumpDirect[8] = ipFindLevLine(1160.9f , 1 , 7 );

	ipT1200 = ipFindLevLine( 1200 , 1 , 7 );

	ipT1085 = ipFindLevLine( 1085 , 2 , 7 );

	ipN3_1749 = ipFindLevLine( 1749 , 3 , 7 );
	ipN3_1747 = ipFindLevLine( 1747 , 3 , 7 );
	ipN3_1754 = ipFindLevLine( 1754 , 3 , 7 );
	ipN3_1752 = ipFindLevLine( 1752 , 3 , 7 );
	ipN3_1751 = ipFindLevLine( 1751 , 3 , 7 );

	ipT990 = ipFindLevLine( 991 , 3 , 7 );

	ipT1486 = ipFindLevLine( 1486 , 4 , 7 );

	ipT765 = ipFindLevLine( 765 , 4 , 7 );

	ipT1243 = ipFindLevLine( 1243 , 5 , 7 );

	ipT1239 = ipFindLevLine( 1239 , 5 , 7 );

	ipT374g = ipFindLevLine( 373 , 3 , 7 );

	/*this is the stronger of the two lines*/
	ipT374x = ipFindLevLine( 374 , 3 , 7 );

	ipT2140 = ipFindLevLine( 2141 , 2 , 7 );

	ipT671 = ipFindLevLine( 671 , 2 , 7 );

	ipT315 = ipFindLevLine( 315 , 3 , 7 );

	ipT324 = ipFindLevLine( 324 , 3 , 7 );

	ipT333 = ipFindLevLine( 333 , 3 , 7 );

	ipT209 = ipFindLevLine( 209 , 5 , 7 );

	/*fine structure lines */
	/*[N II] 121.7*/
	ipT122 = ipFindLevLine( 1217000 , 2 , 7 );

	/*[N II] 205.4*/
	ipT205 = ipFindLevLine( 2054000 , 2 , 7 );

	/*big disagreement in A for this line, other val is 2x larger*/
	/*see review in DEO Seaton 70th birthday*/
	ipT57 = ipFindLevLine( 572100 , 3 , 7 );

	/*oxygen line optical depth data*/
	ipT6300 = ipFindLevLine( 6300 , 1 , 8 );

	ipT6363 = ipFindLevLine( 6363 , 1 , 8 );

	/*A from NISt 96*/
	ipT5577 = ipFindLevLine( 5577 , 1 , 8 );

	ipT834 = ipFindLevLine( 833.8f , 2 , 8 );

	ipT1661 = ipFindLevLine( 1661 , 3 , 8 );

	ipT1666 = ipFindLevLine( 1666 , 3 , 8 );

	ipT835 = ipFindLevLine( 835 , 3 , 8 );

	ipO4_1400 = ipFindLevLine( 1400 , 4 , 8 );
	ipO4_1397 = ipFindLevLine( 1397 , 4 , 8 );
	ipO4_1407 = ipFindLevLine( 1407 , 4 , 8 );
	ipO4_1405 = ipFindLevLine( 1405 , 4 , 8 );
	ipO4_1401 = ipFindLevLine( 1401 , 4 , 8 );

	ipT789 = ipFindLevLine( 789 , 4 , 8 );

	ipT630 = ipFindLevLine( 630 , 5 , 8 );

	/*start OI 6 level atom*/
	ipT1304 = ipFindLevLine( 1304 , 1 , 8 );

	ipT1039 = ipFindLevLine( 1039 , 1 , 8 );

	ipT8446 = ipFindLevLine( 8446 , 1 , 8 );

	ipT4368 = ipFindLevLine( 4368 , 1 , 8 );

	ipTOI13 = ipFindLevLine( 13200 ,  1 , 8 );

	ipTOI11 = ipFindLevLine( 11300 , 1 , 8 );

	ipTOI29 = ipFindLevLine( 29000 , 1 , 8 );

	ipTOI46 = ipFindLevLine( 46000 , 1 , 8 );

	ipTO1025 = ipFindLevLine( 1025 , 1 , 8 );

	/*end of OI 6 level atom*/

	ipT304 = ipFindLevLine( 304 , 3 , 8 );

	ipT1214 = ipFindLevLine( 1218 , 5 , 8 );

	ipT150 = ipFindLevLine( 150 , 6 , 8 );

	/*fine structure lines*/
	/*[O I] 146 microns*/
	ipT146 = ipFindLevLine( 1455300 , 1 , 8 );

	/*[O I] 63 microns*/
	ipT63 = ipFindLevLine( 631700 , 1 , 8 );

	/*[O III] 88.3564 m*/
	ipTO88 = ipFindLevLine( 883300 , 3 , 8 );

	/*[O III] 51.8145*/
	ipT52 = ipFindLevLine( 518000 , 3 , 8 );

	/*[O IV] 25.89mic, A from*/
	ipT26 = ipFindLevLine( 258800 , 4 , 8 );

	ipT1032 = ipFindLevLine( 1031.9f , 6 , 8 );

	ipT1037 = ipFindLevLine( 1037.6f , 6 , 8 );

	/*fluorine*/
	/*f ii 29.33 mic*/
	ipF0229 = ipFindLevLine( 293300 , 2 , 9 );

	/*f ii 67.2 mic*/
	ipF0267 = ipFindLevLine( 672000 , 2 , 9 );

	/*f iv 44.07 mic*/
	ipF444 = ipFindLevLine( 440700 , 4 , 9 );

	/*f iv 25.83 mic*/
	ipF425 = ipFindLevLine( 258300 , 4 , 9 );

	/*neon*/
	ipT770 = ipFindLevLine( 770.4f , 8 , 10 );

	ipT780 = ipFindLevLine( 780.3f , 8 , 10 );

	/*[Ne VI] 7.652 micron*/
	ipxNe0676 = ipFindLevLine( 76520 , 6 , 10 );

	ipT895 = ipFindLevLine( 895 , 7 , 10 );

	ipT88 = ipFindLevLine( 88 , 8 , 10 );

	/*fine structure lines */
	/*[Ne II] 12.8*/
	ipTNe13 = ipFindLevLine( 128100 , 2 , 10 );

	/*[Ne III] 36.013 m*/
	ipTNe36 = ipFindLevLine( 360140 , 3 , 10 );

	/*[Ne III] 15.56 m*/
	ipTNe16 = ipFindLevLine( 155500 , 3 , 10 );

	/*[Ne V] 14.32 m*/
	ipTNe14 = ipFindLevLine( 143200 , 5 , 10 );

	/*[Ne V] 24.318 m*/
	ipTNe24 = ipFindLevLine( 243100 , 5 , 10 );

	/*sodium line optical depth data*/
	ipT5895 = ipFindLevLine( 5891.9f , 1 , 11 );

	/*[Na III] 7.3177mic*/
	ipfsNa373 = ipFindLevLine( 73200 , 3 , 11 );

	/*[Na IV] 9.039 mic*/
	ipfsNa490 = ipFindLevLine( 90390 , 4 , 11 );

	/*[Na IV] 21.29 mic*/
	ipfsNa421 = ipFindLevLine( 212900 , 4 , 11 );

	/*[Na VI] 14.40 mic*/
	ipxNa6143 = ipFindLevLine( 144000 , 6 , 11 );

	/*[Na VI] 8.611 mic*/
	ipxNa6862 = ipFindLevLine( 86110 , 6 , 11 );

	/*[Na VII] 4.68 micron*/
	ipxNa0746 = ipFindLevLine( 46800 , 7 , 11 );

	/*magnesium line optical depth data*/
	ipMgI2853 = ipFindLevLine( 2853 , 1 , 12 );

	ipMgI2026 = ipFindLevLine( 2026 , 1 , 12 );

	ipT2796 = ipFindLevLine( 2795.5f , 2 , 12 );

	ipT2804 = ipFindLevLine( 2802.7f , 2 , 12 );

	ipT705 = ipFindLevLine( 705 , 9 , 12 );

	ipT4561 = ipFindLevLine( 4561 , 1 , 12 );

	/*[Mg V] 1325, 3-1 in three level atom*/
	ipxMg51325 = ipFindLevLine( 1325 , 5 , 12 );

	/*[Mg V] 2417, 3-2 in three level atom*/
	ipxMg52417 = ipFindLevLine( 2417 , 5 , 12 );

	/*[Mg V] 2855, 2-1 in three level atom, really 2928, 2782*/
	ipxMg52855 = ipFindLevLine( 2855 , 5 , 12 );

	/*[Mg VII] 1325, 3-1 in three level atom*/
	ipxMg71190 = ipFindLevLine( 1190 , 7 , 12 );

	/*[Mg VII] 2261, 3-2 in three level atom*/
	ipxMg72261 = ipFindLevLine( 2261 , 7 , 12 );

	/*[Mg VII] 2569, 2-1 in three level atom, really 2509, 2629*/
	ipxMg72569 = ipFindLevLine( 2569 , 7 , 12 );

	/*[Mg VIII] 3.03 micron*/
	ipxMg08303 = ipFindLevLine( 30300 , 8 , 12 );

	/*the Mg X 615 li seq doublet*/
	ipTMg610 = ipFindLevLine( 609.8f , 10 , 12 );

	ipTMg625 = ipFindLevLine( 624.9f , 10 , 12 );

	ipT58 = ipFindLevLine( 57.9f , 10 , 12 );

	/*Mg IV 4.487m*/
	ipTMg4 = ipFindLevLine( 44850 , 4 , 12 );

	/*[Mg V] 13.52*/
	ipTMg14 = ipFindLevLine( 135200 , 5 , 12 );

	/*[Mg V] 5.610*/
	ipTMg6 = ipFindLevLine( 56100 , 5 , 12 );

	/*[Mg VII] 9.033 mic*/
	ipfsMg790 = ipFindLevLine( 90330 , 7 , 12 );

	/*[Mg VII] 5.503 mic7 , 12 );*/
	ipfsMg755 = ipFindLevLine( 55030 , 7 , 12 );

	/* aluminum line optical depth data */
	ipAlI3957 = ipFindLevLine( 3957 , 1 , 13 );

	ipAlI3090 = ipFindLevLine( 3090 , 1 , 13 );

	ipT1855 = ipFindLevLine( 1855 , 3 , 13 );

	ipT1863 = ipFindLevLine( 1863 , 3 , 13 );

	ipT2670 = ipFindLevLine( 2670 , 2 , 13 );

	/*[Al V] 2.905 mic*/
	ipAl529 = ipFindLevLine( 29052 , 5 , 13 );

	/*[Al VI] 3.66 mic*/
	ipAl6366 = ipFindLevLine( 36600 , 6 , 13 );

	/*[Al VI] 9.116 mic*/
	ipAl6912 = ipFindLevLine( 91160 , 6 , 13 );

	/*[Al VIII] 5.848 mic*/
	ipAl8575 = ipFindLevLine( 58480 , 8 , 13 );

	/*[Al VIII] 3.69 mic*/
	ipAl8370 = ipFindLevLine( 36900 , 8 , 13 );

	/*[Al IX] 2.04 micron*/
	ipAl09204 = ipFindLevLine( 20400 , 9 , 13 );

	ipT639 = ipFindLevLine( 639 , 10 , 13 );

	/*Al XI 550, 568 Li seq doublet */
	ipTAl550 = ipFindLevLine( 550 , 11 , 13 );

	ipTAl568 = ipFindLevLine( 568 , 11 , 13 );

	ipTAl48 = ipFindLevLine( 48 , 11 , 13 );

	/*silicon line optical depth data*/
	ipSi1_130m = ipFindLevLine( 1295823.0f , 1 , 14 );
	ipSi1_68m = ipFindLevLine( 683995.25f , 1 , 14 );

	ipSii2518 = ipFindLevLine( 2518 , 1 , 14 );

	ipSii2215 = ipFindLevLine( 2215 , 1 , 14 );

	ipSi2_2334 = ipFindLevLine( 2334 , 2 , 14 );
	ipSi2_2329 = ipFindLevLine( 2329 , 2 , 14 );
	ipSi2_2350 = ipFindLevLine( 2350 , 2 , 14 );
	ipSi2_2344 = ipFindLevLine( 2344 , 2 , 14 );
	ipSi2_2336 = ipFindLevLine( 2336 , 2 , 14 );

	ipT1808 = ipFindLevLine( 1813.99f , 2 , 14 );

	ipT1527 = ipFindLevLine( 1531 , 2 , 14 );

	ipT1305 = ipFindLevLine( 1307.7f , 2 , 14 );

	ipT1260 = ipFindLevLine( 1263.3f , 2 , 14 );

	ipT1207 = ipFindLevLine( 1207 , 3 , 14 );

	ipT1895 = ipFindLevLine( 1892 , 3 , 14 );

	ipT1394 = ipFindLevLine( 1394 , 4 , 14 );

	ipT1403 = ipFindLevLine( 1403 , 4 , 14 );

	/*[Si VI] 1.9641mic*/
	ipSi619 = ipFindLevLine( 19631 , 6 , 14 );

	/*[Si X] 1.43 micron*/
	ipSi10143 = ipFindLevLine( 14300 , 10 , 14 );

	/* Si X] ll 621.1, 611.7, 598.6 */
	ipSi10_606 = ipFindLevLine( 606 , 10 , 14 );

	ipTSi499 = ipFindLevLine( 499 , 12 , 14 );

	ipTSi521 = ipFindLevLine( 521 , 12 , 14 );

	ipTSi41 = ipFindLevLine( 41 , 12 , 14 );

	/*[Si II] 34.8 mic*/
	ipTSi35 = ipFindLevLine( 348140 , 2 , 14 );

	/*[Si VII] 2.481*/
	ipTSi25 = ipFindLevLine( 24810 , 7 , 14 );

	/*[Si VII] 6.4922*/
	ipTSi65 = ipFindLevLine( 64920 , 7 ,14 );

	/*[Si IX] 2.585 mic*/
	ipTSi3 = ipFindLevLine( 25840 , 9 , 14 );

	/*[Si IX] 3.929 mic*/
	ipTSi4 = ipFindLevLine( 39290 , 9 , 14 );

	/*phosphorus line data*/
	/*P II 60.64 mic*/
	ipP0260 = ipFindLevLine( 606400 , 2 , 15 );

	/*P II 32.87 mic*/
	ipP0233 = ipFindLevLine( 328700 , 2 , 15 );

	/*P III 17.885 mic*/
	ipP0318 = ipFindLevLine( 178850 , 3 , 15 );

	/*P VII 1.3745 mic*/
	ipP713 = ipFindLevLine( 13745 , 7 , 15 );

	/*P VIII 4.85 mic*/
	ipP848 = ipFindLevLine( 48500 , 8 , 15 );

	/*P VIII 1.735 mic*/
	ipP817 = ipFindLevLine( 17350 , 8 , 15 );

	/*P X 2.708 mic*/
	ipP1027 = ipFindLevLine( 27080 ,  10 , 15 );

	/*P X 1.868*/
	ipP1018 = ipFindLevLine( 18680 , 10 , 15 );

	/*sulphur line optical depth data*/
	ipS1_25m = ipFindLevLine( 251947.453f , 1 , 16 );
	ipS1_56m = ipFindLevLine( 562909.625f , 1 , 16 );

	ipT1256 = ipFindLevLine( 1256 , 2 , 16 );

	ipT1194 = ipFindLevLine( 1197.55f , 3 , 16 );

	ipTS1720 = ipFindLevLine( 1720 , 3 , 16 );

	ipS4_1405 = ipFindLevLine( 1405 , 4 , 16 );
	ipS4_1398 = ipFindLevLine( 1398 , 4 , 16 );
	ipS4_1424 = ipFindLevLine( 1424 , 4 , 16 );
	ipS4_1417 = ipFindLevLine( 1417 , 4 , 16 );
	ipS4_1407 = ipFindLevLine( 1406 , 4 , 16 );

	ipT1198 = ipFindLevLine( 1198 , 5 , 16 );

	ipT786 = ipFindLevLine( 786.47f , 5 , 16 );

	ipT933 = ipFindLevLine( 933.38f , 6 , 16 );

	ipT944 = ipFindLevLine( 944.52f , 6 , 16 );

	/*[S VIII] 9913*/
	ipfsS810 = ipFindLevLine( 9913 , 8 , 16 );

	/*[S IX] 1.252 mic*/
	ipfsS912 = ipFindLevLine( 12520 , 9 , 16 );

	/*[S IX] 3.755 mic*/
	ipfsS938 = ipFindLevLine( 37550 , 9 , 16 );

	/*[S XI] 1.9201 mic*/
	ipfsS1119 = ipFindLevLine( 19201 , 11 , 16 );

	/*[S XI] 1.3927 mic*/
	ipfsS1114 = ipFindLevLine( 13927 , 11 , 16 );

	/*S XII 7536*/
	/*this is a vacuum wavelength*/
	ipfsS1207 = ipFindLevLine( 7611 , 12 , 16 );

	/* >>chng 03 sep 27,  group of lines moved from expion */
	ipS12_520 = ipFindLevLine( 520 , 12 , 16 );

	ipTSu418 = ipFindLevLine( 418 , 14 , 16 );

	ipTSu446 = ipFindLevLine( 446 , 14 , 16 );

	ipTSu30 = ipFindLevLine( 30 , 14 , 16 );

	/*fine structure lines added in by Jim Kingdon*/
	ipTS19 = ipFindLevLine( 186700 , 3 , 16 );

	/*[S III] 33.48*/
	ipTS34 = ipFindLevLine( 334700 , 3 , 16 );

	ipTS11 = ipFindLevLine( 105100 , 4 , 16 );

	/*chlorine line optical depth data*/

	/* Cl 1 mic */
	ipCl1_11m = ipFindLevLine( 113296.3984f , 1 , 17 );

	/*Cl II 14.3678 mic*/
	ipfsCl214 = ipFindLevLine( 144000 , 2 , 17 );

	/*Cl II 33.281 mic*/
	ipfsCl233 = ipFindLevLine( 333000 , 2 , 17 );

	/*[Cl 4] 20.354 mic*/
	ipCl04203 = ipFindLevLine( 204000 ,  4 , 17 );

	/*[Cl 4] 11.741 mic*/
	ipCl04117 = ipFindLevLine( 117000 ,4 , 17 );

	/*Cl IX 7334A*/
	ipCl973 = ipFindLevLine( 7334 , 9 , 17 );

	/*Cl X 3.05 mic*/
	ipCl1030 = ipFindLevLine( 30500 , 10 , 17 );

	/*Cl X 9332*/
	ipCl1092 = ipFindLevLine( 9332 , 10 , 17 );

	/*argon line optical depth data*/
	/*Ar XVI*/
	ipT354 = ipFindLevLine( 354 , 16 , 18 );

	/*Ar XVI*/
	ipT389 = ipFindLevLine( 389 , 16 , 18 );

	/*Ar XVI*/
	ipT25 = ipFindLevLine( 25 , 16 , 18 );

	/*argon fine structure lines*/
	/*[Ar II]*/
	ipTAr7 = ipFindLevLine( 69800 , 2 , 18 );

	/*[Ar III] 9.0 mic*/
	ipTAr9 = ipFindLevLine( 90000 , 3 , 18 );

	/*[Ar III] 21.83 mic*/
	ipTAr22 = ipFindLevLine( 218300 , 3 , 18 );

	/*[Ar V] 13.1 mic*/
	ipTAr13 = ipFindLevLine( 131000 , 5 , 18 );

	/*[Ar V] 8.0 mic*/
	ipTAr8 = ipFindLevLine( 80000 , 5 , 18 );

	/*[Ar VI] 4.53 micron*/
	ipAr06453 = ipFindLevLine( 45300 , 6 , 18 );

	/*Ar X 5533*/
	ipAr1055 = ipFindLevLine( 5533 , 10 , 18 );

	/*[Ar XI] 2.619 mic*/
	ipAr1126 = ipFindLevLine( 26190 , 11 , 18 );

	/*[Ar XI] 6917.5*/
	ipAr1178 = ipFindLevLine( 6917 , 11 , 18 );

	/*potasium - really should split into two*/
	ipKI7745 = ipFindLevLine( 7676.2f , 1 , 19 );

	/*[K III] 4.62 micron*/
	ipxK03462 = ipFindLevLine( 46200 , 3 , 19 );

	/*[K IV] 5.982 micron*/
	ipxK04598 = ipFindLevLine( 59800 , 4 , 19 );

	/*[K IV] 15.39 micron4 , 19 );*/
	ipxK04154 = ipFindLevLine( 153900 , 4 , 19 );

	/* [K VI] 8.823 micron */
	ipxK06882 = ipFindLevLine( 88200 , 6 , 19 );

	/*[K VI] 5.575 micron6 */
	ipxK06557 = ipFindLevLine( 55700 , 6 , 19 );

	/*[K VII] 3.19 micron7 , 19 );*/
	ipxK07319 = ipFindLevLine( 31905 , 7 , 19 );

	/* [K XI] 4249.99 A  */
	ipxK11425 = ipFindLevLine( 4250 , 11 , 19 );

	/*calcium line optical depth data*/
	ipCaI4228 = ipFindLevLine( 4228 , 1 , 20 );

	ipT3934 = ipFindLevLine( 3934 , 2 , 20 );

	ipT3969 = ipFindLevLine( 3969 , 2 , 20 );

	ipT8498 = ipFindLevLine( 8498 , 2 , 20 );

	ipT8542 = ipFindLevLine( 8542 , 2 , 20 );

	ipT8662 = ipFindLevLine( 8662 , 2 , 20 );

	ipT7291 = ipFindLevLine( 7291 ,  2 , 20 );

	ipT7324 = ipFindLevLine( 7324 ,  2 , 20 );

	ipTCa345 = ipFindLevLine( 345 , 18 , 20 );

	ipTCa19 = ipFindLevLine( 19 , 18 , 20 );

	/*[Ca IV] 3.21 min*/
	ipTCa3 = ipFindLevLine( 32100 , 4 , 20 );

	ipTCa12 = ipFindLevLine( 114000 , 5 , 20 );

	ipTCa4 = ipFindLevLine( 41500 , 5 , 20 );

	/*[Ca VII] 4.09 micron*/
	ipCa0741 = ipFindLevLine( 40900 , 7 , 20 );

	/*[Ca VII] 6.15 micron*/
	ipCa0761 = ipFindLevLine( 61500 , 7 , 20 );

	/*[Ca VIII] 2.32 micron*/
	ipCa08232 = ipFindLevLine( 23200 , 8 , 20 );

	/*[Ca 12] 3328.78 A */
	ipCa12333 = ipFindLevLine( 3329 , 12 , 20 );

	/* Ca 18 302 */
	ipTCa302 = ipFindLevLine( 302 , 18 , 20 );

	/*scandium data*/
	/*[Sc V] 2.31 micron*/
	ipSc05231 = ipFindLevLine( 23100 , 5 , 21 );

	/*[Sc 13] 2637.97 A */
	ipSc13264 = ipFindLevLine( 2638 , 13 , 21 );

	/*Titanium data*/
	/*[Ti VI] 1.72 micron*/
	ipTi06172 = ipFindLevLine( 17200 , 6 , 22 );

	/*[Ti 14] 2117.79 A */
	ipTi14212 = ipFindLevLine( 2118 , 14 , 22 );

	/*Vanadium data*/
	/*[V VII] 1.30 micron*/
	ipVa07130 = ipFindLevLine( 13000 , 7 , 23 );

	/*[V 15] 1721.38 A */
	ipVa15172 = ipFindLevLine( 1721 , 15 , 23 );

	/*Chromium data*/
	/*[Cr VIII] 1.01 micron*/
	ipCr08101 = ipFindLevLine( 10100 , 8 , 24 );

	/*[Cr 16] 1410.60 A*/
	ipCr16141 = ipFindLevLine( 1411 , 16 , 24 );

	/*Manganese data*/
	/*[Mn IX] 7968.5A*/
	ipxMn0979 = ipFindLevLine( 7968 , 9 , 25 );

	/*[Mn 17] 1169.59 A */
	ipxMn1712 = ipFindLevLine( 1170 , 17 , 25 );

	/*iron Fe line optical depth data*/

	/* [Fe I] ground term */
	ipFe1_24m = ipFindLevLine( 240359.546f , 1 , 26 );
	ipFe1_35m = ipFindLevLine( 347043.25f , 1 , 26 );
	ipFe1_54m = ipFindLevLine( 542946.5625f , 1 , 26 );
	ipFe1_111m = ipFindLevLine( 1111549.25f , 1 , 26 );

	ipFeI3884 = ipFindLevLine( 3884 , 1 , 26 );

	ipFeI3729 = ipFindLevLine( 3729 , 1 , 26 );

	ipFeI3457 = ipFindLevLine( 3457 ,  1 , 26 );

	ipFeI3021 = ipFindLevLine( 3021 ,  1 , 26 );

	ipFeI2966 = ipFindLevLine( 2966 ,  1 , 26 );

	/*>>chng 03 oct 07, chng wl from 2360 to 2400, the value
	 * before cnhg of feii lines to stuc */
	ipTuv3 = ipFindLevLine( 2400 ,  2 , 26 );

	ipTr48 = ipFindLevLine( 6200 ,  2 , 26 );

	ipTFe16 = ipFindLevLine( 1080 ,  2 , 26 );

	ipTFe26 = ipFindLevLine( 1500 ,  2 , 26 );

	ipTFe34 = ipFindLevLine( 11500 ,  2 , 26 );

	ipTFe35 = ipFindLevLine( 2500 ,  2 , 26 );

	ipTFe46 = ipFindLevLine( 2300 ,  2 , 26 );

	ipTFe56 = ipFindLevLine( 8900 ,  2 , 26 );

	ipT1122 = ipFindLevLine( 1125.8f ,  3 , 26 );

	/* [Fe VII] 9.53 mic, A */
	ipFe0795 = ipFindLevLine( 95100,  7 , 26 );

	/* [Fe VII] 7.81 mic, A */
	ipFe0778 = ipFindLevLine( 78100, 7 , 26 );

	ipT245 = ipFindLevLine( 245 , 9 , 26 );

	/*this is pump line for upper level of [FeX] 6375*/
	ipT352 = ipFindLevLine( 352.1f , 10 , 26 );

	ipFe106375 = ipFindLevLine( 6373.3f , 10 , 26 );

	/*actually three lines contributing here, important for coronal pumping */
	/*following is whole multiplet average*/
	ipT353 = ipFindLevLine( 353 , 11 , 26 );

	/* >>chng 05 dec 19, go to 5 level Fe 11 */
	/*ipTFe07 = ipFindLevLine( 7892 , 11 , 26 );
	ipTFe61 = ipFindLevLine( 60827 , 11 , 26 );*/

	/*[Fe XIII] 1.07 microns
	ipFe1310 = ipFindLevLine( 10700 , 13 , 26 );*/

	/*[Fe XIII] 1.08 microns
	ipFe1311 = ipFindLevLine( 10800 , 13 , 26 );*/

	// [Fe XVII] 17.1 A M2 line
	ipFe17_17 = ipFindLevLine( (realnum)17.1 , 17 , 26 );

	/* N-seq [Fe20] 721, 578A */
	ipTFe20_721 = ipFindLevLine( 721 , 20 , 26 );
	ipTFe20_578 = ipFindLevLine( 578 , 20 , 26 );

	/* boron sequence fe iron lines */
	/*>>chng 00 oct 26 added these lines */
	ipFe22_247 = ipFindLevLine( (realnum)247.2 , 22 , 26 );
	ipFe22_217 = ipFindLevLine( (realnum)217.3 , 22 , 26 );
	ipFe22_348 = ipFindLevLine( (realnum)349.3 , 22 , 26 );
	ipFe22_292 = ipFindLevLine( (realnum)292.4 , 22 , 26 );
	ipFe22_253 = ipFindLevLine( (realnum)253.2 , 22 , 26 );
	ipFe22_846 = ipFindLevLine( 846 , 22 , 26 );

	ipT347 = ipFindLevLine( 347 , 14 , 26 );

	ipT192 = ipFindLevLine( 192 , 24 , 26 );

	ipT255 = ipFindLevLine( 255 , 24 , 26 );

	ipT11 = ipFindLevLine( 11 , 24 , 26 );

	ipT191 = ipFindLevLine( 1786 , 2 , 26 );


	/* [Fe 18] 974.86 A */
	ipFe18975 = ipFindLevLine( 974 , 18 , 26 );

	ipTFe23 = ipFindLevLine( 2299 , 21 , 26 );

	ipTFe13 = ipFindLevLine( 1354 , 21 , 26 );

	/* Cobalt data */
	/* [Co XI] 5168.A */
	ipCo11527 = ipFindLevLine( 5168 , 11 , 27 );

	/* Nickel data */
	ipNi1_7m = ipFindLevLine( 75046.164f , 1 , 28 );
	ipNi1_11m = ipFindLevLine( 113044.031f , 1 , 28 );

	/* [Ni XII] 4230.8A */
	ipxNi1242 = ipFindLevLine( 4230 , 12 , 28 );

	/* Zinc data */
	/* Zn 04 3.625 microns */
	ipZn04363 = ipFindLevLine( 36250 , 4 , 30 );

	/* flush the line list, freeing the extra storage and checking that all
	 * lines have been claimed */
	endFindLevLine( );

	/** \todo	1	streamline all of this, using TransitionZero
	 * and then setting dangerously large negative numbers. */

	/* only do this one time, and only if number of atom_level2 lines is positive */
	if( lgFirst && nWindLine>0)
	{ 

		lgFirst = false;
		/* these are the massive set of op lines, with g-bar approx cs
		 * confirm that input data are valid */

		for( i=0; i < nWindLine; i++ )
		{
			/* this information was read in in createdata */
			ASSERT( TauLine2[i].Hi->nelem > 0 );
			ASSERT( TauLine2[i].Hi->nelem <= (int)LIMELM );

			ASSERT( TauLine2[i].Hi->IonStg > 0 );
			ASSERT( TauLine2[i].Hi->IonStg <= (int)LIMELM );

			ASSERT( TauLine2[i].Lo->g >0. );

			ASSERT( TauLine2[i].Hi->g > 0. );

			/* check that energy is positive*/
			ASSERT( TauLine2[i].EnergyWN > 0 );

			/* TauLine2[i].Emis->gf this is gf if positive, A if negative */
			/* test whether a or gf entered, convert A to gf */
			if( TauLine2[i].Emis->gf < 0. )
			{
				/* convert A (=-gf) into real gf */
				TauLine2[i].Emis->gf *= (realnum)((-TauLine2[i].Hi->g)/TRANS_PROB_CONST/POW2(TauLine2[i].EnergyWN));
			}

			/*now put into standard format */
			TauLine2[i].WLAng = 1.e8f/TauLine2[i].EnergyWN;
			TauLine2[i].Lo->Pop = 0.;
			TauLine2[i].Hi->Pop = 0.;
			TauLine2[i].Emis->iRedisFun = ipPRD;

			/* these are line optical depth arrays
			 * inward optical depth */
			TauLine2[i].Emis->TauIn = opac.taumin;
			TauLine2[i].Emis->TauCon = opac.taumin;
			TauLine2[i].Emis->ColOvTot = 0.;
			/* outward optical depth */
			TauLine2[i].Emis->TauTot = 1e20f;
			/* escape probability */
			TauLine2[i].Emis->Pesc = 1.;
			/* inward part of line */
			TauLine2[i].Emis->FracInwd = 1.;
			/* destruction probability */
			TauLine2[i].Emis->Pdest = 0.;
			TauLine2[i].Emis->Pelec_esc = 0.;
			/* line pumping rate */
			TauLine2[i].Emis->pump = 0.;
			/* population of lower level */
			TauLine2[i].Lo->Pop = 0.;
			/* population of upper level */
			TauLine2[i].Hi->Pop = 0.;
			/* population of lower level with correction for stim emission */
			TauLine2[i].Emis->PopOpc = 0.;
			/* following two heat exchange excitation, deexcitation */
			TauLine2[i].Coll.cool = 0.;
			TauLine2[i].Coll.heat = 0.;
			/* intensity of line */
			TauLine2[i].Emis->xIntensity = 0.;
			/* number of photons emitted in line */
			TauLine2[i].Emis->phots = 0.;
			/* ots rate */
			TauLine2[i].Emis->ots = 0.;
		}
	}

	for( i=0; i < nUTA; i++ )
	{
		/* this information was read in in createdata */
		ASSERT( UTALines[i].Hi->nelem > 0 );
		ASSERT( UTALines[i].Hi->nelem <= (int)LIMELM );

		ASSERT( UTALines[i].Hi->IonStg > 0 );
		ASSERT( UTALines[i].Hi->IonStg <= (int)LIMELM );

		ASSERT( UTALines[i].Lo->g > 0. );

		ASSERT( UTALines[i].Hi->g > 0. );

		/* check that energy is positive*/
		ASSERT( UTALines[i].EnergyWN > 0 );

		UTALines[i].Lo->Pop = 0.;
		UTALines[i].Hi->Pop = 0.;
		UTALines[i].Emis->iRedisFun = ipPRD;

		/* these are line optical depth arrays
		 * inward optical depth */
		UTALines[i].Emis->TauIn = opac.taumin;
		UTALines[i].Emis->TauCon = opac.taumin;
		UTALines[i].Emis->ColOvTot = 0.;
		/* outward optical depth */
		UTALines[i].Emis->TauTot = 1e20f;
		/* escape probability */
		UTALines[i].Emis->Pesc = 1.;
		/* inward part of line */
		UTALines[i].Emis->FracInwd = 1.;
		/* destruction probability */
		UTALines[i].Emis->Pdest = 0.;
		UTALines[i].Emis->Pelec_esc = 0.;
		/* line pumping rate */
		UTALines[i].Emis->pump = 0.;
		/* population of lower level */
		UTALines[i].Lo->Pop = 0.;
		/* population of upper level */
		UTALines[i].Hi->Pop = 0.;
		/* population of lower level with correction for stim emission */
		UTALines[i].Emis->PopOpc = 0.;
		/* following two heat exchange excitation, deexcitation */
		UTALines[i].Coll.cool = 0.;
		/* heat is the net heat per pump and was set when data read in
		 * this is different from other lines with this structure 
		UTALines[i].Coll.heat = 0.;*/
		/* intensity of line */
		UTALines[i].Emis->xIntensity = 0.;
		/* number of photons emitted in line */
		UTALines[i].Emis->phots = 0.;
		UTALines[i].Emis->ots = 0.;
	}

	for( i=0; i < nHFLines; i++ )
	{
		/* this information was read in in createdata */
		ASSERT( HFLines[i].Hi->nelem > 0 );
		ASSERT( HFLines[i].Hi->nelem <= (int)LIMELM );

		ASSERT( HFLines[i].Hi->IonStg > 0 );
		ASSERT( HFLines[i].Hi->IonStg <= (int)LIMELM );

		ASSERT( HFLines[i].Lo->g > 0. );

		ASSERT( HFLines[i].Hi->g > 0. );

		/* check that energy is positive*/
		ASSERT( HFLines[i].EnergyWN > 0 );
		ASSERT( HFLines[i].Emis->Aul>0 );
		ASSERT( HFLines[i].Emis->damp>0 );

		/* HFLines[i].Emis->gf this is gf if positive, A if negative */
		/* test whether a or gf entered, convert A to gf */
		if( HFLines[i].Emis->gf < 0. )
		{
			/* convert A (=-gf) into real gf */
			HFLines[i].Emis->gf *= (realnum)(-HFLines[i].Hi->g/TRANS_PROB_CONST/POW2(HFLines[i].EnergyWN));
		}

		/*now put into standard format */
		HFLines[i].WLAng = 1.e8f/HFLines[i].EnergyWN;
		HFLines[i].Lo->Pop = 0.;
		HFLines[i].Hi->Pop = 0.;
		/* change from partial to complete redistribution */
		HFLines[i].Emis->iRedisFun = ipCRD;

		/* these are line optical depth arrays
		 * inward optical depth */
		HFLines[i].Emis->TauIn = opac.taumin;
		HFLines[i].Emis->TauCon = opac.taumin;
		HFLines[i].Emis->ColOvTot=0;
		/* outward optical depth */
		HFLines[i].Emis->TauTot = 1e20f;
		/* escape probability */
		HFLines[i].Emis->Pesc = 1.;
		/* inward part of line */
		HFLines[i].Emis->FracInwd = 1.;
		/* destruction probability */
		HFLines[i].Emis->Pdest = 0.;
		HFLines[i].Emis->Pelec_esc = 0.;
		/* line pumping rate */
		HFLines[i].Emis->pump = 0.;
		/* population of lower level */
		HFLines[i].Lo->Pop = 0.;
		/* population of upper level */
		HFLines[i].Hi->Pop = 0.;
		/* population of lower level with correction for stim emission */
		HFLines[i].Emis->PopOpc = 0.;
		/* following two heat exchange excitation, deexcitation */
		HFLines[i].Coll.cool = 0.;
		HFLines[i].Coll.heat = 0.;
		/* intensity of line */
		HFLines[i].Emis->xIntensity = 0.;
		/* number of photons emitted in line */
		HFLines[i].Emis->phots = 0.;
		HFLines[i].Emis->ots = 0.;
	}
	return;
}

/* following used to save whether lines have been claimed by a pointer */
static int *lev2set;

/*generate pointer to level 1 line using wavelengtgh, ion, element */
STATIC long ipFindLevLine( 
		/* realnum ID wavelength, in angstroms */
		realnum wl , 
		/* state of ionization, 1 for neutral atom */
		long ion , 
		/* element number, 1 for H, 26 for Fe */
		long nelem)
{
	long i;/* use for counter in for loop */

	DEBUG_ENTRY( "ipFindLevLine()" );

	ASSERT( wl > 0. );
	ASSERT( ion > 0 );
	ASSERT( ion <= LIMELM );
	ASSERT( nelem > 0 );
	ASSERT( nelem <= LIMELM );

	/* look for the line */
	for( i=1; i<= nLevel1; ++i )
	{
		if( TauLines[i].Hi->nelem == (int)nelem &&
		    TauLines[i].Hi->IonStg == (int)ion &&
			// wl hit within relative accuracy of 0.05A in 1000A
		    fabs(TauLines[i].WLAng - wl)/max(1000.,wl) < 0.05/1000. )
		{
			/* remember that we have hit this line */
			lev2set[i] = true;
			/* and return pointer to the label*/
			return i;
		}
	}
	fprintf(ioQQQ,
		" ipFindLevLine could not find a line with following properties:\n"
		" wavelength=%f\n"
		" ion stage =%li\n"
		" atomic num=%li\n",
		wl , ion, nelem );
	return -1;
}

STATIC void initFindLevLine( void )
{
	long i;

	DEBUG_ENTRY( "initFindLevLine()" );

	/* generate the array of ints to store true and false */
	lev2set = (int*)MALLOC( (size_t)(nLevel1+1)*sizeof(int) );

	/* set them all false, saying that they have not been claimed by
	 * one of the line pointers */
	for( i=1; i<=nLevel1; ++i )
		lev2set[i] = false;
	return;
}

STATIC void endFindLevLine( void )
{
	long i;
	bool lgAbort_loc=false;

	DEBUG_ENTRY( "endFindLevLine()" );

	/* set them all false, saying that they have not been claimed by
	 * one of the line pointers */
	for( i=1; i<=nLevel1; ++i )
	{
		if( !lev2set[i] )
		{
			fprintf(ioQQQ,"PROBLEM endFindLevLine warning; line %li not claimed\n",i);
			fprintf(ioQQQ,
				" line had the following properties:\n"
				" wavelength=%f\n"
				" ion stage =%i\n"
				" atomic num=%i\n",
				TauLines[i].WLAng ,
				TauLines[i].Hi->IonStg ,
				TauLines[i].Hi->nelem );  
			lgAbort_loc = true;
		}
	}

	/* generate the array of ints to store true and false */
	free(lev2set);

	if( lgAbort_loc )
	{
		fprintf(ioQQQ," problems found entering the data.  I live in lines_setup.c\n");
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
