/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtFinal create PrtFinal pages of printout, emission line intensities, etc */
/*StuffComment routine to stuff comments into the stack of comments, def in lines.h */
/*gett2o3 analyze computed [OIII] spectrum to get t^2 */
/*gett2 analyze computed structure to get structural t^2 */
#include "cddefines.h"
#include "iso.h"
#include "cddrive.h"
#include "dynamics.h"
#include "physconst.h"
#include "lines.h"
#include "taulines.h"
#include "warnings.h"
#include "phycon.h"
#include "dense.h"
#include "grainvar.h"
#include "h2.h"
#include "hmi.h"
#include "thermal.h"
#include "hydrogenic.h"
#include "rt.h"
#include "atmdat.h"
#include "timesc.h"
#include "opacity.h"
#include "struc.h"
#include "pressure.h"
#include "conv.h"
#include "geometry.h"
#include "called.h"
#include "iterations.h"
#include "version.h"
#include "colden.h"
#include "input.h"
#include "rfield.h"
#include "radius.h"
#include "peimbt.h"
#include "oxy.h"
#include "ipoint.h"
#include "lines_service.h"
#include "mean.h"
#include "wind.h"
#include "prt.h"

/* DatabasePrintReference print some database references */
void DatabasePrintReference( void )
{
	DEBUG_ENTRY( "DatabasePrintReference()" );

	fprintf(ioQQQ,"\n Some databases used in this calculation:\n");

	fprintf(ioQQQ," Many recombination coefficients are based on Badnell ApJS 167, "
		"334 (2006) and Badnell et al A&A 406, 1151 (2003) and posted on\n"
		"      the web sites http://amdpp.phys.strath.ac.uk/tamoc/RR/ and "
		"http://amdpp.phys.strath.ac.uk/tamoc/DR/.\n");

	if( atmdat.lgLamdaOn )
	{
		fprintf( ioQQQ, " Much of the molecular emission data is from LAMDA "
			"(Schoeier et al 2005, A&A 432, 369-379) as accessed on Feb. 9, 2010.\n");
	}

	if( atmdat.lgChiantiOn )
	{
		fprintf( ioQQQ, " Much of the ionic emission data is from CHIANTI. CHIANTI is a collaborative project involving the NRL (USA),"
			" the Universities of Florence\n      (Italy) and Cambridge (UK), and George Mason University (USA). "
			"(Dere, K. P., Landi, E., Mason, H. E., Monsignori Fossi, B. C., & \n      Young, P. R. 1997, A&AS, 125, 149) "
			"(Dere, K. P., Landi, E., Young, P. R., Del Zanna, G., Landini, M. & Mason, H. E. 2009, A&A, 498,\n      915) using version %s.\n",atmdat.chVersion);
	}

	if( h2.lgH2ON )
	{
		fprintf( ioQQQ, " Much of the H_2 data is from "
			"Wrathmall & Flower 2007, JPhysB, 40 3221\n"
			"      Abgrall et al, 1994, Can. J. Phys., 72, 856 at http://molat.obspm.fr"
			"/index.php?page=pages/Molecules/H2/H2can94.php\n");
	}

	fprintf( ioQQQ, " Please cite those papers that played a significant role in this calculation.\n");
}

// helper routine to center a line in the output
void PrintCenterLine(FILE* io,              // file pointer
		     const char chLine[],   // string to print, should not end with '\n'
		     size_t ArrLen,         // length of chLine
		     size_t LineLen)        // width of the line
{
	unsigned long StrLen = min(strlen(chLine),ArrLen);
	ASSERT( StrLen < LineLen );
	unsigned long pad = (LineLen-StrLen)/2;
	for( unsigned long i=0; i < pad; ++i )
		fprintf( io, " " );
	fprintf( io, "%s\n", chLine );
}

/*gett2o3 analyze computed [OIII] spectrum to get t^2 */
STATIC void gett2o3(realnum *tsqr);

/*gett2 analyze computed structure to get structural t^2 */
STATIC void gett2(realnum *tsqr);

/* helper routine for printing averaged quantities */
inline void PrintRatio(double q1, double q2)
{
	double ratio = ( q2 > SMALLFLOAT ) ? q1/q2 : 0.;
	fprintf( ioQQQ, " " );
	fprintf( ioQQQ, PrintEfmt("%9.2e", ratio) );
	return;
}

/* return is true if calculation ok, false otherwise */
void PrtFinal(void)
{
	short int *lgPrt;
	realnum *wavelength;
	realnum *sclsav , *scaled;
	long int *ipSortLines;
	double *xLog_line_lumin;
	char **chSLab;
	char *chSTyp;
	char chCKey[5], 
	  chGeo[7], 
	  chPlaw[21];

	long int 
	  i, 
	  ipEmType ,
	  ipNegIntensity[33], 
	  ip2500, 
	  ip2kev, 
	  iprnt, 
	  j, 
	  nline, 
	  nNegIntenLines;
	double o4363, 
	  bacobs, 
	  absint, 
	  bacthn, 
	  hbcab, 
	  hbeta, 
	  o5007;

	double a, 
	  ajmass, 
	  ajmin, 
	  alfox, 
	  bot, 
	  bottom, 
	  bremtk, 
	  coleff, 
	  /* N.B. 8 is used for following preset in loop */
	  d[8], 
	  dmean, 
	  ferode, 
	  he4686, 
	  he5876, 
	  heabun, 
	  hescal, 
	  pion, 
	  plow, 
	  powerl, 
	  qion, 
	  qlow, 
	  rate, 
	  ratio, 
	  reclin, 
	  rjeans, 
	  snorm,
	  tmean, 
	  top, 
	  THI,/* HI 21 cm spin temperature */
	  uhel, 
	  uhl, 
	  usp, 
	  wmean, 
	  znit;

	  double bac, 
	  tel, 
	  x;

	DEBUG_ENTRY( "PrtFinal()" );

	/* return if not talking */
	/* >>chng 05 nov 11, also if nzone is zero, sign of abort before model got started */
	if( !called.lgTalk || (lgAbort && nzone< 1)  )
	{ 
		return;
	}

	/* print out header, or parts that were saved */

	/* this would be a major logical error */
	ASSERT( LineSave.nsum > 1 );

	/* print name and version number */
	fprintf( ioQQQ, "\f\n");
	fprintf( ioQQQ, "%23c", ' ' );
	int len = (int)strlen(t_version::Inst().chVersion);
	int repeat = (72-len)/2;
	for( i=0; i < repeat; ++i )
		fprintf( ioQQQ, "*" );
	fprintf( ioQQQ, "> Cloudy %s <", t_version::Inst().chVersion );
	for( i=0; i < 72-repeat-len; ++i )
		fprintf( ioQQQ, "*" );
	fprintf( ioQQQ, "\n" );

	for( i=0; i <= input.nSave; i++ )
	{
		char chCard[INPUT_LINE_LENGTH];
		/* print command line, unless it is a continue line */

		/* copy start of command to key, making it into capitals */
		cap4(chCKey,input.chCardSav[i]);

		/* copy input to CAPS to make sure hide not on it */
		strcpy( chCard , input.chCardSav[i] );
		caps( chCard );

		/* print if not continue or hide */
		/* >>chng 04 jan 21, add hide option */
		if( (strcmp(chCKey,"CONT")!= 0) && !nMatch( "HIDE" , chCard) )
			fprintf( ioQQQ, "%23c* %-80s*\n", ' ', input.chCardSav[i] );
	}

	/* print log of ionization parameter if greater than zero - U==0 for PDR calcs */
	if( rfield.uh > 0. )
	{
		a = log10(rfield.uh);
	}
	else
	{
		a = -37.;
	}

	fprintf( ioQQQ, 
		"                       *********************************> Log(U):%6.2f <*********************************\n", 
	  a );

	if( t_version::Inst().nBetaVer > 0 )
	{
		fprintf( ioQQQ, 
			"\n                        This is a beta test version of the code, and is intended for testing only.\n\n" );
	}

	if( warnings.lgWarngs )
	{
		fprintf( ioQQQ, "  \n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!  Warnings exist, this calculation has serious problems.\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "  \n" );
	}
	else if( warnings.lgCautns )
	{
		fprintf( ioQQQ, 
			"                       >>>>>>>>>> Cautions are present.\n" );
	}

	if( conv.lgBadStop )
	{
		fprintf( ioQQQ, 
			"                       >>>>>>>>>> The calculation stopped for unintended reasons.\n" );
	}

	if( iterations.lgIterAgain )
	{
		fprintf( ioQQQ, 
			"                       >>>>>>>>>> Another iteration is needed.\n" );
	}

	/* open or closed geometry?? */
	if( geometry.lgSphere )
	{
		strcpy( chGeo, "Closed" );
	}
	else
	{
		strcpy( chGeo, "  Open" );
	}

	/* now give description of pressure law */
	if( strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		strcpy( chPlaw, " Constant  Pressure " );
	}

	else if( strcmp(dense.chDenseLaw,"CDEN") == 0 )
	{
		strcpy( chPlaw, "  Constant  Density " );
	}

	else if( (strcmp(dense.chDenseLaw,"POWD") == 0 || strcmp(dense.chDenseLaw
	  ,"POWR") == 0) || strcmp(dense.chDenseLaw,"POWC") == 0 )
	{
		strcpy( chPlaw, "  Power Law Density " );
	}

	else if( strcmp(dense.chDenseLaw,"SINE") == 0 )
	{
		strcpy( chPlaw, " Rapid Fluctuations " );
	}

	else if( strncmp(dense.chDenseLaw , "DLW" , 3) == 0 )
	{
		strcpy( chPlaw, " Special Density Lw " );
	}

	else if( strcmp(dense.chDenseLaw,"HYDR") == 0 )
	{
		strcpy( chPlaw, " Hydrostatic Equlib " );
	}

	else if( strcmp(dense.chDenseLaw,"WIND") == 0 )
	{
		strcpy( chPlaw, "  Radia Driven Wind " );
	}

	else if( strcmp(dense.chDenseLaw,"GLOB") == 0 )
	{
		strcpy( chPlaw, " Globule            " );
	}

	else
	{
		strcpy( chPlaw, " UNRECOGNIZED CPRES " );
	}

	fprintf( ioQQQ, 
		"\n                     Emission Line Spectrum. %20.20sModel.  %6.6s geometry.  Iteration %ld of %ld.\n", 
	  chPlaw, chGeo, iteration, iterations.itermx + 1 );

	/* emission lines as logs of total luminosity */
	if(  radius.distance > 0. && radius.lgRadiusKnown && prt.lgPrintFluxEarth )
	{
		char chLine[INPUT_LINE_LENGTH];
		const char* chUnit;
		if( geometry.iEmissPower == 1 && !geometry.lgSizeSet )
			chUnit = "/arcsec";
		else if( geometry.iEmissPower == 0 && !geometry.lgSizeSet )
			chUnit = "/arcsec^2";
		else
			chUnit = "";

		sprintf( chLine, "Flux observed at the Earth (erg/s/cm^2%s).", chUnit );
		PrintCenterLine( ioQQQ, chLine, sizeof(chLine), 130 );
	}
	else if(  prt.lgSurfaceBrightness )
	{
		char chLine[INPUT_LINE_LENGTH];
		const char* chUnit;
		if( prt.lgSurfaceBrightness_SR )
			chUnit = "sr";
		else
			chUnit = "arcsec^2";

		sprintf( chLine, "Surface brightness (erg/s/cm^2/%s).", chUnit );
		PrintCenterLine( ioQQQ, chLine, sizeof(chLine), 130 );
	}
	else if( radius.lgPredLumin )
	{
		char chLine[INPUT_LINE_LENGTH];
		const char* chUnit;
		if( geometry.iEmissPower == 2 )
			chUnit = "erg/s";
		else if( geometry.iEmissPower == 1 )
			chUnit = "erg/s/cm";
		else if( geometry.iEmissPower == 0 )
			chUnit = "erg/s/cm^2";
		else
			TotalInsanity();

		char chCoverage[INPUT_LINE_LENGTH];
		if( fp_equal( geometry.covgeo, realnum(1.) ) )
			sprintf( chCoverage, "with full coverage" );
		else
			sprintf( chCoverage, "with a covering factor of %.1f%%", geometry.covgeo*100. );

		if( radius.lgCylnOn )
			sprintf( chLine, "Luminosity (%s) emitted by a partial cylinder %s.", chUnit, chCoverage );
		else
			sprintf( chLine, "Luminosity (%s) emitted by a shell %s.", chUnit, chCoverage );
		PrintCenterLine( ioQQQ, chLine, sizeof(chLine), 130 );

		if( geometry.iEmissPower != 2 )
		{
			const char* chAper;
			if( geometry.iEmissPower == 1 )
				chAper = "long slit";
			else if( geometry.iEmissPower == 0 )
				chAper = "pencil beam";
			else
				TotalInsanity();

			sprintf( chLine, "Observed through a %s with aperture covering factor %.1f%%.",
				 chAper, geometry.covaper*100. );
			PrintCenterLine( ioQQQ, chLine, sizeof(chLine), 130 );
		}
	}
	else
	{
		char chLine[INPUT_LINE_LENGTH];
		sprintf( chLine, "Intensity (erg/s/cm^2)." );
		PrintCenterLine( ioQQQ, chLine, sizeof(chLine), 130 );
	}

	fprintf( ioQQQ, "\n" );

	/******************************************************************
	 * kill Hummer & Storey case b predictions if outside their table *
	 ******************************************************************/

	/* >>chng 02 aug 29, from lgHCaseBOK to following - caught by Ryan Porter */
	if( !atmdat.lgHCaseBOK[1][ipHYDROGEN] )
	{
		static const int NWLH = 21;
		/* list of all H case b lines */
		realnum wlh[NWLH] = {6563.e0f, 4861.e0f, 4340.e0f, 4102.e0f, 1.875e4f, 1.282e4f, 
						   1.094e4f, 1.005e4f, 2.625e4f, 2.166e4f, 1.945e4f, 7.458e4f, 
						   4.653e4f, 3.740e4f, 4.051e4f, 7.458e4f, 3.296e4f, 1216.e0f,
						   1026.e0f, 972.5e0f, 949.7e0f };

		/* table exceeded - kill all H case b predictions*/
		for( i=0; i < LineSave.nsum; i++ )
		{
			/* >>chng 04 jun 21, kill both case a and b at same time,
			 * actually lgHCaseBOK has separate A and B flags, but
			 * this is simpler */
			if( (strcmp( LineSv[i].chALab,"Ca B" )==0) || 
				(strcmp( LineSv[i].chALab,"Ca A" )==0) )
			{
				realnum errorwave;
				/* this logic must be kept parallel with which H lines are
				 * added as case B in lines_hydro - any automatic hosing of
				 * case b would kill both H I and He II */
				errorwave = WavlenErrorGet( LineSv[i].wavelength );
				for( j=0; j<NWLH; ++j )
				{
					if( fabs(LineSv[i].wavelength-wlh[j] ) <= errorwave )
					{
						LineSv[i].SumLine[0] = 0.;
						LineSv[i].SumLine[1] = 0.;
						break;
					}
				}
			}
		}
	}

	if( !atmdat.lgHCaseBOK[1][ipHELIUM] )
	{
		/* table exceeded - kill all He case b predictions*/
		static const int NWLHE = 20;
		realnum wlhe[NWLHE] = {1640.e0f, 1215.e0f, 1085.e0f, 1025.e0f, 4686.e0f, 3203.e0f,
						     2733.e0f, 2511.e0f, 1.012e4f, 6560.e0f, 5412.e0f, 4860.e0f,
						     1.864e4f, 1.163e4f, 9345.e0f, 8237.e0f, 303.8e0f, 256.3e0f,
						     243.0e0f, 237.3e0f};
		for( i=0; i < LineSave.nsum; i++ )
		{
			if( (strcmp( LineSv[i].chALab,"Ca B" )==0) || 
				(strcmp( LineSv[i].chALab,"Ca A" )==0) )
			{
				realnum errorwave;
				/* this logic must be kept parallel with which H lines are
				 * added as case B in lines_hydro - any automatic hosing of
				 * case b would kill both H I and He II */
				errorwave = WavlenErrorGet( LineSv[i].wavelength );
				for( j=0; j<NWLHE; ++j )
				{
					if( fabs(LineSv[i].wavelength-wlhe[j] ) <= errorwave )
					{
						LineSv[i].SumLine[0] = 0.;
						LineSv[i].SumLine[1] = 0.;
						break;
					}
				}
			}
		}
	}

	/**********************************************************
	 *analyse line arrays for abundances, temp, etc           *
	 **********************************************************/

	/* find apparent helium abundance, will not find these if helium is turned off */

	if( cdLine("TOTL",4861.36f,&hbeta,&absint)<=0 )
	{
		if( dense.lgElmtOn[ipHYDROGEN] )
		{
			/* this is a major logical error if hydrogen is turned on */
			fprintf( ioQQQ, " PrtFinal could not find TOTL 4861 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* get HeI 5876 */
		/* >>chng 06 may 15, changed this so that it works for up to six sig figs. */
		if( cdLine("He 1",5875.61f,&he5876,&absint)<=0 )
		{
			/* 06 aug 28, from numLevels_max to _local. */
			if( iso.numLevels_local[ipHE_LIKE][ipHELIUM] >= 14 )
			{
				/* this is a major logical error if helium is turned on */
				fprintf( ioQQQ, " PrtFinal could not find He 1 5876 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* get HeII 4686 */
		/* >>chng 06 may 15, changed this so that it works for up to six sig figs. */
		if( cdLine("He 2",4686.01f,&he4686,&absint)<=0 )
		{
			/* 06 aug 28, from numLevels_max to _local. */
			if( iso.numLevels_local[ipH_LIKE][ipHELIUM] > 5 )
			{
				/* this is a major logical error if helium is turned on 
				 * and the model atom has enough levels */
				fprintf( ioQQQ, " PrtFinal could not find He 2 4686 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	else
	{
		he5876 = 0.;
		absint = 0.;
		he4686 = 0.;
	}

	if( hbeta > 0. )
	{
		heabun = (he4686*0.078 + he5876*0.739)/hbeta;
	}
	else
	{
		heabun = 0.;
	}

	if( dense.lgElmtOn[ipHELIUM] )
	{
		hescal = heabun/(dense.gas_phase[ipHELIUM]/dense.gas_phase[ipHYDROGEN]);
	}
	else
	{
		hescal = 0.;
	}

	/* get temperature from [OIII] spectrum, o may be turned off,
	 * but lines still dumped into stack */
	if( cdLine("O  3",5007.,&o5007,&absint)<=0 )
	{
		if( dense.lgElmtOn[ipOXYGEN] )
		{
			/* this is a major logical error if hydrogen is turned on */
			fprintf( ioQQQ, " PrtFinal could not find O  3 5007 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( cdLine("TOTL",4363.,&o4363,&absint)<=0 )
	{
		if( dense.lgElmtOn[ipOXYGEN] )
		{
			/* this is a major logical error if hydrogen is turned on */
			fprintf( ioQQQ, " PrtFinal could not find TOTL 4363 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* first find low density limit OIII zone temp */
	if( (o4363 != 0.) && (o5007 != 0.) )
	{
		/* following will assume coll excitation only, so correct
		 * for 4363's that cascade as 5007 */
		bot = o5007 - o4363/oxy.o3enro;

		if( bot > 0. )
		{
			ratio = o4363/bot;
		}
		else
		{
			ratio = 0.178;
		}

		if( ratio > 0.177 )
		{
			/* ratio was above infinite temperature limit */
			peimbt.toiiir = 1e10;
		}
		else
		{
			/* data for following set in OIII cooling routine
			 * ratio of collision strengths, factor of 4/3 makes 5007+4959
			 * >>chng 96 sep 07, reset cs to values at 10^4K,
			 * had been last temp in model */
			/*>>chng 06 jul 25, update to recent data */
			oxy.o3cs12 = 2.2347f;
			oxy.o3cs13 = 0.29811f;
			ratio = ratio/1.3333/(oxy.o3cs13/oxy.o3cs12);
			/* ratio of energies and branching ratio for 4363 */
			ratio = ratio/oxy.o3enro/oxy.o3br32;
			peimbt.toiiir = (realnum)(-oxy.o3ex23/log(ratio));
		}
	}

	else
	{
		peimbt.toiiir = 0.;
	}

	if( geometry.iEmissPower == 2 )
	{
		/* find temperature from Balmer continuum */
		if( cdLine("Bac ",3646.,&bacobs,&absint)<=0 )
		{
			fprintf( ioQQQ, " PrtFinal could not find Bac 3546 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* we pulled hbeta out of stack with cdLine above */
		if( hbeta > 0. )
		{
			bac = bacobs/hbeta;
		}
		else
		{
			bac = 0.;
		}
	}
	else
	{
		bac = 0.;
	}

	if( bac > 0. )
	{
		/*----------------------------------------------------------*
		 ***** TableCurve c:\tcwin2\balcon.for Sep 6, 1994 11:13:38 AM
		 ***** log bal/Hbet
		 ***** X= log temp
		 ***** Y= 
		 ***** Eqn# 6503  y=a+b/x+c/x^2+d/x^3+e/x^4+f/x^5
		 ***** r2=0.9999987190883581
		 ***** r2adj=0.9999910336185065
		 ***** StdErr=0.001705886369042427
		 ***** Fval=312277.1895753243
		 ***** a= -4.82796940090671
		 ***** b= 33.08493347410885
		 ***** c= -56.08886262205931
		 ***** d= 52.44759454803217
		 ***** e= -25.07958990094209
		 ***** f= 4.815046760060006
		 *----------------------------------------------------------*
		 * bac is double precision!!!! */
		x = 1.e0/log10(bac);
		tel = -4.827969400906710 + x*(33.08493347410885 + x*(-56.08886262205931 + 
		  x*(52.44759454803217 + x*(-25.07958990094209 + x*(4.815046760060006)))));

		if( tel > 1. && tel < 5. )
		{
			peimbt.tbac = (realnum)pow(10.,tel);
		}
		else
		{
			peimbt.tbac = 0.;
		}
	}
	else
	{
		peimbt.tbac = 0.;
	}

	if( geometry.iEmissPower == 2 )
	{
		/* find temperature from optically thin Balmer continuum and case B H-beta */
		if( cdLine("thin",3646.,&bacthn,&absint)<=0 )
		{
			fprintf( ioQQQ, " PrtFinal could not find thin 3646 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* >>chng 06 may 15, changed this so that it works for up to six sig figs. */
		if( cdLine("Ca B",4861.36f,&hbcab,&absint)<=0 )
		{
			fprintf( ioQQQ, " PrtFinal could not find Ca B 4861 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( hbcab > 0. )
		{
			bacthn /= hbcab;
		}
		else
		{
			bacthn = 0.;
		}
	}
	else
	{
		bacthn = 0.;
	}

	if( bacthn > 0. )
	{
		/*----------------------------------------------------------*
		 ***** TableCurve c:\tcwin2\balcon.for Sep 6, 1994 11:13:38 AM
		 ***** log bal/Hbet
		 ***** X= log temp
		 ***** Y= 
		 ***** Eqn# 6503  y=a+b/x+c/x^2+d/x^3+e/x^4+f/x^5
		 ***** r2=0.9999987190883581
		 ***** r2adj=0.9999910336185065
		 ***** StdErr=0.001705886369042427
		 ***** Fval=312277.1895753243
		 ***** a= -4.82796940090671
		 ***** b= 33.08493347410885
		 ***** c= -56.08886262205931
		 ***** d= 52.44759454803217
		 ***** e= -25.07958990094209
		 ***** f= 4.815046760060006
		 *----------------------------------------------------------*
		 * bac is double precision!!!! */
		x = 1.e0/log10(bacthn);
		tel = -4.827969400906710 + x*(33.08493347410885 + x*(-56.08886262205931 + 
		  x*(52.44759454803217 + x*(-25.07958990094209 + x*(4.815046760060006)))));

		if( tel > 1. && tel < 5. )
		{
			peimbt.tbcthn = (realnum)pow(10.,tel);
		}
		else
		{
			peimbt.tbcthn = 0.;
		}
	}
	else
	{
		peimbt.tbcthn = 0.;
	}

	/* we now have temps from OIII ratio and BAC ratio, now
	 * do Peimbert analysis, getting To and t^2 */

	peimbt.tohyox = (realnum)((8.5*peimbt.toiiir - 7.5*peimbt.tbcthn - 228200. + 
	  sqrt(POW2(8.5*peimbt.toiiir-7.5*peimbt.tbcthn-228200.)+9.128e5*
	  peimbt.tbcthn))/2.);

	if( peimbt.tohyox > 0. )
	{
		peimbt.t2hyox = (realnum)((peimbt.tohyox - peimbt.tbcthn)/(1.7*peimbt.tohyox));
	}
	else
	{
		peimbt.t2hyox = 0.;
	}

	/*----------------------------------------------------------
	 *
	 * first set scaled lines */

	/* get space for scaled */
	scaled = (realnum *)MALLOC( sizeof(realnum)*(unsigned long)LineSave.nsum);

	/* get space for xLog_line_lumin */
	xLog_line_lumin = (double *)MALLOC( sizeof(double)*(unsigned long)LineSave.nsum);

	/* this is option to not print certain contributions */
	/* gjf 98 jun 10, get space for array lgPrt */
	lgPrt = (short int *)MALLOC( sizeof(short int)*(unsigned long)LineSave.nsum);

	/* get space for wavelength */
	wavelength = (realnum *)MALLOC( sizeof(realnum)*(unsigned long)LineSave.nsum);

	/* get space for sclsav */
	sclsav = (realnum *)MALLOC( sizeof(realnum)*(unsigned long)LineSave.nsum );

	/* get space for chSTyp */
	chSTyp = (char *)MALLOC( sizeof(char)*(unsigned long)LineSave.nsum );

	/* get space for chSLab,
	 * first define array of pointers*/
	/* char chSLab[NLINES][5];*/
	chSLab = ((char**)MALLOC((unsigned long)LineSave.nsum*sizeof(char*)));

	/* now allocate all the labels for each of the above lines */
	for( i=0; i<LineSave.nsum; ++i)
	{
		chSLab[i] = (char*)MALLOC(5*sizeof(char));
	}

	/* get space for array of indices for lines, for possible sorting */
	ipSortLines = (long *)MALLOC( sizeof(long)*(unsigned long)LineSave.nsum );

	ASSERT( LineSave.ipNormWavL >= 0 );

	/* option to also print usual first two sets of line arrays 
	 * but for two sets of cumulative arrays for time-dependent sims too */
	int nEmType = 2;
	if( prt.lgPrintLineCumulative && iteration > dynamics.n_initial_relax )
		nEmType = 4;

	for( ipEmType=0; ipEmType<nEmType; ++ipEmType )
	{
		/* this is the intensity of the line spectrum will be normalized to */
		snorm = LineSv[LineSave.ipNormWavL].SumLine[ipEmType];

		/* check that this line has positive intensity */
		if( ((snorm <= SMALLDOUBLE ) || (LineSave.ipNormWavL < 0)) || (LineSave.ipNormWavL > LineSave.nsum) )
		{
			fprintf( ioQQQ, "\n\n >>PROBLEM Normalization line has small or zero intensity, its value was %.2e and its intensity was set to 1."
				"\n >>Please consider using another normalization line (this is set with the NORMALIZE command).\n" , snorm);
			fprintf( ioQQQ, " >>The relative intensities will be meaningless, and many lines may appear too faint.\n" );
			snorm = 1.;
		}
		for( i=0; i < LineSave.nsum; i++ )
		{

			/* when normalization line is off-scale small (generally a model
			 * with mis-set parameters) the scaled intensity can be larger than
			 * a realnum - this is not actually a problem since the number will
			 * overflow the format and hence be unreadable */
			double scale = LineSv[i].SumLine[ipEmType]/snorm*LineSave.ScaleNormLine;
			/* this will become a realnum, so limit dynamic range */
			scale = MIN2(BIGFLOAT , scale );
			scale = MAX2( -BIGFLOAT , scale );

			/* find logs of ALL line intensities/luminosities */
			scaled[i] = (realnum)scale;

			if( LineSv[i].SumLine[ipEmType] > 0. )
			{
				xLog_line_lumin[i] = log10(LineSv[i].SumLine[ipEmType]) + radius.Conv2PrtInten;
			}
			else
			{
				xLog_line_lumin[i] = -38.;
			}
		}

		/* now find which lines to print, which to ignore because they are the wrong type */
		for( i=0; i < LineSave.nsum; i++ )
		{
			/* never print unit normalization check, at least in main list */
			if( (strcmp(LineSv[i].chALab,"Unit") == 0) || (strcmp(LineSv[i].chALab,"UntD") == 0) )
				lgPrt[i] = false;
			else if( strcmp(LineSv[i].chALab,"Coll") == 0 && !prt.lgPrnColl )
				lgPrt[i] = false;
			else if( strcmp(LineSv[i].chALab,"Pump") == 0 && !prt.lgPrnPump )
				lgPrt[i] = false;
			else if( strncmp(LineSv[i].chALab,"Inw",3) == 0 && !prt.lgPrnInwd )
				lgPrt[i] = false;
			else if( strcmp(LineSv[i].chALab,"Heat") == 0 && !prt.lgPrnHeat )
				lgPrt[i] = false;
			/* these are diffuse and transmitted continua - normally do not print
			 * nFnu or nInu */
			else if( !prt.lgPrnDiff && 
				( (strcmp(LineSv[i].chALab,"nFnu") == 0) || (strcmp(LineSv[i].chALab,"nInu") == 0) ) )
				lgPrt[i] = false;
			else
				lgPrt[i] = true;
		}

		/* do not print relatively faint lines unless requested */
		nNegIntenLines = 0;

		/* set ipNegIntensity to bomb to make sure set in following */
		for(i=0; i< 32; i++ )
		{
			ipNegIntensity[i] = LONG_MAX;
		}

		for(i=0;i<8;++i)
		{
			d[i] = -DBL_MAX;
		}

		/* create header for blocks of emission line intensities */
		const char chIntensityType[4][100]=
		{"     Intrinsic" , "      Emergent" , "Cumulative intrinsic" , "Cumulative emergent" };
		ASSERT( ipEmType==0 || ipEmType==1 || ipEmType==2 || ipEmType==3 );
		/* if true then printing in 4 columns of lines, this is offset to
		 * center the title */
		fprintf( ioQQQ, "\n" );
		if( prt.lgPrtLineArray )
			fprintf( ioQQQ, "                                              " );
		fprintf( ioQQQ, "%s" , chIntensityType[ipEmType] );
		fprintf( ioQQQ, " line intensities\n" );
		// caution about emergent intensities when outward optical
		// depths are not yet known
		if( ipEmType==1 && iteration==1 )
			fprintf(ioQQQ," Caution: emergent intensities are not reliable on the "
			"first iteration.\n");

		/* option to only print brighter lines */
		if( prt.lgFaintOn )
		{
			iprnt = 0;
			for( i=0; i < LineSave.nsum; i++ )
			{
				/* this insanity can happen when arrays overrun */
				ASSERT( iprnt <= i);
				if( scaled[i] >= prt.TooFaint && lgPrt[i] )
				{
					if( prt.lgPrtLineLog )
					{
						xLog_line_lumin[iprnt] = log10(LineSv[i].SumLine[ipEmType]) + radius.Conv2PrtInten;
					}
					else
					{
						xLog_line_lumin[iprnt] = LineSv[i].SumLine[ipEmType] * pow(10.,radius.Conv2PrtInten);
					}
					sclsav[iprnt] = scaled[i];
					chSTyp[iprnt] = LineSv[i].chSumTyp;
					/* check that null is correct, string overruns have 
					 * been a problem in the past */
					ASSERT( strlen( LineSv[i].chALab )<5 );
					strcpy( chSLab[iprnt], LineSv[i].chALab );
					wavelength[iprnt] = LineSv[i].wavelength;
					++iprnt;
				}
				else if( -scaled[i] > prt.TooFaint && lgPrt[i] )
				{
					/* negative intensities occur if line absorbs continuum */
					ipNegIntensity[nNegIntenLines] = i;
					nNegIntenLines = MIN2(32,nNegIntenLines+1);
				}
				/* special labels to give id for blocks of lines 
				 * do not add these labels when sorting by wavelength since invalid */
				else if( strcmp( LineSv[i].chALab, "####" ) == 0  &&!prt.lgSortLines )
				{
					strcpy( chSLab[iprnt], LineSv[i].chALab );
					xLog_line_lumin[iprnt] = 0.;
					sclsav[iprnt] = 0.;
					chSTyp[iprnt] = LineSv[i].chSumTyp;
					wavelength[iprnt] = LineSv[i].wavelength;
					++iprnt;
				}
			}
		}

		else
		{
			/* print everything */
			iprnt = LineSave.nsum;
			for( i=0; i < LineSave.nsum; i++ )
			{
				if( strcmp( LineSv[i].chALab, "####" ) == 0  )
				{
					strcpy( chSLab[i], LineSv[i].chALab );
					xLog_line_lumin[i] = 0.;
					sclsav[i] = 0.;
					chSTyp[i] = LineSv[i].chSumTyp;
					wavelength[i] = LineSv[i].wavelength;
				}
				else
				{
					sclsav[i] = scaled[i];
					chSTyp[i] = LineSv[i].chSumTyp;
					strcpy( chSLab[i], LineSv[i].chALab );
					wavelength[i] = LineSv[i].wavelength;
				}
				if( scaled[i] < 0. )
				{
					ipNegIntensity[nNegIntenLines] = i;
					nNegIntenLines = MIN2(32,nNegIntenLines+1);
				}
			}
		}

		/* the number of lines to print better be positive */
		ASSERT( iprnt > 0. );

		/* reorder lines according to wavelength for comparison with spectrum
		 * including writing out the results */
		if( prt.lgSortLines )
		{
			int lgFlag;
			if( prt.lgSortLineWavelength )
			{
				/* first check if wavelength range specified */
				if( prt.wlSort1 >-0.1 )
				{
					j = 0;
					/* skip over those lines not desired */
					for( i=0; i<iprnt; ++i )
					{
						if( wavelength[i]>= prt.wlSort1 && wavelength[i]<= prt.wlSort2 )
						{
							if( j!=i )
							{
								sclsav[j] = sclsav[i];
								chSTyp[j] = chSTyp[i];
								strcpy( chSLab[j], chSLab[i] );
								wavelength[j] = wavelength[i];
								/* >>chng 05 feb 03, add this, had been left out, 
								 * thanks to Marcus Copetti for discovering the bug */
								xLog_line_lumin[j] = xLog_line_lumin[i];
							}
							++j;
						}
					}
					iprnt = j;
				}
				/*spsort netlib routine to sort array returning sorted indices */
				spsort(wavelength, 
					   iprnt, 
					  ipSortLines, 
					  /* flag saying what to do - 1 sorts into increasing order, not changing
					   * the original routine */
					  1, 
					  &lgFlag);
				if( lgFlag ) 
					TotalInsanity();
			}
			else if( prt.lgSortLineIntensity )
			{
				/*spsort netlib routine to sort array returning sorted indices */
				spsort(sclsav, 
					   iprnt, 
					  ipSortLines, 
					  /* flag saying what to do - -1 sorts into decreasing order, not changing
					   * the original routine */
					  -1, 
					  &lgFlag);
				if( lgFlag ) 
					TotalInsanity();
			}
			else
			{
				/* only to keep lint happen, or in case vars clobbered */
				TotalInsanity();
			}
		}
		else
		{
			/* do not sort lines (usual case) simply print in incoming order */
			for( i=0; i<iprnt; ++i )
			{
				ipSortLines[i] = i;
			}
		}

		/* print out all lines which made it through the faint filter,
		 * there are iprnt lines to print 
		 * can print in either 4 column (the default ) or one long
		 * column of lines */
		if( prt.lgPrtLineArray )
		{
			/* four or three columns ? - depends on how many sig figs */
			if( LineSave.sig_figs >= 5 )
			{
				nline = (iprnt + 2)/3;
			}
			else
			{
				/* nline will be the number of horizontal lines - 
				* the 4 represents the 4 columns */
				nline = (iprnt + 3)/4;
			}
		}
		else
		{
			/* this option print a single column of emission lines */
			nline = iprnt;
		}

		/* now loop over the spectrum, printing lines */
		for( i=0; i < nline; i++ )
		{
			fprintf( ioQQQ, " " );

			/* this skips over nline per step, to go to the next column in 
			 * the output */
			for( j=i; j<iprnt; j = j + nline)
			{ 
				/* index with possibly reordered set of lines */
				long ipLin = ipSortLines[j];
				/* this is the actual print statement for the emission line
				 * spectrum*/
				if( strcmp( chSLab[ipLin], "####" ) == 0  )
				{
					/* special labels */
					/*fprintf( ioQQQ, "1111222223333333444444444      " ); */
					/* this string was set in */
					fprintf( ioQQQ, "%s",LineSave.chHoldComments[(int)wavelength[ipLin]] ); 
				}
				else
				{
					/* the label for the line */
					fprintf( ioQQQ, "%4.4s ",chSLab[ipLin] ); 

					/* print the wavelength for the line */
					prt_wl( ioQQQ , wavelength[ipLin]);

					/* print the log of the intensity/luminosity of the 
					 * lines, the usual case */
					if( prt.lgPrtLineLog )
					{
						fprintf( ioQQQ, " %7.3f", xLog_line_lumin[ipLin] );
					}
					else
					{
						/* print linear quantity instead */
						fprintf( ioQQQ, "  %.2e ", xLog_line_lumin[ipLin] );
					}

					/* print scaled intensity with various formats,
					 * depending on order of magnitude.  value is
					 * always the same but the format changes. */
					if( sclsav[ipLin] < 9999.5 )
					{
						fprintf( ioQQQ, "%9.4f",sclsav[ipLin] );
					}
					else if( sclsav[ipLin] < 99999.5 )
					{
						fprintf( ioQQQ, "%9.3f",sclsav[ipLin] );
					}
					else if( sclsav[ipLin] < 999999.5 )
					{
						fprintf( ioQQQ, "%9.2f",sclsav[ipLin] );
					}
					else if( sclsav[ipLin] < 9999999.5 )
					{
						fprintf( ioQQQ, "%9.1f",sclsav[ipLin] );
					}
					else
					{
						fprintf( ioQQQ, "*********" );
					}

					/* now end the block with some spaces to set next one off */
					fprintf( ioQQQ, "      " );
				}
			}
			/* now end the lines */
			fprintf( ioQQQ, "      \n" );
		}

		if( nNegIntenLines > 0 )
		{
			fprintf( ioQQQ, " Lines with negative intensities -  Linear intensities relative to normalize line.\n" );
			fprintf( ioQQQ, "          " );

			for( i=0; i < nNegIntenLines; ++i )
			{
				fprintf( ioQQQ, "%ld %s %.0f %.1e, ", 
				  ipNegIntensity[i], 
				  LineSv[ipNegIntensity[i]].chALab, 
				  LineSv[ipNegIntensity[i]].wavelength, 
				  scaled[ipNegIntensity[i]] );
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* now find which were the very strongest, more that 5% of cooling */
	j = 0;
	for( i=0; i < LineSave.nsum; i++ )
	{
		a = (double)LineSv[i].SumLine[0]/(double)thermal.totcol;
		if( (a >= 0.05) && LineSv[i].chSumTyp == 'c' )
		{
			ipNegIntensity[j] = i;
			d[j] = a;
			j = MIN2(j+1,7);
		}
	}

	fprintf( ioQQQ, "\n\n\n %s\n", input.chTitle );
	if( j != 0 )
	{
		fprintf( ioQQQ, " Cooling: " );
		for( i=0; i < j; i++ )
		{
			fprintf( ioQQQ, " %4.4s ", 
				LineSv[ipNegIntensity[i]].chALab);

			prt_wl( ioQQQ, LineSv[ipNegIntensity[i]].wavelength );

			fprintf( ioQQQ, ":%5.3f", 
				d[i] );
		}
		fprintf( ioQQQ, "  \n" );
	}

	/* now find strongest heating, more that 5% of total */
	j = 0;
	for( i=0; i < LineSave.nsum; i++ )
	{
		a = (double)LineSv[i].SumLine[0]/(double)thermal.power;
		if( (a >= 0.05) && LineSv[i].chSumTyp == 'h' )
		{
			ipNegIntensity[j] = i;
			d[j] = a;
			j = MIN2(j+1,7);
		}
	}

	if( j != 0 )
	{
		fprintf( ioQQQ, " Heating: " );
		for( i=0; i < j; i++ )
		{
			fprintf( ioQQQ, " %4.4s ", 
				LineSv[ipNegIntensity[i]].chALab);

			prt_wl(ioQQQ, LineSv[ipNegIntensity[i]].wavelength);

			fprintf( ioQQQ, ":%5.3f", 
				d[i] );
		}
		fprintf( ioQQQ, "  \n" );
	}

	DatabasePrintReference();

	/* print out ionization parameters and radiation making it through */
	qlow = 0.;
	plow = 0.;
	for( i=0; i < (iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s] - 1); i++ )
	{
		/* N.B. in following en1ryd prevents overflow */
		plow += (rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i])*
		  EN1RYD*rfield.anu[i];
		qlow += rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i];
	}

	qlow *= radius.r1r0sq;
	plow *= radius.r1r0sq;
	if( plow > 0. )
	{
		qlow = log10(qlow) + radius.Conv2PrtInten;
		plow = log10(plow) + radius.Conv2PrtInten;
	}
	else
	{
		qlow = 0.;
		plow = 0.;
	}

	qion = 0.;
	pion = 0.;
	for( i=iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1; i < rfield.nflux; i++ )
	{
		/* N.B. in following en1ryd prevents overflow */
		pion += (rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i])*
		  EN1RYD*rfield.anu[i];
		qion += rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i];
	}

	qion *= radius.r1r0sq;
	pion *= radius.r1r0sq;

	if( pion > 0. )
	{
		qion = log10(qion) + radius.Conv2PrtInten;
		pion = log10(pion) + radius.Conv2PrtInten;
	}
	else
	{
		qion = 0.;
		pion = 0.;
	}

	/* derive ionization parameter for spherical geometry */
	if( rfield.qhtot > 0. )
	{
		if( rfield.lgUSphON )
		{
			/* RSTROM is stromgren radius */
			usp = rfield.qhtot/POW2(rfield.rstrom/radius.rinner)/
			  2.998e10/dense.gas_phase[ipHYDROGEN];
			usp = log10(usp);
		}
		else
		{
			/* no stromgren radius found, use outer radius */
			usp = rfield.qhtot/radius.r1r0sq/2.998e10/dense.gas_phase[ipHYDROGEN];
			usp = log10(usp);
		}
	}

	else
	{
		usp = -37.;
	}

	if( rfield.uh > 0. )
	{
		uhl = log10(rfield.uh);
	}
	else
	{
		uhl = -37.;
	}

	if( rfield.uheii > 0. )
	{
		uhel = log10(rfield.uheii);
	}
	else
	{
		uhel = -37.;
	}

	fprintf( ioQQQ, 
		"\n IONIZE PARMET:  U(1-)%8.4f  U(4-):%8.4f  U(sp):%6.2f  "
		"Q(ion):  %7.3f  L(ion)%7.3f    Q(low):%7.3f    L(low)%7.3f\n", 
	  uhl, uhel, usp, qion, pion, qlow, plow );

	a = fabs((thermal.power-thermal.totcol)*100.)/thermal.power;
	/* output power and total cooling; can be neg for shocks, collisional ioniz */
	if( thermal.power > 0. )
	{
		powerl = log10(thermal.power) + radius.Conv2PrtInten;
	}
	else
	{
		powerl = 0.;
	}

	if( thermal.totcol > 0. )
	{
		thermal.totcol = log10(thermal.totcol) + radius.Conv2PrtInten;
	}
	else
	{
		thermal.totcol = 0.;
	}

	if( thermal.FreeFreeTotHeat > 0. )
	{
		thermal.FreeFreeTotHeat = log10(thermal.FreeFreeTotHeat) + radius.Conv2PrtInten;
	}
	else
	{
		thermal.FreeFreeTotHeat = 0.;
	}

	/* following is recombination line intensity */
	reclin = totlin('r');
	if( reclin > 0. )
	{
		reclin = log10(reclin) + radius.Conv2PrtInten;
	}
	else
	{
		reclin = 0.;
	}

	fprintf( ioQQQ, 
		" ENERGY BUDGET:  Heat:%8.3f  Coolg:%8.3f  Error:%5.1f%%  Rec Lin:%8.3f  F-F  H%7.3f    P(rad/tot)max     ", 
	  powerl, 
	  thermal.totcol, 
	  a, 
	  reclin, 
	  thermal.FreeFreeTotHeat );
	PrintE82( ioQQQ, pressure.RadBetaMax );
	fprintf( ioQQQ, "\n");

	/* effective x-ray column density, from 0.5keV attenuation, no scat
	 * IPXRY is pointer to 73.5 Ryd */
	coleff = opac.TauAbsGeo[0][rt.ipxry-1] - rt.tauxry;
	coleff /= 2.14e-22;

	/* find t^2 from H II structure */
	gett2(&peimbt.t2hstr);

	/* find t^2 from OIII structure */
	gett2o3(&peimbt.t2o3str);

	fprintf( ioQQQ, "\n     Col(Heff):      ");
	PrintE93(ioQQQ, coleff);
	fprintf( ioQQQ, "  snd travl time  ");
	PrintE82(ioQQQ, timesc.sound);
	fprintf( ioQQQ, " sec  Te-low: ");
	PrintE82(ioQQQ, thermal.tlowst);
	fprintf( ioQQQ, "  Te-hi: ");
	PrintE82(ioQQQ, thermal.thist);

	/* this is the intensity of the UV continuum at the illuminated face, relative to the Habing value, as
	 * defined by Tielens & Hollenbach 1985 */
	fprintf( ioQQQ, "  G0TH85: ");
	PrintE82( ioQQQ, hmi.UV_Cont_rel2_Habing_TH85_face );
	/* this is the intensity of the UV continuum at the illuminated face, relative to the Habing value, as
	 * defined by Draine & Bertoldi 1985 */
	fprintf( ioQQQ, "  G0DB96:");
	PrintE82( ioQQQ, hmi.UV_Cont_rel2_Draine_DB96_face );

	fprintf( ioQQQ, "\n");

	fprintf( ioQQQ, "  Emiss Measure    n(e)n(p) dl       ");
	PrintE93(ioQQQ, colden.dlnenp);
	fprintf( ioQQQ, "  n(e)n(He+)dl         ");
	PrintE93(ioQQQ, colden.dlnenHep);
	fprintf( ioQQQ, "  En(e)n(He++) dl         ");
	PrintE93(ioQQQ, colden.dlnenHepp);
	fprintf( ioQQQ, "  ne nC+:");
	PrintE82(ioQQQ, colden.dlnenCp);
	fprintf( ioQQQ, "\n");

	/* following is wl where gas goes thick to bremsstrahlung, in cm */
	if( rfield.EnergyBremsThin > 0. )
	{
		bremtk = RYDLAM*1e-8/rfield.EnergyBremsThin;
	}
	else
	{
		bremtk = 1e30;
	}

	/* apparent helium abundance */
	fprintf( ioQQQ, " He/Ha:");
	PrintE82( ioQQQ, heabun);

	/* he/h relative to correct helium abundance */
	fprintf(ioQQQ, "  =%7.2f*true  Lthin:",hescal);

	/* wavelength were structure is optically thick to bremsstrahlung absorption */
	PrintE82( ioQQQ, bremtk);

	/* this is ratio of conv.nTotalIoniz, the number of times ConvBase 
	 * was called during the model, over the number of zones.
	 * This is a measure of the convergence of the model - 
	 * includes initial search so worse for short calculations.
	 * It is a measure of how hard the model was to converge */
	if( nzone > 0 )
	{
		/* >>chng 07 feb 23, subtract n calls to do initial solution
		 * so this is the number of calls needed to converge the zones.
		 * different is to discount careful approach to molecular solutions
		 * in one zone models */
		znit = (double)(conv.nTotalIoniz-conv.nTotalIoniz_start)/(double)(nzone);
	}
	else
	{
		znit = 0.;
	}
	/* >>chng 02 jan 09, format from 5.3f to 5.2f */
	fprintf(ioQQQ, "  itr/zn:%5.2f",znit);

	/* sort-of same thing for large H2 molecule - number is number of level pop solutions per zone */
	fprintf(ioQQQ, "  H2 itr/zn:%6.2f",H2_itrzn());

	/* say whether we used stored opacities (T) or derived them from scratch (F) */
	fprintf(ioQQQ, "  File Opacity: F" );

	/* log of total mass in grams if spherical, or gm/cm2 if plane parallel */
	{
		/* this is mass per unit inner area */
		double xmass = log10( SDIV(dense.xMassTotal) );
		xmass += (realnum)(1.0992099 + 2.*log10(radius.rinner));
		fprintf(ioQQQ,"  Mass tot  %.3f",
			xmass);
	}
	fprintf(ioQQQ,"\n");

	/* this block is a series of prints dealing with 21 cm quantities
	 * first this is the temperature derived from Lya - 21 cm optical depths
	 * get harmonic mean HI temperature, as needed for 21 cm spin temperature */
	if( cdTemp( "opti",0,&THI,"volume" ) )
	{
		fprintf(ioQQQ,"\n>>>>>>>>>>>>>>>>> PrtFinal, impossible problem getting 21cm opt.\n");
		TotalInsanity();
	}
	fprintf( ioQQQ, "   Temps(21 cm)   T(21cm/Ly a)  ");
	PrintE82(ioQQQ, THI );

	/* get harmonic mean HI gas kin temperature, as needed for 21 cm spin temperature 
	 * >>chng 06 jul 21, this was over volume but hazy said radius - change to radius,
	 * bug discovered by Eric Pellegrini & Jack Baldwin  */
	/*if( cdTemp( "21cm",0,&THI,"volume" ) )*/
	if( cdTemp( "21cm",0,&THI,"radius" ) )
	{
		fprintf(ioQQQ,"\n>>>>>>>>>>>>>>>>> PrtFinal, impossible problem getting 21cm temp.\n");
		TotalInsanity();
	}
	fprintf(ioQQQ, "        T(<nH/Tkin>)  ");
	PrintE82(ioQQQ, THI);

	/* get harmonic mean HI 21 21 spin temperature, as needed for 21 cm spin temperature 
	 * for this, always weighted by radius, volume would be ignored were it present */
	if( cdTemp( "spin",0,&THI,"radius" ) )
	{
		fprintf(ioQQQ,"\n>>>>>>>>>>>>>>>>> PrtFinal, impossible problem getting 21cm spin.\n");
		TotalInsanity();
	}
	fprintf(ioQQQ, "          T(<nH/Tspin>)    ");
	PrintE82(ioQQQ, THI);

	/* now convert this HI spin temperature into a brightness temperature */
	THI *= (1. - sexp( HFLines[0].Emis->TauIn ) );
	fprintf( ioQQQ, "          TB21cm:");
	PrintE82(ioQQQ, THI);
	fprintf( ioQQQ, "\n");

	fprintf( ioQQQ, "                   N(H0/Tspin)  ");
	PrintE82(ioQQQ, colden.H0_ov_Tspin );
	fprintf( ioQQQ, "        N(OH/Tkin)    ");
	PrintE82(ioQQQ, colden.OH_ov_Tspin );

	/* mean B weighted wrt 21 cm absorption */
	bot = cdB21cm();
	fprintf( ioQQQ, "          B(21cm)          ");
	PrintE82(ioQQQ, bot );

	/* end prints for 21 cm */
	fprintf(ioQQQ, "\n");

	/* timescale (sec here) for photoerosion of Fe-Ni */
	rate = timesc.TimeErode*2e-26;
	if( rate > SMALLFLOAT )
	{
		ferode = 1./rate;
	}
	else
	{
		ferode = 0.;
	}

	/* mean acceleration */
	if( wind.acldr > 0. )
	{
		wind.AccelAver /= wind.acldr;
	}
	else
	{
		wind.AccelAver = 0.;
	}

	/* DMEAN is mean density (gm per cc); mean temp is weighted by mass density */
	wmean = colden.wmas/SDIV(colden.TotMassColl);
	/* >>chng 02 aug 21, from radius.depth_x_fillfac to integral of radius times fillfac */
	dmean = colden.TotMassColl/SDIV(radius.depth_x_fillfac);
	tmean = colden.tmas/SDIV(colden.TotMassColl);
	/* mean mass per particle */
	wmean = colden.wmas/SDIV(colden.TotMassColl);

	fprintf( ioQQQ, "   <a>:");
	PrintE82(ioQQQ , wind.AccelAver);
	fprintf( ioQQQ, "  erdeFe");
	PrintE71(ioQQQ , ferode);
	fprintf( ioQQQ, "  Tcompt");
	PrintE82(ioQQQ , timesc.tcmptn);
	fprintf( ioQQQ, "  Tthr");
	PrintE82(ioQQQ , timesc.time_therm_long);
	fprintf( ioQQQ, "  <Tden>: ");
	PrintE82(ioQQQ , tmean);
	fprintf( ioQQQ, "  <dens>:");
	PrintE82(ioQQQ , dmean);
	fprintf( ioQQQ, "  <MolWgt>");
	PrintE82(ioQQQ , wmean);
	fprintf(ioQQQ,"\n");

	/* log of Jeans length and mass - this is mean over model */
	if( tmean > 0. )
	{
		rjeans = 7.79637 + (log10(tmean) - log10(dense.wmole) - log10(dense.xMassDensity*
			geometry.FillFac))/2.;
	}
	else
	{
		/* tmean undefined, set rjeans to large value so 0 printed below */
		rjeans = 40.f;
	}

	if( rjeans < 36. )
	{
		rjeans = (double)pow(10.,rjeans);
		/* AJMASS is Jeans mass in solar units */
		ajmass = 3.*log10(rjeans/2.) + log10(4.*PI/3.*dense.xMassDensity*
		  geometry.FillFac) - log10(SOLAR_MASS);
		if( ajmass < 37 )
		{
			ajmass = pow(10.,ajmass);
		}
		else
		{
			ajmass = 0.;
		}
	}
	else
	{
		rjeans = 0.;
		ajmass = 0.;
	}

	/* Jeans length and mass - this is smallest over model */
	ajmin = colden.ajmmin - log10(SOLAR_MASS);
	if( ajmin < 37 )
	{
		ajmin = pow(10.,ajmin);
	}
	else
	{
		ajmin = 0.;
	}

	/* estimate alpha (o-x) */
	if( rfield.anu[rfield.nflux-1] > 150. )
	{
		/* generate pointers to energies where alpha ox will be evaluated */
		ip2kev = ipoint(147.);
		ip2500 = ipoint(0.3648);

		/* now get fluxes there */
		bottom = rfield.flux[0][ip2500-1]*
		  rfield.anu[ip2500-1]/rfield.widflx[ip2500-1];

		top = rfield.flux[0][ip2kev-1]*
		  rfield.anu[ip2kev-1]/rfield.widflx[ip2kev-1];

		/* generate alpha ox if denominator is non-zero */
		if( bottom > 1e-30 && top > 1e-30 )
		{
			ratio = log10(top) - log10(bottom);
			if( ratio < 36. && ratio > -36. )
			{
				ratio = pow(10.,ratio);
			}
			else
			{
				ratio = 0.;
			}
		}

		else
		{
			ratio = 0.;
		}

		if( ratio > 0. )
		{
			// the separate variable freq_ratio is needed to work around a bug in icc 10.0
			double freq_ratio = rfield.anu[ip2kev-1]/rfield.anu[ip2500-1];
			alfox = log(ratio)/log(freq_ratio);
		}
		else
		{
			alfox = 0.;
		}
	}
	else
	{
		alfox = 0.;
	}

	if( colden.rjnmin < 37 )
	{
		colden.rjnmin = (realnum)pow((realnum)10.f,colden.rjnmin);
	}
	else
	{
		colden.rjnmin = FLT_MAX;
	}

	fprintf( ioQQQ, "     Mean Jeans  l(cm)");
	PrintE82(ioQQQ, rjeans);
	fprintf( ioQQQ, "  M(sun)");
	PrintE82(ioQQQ, ajmass);
	fprintf( ioQQQ, "  smallest:     len(cm):");
	PrintE82(ioQQQ, colden.rjnmin);
	fprintf( ioQQQ, "  M(sun):");
	PrintE82(ioQQQ, ajmin);
	fprintf( ioQQQ, "  a_ox tran:%6.2f\n" , alfox);

	fprintf( ioQQQ, " Rion:");
	double R_ion;
	if( rfield.lgUSphON )
		R_ion = rfield.rstrom;
	else
		R_ion = radius.Radius;
	PrintE93(ioQQQ, R_ion);
	fprintf( ioQQQ, "  Dist:");
	PrintE82(ioQQQ, radius.distance);
	fprintf( ioQQQ, "  Diam:");
	if( radius.distance > 0. )
		PrintE93(ioQQQ, 2.*R_ion*AS1RAD/radius.distance);
	else
		PrintE93(ioQQQ, 0.);
	fprintf( ioQQQ, "\n");

	if( prt.lgPrintTime )
	{
		/* print execution time by default */
		alfox = cdExecTime();
	}
	else
	{
		/* flag set false with no time command, so that different runs can
		 * compare exactly */
		alfox = 0.;
	}

	/* some details about the hydrogen and helium ions */
	fprintf( ioQQQ, 
		" Hatom level%3ld  HatomType:%4.4s  HInducImp %2c"
		"  He tot level:%3ld  He2 level:  %3ld  ExecTime%8.1f\n", 
		/* 06 aug 28, from numLevels_max to _local. */
	  iso.numLevels_local[ipH_LIKE][ipHYDROGEN], 
	  hydro.chHTopType, 
	  TorF(hydro.lgHInducImp), 
	  /* 06 aug 28, from numLevels_max to _local. */
	  iso.numLevels_local[ipHE_LIKE][ipHELIUM],
	  /* 06 aug 28, from numLevels_max to _local. */
	  iso.numLevels_local[ipH_LIKE][ipHELIUM] , 
	  alfox );

	/* now give an indication of the convergence error budget */
	fprintf( ioQQQ, 
		" ConvrgError(%%)  <eden>%7.3f  MaxEden%7.3f  <H-C>%7.2f  Max(H-C)%8.2f  <Press>%8.3f  MaxPrs er%7.3f\n", 
		conv.AverEdenError/nzone*100. , 
		conv.BigEdenError*100. ,
		conv.AverHeatCoolError/nzone*100. , 
		conv.BigHeatCoolError*100. ,
		conv.AverPressError/nzone*100. , 
		conv.BigPressError*100. );

	fprintf(ioQQQ,
		"  Continuity(%%)  chng Te%6.1f  elec den%6.1f  n(H2)%7.1f  n(CO)    %7.1f\n",
		phycon.BigJumpTe*100. ,
		phycon.BigJumpne*100. ,
		phycon.BigJumpH2*100. ,
		phycon.BigJumpCO*100. );

	/* print out some average quantities */
	fprintf( ioQQQ, "\n                                                        Averaged Quantities\n" );
	fprintf( ioQQQ, "             Te      Te(Ne)   Te(NeNp)  Te(NeHe+)Te(NeHe2+) Te(NeO+)  Te(NeO2+)"
		 "  Te(H2)     N(H)     Ne(O2+)   Ne(Np)\n" );
	static const char* weight[3] = { "Radius", "Area", "Volume" };
	int dmax = geometry.lgGeoPP ? 1 : 3;
	for( int dim=0; dim < dmax; ++dim )
	{
		fprintf( ioQQQ, " %6s:", weight[dim] );
		// <Te>/<1>
		PrintRatio( mean.TempMean[dim][0], mean.TempMean[dim][1] );
		// <Te*ne>/<ne>
		PrintRatio( mean.TempEdenMean[dim][0], mean.TempEdenMean[dim][1] );
		// <Te*ne*nion>/<ne*nion>
		PrintRatio( mean.TempIonEdenMean[dim][ipHYDROGEN][1][0], mean.TempIonEdenMean[dim][ipHYDROGEN][1][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipHELIUM][1][0], mean.TempIonEdenMean[dim][ipHELIUM][1][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipHELIUM][2][0], mean.TempIonEdenMean[dim][ipHELIUM][2][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipOXYGEN][1][0], mean.TempIonEdenMean[dim][ipOXYGEN][1][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipOXYGEN][2][0], mean.TempIonEdenMean[dim][ipOXYGEN][2][1] );
		// <Te*nH2>/<nH2>
		PrintRatio( mean.TempIonMean[dim][ipHYDROGEN][2][0], mean.TempIonMean[dim][ipHYDROGEN][2][1] );
		// <nH>/<1>
		PrintRatio( mean.xIonMean[dim][ipHYDROGEN][0][1], mean.TempMean[dim][1] );
		// <ne*nion>/<nion>
		PrintRatio( mean.TempIonEdenMean[dim][ipOXYGEN][2][1], mean.TempIonMean[dim][ipOXYGEN][2][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipHYDROGEN][1][1], mean.TempIonMean[dim][ipHYDROGEN][1][1] );
		fprintf( ioQQQ, "\n" );
	}

	/* print out Peimbert analysis, tsqden default 1e7, changed
	 * with set tsqden command */
	if( dense.gas_phase[ipHYDROGEN] < peimbt.tsqden )
	{
		fprintf(  ioQQQ, " \n" );

		/* temperature from the [OIII] 5007/4363 ratio */
		fprintf(  ioQQQ, " Peimbert T(OIIIr)");
		PrintE82( ioQQQ , peimbt.toiiir);

		/* temperature from predicted Balmer jump relative to Hbeta */
		fprintf(  ioQQQ, " T(Bac)");
		PrintE82( ioQQQ , peimbt.tbac);

		/* temperature predicted from optically thin Balmer jump rel to Hbeta */
		fprintf(  ioQQQ, " T(Hth)");
		PrintE82( ioQQQ , peimbt.tbcthn);

		/* t^2 predicted from the structure, weighted by H */
		fprintf(  ioQQQ, " t2(Hstrc)");
		fprintf( ioQQQ,PrintEfmt("%9.2e", peimbt.t2hstr));

		/* temperature from both [OIII] and the Balmer jump rel to Hbeta */
		fprintf(  ioQQQ, " T(O3-BAC)");
		PrintE82( ioQQQ , peimbt.tohyox);

		/* t2 from both [OIII] and the Balmer jump rel to Hbeta */
		fprintf(  ioQQQ, " t2(O3-BC)");
		fprintf( ioQQQ,PrintEfmt("%9.2e", peimbt.t2hyox));

		/* structural t2 from the O+2 predicted radial dependence */
		fprintf(  ioQQQ, " t2(O3str)");
		fprintf( ioQQQ,PrintEfmt("%9.2e", peimbt.t2o3str));

		fprintf(  ioQQQ, "\n");

		if( gv.lgDustOn() )
		{
			fprintf( ioQQQ, " Be careful: grains exist.  This spectrum was not corrected for reddening before analysis.\n" );
		}
	}

	/* print average (over radius) grain dust parameters if lgDustOn() is true,
	 * average quantities are incremented in radius_increment, zeroed in IterStart */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		char chQHeat;
		double AV , AB;
		double total_dust2gas = 0.;

		fprintf( ioQQQ, "\n Average Grain Properties (over radius):\n" );

		for( size_t i0=0; i0 < gv.bin.size(); i0 += 10 ) 
		{
			if( i0 > 0 )
				fprintf( ioQQQ, "\n" );

			/* this is upper limit to how many grain species we will print across line */
			size_t i1 = min(i0+10,gv.bin.size());

			fprintf( ioQQQ, "       " );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				chQHeat = (char)(gv.bin[nd]->lgEverQHeat ? '*' : ' ');
				fprintf( ioQQQ, "%-12.12s%c", gv.bin[nd]->chDstLab, chQHeat );
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "    nd:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ, "%7ld   ", (unsigned long)nd );
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <Tgr>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd]->avdust/radius.depth_x_fillfac));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <Vel>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd]->avdft/radius.depth_x_fillfac));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <Pot>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd]->avdpot/radius.depth_x_fillfac));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <D/G>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd]->avDGRatio/radius.depth_x_fillfac));
				/* add up total dust to gas mass ratio */
				total_dust2gas += gv.bin[nd]->avDGRatio/radius.depth_x_fillfac;
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf(ioQQQ," Dust to gas ratio (by mass):");
		fprintf(ioQQQ,PrintEfmt("%10.3e", total_dust2gas));

		/* total extinction (conv to mags) at V and B per hydrogen, this includes
		 * forward scattering as an extinction process, so is what would be measured
		 * for a star, but not for an extended source where forward scattering
		 * should be discounted */
		AV = rfield.extin_mag_V_point/SDIV(colden.colden[ipCOL_HTOT]);
		AB = rfield.extin_mag_B_point/SDIV(colden.colden[ipCOL_HTOT]);
		/* print A_V/N(Htot) for point and extended sources */
		fprintf(ioQQQ,", A(V)/N(H)(pnt):%.3e, (ext):%.3e", 
			AV,
			rfield.extin_mag_V_extended/SDIV(colden.colden[ipCOL_HTOT]) );

		/* ratio of total to selective extinction */
		fprintf(ioQQQ,", R:");

		/* gray grains have AB - AV == 0 */
		if( fabs(AB-AV)>SMALLFLOAT )
		{
			fprintf(ioQQQ,"%.3e", AV/(AB-AV) );
		}
		else
		{
			fprintf(ioQQQ,"%.3e", 0. );
		}
		fprintf(ioQQQ," AV(ext):%.3e (pnt):%.3e\n",
			rfield.extin_mag_V_extended, 
			rfield.extin_mag_V_point);
	}

	/* now release saved arrays */
	free( wavelength );
	free( ipSortLines );
	free( sclsav );
	free( lgPrt );
	free( chSTyp );

	/* now return the space claimed for the chSLab array */
	for( i=0; i < LineSave.nsum; ++i )
	{
		free( chSLab[i] );
	}
	free( chSLab );

	free( scaled );
	free( xLog_line_lumin );

	/* option to make short printout */
	if( !prt.lgPrtShort && called.lgTalk )
	{
		/* print log of optical depths, 
		 * calls prtmet if print line optical depths command entered */
		PrtAllTau();

		/* only print mean ionization and emergent continuum on last iteration */
		if( iterations.lgLastIt )
		{
			/* option to print column densities, set with print column density command */
			if( prt.lgPrintColumns )
				PrtColumns(ioQQQ);
			/* print mean ionization fractions for all elements */
			PrtMeanIon('i', false, ioQQQ);
			/* print mean ionization fractions for all elements with density weighting*/
			PrtMeanIon('i', true , ioQQQ);
			/* print mean temperature fractions for all elements */
			PrtMeanIon('t', false , ioQQQ);
			/* print mean temperature fractions for all elements with density weighting */
			PrtMeanIon('t', true , ioQQQ);
			/* print emergent continuum */
			PrtContinuum();
		}
	}

	/* print input title for model */
	fprintf( ioQQQ, "%s\n\n", input.chTitle );
	return;
}

/* routine to stuff comments into the stack of comments,
 * return is index to this comment */
long int StuffComment( const char * chComment )
{
	long int n , i;

	DEBUG_ENTRY( "StuffComment()" );

	/* only do this one time per core load */
	if( LineSave.ipass == 0 )
	{
		if( LineSave.nComment >= NHOLDCOMMENTS )
		{
			fprintf( ioQQQ, " Too many comments have been entered into line array; increase the value of NHOLDCOMMENTS.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* want this to finally be 33 char long to match format */
		static const int NWIDTH = 33;
		strcpy( LineSave.chHoldComments[LineSave.nComment], chComment );

		/* current length of this string */
		n = (long)strlen( LineSave.chHoldComments[LineSave.nComment] );
		for( i=0; i<NWIDTH-n-1-6; ++i )
		{
			strcat( LineSave.chHoldComments[LineSave.nComment], ".");
		}

		strcat( LineSave.chHoldComments[LineSave.nComment], "..");

		for( i=0; i<6; ++i )
		{
			strcat( LineSave.chHoldComments[LineSave.nComment], " ");
		}
	}

	++LineSave.nComment;
	return( LineSave.nComment-1 );
}

/*gett2 analyze computed structure to get structural t^2 */
STATIC void gett2(realnum *tsqr)
{
	long int i;

	double tmean;
	double a, 
	  as, 
	  b;

	DEBUG_ENTRY( "gett2()" );

	/* get T, t^2 */
	a = 0.;
	b = 0.;

	ASSERT( nzone < struc.nzlim);
	// struc.volstr[] is radius.dVeffAper saved as a function of nzone
	// we need this version of radius.dVeff since we want to compare to
	// emission lines that react to the APERTURE command
	for( i=0; i < nzone; i++ )
	{
		as = (double)(struc.volstr[i])*(double)(struc.hiistr[i])*
		  (double)(struc.ednstr[i]);
		a += (double)(struc.testr[i])*as;
		/* B is used twice below */
		b += as;
	}

	if( b <= 0. )
	{
		*tsqr = 0.;
	}
	else
	{
		/* following is H+ weighted mean temp over vol */
		tmean = a/b;
		a = 0.;

		ASSERT( nzone < struc.nzlim );
		for( i=0; i < nzone; i++ )
		{
			as = (double)(struc.volstr[i])*(double)(struc.hiistr[i])*
			  struc.ednstr[i];
			a += (POW2((double)(struc.testr[i]-tmean)))*as;
		}
		*tsqr = (realnum)(a/(b*(POW2(tmean))));
	}

	return;
}

/*gett2o3 analyze computed [OIII] spectrum to get t^2 */
STATIC void gett2o3(realnum *tsqr)
{
	long int i;
	double tmean;
	double a, 
	  as, 
	  b;

	DEBUG_ENTRY( "gett2o3()" );

	/* get T, t^2 */	a = 0.;
	b = 0.;
	ASSERT(nzone<struc.nzlim);
	// struc.volstr[] is radius.dVeffAper saved as a function of nzone
	// we need this version of radius.dVeff since we want to compare to
	// emission lines that react to the APERTURE command
	for( i=0; i < nzone; i++ )
	{
		as = (double)(struc.volstr[i])*(double)(struc.o3str[i])*
		  (double)(struc.ednstr[i]);
		a += (double)(struc.testr[i])*as;

		/* B is used twice below */
		b += as;
	}

	if( b <= 0. )
	{
		*tsqr = 0.;
	}

	else
	{
		/* following is H+ weighted mean temp over vol */
		tmean = a/b;
		a = 0.;
		ASSERT( nzone < struc.nzlim );
		for( i=0; i < nzone; i++ )
		{
			as = (double)(struc.volstr[i])*(double)(struc.o3str[i])*
			  struc.ednstr[i];
			a += (POW2((double)(struc.testr[i]-tmean)))*as;
		}
		*tsqr = (realnum)(a/(b*(POW2(tmean))));
	}
	return;
}
