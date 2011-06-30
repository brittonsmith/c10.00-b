/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdDrive main routine to call cloudy under all circumstances) */
/*cdReasonGeo write why the model stopped and type of geometry on io file */
/*cdWarnings write all warnings entered into comment stack */
/*cdEms obtain the local emissivity for a line, for the last computed zone */
/*cdColm get the column density for a constituent  */
/*cdLine get the predicted line intensity, also index for line in stack */
/*cdLine_ip get the predicted line intensity, using index for line in stack */
/*cdCautions print out all cautions after calculation, on arbitrary io unit */
/*cdTemp_last routine to query results and return temperature of last zone */
/*cdDepth_depth get depth structure from previous iteration */
/*cdTimescales returns thermal, recombination, and H2 formation timescales */
/*cdSurprises print out all surprises on arbitrary unit number */
/*cdNotes print stack of notes about current calculation */
/*cdPressure_last routine to query results and return pressure of last zone */
/*cdTalk tells the code whether to print results or be silent */
/*cdOutput redirect output to arbitrary file */
/*cdInput redirect input from arbitrary file */
/*cdRead routine to read in command lines when cloudy used as subroutine */
/*cdErrors produce summary of all warnings, cautions, etc, on arbitrary io unit */
/*cdIonFrac get ionization fractions for a constituent */
/*cdTemp get mean electron temperature for any element */
/*cdCooling_last routine to query results and return cooling of last zone */
/*cdHeating_last routine to query results and return heating of last zone */
/*cdEDEN_last return electron density of last zone */
/*cdNoExec call this routine to tell code not to actually execute */
/*cdDate - puts date of code into string */
/*cdVersion produces string that gives version number of the code */
/*cdExecTime any routine can call this, find the time [s] since cdInit was called */
/*cdPrintCommands( FILE *) prints all input commands into file */
/*cdDrive main routine to call cloudy under all circumstances) */
/*cdNwcns get the number of cautions and warnings, to tell if calculation is ok */
/*debugLine provides a debugging hook into the main line array  */
/*cdEms_ip obtain the local emissivity for a line with known index */
/*cdnZone gets number of zones */
/*cdClosePunchFiles closes all the save files that have been used */
/*cdB21cm - returns B as measured by 21 cm */
/*cdPrtWL print line wavelengths in Angstroms in the standard format */
#include "cddefines.h"
#include "trace.h"
#include "conv.h"
#include "init.h"
#include "lines.h"
#include "pressure.h"
#include "prt.h"
#include "colden.h"
#include "dense.h"
#include "radius.h"
#include "struc.h"
#include "mole.h"
#include "elementnames.h"
#include "mean.h"
#include "phycon.h"
#include "called.h"
#include "parse.h"
#include "input.h"
#include "taulines.h"
#include "version.h"
#include "thermal.h"
#include "optimize.h"
#include "grid.h"
#include "timesc.h"
#include "cloudy.h"
#include "warnings.h"
#include "lines_service.h"
#include "cddrive.h"

/*************************************************************************
 *
 * cdDrive - main routine to call cloudy - returns 0 if all ok, 1 if problems
 *
 ************************************************************************/

int cdDrive()
{
	bool lgBAD;

	DEBUG_ENTRY( "cdDrive()" );
	/*********************************
	 * main routine to call cloudy   *
	 *********************************/

	/* this is set false when code loaded, set true when cdInit called,
	 * this is check that cdInit was called first */
	if( !lgcdInitCalled )
	{
		printf(" cdInit was not called first - this must be the first call.\n");
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"cdDrive: lgOptimr=%1i lgVaryOn=%1i lgNoVary=%1i input.nSave:%li\n",
			optimize.lgOptimr , optimize.lgVaryOn , optimize.lgNoVary, input.nSave );
	}

	/* should we call cloudy, or the optimization driver? 
	 * possible to have VARY on line without OPTIMIZE being set 
	 * lgNoVary set with "no optimize" command */
	if( optimize.lgOptimr && optimize.lgVaryOn && !optimize.lgNoVary )
		optimize.lgVaryOn = true;
	else
		optimize.lgVaryOn = false;

	/* one time initialization of core load - returns if already called 
	 * called here rather than in cdInit since at this point we know if
	 * single sim or grid */
	InitCoreload();

	if( optimize.lgVaryOn )
	{
		/* this branch called if optimizing or grid calculation */
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: calling grid_do\n");
		/* option to drive optimizer set if OPTIMIZE was in input stream */
		lgBAD = grid_do();
	}
	else
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: calling cloudy\n");

		/* optimize did not occur, only compute one model, call cloudy */
		lgBAD = cloudy();
	}

	/* reset flag saying that cdInit has not been called */
	lgcdInitCalled = false;

	if( lgAbort || lgBAD )
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: returning failure during call. \n");
		/* lgAbort set true if something wrong, so return lgBAD false. */
		return 1;
	}
	else
	{
		/* everything is ok, return 0 */
		return 0;
	}
}

/*************************************************************************
*
* cdPrtWL write emission line wavelength
*
************************************************************************/

/*cdPrtWL print line wavelengths in Angstroms in the standard format - 
 * a wrapper for prt_wl which must be kept parallel with sprt_wl
 * both of those live in pdt.c */
void cdPrtWL( FILE *io , realnum wl )
{

	DEBUG_ENTRY( "cdPrtWL()" );

	prt_wl( io , wl );
	return;
}


/*************************************************************************
 *
 * cdReasonGeo write why the model stopped and type of geometry on io file 
 *
 ************************************************************************/


/*cdReasonGeo write why the model stopped and type of geometry on io file */
void cdReasonGeo(FILE * ioOUT)
{

	DEBUG_ENTRY( "cdReasonGeo()" );

	/*this is the reason the calculation stopped*/
	fprintf( ioOUT, "%s", warnings.chRgcln[0] );
	fprintf( ioOUT , "\n" );
	/* this is the geometry */
	fprintf( ioOUT, "%s", warnings.chRgcln[1] );
	fprintf( ioOUT , "\n" );
	return;
}


/*************************************************************************
 *
 * cdWarnings write all warnings entered into comment stack 
 *
 ************************************************************************/

/*cdWarnings write all warnings entered into comment stack */

void cdWarnings(FILE *ioPNT )
{
	long int i;

	DEBUG_ENTRY( "cdWarnings()" );

	if( warnings.nwarn > 0 )
		{
			for( i=0; i < warnings.nwarn; i++ )
			{
				/* these are all warnings that were entered in comment */
				fprintf( ioPNT, "%s", warnings.chWarnln[i] );
				fprintf( ioPNT, "\n" );
			}
		}

	return;
}


/*************************************************************************
 *
 * cdCautions print out all cautions after calculation, on arbitrary io unit 
 *
 ************************************************************************/

/*cdCautions print out all cautions after calculation, on arbitrary io unit */

void cdCautions(FILE * ioOUT)
{
	long int i;

	DEBUG_ENTRY( "cdCautions()" );

	if( warnings.ncaun > 0 )
	{
		for( i=0; i < warnings.ncaun; i++ )
		{
			fprintf( ioOUT, "%s", warnings.chCaunln[i] );
			fprintf( ioOUT, "\n" );
		}
	}
	return;
}

/*************************************************************************
 *
 * cdTimescales returns thermal, recombination, and H2 formation timescales 
 *
 ************************************************************************/

void cdTimescales(
	/* the thermal cooling timescale */
	double *TTherm , 
	/* the hydrogen recombination timescale */
	double *THRecom , 
	/* the H2 formation timescale */
	double *TH2 )
{

	DEBUG_ENTRY( "cdTimescales()" );

	/* these were all evaluated in AgeCheck, which was called by PrtComment */

	/* thermal or cooling timescale */
	*TTherm = timesc.time_therm_long;

	/* the hydrogen recombination timescale */
	*THRecom = timesc.time_Hrecom_long;

	/* longer of the the H2 formation and destruction timescales */
	*TH2 = MAX2( timesc.time_H2_Dest_longest , timesc.time_H2_Form_longest );
	return;
}


/*************************************************************************
 *
 * cdSurprises print out all surprises on arbitrary unit number 
 *
 ************************************************************************/

/*cdSurprises print out all surprises on arbitrary unit number */

void cdSurprises(FILE * ioOUT)
{
	long int i;

	DEBUG_ENTRY( "cdSurprises()" );

	if( warnings.nbang > 0 )
	{
		for( i=0; i < warnings.nbang; i++ )
		{
			fprintf( ioOUT, "%s", warnings.chBangln[i] );
			fprintf( ioOUT, "\n" );
		}
	}

	return;
}


/*************************************************************************
 *
 * cdNotes print stack of notes about current calculation 
 *
 ************************************************************************/

/*cdNotes print stack of notes about current calculation */

void cdNotes(FILE * ioOUT)
{
	long int i;

	DEBUG_ENTRY( "cdNotes()" );

	if( warnings.nnote > 0 )
	{
		for( i=0; i < warnings.nnote; i++ )
		{
			fprintf( ioOUT, "%s", warnings.chNoteln[i] );
			fprintf( ioOUT, "\n" );
		}
	}
	return;
}

/*************************************************************************
 *
 * cdCooling_last routine to query results and return cooling of last zone 
 *
 ************************************************************************/

/*cdCooling_last routine to query results and return cooling of last zone */
double cdCooling_last() /* return cooling for last zone */
{
	return thermal.ctot;
}

/*************************************************************************
 *
 * cdVersion - puts version number of code into string 
 * incoming string must have sufficient length and will become null
 * terminated string
 *
 ************************************************************************/

void cdVersion(char chString[])
{
	strcpy( chString , t_version::Inst().chVersion );
	return;
}

/*************************************************************************
 *
 * cdDate - puts date of code into string 
 * incoming string must have at least 8 char and will become null
 * terminated string
 *
 ************************************************************************/

/* cdDate - puts date of code into string  */
void cdDate(char chString[])
{
	strcpy( chString , t_version::Inst().chDate );
	return;
}


/*************************************************************************
 *
 * cdHeating_last routine to query results and return heating of last zone
 *
 ************************************************************************/

/*cdHeating_last routine to query results and return heating of last zone */

double cdHeating_last() /* return heating for last zone */
{
	return thermal.htot;
}


/*************************************************************************
 *
 * cdEDEN_last return electron density of last zone
 *
 ************************************************************************/

/*cdEDEN_last return electron density of last zone */

double cdEDEN_last() /* return electron density for last zone */
{
	return dense.eden;
}

/*************************************************************************
 *
 * cdNoExec call this routine to tell code not to actually execute
 *
 ************************************************************************/

/*cdNoExec call this routine to tell code not to actually execute */
#include "noexec.h"

void cdNoExec()
{

	DEBUG_ENTRY( "cdNoExec()" );

	/* option to read in all input quantiles and NOT execute the actual model
	 * only check on input parameters - set by calling cdNoExec */
	noexec.lgNoExec = true;

	return;
}


/*************************************************************************
 *
 * cdSetExecTime routine to initialize variables keeping track of time at start of calculation
 *
 ************************************************************************/

/* set to false initially, then to true when cdSetExecTime() is called to
 * initialize the clock */
static bool lgCalled=false;

/* >>chng 06 dec 19, RP rm "|| defined(__HP_aCC)" to run on hp */
#if defined(_MSC_VER) 
/* _MSC_VER branches assume getrusage isn't implemented by MS 
 * also is not implemented on our HP superdome */
struct timeval {
	long    tv_sec;         /* seconds */
	long    tv_usec;        /* microseconds */
};
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

/* will be used to save initial time */
static struct timeval before;

/* cdClock stores time since arbitrary datum in clock_dat           */
STATIC void cdClock(struct timeval *clock_dat)
{
	DEBUG_ENTRY( "cdClock()" );

/* >>chng 06 sep 2 rjrw: use long recurrence, fine grain UNIX clock *
 * -- maintain system dependences in a single routine               */
#if defined(_MSC_VER) || defined(__HP_aCC)
	clock_t clock_val;
	/* >>chng 05 dec 21, from above to below due to negative exec times */
	/*return (double)(clock() - before) / (double)CLOCKS_PER_SEC;*/
	clock_val = clock();
	clock_dat->tv_sec = clock_val/CLOCKS_PER_SEC;
	clock_dat->tv_usec = 1000000*(clock_val-(clock_dat->tv_sec*CLOCKS_PER_SEC))/CLOCKS_PER_SEC;
	/*>>chng 06 oct 05, this produced very large number, time typically 50% too long 
	clock_dat->tv_usec = 0;*/
#else
	struct rusage rusage;
	if(getrusage(RUSAGE_SELF,&rusage) != 0)
	{ 
		fprintf( ioQQQ, "DISASTER cdClock called getrusage with invalid arguments.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	clock_dat->tv_sec = rusage.ru_utime.tv_sec;
	clock_dat->tv_usec = rusage.ru_utime.tv_usec;
#endif

}
/* cdSetExecTime is called by cdInit when everything is initialized,
 * so that every time cdExecTime is called the elapsed time is returned */
void cdSetExecTime()
{
	cdClock(&before);
	lgCalled = true;
}
/* cdExecTime returns the elapsed time cpu time (sec) that has elapsed 
 * since cdInit called cdSetExecTime.*/
double cdExecTime()
{
	DEBUG_ENTRY( "cdExecTime()" );

	struct timeval clock_dat;
	/* check that we were properly initialized */
	if( lgCalled)
	{
		cdClock(&clock_dat);
		/*fprintf(ioQQQ,"\n DEBUG sec %.2e usec %.2e\n",
			(double)(clock_dat.tv_sec-before.tv_sec),
			1e-6*(double)(clock_dat.tv_usec-before.tv_usec));*/
		return (double)(clock_dat.tv_sec-before.tv_sec)+1e-6*(double)(clock_dat.tv_usec-before.tv_usec);
	}
	else
	{
		/* this is a big problem, we were asked for the elapsed time but
		 * the timer was not initialized by calling SetExecTime */
		fprintf( ioQQQ, "DISASTER cdExecTime was called before SetExecTime, impossible.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}

/*************************************************************************
 *
 * cdPrintCommands prints all input commands into file
 *
 ************************************************************************/

/* cdPrintCommands( FILE *)
 * prints all input commands into file */
void cdPrintCommands( FILE * ioOUT )
{
	long int i;
	fprintf( ioOUT, " Input commands follow:\n" );
	fprintf( ioOUT, "c ======================\n" );

	for( i=0; i <= input.nSave; i++ )
	{
		fprintf( ioOUT, "%s\n", input.chCardSav[i] );
	}
	fprintf( ioOUT, "c ======================\n" );
}


/*************************************************************************
 *
 * cdEms obtain the local emissivity for a line, for the last computed zone
 *
 ************************************************************************/
// if called without last parameter, return intrinsic emission
long int cdEmis(
	/* return value will be index of line within stack,
	 * negative of number of lines in the stack if the line could not be found*/
	/* 4 char null terminated string label */
	char *chLabel,
	/* line wavelength */
	realnum wavelength, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss  )
{
	DEBUG_ENTRY( "cdEms()" );
	long int i = cdEmis( chLabel , wavelength , emiss, false );
	return i;
}

/* \todo 2 This routine, cdLine, cdEmis_ip, and cdLine_ip should be consolidated somewhat.
 * in particular so that this routine has the same "closest line" reporting that cdLine has. */
long int cdEmis(
	/* return value will be index of line within stack,
	 * negative of number of lines in the stack if the line could not be found*/
	/* 4 char null terminated string label */
	char *chLabel,
	/* line wavelength */
	realnum wavelength, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// false to return intrinsic emission, true for emergent
	bool lgEmergent )
{
	/* use following to store local image of character strings */
	char chCARD[INPUT_LINE_LENGTH];
	char chCaps[5];
	long int j;
	realnum errorwave;

	DEBUG_ENTRY( "cdEms()" );

	/* returns the emissivity in the desired line
	 * used internally by the code, to do save lines structure */

	// can return intrinsic or emergent emission, default is intrinsic
	long int ipEmisType = 0;
	if( lgEmergent )
		ipEmisType = 1;

	strcpy( chCARD, chLabel );

	/* make sure chLabel is all caps */
	caps(chCARD);/* convert to caps */

	/* get the error associated with 4 significant figures */
	errorwave = WavlenErrorGet( wavelength );

	for( j=0; j < LineSave.nsum; j++ )
	{
		/* change chLabel to all caps to be like input chLineLabel */
		cap4(chCaps, LineSv[j].chALab);

		/* check wavelength and chLabel for a match */
		/*if( fabs(LineSv[j].wavelength- wavelength)/MAX2(DELTA,wavelength)<errorwave && 
			strcmp(chCaps,chCARD) == 0 ) */
		if( fabs(LineSv[j].wavelength-wavelength) < errorwave && strcmp(chCaps,chCARD) == 0 )
		{
			/* match, so set emiss to emissivity in line */
			*emiss = LineSv[j].emslin[ipEmisType];
			/* and announce success by returning line index within stack */
			return j;
		}
	}

	/* we fell through without finding the line - return false */
	return -LineSave.nsum;
}

/*************************************************************************
 *
 * cdEms_ip obtain the local emissivity for a line with known index
 *
 ************************************************************************/

void cdEmis_ip(
	/* index of the line in the stack */
	long int ipLine, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis_ip()" );

	/* returns the emissivity in a line - implements save lines structure 
	 * this version uses previously stored line index to save speed */
	ASSERT( ipLine >= 0 && ipLine < LineSave.nsum );
	*emiss = LineSv[ipLine].emslin[lgEmergent];
	return;
}

/*************************************************************************
 *
 * cdColm get the column density for a constituent - 0 return if ok, 1 if problems 
 *
 ************************************************************************/

int cdColm(
	/* return value is zero if all ok, 1 if errors happened */
	/* 4-char + eol string that is first 
	 * 4 char of element name as spelled by Cloudy, upper or lower case */
	const char *chLabel,	

	/* integer stage of ionization, 1 for atom, 2 for A+, etc, 
	 * 0 is special flag for CO, H2, OH, or excited state */
	long int ion,

	/* the column density derived by the code [cm-2] */
	double *theocl )	
{
	long int nelem;
	/* use following to store local image of character strings */
	char chLABEL_CAPS[20];

	DEBUG_ENTRY( "cdColm()" );

	/* check that chLabel[4] is null - supposed to be 4 char + end */
	if( chLabel[4] != '\0' || strlen(chLabel) != 4 )
	{
		fprintf( ioQQQ, " cdColm called with insane chLabel (between quotes) \"%s\", must be exactly 4 characters long.\n",
		  chLabel );
		return 1;
	}

	strcpy( chLABEL_CAPS, chLabel );
	/* convert element label to all caps */
	caps(chLABEL_CAPS);

	/* zero ionization stage has special meaning.  The quantities recognized are
	 * the molecules, "H2  ", "OH  ", "CO  ", etc 
	 * "CII*" excited state C+ */
	if( ion < 0 )
	{
		fprintf( ioQQQ, " cdColm called with insane ion, =%li\n", 
		  ion );
		return 1;
	}

	else if( ion == 0 )
	{
		struct molecule *sp;
		/* this case molecular column density */
		/* want the molecular hydrogen H2 column density */
		if( strcmp( chLABEL_CAPS , "H2  " )==0 )
		{
			*theocl = colden.colden[ipCOL_H2g] + colden.colden[ipCOL_H2s];
		}

		/* H- column density  */
		else if( strcmp( chLABEL_CAPS , "H-  " )==0 )
		{
			*theocl = colden.colden[ipCOL_HMIN];
		}

		/* H2+ column density ipCOL_H2p is 4 */
		else if( strcmp( chLABEL_CAPS , "H2+ " )==0 )
		{
			*theocl = colden.colden[ipCOL_H2p];
		}

		/* H3+ column density */
		else if( strcmp( chLABEL_CAPS , "H3+ " )==0 )
		{
			*theocl = colden.colden[ipCOL_H3p];
		}

		/* H2g - ground H2 column density */
		else if( strcmp( chLABEL_CAPS , "H2G " )==0 )
		{
			*theocl = colden.colden[ipCOL_H2g];
		}

		/* H2* - excited H2 - column density */
		else if( strcmp( chLABEL_CAPS , "H2* " )==0 )
		{
			*theocl = colden.colden[ipCOL_H2s];
		}

		/* HeH+ column density */
		else if( strcmp( chLABEL_CAPS , "HEH+" )==0 )
		{
			*theocl = colden.colden[ipCOL_HeHp];
		}

		/* molecular species column density */
		else if( (sp = findspecies(chLabel)) != &null_mole )
		{
			*theocl = sp->hevcol;
		}

		/* ===========================================================*/
		/* end special case molecular column densities, start special cases
		 * excited state column densities */
		/* CII^* column density, population of J=3/2 upper level of split ground term */
		else if( strcmp( chLABEL_CAPS , "CII*" )==0 )
		{
			*theocl = colden.C2Colden[1];
		}
		else if( strcmp( chLABEL_CAPS , "C11*" )==0 )
		{
			*theocl = colden.C1Colden[0];
		}
		else if( strcmp( chLABEL_CAPS , "C12*" )==0 )
		{
			*theocl = colden.C1Colden[1];
		}
		else if( strcmp( chLABEL_CAPS , "C13*" )==0 )
		{
			*theocl = colden.C1Colden[2];
		}
		else if( strcmp( chLABEL_CAPS , "O11*" )==0 )
		{
			*theocl = colden.O1Colden[0];
		}
		else if( strcmp( chLABEL_CAPS , "O12*" )==0 )
		{
			*theocl = colden.O1Colden[1];
		}
		else if( strcmp( chLABEL_CAPS , "O13*" )==0 )
		{
			*theocl = colden.O1Colden[2];
		}
		/* CIII excited states, upper level of 1909 */
		else if( strcmp( chLABEL_CAPS , "C30*" )==0 )
		{
			*theocl = colden.C3Colden[1];
		}
		else if( strcmp( chLABEL_CAPS , "C31*" )==0 )
		{
			*theocl = colden.C3Colden[2];
		}
		else if( strcmp( chLABEL_CAPS , "C32*" )==0 )
		{
			*theocl = colden.C3Colden[3];
		}
		else if( strcmp( chLABEL_CAPS , "SI2*" )==0 )
		{
			*theocl = colden.Si2Colden[1];
		}
		else if( strcmp( chLABEL_CAPS , "HE1*" )==0 )
		{
			*theocl = colden.He123S;
		}
		/* special option, "H2vJ" */
		else if( strncmp(chLABEL_CAPS , "H2" , 2 ) == 0 )
		{
			long int iVib = chLABEL_CAPS[2] - '0';
			long int iRot = chLABEL_CAPS[3] - '0';
			if( iVib<0 || iRot < 0 )
			{
				fprintf( ioQQQ, " cdColm called with insane v,J for H2=\"%4.4s\" caps=\"%4.4s\"\n", 
				  chLabel , chLABEL_CAPS );
				return 1;
			}
			*theocl = cdH2_colden( iVib , iRot );
		}

		/* clueless as to what was meant - bail */
		else
		{
			fprintf( ioQQQ, " cdColm called with unknown element chLabel, org=\"%4.4s \" 0 caps=\"%4.4s\" 0\n", 
			  chLabel , chLABEL_CAPS );
			return 1;
		}
	}
	else
	{
		/* this case, ionization stage of some element */
		/* find which element this is */
		nelem = 0;
		while( nelem < LIMELM && 
			strncmp(chLABEL_CAPS,elementnames.chElementNameShort[nelem],4) != 0 )
		{
			++nelem;
		}

		/* this is true if we have one of the first 30 elements in the label,
		 * nelem is on C scale */
		if( nelem < LIMELM )
		{

			/* sanity check - does this ionization stage exist?
			 * max2 is to pick up H2 as H 3 */
			if( ion > MAX2(3,nelem + 2) )
			{
				fprintf( ioQQQ, 
				  " cdColm asked to return ionization stage %ld for element %s but this is not physical.\n", 
				  ion, chLabel );
				return 1;
			}

			/* the column density, ion is on physics scale, but means are on C scale */
			*theocl = mean.xIonMean[0][nelem][ion-1][0];
			/*>>chng 06 jan 23, div by factor of two
			 * special case of H2 when being tricked as H 3 - this stores 2H_2 so that
			 * the fraction of H in H0 and H+ is correct - need to remove this extra
			 * factor of two here */
			if( nelem==ipHYDROGEN && ion==3 )
				*theocl /= 2.;
		}
		else
		{
			fprintf( ioQQQ, 
			  " cdColm did not understand this combination of ion %4ld and element %4.4s.\n", 
			  ion, chLabel );
			return 1;
		}
	}
	return 0;
}


/*************************************************************************
 *
 *cdErrors produce summary of all warnings, cautions, etc, on arbitrary io unit 
 *
 ************************************************************************/

void cdErrors(FILE *ioOUT)
{
	long int nc, 
	  nn, 
	  npe, 
	  ns, 
	  nte, 
	  nw ,
	  nIone,
	  nEdene;
	bool lgAbort_loc;

	DEBUG_ENTRY( "cdErrors()" );

	/* first check for number of warnings, cautions, etc */
	cdNwcns(&lgAbort_loc,&nw,&nc,&nn,&ns,&nte,&npe, &nIone, &nEdene );

	/* only say something is one of these problems is nonzero */
	if( nw || nc || nte || npe ||	nIone || nEdene || lgAbort_loc )
	{
		/* say the title of the model */
		fprintf( ioOUT, "%75.75s\n", input.chTitle );

		if( lgAbort_loc )
			fprintf(ioOUT," Calculation ended with abort!\n");

		/* print warnings on the io unit */
		if( nw != 0 )
		{
			cdWarnings(ioOUT);
		}

		/* print cautions on the io unit */
		if( nc != 0 )
		{
			cdCautions(ioOUT);
		}

		if( nte != 0 )
		{
			fprintf( ioOUT , "Te failures=%4ld\n", nte );
		}

		if( npe != 0 )
		{
			fprintf( ioOUT , "Pressure failures=%4ld\n", npe );
		}

		if( nIone != 0 )
		{
			fprintf( ioOUT , "Ionization failures=%4ld\n", nte );
		}

		if( nEdene != 0 )
		{
			fprintf( ioOUT , "Electron density failures=%4ld\n", npe );
		}
	}
	return;
}

/*************************************************************************
 *
 *cdDepth_depth get depth structure from previous iteration 
 *
 ************************************************************************/
void cdDepth_depth( double cdDepth[] )
{
	long int nz;

	DEBUG_ENTRY( "cdDepth_depth()" );

	for( nz = 0; nz<nzone; ++nz )
	{
		cdDepth[nz] = struc.depth[nz];
	}
	return;
}

/*************************************************************************
 *
 *cdPressure_depth routine to query results and return pressure of last iteration 
 *
 ************************************************************************/

/*
 * cdPressure_depth
 * This returns the pressure and its constituents for the last iteration. 
 * space was allocated in the calling routine for the vectors - 
 * before calling this, cdnZone should have been called to get the number of
 * zones, then space allocated for the arrays */
void cdPressure_depth(
	/* total pressure, all forms*/
	double TotalPressure[],			
	/* gas pressure */
	double GasPressure[],				
	/* line radiation pressure */
	double RadiationPressure[])
{
	long int nz;

	DEBUG_ENTRY( "cdPressure_depth()" );

	for( nz = 0; nz<nzone; ++nz )
	{
		TotalPressure[nz] = struc.pressure[nz];
		GasPressure[nz] = struc.GasPressure[nz];
		RadiationPressure[nz] = struc.pres_radiation_lines_curr[nz];
	}
	return;
}

/*************************************************************************
 *
 *cdPressure_last routine to query results and return pressure of last zone 
 *
 ************************************************************************/

void cdPressure_last(
		double *PresTotal,  /* total pressure, all forms, for the last computed zone*/
		double *PresGas,    /* gas pressure */
		double *PresRad)    /* line radiation pressure */
{

	DEBUG_ENTRY( "cdPressure_last()" );

	*PresGas = pressure.PresGasCurr;
	*PresRad = pressure.pres_radiation_lines_curr;
	*PresTotal = pressure.PresTotlCurr;
	return;
}

/*************************************************************************
 *
 *cdnZone gets number of zones
 *
 ************************************************************************/

/* returns number of zones */
long int cdnZone() 
{
	return nzone;
}

/*************************************************************************
 *
 *cdTemp_last routine to query results and return temperature of last zone 
 *
 ************************************************************************/


double cdTemp_last()
{
	return phycon.te;
}

/*************************************************************************
 *
 *cdIonFrac get ionization fractions for a constituent
 *
 ************************************************************************/

int cdIonFrac(
	/* four char string, null terminated, giving the element name */
	const char *chLabel, 
	/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
	 * 0 says special case */
	long int IonStage, 
	/* will be fractional ionization */
	double *fracin, 
	/* how to weight the average, must be "VOLUME" or "RADIUS" */
	const char *chWeight ,
	/* if true then weighting also has electron density, if false then only volume or radius */
	bool lgDensity ) 
	/* return value is 0 if element was found, non-zero if failed */
{
	long int ip, 
		ion, /* used as index within aaa vector*/
		nelem;
	realnum aaa[LIMELM + 1];
	/* use following to store local image of character strings */
	char chCARD[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "cdIonFrac()" );

	strcpy( chCARD, chWeight );
	/* make sure chWeight is all caps */
	caps(chCARD);/* convert to caps */

	int dim;
	if( strcmp(chCARD,"RADIUS") == 0 )
		dim = 0;
	else if( strcmp(chCARD,"AREA") == 0 )
		dim = 1;
	else if( strcmp(chCARD,"VOLUME") == 0 )
		dim = 2;
	else
	{
		fprintf( ioQQQ, " cdIonFrac: chWeight=%6.6s makes no sense to me, valid options are RADIUS, AREA, and VOLUME\n", 
		  chWeight );
		*fracin = 0.;
		return 1;
	}

	/* first ensure that chLabel is all caps */
	strcpy( chCARD, chLabel );
	/* make sure chLabel is all caps */
	caps(chCARD);/* convert to caps */

	if( IonStage==0 )
	{
		/* special case */
		if( strcmp(chCARD,"H2  " ) == 0 )
		{
			/* this will be request for H2, treated as third stage of hydrogen */
			nelem = 0;
			IonStage = 3;
		}
		else
		{
			fprintf( ioQQQ, " cdIonFrac: ion stage of zero and element %s makes no sense to me\n", 
			chCARD );
			*fracin = 0.;
			return 1;
		}
	}

	else 
	{
		/* find which element this is */
		nelem = 0;
		while( nelem < LIMELM && 
		       strcmp(chCARD,elementnames.chElementNameShort[nelem]) != 0 )
		{
			++nelem;
		}
	}

	/* if element not recognized and above loop falls through, nelem is LIMELM */
	if( nelem >= LIMELM )
	{
		fprintf( ioQQQ, " cdIonFrac called with unknown element chLabel, =%4.4s\n", 
		  chLabel );
		return 1;
	}

	/* sanity check - does this ionization stage exist? 
	 * IonStage is on spectroscopic scale and nelem is on C scale */

	/* ion will be used as pointer within the aaa array that contains average values,
	 * convert to C scale */
	ion = IonStage - 1;

	if( (ion > nelem+1 || ion < 0 ) && !(nelem==ipHYDROGEN&&ion==2))
	{
		fprintf( ioQQQ, " cdIonFrac asked to return ionization stage %ld for element %4.4s but this is not physical.\n", 
		  IonStage, chLabel );
		*fracin = -1.;
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest argument is LIMELM
	 * 'i' means ionization, not temperature */
	/* last argument says whether to include electron density */
	/* MeanIon uses c scale for nelem */
	mean.MeanIon('i',nelem,dim,&ip,aaa,lgDensity);
	*fracin = pow((realnum)10.f,aaa[ion]);

	/* we succeeded - say so */
	return 0;
}

/*************************************************************************
 *
 * debugLine provides a debugging hook into the main line array 
 *
 ************************************************************************/

 /* debugLine provides a debugging hook into the main line array 
  * loops over whole array and finds every line that matches length,
  * the wavelength, the argument to the function
  * put breakpoint inside if test 
  * return value is number of matches, also prints all matches*/
long debugLine( realnum wavelength )
{
	long j, kount;
	realnum errorwave;

	kount = 0;

	/* get the error associated with 4 significant figures */
	errorwave = WavlenErrorGet( wavelength );

	for( j=0; j < LineSave.nsum; j++ )
	{
		/* check wavelength and chLabel for a match */
		/* if( fabs(LineSv[j].wavelength- wavelength)/MAX2(DELTA,wavelength) < errorwave ) */
		if( fabs(LineSv[j].wavelength-wavelength) < errorwave )
		{
			printf("%s\n", LineSv[j].chALab);
			++kount;
		}
	}
	printf(" hits = %li\n", kount );
	return kount;
}

/*************************************************************************
 *
 *cdLine get the predicted line intensity, also index for line in stack 
 *
 ************************************************************************/

 /* returns array index for line in array stack if we found the line, 
  * return negative of total number of lines as debugging aid if line not found */
// emergent or intrinsic not specified - use intrinsic
long int cdLine(
		const char *chLabel, 
		/* wavelength of line in angstroms, not format printed by code */
		realnum wavelength, 
		/* linear intensity relative to normalization line*/
		double *relint, 
		/* log of luminosity or intensity of line */
		double *absint )
{
	DEBUG_ENTRY( "cdLine()" );
	long int i = cdLine( chLabel , wavelength , relint , absint, false );
	return i;
}
long int cdLine(
	  const char *chLabel, 
		/* wavelength of line in angstroms, not format printed by code */
	  realnum wavelength, 
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint ,
	  bool lgEmergent )
{
	char chCaps[5], 
	  chFind[5];
	long int ipobs, 
	  j, index_of_closest=LONG_MIN,
	  index_of_closest_w_correct_label=-1;
	realnum errorwave, smallest_error=BIGFLOAT,
		smallest_error_w_correct_label=BIGFLOAT;

	DEBUG_ENTRY( "cdLine()" );

	/* this is zero when cdLine called with cdNoExec called too */
	if( LineSave.nsum == 0 )
	{
		*relint =  0.;
		*absint = 0.;
		return 0;
	}
	ASSERT( LineSave.ipNormWavL >= 0);
	ASSERT( LineSave.nsum > 0);

	/* check that chLabel[4] is null - supposed to be 4 char + end */
	if( chLabel[4] != '\0' || strlen(chLabel) != 4 )
	{
		fprintf( ioQQQ, " cdLine called with insane chLabel (between quotes) \"%s\", must be exactly 4 characters long.\n",
		  chLabel );
		return 1;
	}

	/* change chLabel to all caps */
	cap4(chFind, chLabel);

	/* this replaces tabs with spaces. */
	/* \todo	2 keep this in, do it elsewhere, just warn and bail? */
	for( j=0; j<=3; j++ )
	{
		if( chFind[j] == '\t' )
		{
			chFind[j] = ' ';
		}
	}

	/* get the error associated with 4 significant figures */
	errorwave = WavlenErrorGet( wavelength );

	{
		/*@-redef@*/
		enum{DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && fabs(wavelength-1000.)<0.01 )
		{
			fprintf(ioQQQ,"cdDrive wl %.4e error %.3e\n",
				wavelength, errorwave );
		}
	}

	/* now go through entire line stack, do not do 0, which is unity integration  */
	for( j=1; j < LineSave.nsum; j++ )
	{
		/* find closest line to requested wavelength to
		 * report when we don't get exact match */
		realnum current_error;
		current_error = (realnum)fabs(LineSv[j].wavelength-wavelength);

		/* change chLabel to all caps to be like input chALab */
		cap4(chCaps, LineSv[j].chALab);

		if( current_error < smallest_error )
		{
			index_of_closest = j;
			smallest_error = current_error;
		}

		if( current_error < smallest_error_w_correct_label &&
			(strcmp(chCaps,chFind) == 0) )
		{
			index_of_closest_w_correct_label = j;
			smallest_error_w_correct_label = current_error;
		}

		/* check wavelength and chLabel for a match */
		/* DELTA since often want wavelength of zero */
		if( current_error <= errorwave || 
			fp_equal( wavelength + errorwave, LineSv[j].wavelength ) ||
			fp_equal( wavelength - errorwave, LineSv[j].wavelength ) )
		{
			/* now see if labels agree */
			if( strcmp(chCaps,chFind) == 0 )
			{
				/* match, so set pointer */
				ipobs = j;

				/* does the normalization line have a positive intensity*/
				if( LineSv[LineSave.ipNormWavL].SumLine[lgEmergent] > 0. )
				{
					*relint = LineSv[ipobs].SumLine[lgEmergent]/
						LineSv[LineSave.ipNormWavL].SumLine[lgEmergent]*
						LineSave.ScaleNormLine;
				}
				else
				{
					*relint = 0.;
				}

				/* return log of current line intensity if it is positive */
				if( LineSv[ipobs].SumLine[lgEmergent] > 0. )
				{
					*absint = log10(LineSv[ipobs].SumLine[lgEmergent]) + 
						radius.Conv2PrtInten;
				}
				else
				{
					/* line intensity is actually zero, return small number */
					*absint = -37.;
				}
				/* we found the line, return pointer to its location */
				return ipobs;
			}
		}
	}

	/* >>chng 05 dec 21, report closest line if we did not find exact match, note that
	 * exact match returns above, where we will return negative number of lines in stack */
	fprintf( ioQQQ," PROBLEM cdLine did not find line with label (between quotes) \"%4s\" and wavelength ", chFind );
	prt_wl(ioQQQ,wavelength);
	if( index_of_closest >= 0 )
	{
		fprintf( ioQQQ,".\n  The closest line (any label) was   \"%4s\"\t", 
			LineSv[index_of_closest].chALab );
		prt_wl(ioQQQ,LineSv[index_of_closest].wavelength );
		if( index_of_closest_w_correct_label >= 0 )
		{
			fprintf( ioQQQ,"\n  The closest with correct label was \"%4s\"\t", chFind );
			prt_wl(ioQQQ,LineSv[index_of_closest_w_correct_label].wavelength );
			fprintf( ioQQQ,"\n" );
		}
		else
			fprintf( ioQQQ,"\n  No line found with label \"%s\".\n", chFind );
		fprintf( ioQQQ,"\n" );
	}
	else
	{
		fprintf( ioQQQ,".\n PROBLEM No close line was found\n" );
		TotalInsanity();
	}

	*absint = 0.;
	*relint = 0.;

	/* if we fell down to here we did not find the line 
	 * return negative of total number of lines as debugging aid */
	return -LineSave.nsum;
}

/*cdLine_ip get the predicted line intensity, using index for line in stack */
void cdLine_ip(long int ipLine, 
			   /* linear intensity relative to normalization line*/
			   double *relint, 
			   /* log of luminosity or intensity of line */
			   double *absint )
{

	DEBUG_ENTRY( "cdLine_ip()" );
	cdLine_ip( ipLine , relint , absint , false );
}
void cdLine_ip(long int ipLine, 
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint ,
	  bool lgEmergent )
{

	DEBUG_ENTRY( "cdLine_ip()" );

	/* this is zero when cdLine called with cdNoExec called too */
	if( LineSave.nsum == 0 )
	{
		*relint =  0.;
		*absint = 0.;
		return;
	}
	ASSERT( LineSave.ipNormWavL >= 0);
	ASSERT( LineSave.nsum > 0);

	/* does the normalization line have a positive intensity*/
	if( LineSv[LineSave.ipNormWavL].SumLine[lgEmergent] > 0. )
		*relint = LineSv[ipLine].SumLine[lgEmergent]/
			LineSv[LineSave.ipNormWavL].SumLine[lgEmergent]*
			LineSave.ScaleNormLine;
	else
		*relint = 0.;

	/* return log of current line intensity if it is positive */
	if( LineSv[ipLine].SumLine[lgEmergent] > 0. )
		*absint = log10(LineSv[ipLine].SumLine[lgEmergent]) + 
			radius.Conv2PrtInten;
	else
		/* line intensity is actually zero, return small number */
		*absint = -37.;

	return;
}

/*************************************************************************
 *
 *cdNwcns get the number of cautions and warnings, to tell if calculation is ok
 *
 ************************************************************************/

void cdNwcns(
  /* abort status, this better be false */
  bool *lgAbort_ret ,
  /* the number of warnings, cautions, notes, and surprises */
  long int *NumberWarnings, 
  long int *NumberCautions, 
  long int *NumberNotes, 
  long int *NumberSurprises, 
  /* the number of temperature convergence failures */
  long int *NumberTempFailures, 
  /* the number of pressure convergence failures */
  long int *NumberPresFailures,
  /* the number of ionization convergence failures */
  long int *NumberIonFailures, 
  /* the number of electron density convergence failures */
  long int *NumberNeFailures )
{

	DEBUG_ENTRY( "cdNwcns()" );

	/* this would be set true if code ended with abort - very very bad */
	*lgAbort_ret = lgAbort;
	/* this sub is called after comment, to find the number of various comments */
	*NumberWarnings = warnings.nwarn;
	*NumberCautions = warnings.ncaun;
	*NumberNotes = warnings.nnote;
	*NumberSurprises = warnings.nbang;

	/* these are counters that were incremented during convergence failures */
	*NumberTempFailures = conv.nTeFail;
	*NumberPresFailures = conv.nPreFail;
	*NumberIonFailures = conv.nIonFail;
	*NumberNeFailures = conv.nNeFail;
	return;
}

void cdOutput( const char* filename, const char *mode )
{
	DEBUG_ENTRY( "cdOutput()" );

	if( ioQQQ != stdout && ioQQQ != NULL )
		fclose(ioQQQ);
	FILE *fp = stdout;
	if( *filename != '\0' )
		fp = open_data( filename, mode, AS_LOCAL_ONLY );
	ioQQQ = fp;
}

void cdInput( const char* filename, const char *mode )
{
	DEBUG_ENTRY( "cdInput()" );

	if( ioStdin != stdin && ioStdin != NULL )
		fclose(ioStdin);
	FILE *fp = stdin;
	if( *filename != '\0' )
		fp = open_data( filename, mode, AS_LOCAL_ONLY );
	ioStdin = fp;
}

/*************************************************************************
 *
 * cdTalk tells the code whether to print results or be silent
 *
 ************************************************************************/

void cdTalk(bool lgTOn)
{

	DEBUG_ENTRY( "cdTalk()" );

	/* MPI talking rules must be obeyed, otherwise loss of output may result */
	called.lgTalk = lgTOn && cpu.lgMPI_talk();
	/* has talk been forced off? */
	called.lgTalkForcedOff = !lgTOn;
	return;
}

/*cdB21cm - returns B as measured by 21 cm */
double cdB21cm()
{
	double ret;

	DEBUG_ENTRY( "cdB21cm()" );

	// return average over radius
	if( mean.TempB_HarMean[0][1] > SMALLFLOAT )
	{
		ret = mean.TempB_HarMean[0][0]/mean.TempB_HarMean[0][1];
	}
	else
	{
		ret = 0.;
	}
	return ret;
}

/*************************************************************************
 *
 * cdTemp get mean electron temperature for any element 
 *
 ************************************************************************/

/* This routine finds the mean electron temperature for any ionization stage 
 * It returns 0 if it could find the species, 1 if it could not find the species.
 * The first argument is a null terminated 4 char string that gives the element
 * name as spelled by Cloudy.  
 * The second argument is ion stage, 1 for atom, 2 for first ion, etc
 * This third argument will be returned as result,
 * Last parameter is either "RADIUS", "AREA", or "VOLUME" to give weighting 
 *
 * if the ion stage is zero then the element label will have a special meaning.
 * The string "21CM" is will return the 21 cm temperature.
 * The string "H2  " will return the temperature weighted wrt the H2 density
 * There are several other options as listed below */
/** \todo	2	this should have a last argument like cdIonFrac for whether or not weighting
 * is wrt electron density */

/* return value is o if things are ok and element was found, 
 * non-zero if element not found or there are problems */
int cdTemp(
	/* four char string, null terminated, giving the element name */
	const char *chLabel, 
	/* IonStage is ionization stage, 1 for atom, up to Z+1 where Z is atomic number,
	 * 0 means that chLabel is a special case */
	long int IonStage, 
	/* will be temperature */
	double *TeMean, 
	/* how to weight the average, must be "RADIUS", "AREA, or "VOLUME" */
	const char *chWeight ) 
{
	long int ip, 
		ion, /* used as pointer within aaa vector*/
		nelem;
	realnum aaa[LIMELM + 1];
	/* use following to store local image of character strings */
	char chWGHT[INPUT_LINE_LENGTH] , chELEM[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "cdTemp()" );

	/* make sure strings are all caps */
	strcpy( chWGHT, chWeight );
	caps(chWGHT);
	strcpy( chELEM, chLabel );
	caps(chELEM);

	/* now see which weighting */
	int dim;
	if( strcmp(chWGHT,"RADIUS") == 0 )
		dim = 0;
	else if( strcmp(chWGHT,"AREA") == 0 )
		dim = 1;
	else if( strcmp(chWGHT,"VOLUME") == 0 )
		dim = 2;
	else
	{
		fprintf( ioQQQ, " cdTemp: chWeight=%6.6s makes no sense to me, the options are RADIUS, AREA, and VOLUME.\n", 
		  chWeight );
		*TeMean = 0.;
		return 1;
	}

	if( IonStage == 0 )
	{
		/* return atomic hydrogen weighted harmonic mean of gas kinetic temperature */
		if( strcmp(chELEM,"21CM") == 0 )
		{
			if( mean.TempHarMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempHarMean[dim][0]/mean.TempHarMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return atomic hydrogen weighted harmonic mean 21 cm spin temperature,
		 * using actual level populations with 1s of H0 */
		else if( strcmp(chELEM,"SPIN") == 0 )
		{
			*TeMean = mean.TempH_21cmSpinMean[dim][0] / SDIV(mean.TempH_21cmSpinMean[dim][1]);
		}
		/* return temperature deduced from ratio of 21 cm and Lya optical depths */
		else if( strcmp(chELEM,"OPTI") == 0 )
		{
			*TeMean = 
				3.84e-7 * Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->TauCon /
				SDIV( HFLines[0].Emis->TauCon );
		}
		/* mean temp weighted by H2 density */
		else if( strcmp(chELEM,"H2  ") == 0 )
		{
			if( mean.TempIonMean[dim][ipHYDROGEN][2][1] > SMALLFLOAT )
				*TeMean = mean.TempIonMean[dim][ipHYDROGEN][2][0] /
					mean.TempIonMean[dim][ipHYDROGEN][2][1];
					
			else
				*TeMean = 0.;
		}
		/* this is temperature weighted by eden */
		else if( strcmp(chELEM,"TENE") == 0 )
		{
			if( mean.TempEdenMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempEdenMean[dim][0]/mean.TempEdenMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* four spaces mean to return simple mean of Te */
		else if( strcmp(chELEM,"    ") == 0 )
		{
			if( mean.TempMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempMean[dim][0]/mean.TempMean[dim][1];
			else
				*TeMean = 0.;
		}
		else
		{
			fprintf( ioQQQ, " cdTemp called with ion=0 and unknown quantity, =%4.4s\n", 
			  chLabel );
			return 1;
		}

		/* say things are ok */
		return 0;
	}

	/* find which element this is */
	nelem = 0;
	while( nelem < LIMELM && 
	       strcmp(chELEM,elementnames.chElementNameShort[nelem]) != 0 )
	{
		++nelem;
	}

	/* if element not recognized and above loop falls through, nelem is LIMELM */
	if( nelem >= LIMELM )
	{
		fprintf( ioQQQ, " cdTemp called with unknown element chLabel, =%4.4s\n", 
		  chLabel );
		return 1;
	}

	/* sanity check - does this ionization stage exist? 
	 * IonStage on spectroscopic scale, nelem on c */

	/* ion will be used as pointer within the aaa array that contains average values,
	 * done this way to prevent lint from false problem in access to aaa array */
	ion = IonStage - 1;

	if( ion > nelem+1 || ion < 0 )
	{
		fprintf( ioQQQ, " cdTemp asked to return ionization stage %ld for element %4.4s but this is not physical.\n", 
		  IonStage, chLabel );
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest arg is LIMELM */
	/* MeanIon uses C scale for nelem */
	mean.MeanIon('t', nelem,dim,&ip,aaa,false);
	*TeMean = pow((realnum)10.f,aaa[ion]);
	return 0;
}

/*************************************************************************
 *
 * cdRead routine to read in command lines 
 * called by user when cloudy used as subroutine 
 * called by maincl when used as a routine
 *
 ************************************************************************/

/* returns the number of lines available in command stack 
 * this is limit to how many more commands can be read */
int cdRead(
	/* the string containing the command */
	const char *chInputLine )	
{
	char *chEOL , /* will be used to search for end of line symbols */
		chCARD[INPUT_LINE_LENGTH],
		chLocal[INPUT_LINE_LENGTH];
	bool lgComment;

	DEBUG_ENTRY( "cdRead()" );

	/* this is set false when code loaded, set true when cdInit called,
	 * this is check that cdInit was called first */
	if( !lgcdInitCalled )
	{
		printf(" cdInit was not called first - this must be the first call.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* totally ignore any card starting with a #, *, space, //, or % 
	 * but want to include special "c " type of comment 
	 * >>chng 06 sep 04 use routine to check for comments */
	if( ( lgInputComment( chInputLine ) ||
	      /* these two are end-of-input-stream sentinels */
	      chInputLine[0] == '\n' || chInputLine[0] == ' ' ) &&
	    /* option to allow "C" lines through */
	    ! ( chInputLine[0] == 'C' || chInputLine[0] == 'c' ) )
	{
		/* return value is number of lines that can still be stuffed in */
		return NKRD - input.nSave;
	}

	/***************************************************************************
	* validate a location to store this line image, then store the version     *
	* that has been truncated from special end of line characters              *
	* stored image will have < INPUT_LINE_LENGTH valid characters              *
	***************************************************************************/

	/* this will now point to the next free slot in the line image save array 
	 * this is where we will stuff this line image */
	++input.nSave;

	if( input.nSave >= NKRD )
	{
		/* too many input commands were entered - bail */
		fprintf( ioQQQ, 
			" Too many line images entered to cdRead.  The limit is %d\n", 
		  NKRD );
		fprintf( ioQQQ, 
			" The limit to the number of allowed input lines is %d.  This limit was exceeded.  Sorry.\n", 
		  NKRD );
		fprintf( ioQQQ, 
			" This limit is set by the variable NKRD which appears in input.h \n" );
		fprintf( ioQQQ, 
			" Increase it everywhere it appears.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	strncpy( chLocal, chInputLine, INPUT_LINE_LENGTH );
	// strncpy will pad chLocal with null bytes if chInputLine is shorter than
	// INPUT_LINE_LENGTH characters, so this indicates an overlong input string
	if( chLocal[INPUT_LINE_LENGTH-1] != '\0' )
	{
		chLocal[INPUT_LINE_LENGTH-1] = '\0';
		fprintf(ioQQQ," PROBLEM cdRead, while parsing the following input line:\n %s\n",
			chInputLine);
		fprintf(ioQQQ," found that the line is longer than %i characters, the longest possible line.\n",
			INPUT_LINE_LENGTH-1);
		fprintf(ioQQQ," Please make the command line no longer than this limit.\n");
	}

	/* now kill any part of line image after special end of line character,
	 * this stops info user wants ignored from entering after here */
	if( (chEOL = strchr_s(chLocal, '\n' ) ) != NULL )
	{
		*chEOL = '\0';
	}
	if( (chEOL = strchr_s(chLocal, '%' ) ) != NULL )
	{
		*chEOL = '\0';
	}
	/* >>chng 02 apr 10, add this char */
	if( (chEOL = strchr_s(chLocal, '#' ) ) != NULL )
	{
		*chEOL = '\0';
	}
	if( (chEOL = strchr_s(chLocal, ';' ) ) != NULL )
	{
		*chEOL = '\0';
	}
	if( (chEOL = strstr_s(chLocal, "//" ) ) != NULL )
	{
		*chEOL = '\0';
	}

	/* now do it again, since we now want to make sure that there is a trailing space
	 * if the line is shorter than 80 char, test on null is to keep lint happy */
	if( (chEOL = strchr_s( chLocal, '\0' )) == NULL )
		TotalInsanity();

	/* pad with two spaces if short enough,
	 * if not short enough for this to be done, then up to user to create correct input */
	if( chEOL-chLocal < INPUT_LINE_LENGTH-2 )
	{
		strcat( chLocal, "  " );
	}

	/* save string in master array for later use in readr */
	strcpy( input.chCardSav[input.nSave], chLocal );

	/* copy string to chCARD, then convert this to caps to check for keywords*/
	strcpy( chCARD, chLocal );
	caps(chCARD);/* convert to caps */

	// is this a comment?  if so do not check for keywords
	lgComment = false;
	if( strncmp(chCARD , "C ", 2 )==0 )
		lgComment = true;

	/* check whether this is a trace command, turn on printout if so */
	if( strncmp(chCARD,"TRACE",5) == 0 )
		trace.lgTrace = true;

	/* print input lines if trace specified */
	if( trace.lgTrace )
		fprintf( ioQQQ,"cdRead=%s=\n",input.chCardSav[input.nSave] );

	{
		// remove string between double quotes
		char chFilename[INPUT_LINE_LENGTH],
			chTemp[INPUT_LINE_LENGTH];
		// this has to have copy of line for GetQuote to function
		strcpy( chTemp , input.chCardSav[input.nSave] );
		GetQuote( chFilename , chCARD , chTemp , false );
	}

	/* now check whether VARY is specified */
	if( !lgComment && nMatch("VARY",chCARD) )
		/* a optimize flag was on the line image */
		optimize.lgVaryOn = true;

	/* now check whether line is "no vary" command */
	if( strncmp(chCARD,"NO VARY",7) == 0 )
		optimize.lgNoVary = true;

	/* now check whether line is "grid" command */
	if( strncmp(chCARD,"GRID",4) == 0 )
	{
		grid.lgGrid = true;
		++grid.nGridCommands;
	}

	/* NO BUFFERING turn off buffered io for standard output, 
	 * used to get complete output when code crashes */
	if( strncmp(chCARD,"NO BUFF",7) == 0 )
	{
		/* if output has already been redirected (e.g. using cdOutput()) then
		 * ignore this command, a warning will be printed in ParseDont() */
		/* stdout may be a preprocessor macro, so lets be really careful here */
		FILE *test = stdout;
		if( ioQQQ == test )
		{
			fprintf( ioQQQ, " cdRead found NO BUFFERING command, redirecting output to stderr now.\n" );
			/* make sure output is not lost */
			fflush( ioQQQ );
			/* stderr is always unbuffered */
			ioQQQ = stderr;
			/* will be used to generate comment at end */
			input.lgSetNoBuffering = false;
		}
	}

	/* use grid command as substitute for optimize command */
	if( strncmp(chCARD,"OPTI",4) == 0 || strncmp(chCARD,"GRID",4) == 0 )
		/* optimize command read in */
		optimize.lgOptimr = true;

	return NKRD - input.nSave;
}

/* wrapper to close all save files */
void cdClosePunchFiles()
{

	DEBUG_ENTRY( "cdClosePunchFiles()" );

	CloseSaveFiles( true );
	return;
}
