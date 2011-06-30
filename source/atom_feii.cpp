/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*FeII_Colden maintain H2 column densities within X */
/*FeIILevelPops  main FeII routine, called by CoolIron to evaluate iron cooling */
/*FeIICreate read in needed data from file 
 * convert form of FeII data, as read in from file within routine FeIICreate
 * into physical form.  called by atmdat_readin  */
/*FeIIPunPop - save level populations */
/*AssertFeIIDep called by assert FeII depart coef command */
/*FeIIPrint print FeII information */
/*FeIICollRatesBoltzmann evaluate collision strenths, 
 * both interpolating on r-mat and creating g-bar 
 * find Boltzmann factors, evaluate collisional rate coefficients */
/*FeIIPrint print output from large FeII atom, called by prtzone */
/*FeIISumBand sum up large FeII emission over certain bands, called in lineset4 */
/*FeII_RT_TauInc called once per zone in RT_tau_inc to increment large FeII atom line optical depths */
/*FeII_RT_tau_reset reset optical depths for large FeII atom, called by update after each iteration */
/*FeIIPoint called by ContCreatePointers to create pointers for lines in large FeII atom */
/*FeIIAccel called by rt_line_driving to compute radiative acceleration due to FeII lines */
/*FeIIAddLines save accumulated FeII intensities called by lineset4 */
/*FeIISaveLines save FeII lines at end of calculation, if save verner set, called by punch_do*/
/*FeIIPunchOpticalDepth save FeII line optical depth at end of calculation, called by punch_do*/
/*FeIIPunchLevels save FeII levels and energies */
/*FeII_LineZero zero out storage for large FeII atom, called by tauout */
/*FeIIIntenZero zero out intensity of FeII atom */
/*FeIIZero initialize some variables, called by zero one time before commands parsed */
/*FeIIReset reset some variables, called by zero */
/*FeIIPunData save line data */ 
/*FeIIPunDepart save some departure coef for large atom, 
 * set with save FeII departure command*/
/*FeIIPun1Depart send the departure coef for physical level nPUN to unit ioPUN */
/*FeIIContCreate create FeII continuum bins to add lines into ncell cells between wavelengths lambda low and high,
 * returns number of cells used */
/*FeIIBandsCreate returns number of FeII bands */
/*FeII_RT_Out do outward rates for FeII, called by RT_diffuse */
/*FeII_OTS do ots rates for FeII, called by RT_OTS */
/*FeII_RT_Make called by RT_line_all, does large FeII atom radiative transfer */
/*FeIILyaPump find rate of Lya excitation of the FeII atom */
/*ParseAtomFeII parse the atom FeII command */
/*FeIIPunchLineStuff include FeII lines in punched optical depths, etc, called from SaveLineStuff */
#include "cddefines.h"
#include "cddrive.h"
#include "thermal.h"
#include "physconst.h"
#include "doppvel.h"
#include "taulines.h"
#include "dense.h"
#include "rfield.h"
#include "radius.h"
#include "lines_service.h"
#include "ipoint.h"
#include "thirdparty.h"
#include "hydrogenic.h"
#include "lines.h"
#include "rt.h"
#include "trace.h"
#include "save.h"
#include "phycon.h"
#include "atomfeii.h"
#include "iso.h"
#include "pressure.h"
#include "parser.h"

/* add FeII lines into ncell cells between wavelengths lambda low and high */
STATIC void FeIIContCreate(double xLamLow , double xLamHigh , long int ncell );

/* read in the FeII Bands file, and sets nFeIIBands, the number of bands,
 * if argument is "" then use default file with bands, 
 * if filename within "" then use it instead,
 * return value is 0 if success, 1 if failed */
STATIC int FeIIBandsCreate( const char chFile[] );

/* this will be the smallest collision strength we will permit with the gbar
 * collision strengths, or for the data that have zeroes */
/* >>chng 99 jul 24, this was 1e-9 in the old fortran version and in baldwin et al. 96,
 * and has been changed several times, and affects results.  this is the smallest
 * non-zero collision strength in the r-matrix calculations */
realnum CS2SMALL = (realnum)1e-5;
/* routines used only within this file */

/*FeIICollRatesBoltzmann evaluate collision strenths, 
 * both interpolating on r-mat and creating g-bar 
 * find Boltzmann factors, evaluate collisional rate coefficients */
STATIC void FeIICollRatesBoltzmann(void);
/* find rate of Lya excitation of the FeII atom */
STATIC void FeIILyaPump(void);

/*extern realnum Fe2LevN[NFE2LEVN][NFE2LEVN][NTA];*/
/*extern realnum Fe2LevN[ipHi][ipLo].t[NTA];*/
/*realnum ***Fe2LevN;*/
/* >>chng 06 mar 01, boost to global namespace */
/*transition **Fe2LevN;*/
/* flag for the collision strength */
int **ncs1;

/* all following variables have file scope
#define	NFEIILINES	68635 */

/* this is size of nnPradDat array */
#define NPRADDAT 159

/* band wavelength, lower and upper bounds, in vacuum Angstroms */
/* FeII_Bands[n][3], where n is the number of bands in fe2Bands.dat
 * these bands are defined in fe2Bands.dat and read in at startup
 * of calculation */
realnum **FeII_Bands; 

// FeII_Cont[n][0] gives lower and upper bounds continuum wavelengths
// in vacuum Angstroms,
// [1] and [2] are inward and outward integrated intensity
// n is the number of bands, one less that number of cells needed
// these bands are defined in FeIIContCreate
realnum **FeII_Cont; 

/* this is the number of bands read in from FeII_bands.ini */
long int nFeIIBands;

/* number of bands in continuum array */
long int nFeIIConBins;

/* the dim of this vector this needs one extra since there are 20 numbers per line,
 * highest not ever used for anything */
/*long int nnPradDat[NPRADDAT+1];*/
static long int *nnPradDat;

/* malloced in feiidata */
/* realnum sPradDat[NPRADDAT][NPRADDAT][8];*/
/* realnum sPradDat[ipHi][ipLo].t[8];*/
static realnum ***sPradDat;

/* array used to integrate FeII line intensities over model,
 * Fe2SavN[upper][lower],
 *static realnum Fe2SavN[NFE2LEVN][NFE2LEVN];*/
static double **Fe2SavN;

/* save effective transition rates */
static double **Fe2A;

/* induced transition rates */
static double **Fe2LPump, **Fe2CPump;

/* energies read in from fe2energies.dat data file */
realnum *Fe2Energies;

/* collision rates */
static realnum **Fe2Coll;

/* Fe2DepCoef[NFE2LEVN];
realnum cli[NFEIILINES], 
  cfe[NFEIILINES];*/
static double 
	/* departure coefficients */
	*Fe2DepCoef , 
	/* level populations */
	*Fe2LevelPop ,
	/* column densities */
	*Fe2ColDen ,
	/* this will become array of Boltzmann factors */
	*FeIIBoltzmann;
	/*FeIIBoltzmann[NFE2LEVN] ,*/

static double EnerLyaProf1, 
  EnerLyaProf4, 
  PhotOccNumLyaCenter;
static double 
		/* the yVector - will become level populations after matrix inversion */
		*yVector,
	  /* this is used to call matrix routines */
	  /*xMatrix[NFE2LEVN][NFE2LEVN] ,*/
	  **xMatrix , 
	  /* this will become the very large 1-D array that 
	   * is passed to the matrix inversion routine*/
	  *amat;


/*FeII_Colden maintain H2 column densities within X */
void FeII_Colden( const char *chLabel )
{
	long int n;

	DEBUG_ENTRY( "FeII_Colden()" );

	if( strcmp(chLabel,"ZERO") == 0 )
	{
		/* zero out column density */
		for( n=0; n < FeII.nFeIILevel_malloc; ++n )
		{
			/* space for the rotation quantum number */
			Fe2ColDen[n] = 0.;
		}
	}

	else if( strcmp(chLabel,"ADD ") == 0 )
	{
		/*  add together column densities */
		for( n=0; n < FeII.nFeIILevel_local; ++n )
		{
			/* state-specific FeII column density */
			Fe2ColDen[n] += Fe2LevelPop[n]*radius.drad_x_fillfac;
		}
	}

	/* check for the print option, which we can't do, else  we have a problem */
	else if( strcmp(chLabel,"PRIN") != 0 )
	{
		fprintf( ioQQQ, " FeII_Colden does not understand the label %s\n", 
		  chLabel );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}

/*
 *=====================================================================
 */
/* FeIICreate read in FeII data from files on disk.  called by atmdat_readin 
 * but only if FeII. lgFeIION is true, set with atom FeII verner command */
void FeIICreate(void)
{
	FILE *ioDATA;

	char chLine[FILENAME_PATH_LENGTH_2];

	long int i, 
	  ipHi ,
	  ipLo,
	  lo,
	  ihi,
	  k, 
	  m1, 
	  m2, 
	  m3;

	DEBUG_ENTRY( "FeIICreate()" );

	if( lgFeIIMalloc )
	{
		/* we have already been called one time, just bail out */

		return;
	}

	/* now set flag so never done again - this is set false when initi
	 * when this is true it is no longer possible to change the number of levels
	 * in the model atom with the atom FeII levels command */
	lgFeIIMalloc = true;

	/* remember how many levels this was, so that in future calculations
	 * we can reset the atom to full value */
	FeII.nFeIILevel_malloc = FeII.nFeIILevel_local;

	/* set up array to save FeII transition probabilities */
	Fe2A = (double **)MALLOC(sizeof(double *)*(unsigned long)FeII.nFeIILevel_malloc );

	/* second dimension, lower level, for line save array */
	for( ipHi=0; ipHi < FeII.nFeIILevel_malloc; ++ipHi )
	{
		Fe2A[ipHi]=(double*)MALLOC(sizeof(double )*(unsigned long)FeII.nFeIILevel_malloc);
	}

	/* set up array to save FeII pumping rates */
	Fe2CPump = (double **)MALLOC(sizeof(double *)*(unsigned long)FeII.nFeIILevel_malloc );

	/* set up array to save FeII pumping rates */
	Fe2LPump = (double **)MALLOC(sizeof(double *)*(unsigned long)FeII.nFeIILevel_malloc );

	/* second dimension, lower level, for line save array */
	for( ipHi=0; ipHi < FeII.nFeIILevel_malloc; ++ipHi )
	{
		Fe2CPump[ipHi]=(double*)MALLOC(sizeof(double )*(unsigned long)FeII.nFeIILevel_malloc);

		Fe2LPump[ipHi]=(double*)MALLOC(sizeof(double )*(unsigned long)FeII.nFeIILevel_malloc);
	}

	/* set up array to save FeII collision rates */
	Fe2Energies = (realnum *)MALLOC(sizeof(realnum)*(unsigned long)FeII.nFeIILevel_malloc );

	/* set up array to save FeII collision rates */
	Fe2Coll = (realnum **)MALLOC(sizeof(realnum *)*(unsigned long)FeII.nFeIILevel_malloc );

	/* second dimension, lower level, for line save array */
	for( ipHi=0; ipHi < FeII.nFeIILevel_malloc; ++ipHi )
	{
		Fe2Coll[ipHi]=(realnum*)MALLOC(sizeof(realnum )*(unsigned long)FeII.nFeIILevel_malloc);
	}

	/* MALLOC space for the 2D matrix array */
	xMatrix = (double **)MALLOC(sizeof(double *)*(unsigned long)FeII.nFeIILevel_malloc );

	/* now do the second dimension */
	for( i=0; i<FeII.nFeIILevel_malloc; ++i )
	{
		xMatrix[i] = (double *)MALLOC(sizeof(double)*(unsigned long)FeII.nFeIILevel_malloc );
	}
	/* MALLOC space for the  1-yVector array */
	amat=(double*)MALLOC( (sizeof(double)*(unsigned long)(FeII.nFeIILevel_malloc*FeII.nFeIILevel_malloc) ));

	/* MALLOC space for the  1-yVector array */
	yVector=(double*)MALLOC( (sizeof(double)*(unsigned long)(FeII.nFeIILevel_malloc) ));

	/* set up array to save FeII line intensities */
	Fe2SavN = (double **)MALLOC(sizeof(double *)*(unsigned long)FeII.nFeIILevel_malloc );

	/* second dimension, lower level, for line save array */
	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ++ipHi )
	{
		Fe2SavN[ipHi]=(double*)MALLOC(sizeof(double )*(unsigned long)ipHi);
	}

	/* now MALLOC space for energy level table*/
	nnPradDat = (long*)MALLOC( (NPRADDAT+1)*sizeof(long) );

	/*Fe2DepCoef[NFE2LEVN];*/
	Fe2DepCoef = (double*)MALLOC( (unsigned long)FeII.nFeIILevel_malloc*sizeof(double) );

	/*Fe2LevelPop[NFE2LEVN];*/
	Fe2LevelPop = (double*)MALLOC( (unsigned long)FeII.nFeIILevel_malloc*sizeof(double) );

	/*Fe2ColDen[NFE2LEVN];*/
	Fe2ColDen = (double*)MALLOC( (unsigned long)FeII.nFeIILevel_malloc*sizeof(double) );

	/*FeIIBoltzmann[NFE2LEVN];*/
	FeIIBoltzmann = (double*)MALLOC( (unsigned long)FeII.nFeIILevel_malloc*sizeof(double) );


	/* MALLOC the realnum sPradDat[NPRADDAT][NPRADDAT][8] array */
	/* MALLOC the realnum sPradDat[ipHi][ipLo].t[8] array */
	sPradDat = ((realnum ***)MALLOC(NPRADDAT*sizeof(realnum **)));

	/* >>chng 00 dec 06, changed lower limit of loop to 1, Tru64 alpha's will not allocate 0 bytes!, PvH */
	sPradDat[0] = NULL;
	for(ipHi=1; ipHi < NPRADDAT; ipHi++) 
	{
		/* >>chng 00 dec 06, changed sizeof(realnum) to sizeof(realnum*), PvH */
		sPradDat[ipHi] = (realnum **)MALLOC((unsigned long)ipHi*sizeof(realnum *));

		/* now make space for the second dimension */
		for( ipLo=0; ipLo< ipHi; ipLo++ )
		{ 
			sPradDat[ipHi][ipLo] = (realnum *)MALLOC(8*sizeof(realnum ));
		}
	}

	/* now set junk to make sure reset before used */
	for(ipHi=0; ipHi < NPRADDAT; ipHi++) 
	{
		for( ipLo=0; ipLo< ipHi; ipLo++ )
		{ 
			for( k=0; k<8; ++k )
			{
				sPradDat[ipHi][ipLo][k] = -FLT_MAX;
			}
		}
	}

	/* now create the main emission line array and a helper array for the cs flag  */
	Fe2LevN=(transition**)MALLOC(sizeof(transition *)*(unsigned long)FeII.nFeIILevel_malloc);
	ncs1=(int**)MALLOC(sizeof(int *)*(unsigned long)FeII.nFeIILevel_malloc);

	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ++ipHi )
	{
		Fe2LevN[ipHi]=(transition*)MALLOC(sizeof(transition)*(unsigned long)ipHi);

		ncs1[ipHi]=(int*)MALLOC(sizeof(int)*(unsigned long)ipHi);
	}

#if	0
	/* now that its created, set to junk */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			Fe2LevN[ipHi][ipLo].Junk();
		}
	}
#endif

	/* now assign state and Emis pointers */
	Fe2LevN[1][0].Lo = AddState2Stack();
	/* now that its created, set to junk */
	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
	{
		/* add this upper level */
		Fe2LevN[ipHi][0].Hi = AddState2Stack();
		for( ipLo=0; ipLo < ipHi; ipLo++ )
		{
			if( ipLo == 0 )
			{
				Fe2LevN[ipHi][ipLo].Lo = Fe2LevN[1][0].Lo;
			}
			else
			{
				/* lower state is same as a previous upper state. */
				Fe2LevN[ipHi][ipLo].Lo = Fe2LevN[ipLo][0].Hi;
			}

			Fe2LevN[ipHi][ipLo].Hi = Fe2LevN[ipHi][0].Hi;
			Fe2LevN[ipHi][ipLo].Emis = AddLine2Stack( true );
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ," FeIICreate opening fe2nn.dat:");
	}

	ioDATA = open_data( "fe2nn.dat", "r" );

	ASSERT( ioDATA !=NULL );
	/* read in the fe2nn.dat file - this gives the Zheng and Pradhan number of level
	 * for every cloudy level.  So nnPradDat[1] is the index in the cloudy level
	 * counting for level 1 of Zheng & Pradan
	 * note that the order of some levels is different, the nnPradDat file goes 
	 * 21 23 22 - also that many levels are missing, the file goes 95 99 94 93 116
	 */
	for( i=0; i < 8; i++ )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		{
			fprintf( ioQQQ, " fe2nn.dat error reading data\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* get these integers from fe2nn.dat */
		k = 20*i;
		/* NPRADDAT is size of nnPradDat array, 159+1, make sure we do not exceed it */
		ASSERT( k+19 < NPRADDAT+1 );
		sscanf( chLine ,
			"%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld%4ld",
			&nnPradDat[k+0], &nnPradDat[k+1],  &nnPradDat[k+2], &nnPradDat[k+3], &nnPradDat[k+4],
			&nnPradDat[k+5], &nnPradDat[k+6],  &nnPradDat[k+7], &nnPradDat[k+8], &nnPradDat[k+9],
			&nnPradDat[k+10],&nnPradDat[k+11], &nnPradDat[k+12],&nnPradDat[k+13],&nnPradDat[k+14],
			&nnPradDat[k+15],&nnPradDat[k+16], &nnPradDat[k+17],&nnPradDat[k+18],&nnPradDat[k+19]
			);
#		if !defined(NDEBUG)
		for( m1=0; m1<20; ++m1 )
		{
			ASSERT( nnPradDat[k+m1] >= 0 && nnPradDat[k+m1] <= NFE2LEVN );
		}
#		endif
	}
	fclose(ioDATA);

	/* now get energies */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ," FeIICreate opening fe2energies.dat:");
	}

	ioDATA = open_data( "fe2energies.dat", "r" );

	/* file now open, read the data */
	for( ipHi=0; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
	{
		/* keep reading until have non-comment line, one that does not
		 * start with # */
		chLine[0] = '#';
		while( chLine[0] == '#' )
		{
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			{
				fprintf( ioQQQ, " fe2energies.dat error reading data\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* first and last number on this line */
		double help;
		sscanf( chLine, "%lf", &help );
		Fe2Energies[ipHi] = (realnum)help;
	}
	fclose(ioDATA);

	/* transition probabilities */

	if( trace.lgTrace )
		fprintf( ioQQQ," FeIICreate opening fe2rad.dat:");

	ioDATA = open_data( "fe2rad.dat", "r" );

	/* get the first line, this is a version number */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " fe2rad.dat error reading data\n" );
		cdEXIT(EXIT_FAILURE);
	}
	/* scan off three integers */
	sscanf( chLine ,"%ld%ld%ld",&lo, &ihi, &m1);
	const int nYrFe2Rad=8, nMoFe2Rad=8, nDyFe2Rad=24;
	if( lo!=nYrFe2Rad || ihi!=nMoFe2Rad || m1!=nDyFe2Rad )
	{
		fprintf( ioQQQ, "DISASTER fe2rad.dat has the wrong magic numbers, expected "
			"%2i %2i %2i and got %2li %2li %2li\n" ,
			nYrFe2Rad, nMoFe2Rad, nDyFe2Rad,
			lo, ihi, m1);
		cdEXIT(EXIT_FAILURE);
	}

	/* file now open, read the data */
	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
	{
		for( ipLo=0; ipLo < ipHi; ipLo++ )
		{
			/* following double since used in sscanf */
			double save[2];
			/* keep reading until have non-comment line, one that does not
			 * start with # */
			chLine[0] = '#';
			while( chLine[0] == '#' )
			{
				if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
				{
					fprintf( ioQQQ, " fe2nn.dat error reading data\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}

			/* first and last number on this line */
			sscanf( chLine ,
				"%ld%ld%ld%ld%lf%lf%ld",
				&lo, &ihi, &m1, &m2 ,
				&save[0], &save[1] , &m3);
			/* the indices ihi and lo within this array were 
			 * read in from the line above */
			Fe2LevN[ihi-1][lo-1].Lo->g = (realnum)m1;
			Fe2LevN[ihi-1][lo-1].Hi->g = (realnum)m2;
			Fe2LevN[ihi-1][lo-1].Emis->Aul = (realnum)save[0];
			/*>>chng 06 apr 10, option to use table of energies */
#			define	USE_OLD	true
			if( USE_OLD )
				Fe2LevN[ihi-1][lo-1].EnergyWN = (realnum)save[1];
			else
				Fe2LevN[ihi-1][lo-1].EnergyWN = Fe2Energies[ihi-1]-Fe2Energies[lo-1];

			/* Aul == -1 indicates intercombination line with no real Aul */
			if( fp_equal( Fe2LevN[ihi-1][lo-1].Emis->Aul , (realnum)-1.0f ) )
			{
				/* these are made-up intercombination lines, set gf to 1e-5 */
				Fe2LevN[ihi-1][lo-1].Emis->gf = 1e-5f;
				
				/* get corresponding A */
				Fe2LevN[ihi-1][lo-1].Emis->Aul = Fe2LevN[ihi-1][lo-1].Emis->gf*(realnum)TRANS_PROB_CONST*
					POW2(Fe2LevN[ihi-1][lo-1].EnergyWN)*Fe2LevN[ihi-1][lo-1].Lo->g/Fe2LevN[ihi-1][lo-1].Hi->g;
			}

			/* the last column of fe2rad.dat, and is 1, 2, or 3.  
			 * 1 signifies that transition is permitted,
			 * 2 is semi-forbidden
			 * 3 forbidden, within lowest 63 levels are forbidden, first permitted
			 * transition is from 64 */
			ncs1[ihi-1][lo-1] = (int)m3;
		}
	}
	fclose(ioDATA);

	/* now read collision data in fe2col.dat 
	 * These are from the following sources
	 >>refer	fe2	cs	Zhang, H.L., & Pradhan, A.K., 1995, A&A 293, 953 
	 >>refer	fe2	cs	Bautista, M., (private communication), 
	 >>refer	fe2	cs	Mewe, R., 1972, A&AS 20, 215 (the g-bar approximation)
	 */

	if( trace.lgTrace )
	{
		fprintf( ioQQQ," FeIICreate opening fe2col.dat:");
	}

	ioDATA = open_data( "fe2col.dat", "r" );

	ASSERT( ioDATA !=NULL);
	for( ipHi=1; ipHi<NPRADDAT; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			/* double since used in sscanf */
			double save[8];
 			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			{
				fprintf( ioQQQ, " fe2col.dat error reading data\n" );
				cdEXIT(EXIT_FAILURE);
			}
			sscanf( chLine ,
				"%ld%ld%lf%lf%lf%lf%lf%lf%lf%lf",
				&m1, &m2,
				&save[0], &save[1] , &save[2],&save[3], &save[4] , &save[5],
				&save[6], &save[7] 
				);
			for( k=0; k<8; ++k )
			{
				/* the max is here because there are some zeroes in the data file.
				 * this is unphysical but is part of their distribution.  As a result
				 * don't let the cs fall below 0.01 */
				sPradDat[m2-1][m1-1][k] = max(CS2SMALL , (realnum)save[k] );
			}
		}
	}

	fclose( ioDATA );

	/*generate needed opacity data for the large FeII atom */

	/* this routine is called only one time when cloudy is init
	 * for the very first time.  It converts the FeII data stored
	 * in the large FeII arrays into the array storage needed by cloudy
	 * for its line optical depth arrays
	 */

	/* convert large FeII line arrays into standard heavy el ar */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* pull information out of existing data arrays */

			ASSERT( Fe2LevN[ipHi][ipLo].EnergyWN > 0. );
			ASSERT( Fe2LevN[ipHi][ipLo].Emis->Aul> 0. );

			/* now put into standard internal line format
			Fe2LevN[ipHi][ipLo].WLAng = (realnum)(1.e8/Fe2LevN[ipHi][ipLo].EnergyWN); */
			/* >>chng 03 oct 28, above neglected index of refraction of air -
			 * convert to below */
			Fe2LevN[ipHi][ipLo].WLAng = 
				(realnum)(1.0e8/
			Fe2LevN[ipHi][ipLo].EnergyWN/
			RefIndex( Fe2LevN[ipHi][ipLo].EnergyWN ));

			/* generate gf from A */
			Fe2LevN[ipHi][ipLo].Emis->gf = 
				(realnum)(Fe2LevN[ipHi][ipLo].Emis->Aul*Fe2LevN[ipHi][ipLo].Hi->g/
				TRANS_PROB_CONST/POW2(Fe2LevN[ipHi][ipLo].EnergyWN));

			Fe2LevN[ipHi][ipLo].Hi->IonStg = 2;
			Fe2LevN[ipHi][ipLo].Hi->nelem = 26;
			/* which redistribution function??  
			 * For resonance line use K2 (-1),
			 * for subordinate lines use CRD with wings */
			/* >>chng 01 mar 09, all had been 1, inc redis with wings */
			/* to reset this, so that the code works as it did pre 01 mar 29,
			 * use command 
			 * atom FeII redistribution resonance pdr 
			 * atom FeII redistribution subordinate pdr */
			if( ipLo == 0 )
			{
				Fe2LevN[ipHi][ipLo].Emis->iRedisFun = FeII.ipRedisFcnResonance;
			}
			else
			{
				/* >>chng 01 feb 27, had been -1, crd with core only,
				 * change to crd with wings as per discussion with Ivan Hubeny */
				Fe2LevN[ipHi][ipLo].Emis->iRedisFun = FeII.ipRedisFcnSubordinate;
			}
			Fe2LevN[ipHi][ipLo].Emis->phots = 0.;
			Fe2LevN[ipHi][ipLo].Emis->ots = 0.;
			Fe2LevN[ipHi][ipLo].Emis->FracInwd = 1.;
		}
	}

	/* finally get FeII bands, this sets  */
	if( FeIIBandsCreate("") )
	{
		fprintf( ioQQQ," failed to open FeII bands file \n");
		cdEXIT(EXIT_FAILURE);
	}

	/* now establish the FeII continuum, these are set to
	 * 1000, 7000, and 1000 in FeIIZero in this file, and
	 * are reset with the atom FeII continuum command */
	FeIIContCreate( FeII.feconwlLo , FeII.feconwlHi , FeII.nfe2con );

	/* this must be true */
	ASSERT( FeII.nFeIILevel_malloc == FeII.nFeIILevel_local );

	return;
}

/*
 *=====================================================================
 */
/***********************************************************************
 *** version of FeIILevelPops.f with overlap in procces 05.28.97 ooooooo
 ******change in common block *te* sqrte 05.28.97
 *******texc is fixed 03.28.97
 *********this version has work on overlap 
 *********this version has # of zones (ZoneCnt) 07.03.97
 *********taux - optical depth depends on iter correction 03.03.97
 *********calculate cooling (Fe2_large_cool) and output fecool from Cloudy 01.29.97god
 *********and cooling derivative (ddT_Fe2_large_cool)
 ************ this program for 371 level iron model 12/14/1995
 ************ this program for 371 level iron model 1/11/1996
 ************ this program for 371 level iron model 3/21/1996
 ************ this program without La 3/27/1996
 ************ this program for 371 level iron model 4/9/1996
 ************ includes:FeIICollRatesBoltzmann;cooling;overlapping in lines */
STATIC bool lgFeIIEverCalled=false;
void FeIILevelPops( void )
{
	long int  i, 
		ipHi ,
		ipLo ,
		n;
	/* used in test for non-positive level populations */
	bool lgPopNeg;

	double 
	  EnergyWN,
	  pop ,
	  sum;

	int32 info;
	int32 ipiv[NFE2LEVN];

	DEBUG_ENTRY( "FeIILevelPops()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ,"   FeIILevelPops fe2 pops called\n");
	}

	/* FeII.lgSimulate was set true with simulate flag on atom FeII command,
	 * for bebugging without actually calling the routine */
	if( FeII.lgSimulate )
		return;

	lgFeIIEverCalled = true;
	/* zero out some arrays */
	for( n=0; n<FeII.nFeIILevel_local; ++n)
	{
		for( ipLo=0; ipLo<FeII.nFeIILevel_local; ++ipLo )
		{
			Fe2CPump[ipLo][n] = 0.;
			Fe2LPump[ipLo][n] = 0.;
			Fe2A[ipLo][n] = 0.;
			xMatrix[ipLo][n] = 0.;
		}
	}

	/* make the g-bar collision strengths and do linear interpolation on r-matrix data.
	 * this also sets Boltzmann factors for all levels,
	 * sets values of FeColl used below, but only if temp has changed */
	FeIICollRatesBoltzmann();

	/* pumping and spontantous decays */
	for( n=0; n<FeII.nFeIILevel_local; ++n)
	{
		for( ipHi=n+1; ipHi<FeII.nFeIILevel_local; ++ipHi )
		{
			/* continuum pumping rate from n to upper ipHi */
			Fe2CPump[n][ipHi] = Fe2LevN[ipHi][n].Emis->pump;

			/* continuum pumping rate from ipHi to lower n */
			Fe2CPump[ipHi][n] = Fe2LevN[ipHi][n].Emis->pump*
				Fe2LevN[ipHi][n].Hi->g/Fe2LevN[ipHi][n].Lo->g;

			/* spontaneous decays */
			Fe2A[ipHi][n] = Fe2LevN[ipHi][n].Emis->Aul*(Fe2LevN[ipHi][n].Emis->Pesc + Fe2LevN[ipHi][n].Emis->Pelec_esc +
			  Fe2LevN[ipHi][n].Emis->Pdest);
		}
	}

	/* now do pumping of atom by Lya */
	FeIILyaPump();

	/* ************************************************************************** 
	 *
	 * final rates into matrix 
	 *
	 ***************************************************************************/

	/* fill in xMatrix with matrix elements */
	for( n=0; n<FeII.nFeIILevel_local; ++n)
	{
		/* all processes leaving level n going down*/
		for( ipLo=0; ipLo<n; ++ipLo )
		{
			xMatrix[n][n] = xMatrix[n][n] + Fe2CPump[n][ipLo] + Fe2LPump[n][ipLo]+ Fe2A[n][ipLo] + 
				Fe2Coll[n][ipLo]*dense.eden;
		}
		/* all processes leaving level n going up*/
		for( ipHi=n+1; ipHi<FeII.nFeIILevel_local; ++ipHi )
		{
			xMatrix[n][n] = xMatrix[n][n] + Fe2CPump[n][ipHi] + Fe2LPump[n][ipHi] + Fe2Coll[n][ipHi]*dense.eden;
		}
		/* all processes entering level n from below*/
		for( ipLo=0; ipLo<n; ++ipLo )
		{
			xMatrix[ipLo][n] = xMatrix[ipLo][n] - Fe2CPump[ipLo][n] - Fe2LPump[ipLo][n] - Fe2Coll[ipLo][n]*dense.eden;
		}
		/* all processes entering level n from above*/
		for( ipHi=n+1; ipHi<FeII.nFeIILevel_local; ++ipHi )
		{
			xMatrix[ipHi][n] = xMatrix[ipHi][n] - Fe2CPump[ipHi][n] - Fe2LPump[ipHi][n] - Fe2Coll[ipHi][n]*dense.eden - 
				Fe2A[ipHi][n];
		}

		/* the top row of the matrix is the sum of level populations */
		xMatrix[n][0] = 1.0;
	}

	{
		/* option to print out entire matrix */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			/* print the matrices */
			for( n=0; n<FeII.nFeIILevel_local; ++n)
			{
				/*fprintf(ioQQQ,"\n");*/
				/* now print the matrix*/
				for( ipLo=0; ipLo<FeII.nFeIILevel_local; ++ipLo )
				{
					fprintf(ioQQQ," %.1e", xMatrix[n][ipLo]);
				}
				fprintf(ioQQQ,"\n");
			}
		}
	}

	/* define the Y Vector.  The oth element is the sum of all level populations
	 * adding up to the total population.  The remaining elements are the level
	 * balance equations adding up to zero */
	yVector[0] = 1.0;
	for( n=1; n < FeII.nFeIILevel_local; n++ )
	{
		yVector[n] = 0.0;
	}

	/* create the 1-yVector array that will save vector,
	 * this is the macro trick */
#	ifdef AMAT
#		undef AMAT
#	endif
#	define AMAT(I_,J_)	(*(amat+(I_)*FeII.nFeIILevel_local+(J_)))

	/* copy current contents of xMatrix array over to special amat array*/
	for( ipHi=0; ipHi < FeII.nFeIILevel_local; ipHi++ )
	{
		for( i=0; i < FeII.nFeIILevel_local; i++ )
		{
			AMAT(i,ipHi) = xMatrix[i][ipHi];
		}
	}

	info = 0;

	/* do the linear algebra to find the level populations */
	getrf_wrapper(FeII.nFeIILevel_local, FeII.nFeIILevel_local, amat, FeII.nFeIILevel_local, ipiv, &info);
	getrs_wrapper('N', FeII.nFeIILevel_local, 1, amat, FeII.nFeIILevel_local, ipiv, yVector, FeII.nFeIILevel_local, &info);

	if( info != 0 )
	{
		fprintf( ioQQQ, "DISASTER FeIILevelPops: dgetrs finds singular or ill-conditioned matrix\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* yVector now contains the level populations */

	/* this better be false after this loop - if not then non-positive level pops */
	lgPopNeg = false;
	/* copy all level pops over to Fe2LevelPop */
	for( ipLo=0; ipLo < FeII.nFeIILevel_local; ipLo++ )
	{
		if(yVector[ipLo] < 0. )
		{
			lgPopNeg = true;
			fprintf(ioQQQ,"PROBLEM FeIILevelPops finds non-positive level population, level is %ld pop is %g\n",
				ipLo , yVector[ipLo] );
		}
		/* this is now correct level population, cm^-3 */
		Fe2LevelPop[ipLo] = yVector[ipLo] * dense.xIonDense[ipIRON][1];
	}
	if( lgPopNeg )
	{
		/* this is here to use the lgPopNeg value for something, if not here then
		 * lint and some compilers will note that var is never used */
		fprintf(ioQQQ , "PROBLEM FeIILevelPops exits with negative level populations.\n");
	}

	/* >>chng 06 jun 24, make sure remainder of populations up through max
	 * limit are zero - this makes safe the case where the number
	 * of levels actually computed has been trimmed down from largest
	 * possible number of levels, for instance, in cool gas */
	for( ipLo=FeII.nFeIILevel_local; ipLo < FeII.nFeIILevel_malloc; ++ipLo )
	{
		Fe2LevelPop[ipLo] = 0.;
	}

	/* now set line opacities, intensities, and level populations 
	 * >>chng 06 jun ipLo should go up to FeII.nFeIILevel_local-1 since this
	 * is the largest lower level with non-zero population */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		/* >>chng 06 jun 24, upper level should go to limit
		 * of all that were allocated */
		for( ipHi=ipLo+1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* >>chng 06 jun 24, in all of these the product 
			 * yVector[ipHi]*dense.xIonDense[ipIRON][1] has been replaced
			 * with Fe2LevelPop[ipLo] - should always have been this way,
			 * and saves a mult */
			Fe2LevN[ipHi][ipLo].Emis->PopOpc = (Fe2LevelPop[ipLo] - 
				Fe2LevelPop[ipHi]*Fe2LevN[ipHi][ipLo].Lo->g/Fe2LevN[ipHi][ipLo].Hi->g);

			Fe2LevN[ipHi][ipLo].Lo->Pop = Fe2LevelPop[ipLo];

			Fe2LevN[ipHi][ipLo].Hi->Pop = Fe2LevelPop[ipHi];

			Fe2LevN[ipHi][ipLo].Emis->phots = Fe2LevelPop[ipHi]*
			  Fe2LevN[ipHi][ipLo].Emis->Aul*(Fe2LevN[ipHi][ipLo].Emis->Pesc + Fe2LevN[ipHi][ipLo].Emis->Pelec_esc );

			Fe2LevN[ipHi][ipLo].Emis->xIntensity = Fe2LevN[ipHi][ipLo].Emis->phots*
			  Fe2LevN[ipHi][ipLo].EnergyErg;

			/* ratio of collisional (new) to pumped excitations */
			/* >>chng 02 mar 04, add test MAX2 to prevent div by zero */
			Fe2LevN[ipHi][ipLo].Emis->ColOvTot = (realnum)(Fe2Coll[ipLo][ipHi]*dense.eden /
				SDIV( Fe2Coll[ipLo][ipHi]*dense.eden + Fe2CPump[ipLo][ipHi] + Fe2LPump[ipLo][ipHi] ) );
		}
	}

	/* only do this if large atom is on and more than ground terms computed - 
	 * do not if only lowest levels are computed */
	if( FeII.lgFeIILargeOn )
	{
		/* the hydrogen Lya destruction rate, then probability */
		hydro.dstfe2lya = 0.;
		EnergyWN = 0.;
		/* count how many photons were removed-added */
		for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
		{
			for( ipHi=ipLo+1; ipHi < FeII.nFeIILevel_local; ipHi++ )
			{
				EnergyWN += Fe2LPump[ipLo][ipHi] + Fe2LPump[ipHi][ipLo];
				hydro.dstfe2lya += (realnum)(
					Fe2LevN[ipHi][ipLo].Lo->Pop*Fe2LPump[ipLo][ipHi] -
					Fe2LevN[ipHi][ipLo].Hi->Pop*Fe2LPump[ipHi][ipLo]); 
			}
		}
		/* the destruction prob comes from
		 * dest rate = n(2p) * A(21) * PDest */
		pop = StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop;
		if( pop > SMALLFLOAT )
		{
			hydro.dstfe2lya /= (realnum)(pop * Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Aul);
		}
		else
		{
			hydro.dstfe2lya = 0.;
		}
		/* NB - do not update hydro.dstfe2lya here if large FeII not on since
		 * done in FeII overlap */
	}

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			/*lint -e644 EnergyWN not init */
			fprintf(ioQQQ," sum all %.1e dest rate%.1e escR= %.1e\n", 
				EnergyWN,hydro.dstfe2lya, 
				Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Pesc);
			/*lint +e644 EnergyWN not init */
		}
	}

	/* next two blocks determine departure coefficients for the atom */

	/* first sum up partition function for the model atom  */
	Fe2DepCoef[0] = 1.0;
	sum = 1.0;
	for( i=1; i < FeII.nFeIILevel_local; i++ )
	{
		/* Boltzmann factor relative to ground times ratio of statistical weights */
		Fe2DepCoef[i] = Fe2DepCoef[0]*FeIIBoltzmann[i]*
			Fe2LevN[i][0].Hi->g/ Fe2LevN[1][0].Lo->g;
		/* this sum is the partition function for the atom */
		sum += Fe2DepCoef[i];
	}

	/* now renormalize departure coefficients with partition function */
	for( i=0; i < FeII.nFeIILevel_local; i++ )
	{
		/* divide by total partition function, Fe2DepCoef is now the fraction of the
		 * population that is in this level in TE */
		Fe2DepCoef[i] /= sum;

		/* convert to true departure coefficient */
		if( Fe2DepCoef[i]>SMALLFLOAT )
		{
			Fe2DepCoef[i] = yVector[i]/Fe2DepCoef[i];
		}
		else
		{
			Fe2DepCoef[i] = 0.;
		}
	}

	/*calculate cooling, heating, and cooling derivative */
	/* this is net cooling, cooling minus heating */
	FeII.Fe2_large_cool = 0.0f;
	/* this is be heating, not heating minus cooling */
	FeII.Fe2_large_heat = 0.f;
	/* derivative of cooling */
	FeII.ddT_Fe2_large_cool = 0.0f;

	/* compute heating and cooling due to model atom */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_local; ipHi++ )
		{
			double OneLine;

			/* net cooling due to single line */
			OneLine = (Fe2Coll[ipLo][ipHi]*Fe2LevelPop[ipLo] - Fe2Coll[ipHi][ipLo]*Fe2LevelPop[ipHi])*
				dense.eden*Fe2LevN[ipHi][ipLo].EnergyErg;

			/* net cooling due to this line */
			FeII.Fe2_large_cool += (realnum)MAX2(0., OneLine);

			/* net heating due to line */
			FeII.Fe2_large_heat += (realnum)MAX2(0., -OneLine);

			/* derivative of FeII cooling */
			if( OneLine >= 0. )
			{
				/* net coolant, exponential dominates deriv */
				FeII.ddT_Fe2_large_cool += (realnum)OneLine*
					(Fe2LevN[ipHi][0].EnergyK*thermal.tsq1 - thermal.halfte);
			}
			else
			{
				/* net heating, te^-1/2 dominates */
				FeII.ddT_Fe2_large_cool -= (realnum)OneLine*thermal.halfte;
			}
		}
	}

	return;
}

/*FeIICollRatesBoltzmann evaluate collision strenths, 
 * both interpolating on r-mat and creating g-bar 
 * find Boltzmann factors, evaluate collisional rate coefficients */
/* >>05 dec 03, verified that this rountine only changes temperature dependent
 * quantities - nothing with temperature dependence is done */
STATIC void FeIICollRatesBoltzmann(void)
{
	/* this will be used to only reevaluate cs when temp changes */
	static double OldTemp = -1.;
	long int i,
		ipLo ,
		ipHi,
		ipTrim;
	realnum ag, 
	  cg, 
	  dg, 
	  gb, 
	  y;
	realnum FracLowTe , FracHighTe;
	static realnum tt[8]={1e3f,3e3f,5e3f,7e3f,1e4f,12e3f,15e3f,2e4f};
	long int ipTemp,
		ipTempp1 , 
		ipLoFe2, 
		ipHiFe2;
	static long int nFeII_old = -1;

	DEBUG_ENTRY( "FeIICollRatesBoltzmann()" );

	/* return if temperature has not changed */
	/* >>chng 06 feb 14, add test on number of levels changing */
	if( fp_equal( phycon.te, OldTemp ) && FeII.nFeIILevel_local == nFeII_old )
	{

		return;
	}
	OldTemp = phycon.te;
	nFeII_old = FeII.nFeIILevel_local;

	/* find g-bar collision strength for all levels 
	 * expression for g-bar taken from 
	 *>>refer	Fe2	g-bar	Mewe, R. 1972, A&A, 20, 215 */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_local; ipHi++ )
		{
			/* these coefficients are from page 221 of Mewe 1972 */
			if( ncs1[ipHi][ipLo] == 3 ) 
			{
				/************************forbidden tr**************************/
				ag = 0.15f;
				cg = 0.f;
				dg = 0.f;
			}
			/************************allowed tr*******************************/
			else if( ncs1[ipHi][ipLo] == 1 )
			{
				ag = 0.6f;
				cg = 0.f;
				dg = 0.28f;
			}
			/************************semiforbidden****************************/
			else if( ncs1[ipHi][ipLo] == 2 )
			{
				ag = 0.0f;
				cg = 0.1f;
				dg = 0.0f;
			}
			else
			{
				/* this branch is impossible, since cs must be 1, 2, or 3 */
				ag = 0.0f;
				cg = 0.1f;
				dg = 0.0f;
				fprintf(ioQQQ,">>>PROBLEM impossible ncs1 in FeIILevelPops\n");
			}

			/*y = 1.438826f*Fe2LevN[ipHi][ipLo].EnergyWN/ phycon.te;*/
			y = Fe2LevN[ipHi][ipLo].EnergyWN/ (realnum)phycon.te_wn;

			gb = (realnum)(ag + (-cg*POW2(y) + dg)*(log(1.0+1.0/y) - 0.4/
			  POW2(y + 1.0)) + cg*y);

			Fe2LevN[ipHi][ipLo].Coll.col_str = 23.861f*1e5f*gb*
			  Fe2LevN[ipHi][ipLo].Emis->Aul*Fe2LevN[ipHi][ipLo].Hi->g/
			  POW3(Fe2LevN[ipHi][ipLo].EnergyWN);

			/* confirm that collision strength is positive */
			ASSERT( Fe2LevN[ipHi][ipLo].Coll.col_str > 0.);

			/* g-bar cs becomes unphysically small for forbidden transitions -
			 * this sets a lower limit to the final cs - CS2SMALL is defined above */
			Fe2LevN[ipHi][ipLo].Coll.col_str = MAX2( CS2SMALL, Fe2LevN[ipHi][ipLo].Coll.col_str);
			/* this was the logic used in the old fortran version,
			 * and reproduces the results in Baldwin et al '96
			 if( Fe2LevN[ipHi][ipLo].cs < 1e-10 )
			{
				Fe2LevN[ipHi][ipLo].cs = 0.01f;
			}
			*/
		}
	}
	/* end loop setting Mewe 1972 g-bar approximation */

	/* we will interpolate on the set of listed collision strengths -
	 * where in this set are we? */
	if( phycon.te <= tt[0] )
	{
		/* temperature is lower than lowest tabulated, use the
		 * lowest tabulated point */
		/* ipTemp usually points to the cell cooler than te, ipTemp+1 to one higher,
		 * here both are lowest */
		ipTemp = 0;
		ipTempp1 = 0;
		FracHighTe = 0.;
	}
	else if( phycon.te > tt[7] )
	{
		/* temperature is higher than highest tabulated, use the
		 * highest tabulated point */
		/* ipTemp usually points to the cell cooler than te, ipTemp+1 to one higher,
		 * here both are highest */
		ipTemp = 7;
		ipTempp1 = 7;
		FracHighTe = 0.;
	}
	else
	{
		/* where in the array is the temperature we need? */
		ipTemp = -1;
		for( i=0; i < 8-1; i++ )
		{
			if( phycon.te <= tt[i+1] )
			{
				ipTemp = i;
				break;
			}

		}
		/* this cannot possibly happen */
		if( ipTemp < 0 )
		{
			fprintf( ioQQQ, " Insanity while looking for temperature in coll str array, te=%g.\n", 
			  phycon.te );
			cdEXIT(EXIT_FAILURE);
		}
		/* ipTemp points to the cell cooler than te, ipTemp+1 to one higher,
		 * do linear interpolation between */
		ipTemp = i;
		ipTempp1 = i+1;
		FracHighTe = ((realnum)phycon.te - tt[ipTemp])/(tt[ipTempp1] - tt[ipTemp]);
	}

	/* this is fraction of final linear interpolated collision strength that
	 * is weighted by the lower bound cs */
	FracLowTe = 1.f - FracHighTe;

	for( ipHi=1; ipHi < NPRADDAT; ipHi++ )
	{
		for( ipLo=0; ipLo < ipHi; ipLo++ )
		{
			/* ipHiFe2 should point to upper level of this pair, and
			 * ipLoFe2 should point to lower level */
			ipHiFe2 = MAX2( nnPradDat[ipHi] , nnPradDat[ipLo] );
			ipLoFe2 = MIN2( nnPradDat[ipHi] , nnPradDat[ipLo] );
			ASSERT( ipHiFe2-1 < NFE2LEVN );
			ASSERT( ipHiFe2-1 >= 0 );
			ASSERT( ipLoFe2-1 < NFE2LEVN );
			ASSERT( ipLoFe2-1 >= 0 );

			/* do linear interpolation for CS, these are dimensioned NPRADDAT = 159 */
			if( ipHiFe2-1 < FeII.nFeIILevel_local )
			{
				/*fprintf(ioQQQ,"DEBUG\t%.2e", Fe2LevN[ipHiFe2-1][ipLoFe2-1].cs);*/
				/* do linear interpolation */
				Fe2LevN[ipHiFe2-1][ipLoFe2-1].Coll.col_str = 
					sPradDat[ipHi][ipLo][ipTemp] * FracLowTe + 
					sPradDat[ipHi][ipLo][ipTempp1] * FracHighTe;
				/*fprintf(ioQQQ,"\t%.2e\n", Fe2LevN[ipHiFe2-1][ipLoFe2-1].cs);*/

				/* confirm that this is positive */
				ASSERT( Fe2LevN[ipHiFe2-1][ipLoFe2-1].Coll.col_str > 0. );
			}
		}
	}

	/* create Boltzmann factors for all levels */
	FeIIBoltzmann[0] = 1.0;
	for( ipHi=1; ipHi < FeII.nFeIILevel_local; ipHi++ )
	{
		/*FeIIBoltzmann[ipHi] = sexp( 1.438826f*Fe2LevN[ipHi][0].EnergyWN/phycon.te );*/
		/* >>chng 99 may 21, from above to following slight different number from phyconst.h 
		 * >>chng 05 dec 01, from sexp to dsexp to avoid all 0 Boltzmann factors in low temps
		 * now that FeII is on all the time */
		FeIIBoltzmann[ipHi] = dsexp( Fe2LevN[ipHi][0].EnergyWN/phycon.te_wn );
	}

	/* now possibly trim down atom if Boltzmann factors for upper levels are zero */
	ipTrim = FeII.nFeIILevel_local - 1;
	while( FeIIBoltzmann[ipTrim] == 0. && ipTrim > 0 )
	{
		--ipTrim;
	}
	/* ipTrim now is the highest level with finite Boltzmann factor - 
	 * use this as the number of levels 
	 *>>chng 05 dec 01, from <=1 to <=0 - func_map does 10K, and large FeII had never
	 * been tested in that limit.  1 is ok. */
	ASSERT( ipTrim > 0 );
	/* trim FeII.nFeIILevel_local so that FeIIBoltzmann is positive
	 * in nearly all cases this does nothing since ipTrim and FeII.nFeIILevel_local
	 * are equal . . . */
	if( ipTrim+1 < FeII.nFeIILevel_local )
	{
		/* >>chng 06 jun 27, zero out collision data */
		for( ipLo=0; ipLo<FeII.nFeIILevel_local; ++ipLo )
		{
			for( ipHi=ipTrim; ipHi<FeII.nFeIILevel_local; ++ipHi )
			{
				Fe2Coll[ipLo][ipHi] = 0.;
				Fe2Coll[ipHi][ipLo] = 0.;
			}
		}
	}
	FeII.nFeIILevel_local = ipTrim+1;
	/*fprintf(ioQQQ," levels reset to %li\n",FeII.nFeIILevel_local);*/

	/* debug code to print out the collision strengths for some levels */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			for( ipLo=0; ipLo < 52; ipLo++ )
			{
				fprintf(ioQQQ,"%e %e\n", 
					Fe2LevN[51][ipLo].Coll.col_str,Fe2LevN[52][ipLo].Coll.col_str);
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* collisional excitation and deexcitation */
	for( ipLo=0; ipLo<FeII.nFeIILevel_local; ++ipLo)
	{
		for( ipHi=ipLo+1; ipHi<FeII.nFeIILevel_local; ++ipHi )
		{
			/* collisional deexcitation rate coefficient from ipHi to lower ipLo
			 * note that it needs eden to become rate
			 * units cm3 s-1 */
			Fe2Coll[ipHi][ipLo] = (realnum)(COLL_CONST/phycon.sqrte*Fe2LevN[ipHi][ipLo].Coll.col_str/
			  Fe2LevN[ipHi][ipLo].Hi->g);
			/*fprintf(ioQQQ,"DEBUG fe2coll %li %li %.2e", ipHi , ipLo , Fe2Coll[ipHi][ipLo] );*/

			/* collisional excitation rate coefficient from ipLo to upper ipHi 
			 * units cm3 s-1 
			 * FeIIBoltzmann guaranteed to be > 0 by nFeIILevel_local trimming */
			Fe2Coll[ipLo][ipHi] = (realnum)(Fe2Coll[ipHi][ipLo]*FeIIBoltzmann[ipHi]/FeIIBoltzmann[ipLo]*
				Fe2LevN[ipHi][ipLo].Hi->g/Fe2LevN[ipHi][ipLo].Lo->g);
			/*fprintf(ioQQQ,"DEBUG fe2coll %.2e\n", Fe2Coll[ipLo][ipHi] );*/
		}
	}

	return;
}

/*
 *=====================================================================
 */
/*subroutine FeIIPrint PhotOccNum_at_nu raspechatki naselennostej v cloudy.out ili v
 * otdel'nom file unit=33
 *!!nado takzhe vklyuchit raspechatku iz perekrytiya linii */
/*FeIIPrint - print output from large FeII atom, called by prtzone */
void FeIIPrint(void)
{

	DEBUG_ENTRY( "FeIIPrint()" );

	return;
}

/*
 *=====================================================================
 */
/*FeIISumBand, sum up large FeII emission over certain bands
 * units are erg s-1 cm-3, same units as xIntensity in line structure */
/* >>chng 06 apr 11, from using erg as energy unit to vacuum angstroms,
 * this fixes bug in physical constant and also uses air wavelengths for wl > 2000A */
double FeIISumBand(
	/* lower and upper range to wavelength in Angstroms */
	realnum wl1, 
	realnum wl2,
	double *SumBandInward )
{
	long int ipHi, 
	  ipLo;
	double SumBandFe2_v;

	DEBUG_ENTRY( "FeIISumBand()" );
	/*sum up large FeII emission over certain bands */

	SumBandFe2_v = 0.;
	*SumBandInward = 0.;
	if( dense.xIonDense[ipIRON][1] > SMALLFLOAT )
	{
		/* line energy in wavenumber */
		ASSERT( wl2 > wl1 );
		for( ipHi=1; ipHi < FeII.nFeIILevel_local; ++ipHi )
		{
			for( ipLo=0; ipLo < ipHi; ++ipLo )
			{
				if( Fe2LevN[ipHi][ipLo].WLAng >= wl1 && 
				  Fe2LevN[ipHi][ipLo].WLAng < wl2 )
				{
					SumBandFe2_v += Fe2LevN[ipHi][ipLo].Emis->xIntensity;
					*SumBandInward += Fe2LevN[ipHi][ipLo].Emis->xIntensity *
						Fe2LevN[ipHi][ipLo].Emis->FracInwd;
				}
			}
		}
	}
	return( SumBandFe2_v );
}

/*
 *=====================================================================
 */
/*FeII_RT_TauInc called once per zone in RT_tau_inc to increment large FeII atom line optical depths */
void FeII_RT_TauInc(void)
{
	long int ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeII_RT_TauInc()" );

	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		/* >>chng 06 jun 24, for upper level go all the way to the
		 * largest possible number of levels so that optical depths
		 * of UV transitions are correct even for very cold gas where
		 * the high level populations are not computed */
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* >>chng 04 aug 28, add test on bogus line */
			if( Fe2LevN[ipHi][ipLo].ipCont > 0 )
				RT_line_one_tauinc( &Fe2LevN[ipHi][ipLo], -8, -8, ipHi, ipLo, GetDopplerWidth(dense.AtomicWeight[ipIRON]) );
		}
	}

	return;
}

/*
 *=====================================================================
 */
/*FeII_RT_tau_reset reset optical depths for large FeII atom, called by update after each iteration */
void FeII_RT_tau_reset(void)
{
	long int ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeII_RT_tau_reset()" );

	/*fprintf(ioQQQ,"DEBUG FeIITauAver1 %li %.2e %.2e nFeIILevel_local %li \n",
		iteration,
		Fe2LevN[9][0].Emis->TauIn ,
		Fe2LevN[9][0].Emis->TauTot, 
		FeII.nFeIILevel_local);*/

	/* called at end of iteration */
	/* >>chng 05 dec 07, had been FeII.nFeIILevel_local but this may have been trimmed down
	 * if previous iteration went very deep - not reset until FeIIReset called */
	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
	{
		for( ipLo=0; ipLo < ipHi; ipLo++ )
		{
			RT_line_one_tau_reset( &Fe2LevN[ipHi][ipLo] );
		}
	}

	return;
}

/*
 *=====================================================================
 */
/*FeIIPoint called by ContCreatePointers to create pointers for lines in large FeII atom */
void FeIIPoint(void)
{
	long int ipHi, 
	  ip, 
	  ipLo;

	DEBUG_ENTRY( "FeIIPoint()" );

	/* routine called when cloudy sets continuum array indices for Fe2 lines, 
	 * once per coreload */
	for( ipLo=0; ipLo < FeII.nFeIILevel_malloc-1; ipLo++ )
	{
		for( ipHi=ipLo+1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{

			/* >>chng 02 feb 11, set continuum index to negative value for fake transition */
			if( fabs(Fe2LevN[ipHi][ipLo].Emis->Aul- 1e-5 ) > 1e-8 ) 
			{
				ip = ipoint(Fe2LevN[ipHi][ipLo].EnergyWN * WAVNRYD);
				Fe2LevN[ipHi][ipLo].ipCont = ip;

				/* do not over write other pointers with FeII since FeII is everywhere */
				if( strcmp( rfield.chLineLabel[ip-1], "    " ) == 0 )
					strcpy( rfield.chLineLabel[ip-1], "FeII" );

				Fe2LevN[ipHi][ipLo].Emis->ipFine = ipFineCont( Fe2LevN[ipHi][ipLo].EnergyWN * WAVNRYD);
			}
			else
			{
				Fe2LevN[ipHi][ipLo].ipCont = -1;
				Fe2LevN[ipHi][ipLo].Emis->ipFine = -1;
			}

			Fe2LevN[ipHi][ipLo].Emis->dampXvel = 
				(realnum)(Fe2LevN[ipHi][ipLo].Emis->Aul/
			  Fe2LevN[ipHi][ipLo].EnergyWN/PI4);

			/* derive the abs coef, call to function is gf, wl (A), g_low */
			Fe2LevN[ipHi][ipLo].Emis->opacity = 
				(realnum)abscf(Fe2LevN[ipHi][ipLo].Emis->gf,
			  Fe2LevN[ipHi][ipLo].EnergyWN,
			  Fe2LevN[ipHi][ipLo].Lo->g);

			/* excitation energy of transition in degrees kelvin */
			Fe2LevN[ipHi][ipLo].EnergyK = 
				(realnum)(T1CM*Fe2LevN[ipHi][ipLo].EnergyWN);

			/* energy of photon in ergs */
			Fe2LevN[ipHi][ipLo].EnergyErg = 
				(realnum)(ERG1CM*Fe2LevN[ipHi][ipLo].EnergyWN);
		}
	}

	return;
}

/*
 *=====================================================================
 */
/*FeIIAccel called by rt_line_driving to compute radiative acceleration due to FeII lines */
void FeIIAccel(double *fe2drive)
{
	long int ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeIIAccel()" );
	/*compute acceleration due to large FeII atom */

	/* this routine computes the line driven radiative acceleration
	 * due to large FeII atom*/

	*fe2drive = 0.;
	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		for( ipHi=ipLo+1; ipHi < FeII.nFeIILevel_local; ipHi++ )
		{
			*fe2drive += Fe2LevN[ipHi][ipLo].Emis->pump*
			  Fe2LevN[ipHi][ipLo].EnergyErg*Fe2LevN[ipHi][ipLo].Emis->PopOpc;
		}
	}

	return;
}

/*
 *=====================================================================
 */
/*FeII_RT_Make called by RT_line_all, does large FeII atom radiative transfer */
void FeII_RT_Make( void )
{
	long int ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeII_RT_Make()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ,"   FeII_RT_Make called\n");
	}

	/* this routine drives calls to make RT relations with large FeII atom */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		/* >>chng 06 jun 24, go up to number allocated to keep UV resonance lines */
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* only evaluate real transitions */
			if( Fe2LevN[ipHi][ipLo].ipCont > 0 ) 
			{
				/*RT_line_one do rt for emission line structure - 
				 * calls RT_line_escape or RT_line_wind */
				RT_line_one( &Fe2LevN[ipHi][ipLo], true, 0.f, GetDopplerWidth(dense.AtomicWeight[ipIRON]) );
			}
		}
	}

	return;
}

/*
 *=====================================================================
 */
/* called in LineSet4 to add FeII lines to save array */
void FeIIAddLines( void )
{

	DEBUG_ENTRY( "FeIIAddLines()" );

	/* this routine puts the emission from the large FeII atom
	 * into an array to save and integrate them*/

	/* add lines option called from lines, add intensities into storage array */

	/* routine is called three different ways, ipass < 0 before space allocated,
	 * =0 when time to generate labels (and we zero out array here), and ipass>0
	 * when time to save intensities */
	if( LineSave.ipass == 0 )
	{
		/* we were called by lines, and we want to zero out Fe2SavN */
		for( long ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
		{
			for( long ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
			{
				Fe2SavN[ipHi][ipLo] = 0.;
			}
		}
	}

	/* this call during calculations, save intensities */
	else if( LineSave.ipass == 1 )
	{
		/* total emission from vol */
		for( long ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
		{
			for( long ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_local; ipHi++ )
			{
				Fe2SavN[ipHi][ipLo] += 
					radius.dVeffAper*Fe2LevN[ipHi][ipLo].Emis->xIntensity;
			}
		}
	}

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /*&& (iteration==2)*/ )
		{
			fprintf(ioQQQ," 69-1\t%li\t%e\n", nzone , Fe2SavN[68][0] );
		}
	}

	return;
}

/*FeIIPunchLevels save level energies and stat weights */
void FeIIPunchLevels(
  /* file we will save to */
  FILE *ioPUN )
{

	long int ipHi;

	DEBUG_ENTRY( "FeIIPunchLevels()" );

	fprintf(ioPUN , "%.2f\t%li\n", 0., (long)Fe2LevN[1][0].Lo->g );
	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ++ipHi )
	{
		fprintf(ioPUN , "%.2f\t%li\n", 
			Fe2LevN[ipHi][0].EnergyWN ,
			(long)Fe2LevN[ipHi][0].Hi->g);
	}

	return;
}

/*FeIIPunchColden save level energies, stat weights, column density */
void FeIIPunchColden(
  /* file we will save to */
  FILE *ioPUN )
{

	long int n;

	DEBUG_ENTRY( "FeIIPunchColden()" );

	fprintf(ioPUN , "%.2f\t%li\t%.3e\n", 0., (long)Fe2LevN[1][0].Lo->g , Fe2ColDen[0]);
	for( n=1; n < FeII.nFeIILevel_malloc; ++n )
	{
		fprintf(ioPUN , "%.2f\t%li\t%.3e\n", 
			Fe2LevN[n][0].EnergyWN ,
			(long)Fe2LevN[n][0].Hi->g,
			Fe2ColDen[n]);
	}

	return;
}


/*FeIIPunchOpticalDepth save FeII optical depths,
 * called by punch_do to save them out,
 * save turned on with save lines command */
void FeIIPunchOpticalDepth(
  /* file we will save to */
  FILE *ioPUN )
{
	long int 
	  ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeIIPunchOpticalDepth()" );

	/* this routine punches the optical depths predicted by the large FeII atom */
	/*>>chng 06 may 19, must use total number of levels allocated not current
	 * number since this decreases as gas grows colder with depth nFeIILevel_local->nFeIILevel_malloc*/
	for( ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* fe2ener(1) and (2) are lower and upper limits to range of
			 * wavelengths of interest.  entered in ryd with
			 * save FeII verner command, where they are converted
			 * to wavenumbers, */
			fprintf( ioPUN, "%ld\t%ld\t%.2f\t%.2e\n", 
				ipHi+1, 
				ipLo+1, 
				Fe2LevN[ipHi][ipLo].WLAng ,
				Fe2LevN[ipHi][ipLo].Emis->TauIn );
		}
	}

	return;
}

/*FeIISaveLines save accumulated FeII intensities, at end of calculation,
 * called by punch_do to save them out,
 * save turned on with save lines command */
void FeIISaveLines(
  /* file we will save to */
  FILE *ioPUN )
{
	long int MaseHi, 
	  MaseLow, 
	  ipHi, 
	  ipLo;
	double hbeta, absint , renorm;
	/* >>chng 00 jun 02, demoted next two to realnum, PvH */
	realnum TauMase, 
	  thresh;

	DEBUG_ENTRY( "FeIISaveLines()" );

	/* this routine puts the emission from the large FeII atom
	 * into a line array, and eventually will save it out */

	/* get the normalization line */
	if( LineSv[LineSave.ipNormWavL].SumLine[0] > 0. )
		renorm = LineSave.ScaleNormLine/LineSv[LineSave.ipNormWavL].SumLine[0];
	else
		renorm = 1.;

	fprintf( ioPUN, " up low log I, I/I(LineSave), Tau\n" );

	/* first look for any masing lines */
	MaseLow = -1;
	MaseHi = -1;
	TauMase = 0.f;
	for( ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			if( Fe2LevN[ipHi][ipLo].Emis->TauIn < TauMase )
			{
				TauMase = Fe2LevN[ipHi][ipLo].Emis->TauIn;
				MaseLow = ipLo;
				MaseHi = ipHi;
			}
		}
	}

	if( TauMase < 0.f )
	{
		fprintf( ioPUN, " Most negative optical depth was %4ld%4ld%10.2e\n", 
		  MaseHi, MaseLow, TauMase );
	}

	/* now print actual line intensities, Hbeta first */
	if( cdLine("TOTL", 4861.36f , &hbeta , &absint)<=0 )
	{
		fprintf( ioQQQ, " FeIILevelPops could not find Hbeta with cdLine.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	fprintf( ioPUN, "Hbeta=%7.3f %10.2e\n", 
	  absint , 
	  hbeta );

	if( renorm > SMALLFLOAT )
		/* this is threshold for faintest line, normally 0, set with 
		 * number on save feii lines command */
		thresh = FeII.fe2thresh/(realnum)renorm;
	else
		thresh = 0.f;

	/* must use total number of levels allocated not current
	 * number since this decreases as gas grows colder with depth nFeIILevel_malloc*/
	for( ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* fe2ener(1) and (2) are lower and upper limits to range of
			 * wavelengths of interest.  entered in ryd with
			 * save FeII verner command, where they are converted
			 * to wavenumbers, */
			if( (Fe2SavN[ipHi][ipLo] > thresh && 
				Fe2LevN[ipHi][ipLo].EnergyWN > FeII.fe2ener[0]) &&
				Fe2LevN[ipHi][ipLo].EnergyWN < FeII.fe2ener[1] )
			{
				if( FeII.lgShortFe2 )
				{
					/* short option on save FeII line command 
					 * does not include rel intensity or optical dep */
					fprintf( ioPUN, "%ld\t%ld\t%.2f\t%.3f\n", 
					  ipHi+1, 
					  ipLo+1, 
					  Fe2LevN[ipHi][ipLo].WLAng ,
					  log10(MAX2(1e-37,Fe2SavN[ipHi][ipLo])) + radius.Conv2PrtInten );
				}
				else
				{
					/* long printout does */
					fprintf( ioPUN, "%ld\t%ld\t%.2f\t%.3f\t%.2e\t%.2e\n", 
					  ipHi+1, 
					  ipLo+1, 
					  Fe2LevN[ipHi][ipLo].WLAng ,
					  log10(MAX2(1e-37,Fe2SavN[ipHi][ipLo])) + radius.Conv2PrtInten, 
					  Fe2SavN[ipHi][ipLo]*renorm ,
					  Fe2LevN[ipHi][ipLo].Emis->TauIn );
				}
			}
		}
	}

	return;
}


/*
 *=====================================================================
 */
/*FeII_LineZero zero out storage for large FeII atom, called in tauout */
void FeII_LineZero(void)
{
	long int ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeII_LineZero()" );

	/* this routine is called in routine zero and it
	* zero's out various elements of the FeII arrays
	* it is called on every iteration
	* */
	for( ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* >>chng 03 feb 14, use EmLineZero rather than explicit sets */
			Fe2LevN[ipHi][ipLo].Zero();
		}
	}

	return;
}
/*
 *=====================================================================
 */
/*FeIIIntenZero zero out storage for large FeII atom, called in tauout */
void FeIIIntenZero(void)
{
	long int ipHi, 
	  ipLo;

	DEBUG_ENTRY( "FeIIIntenZero()" );

	/* this routine is called by routine cool_iron and it
	 * zero's out various elements of the FeII arrays
	 * */
	Fe2LevelPop[0] = 0.;
	for( ipHi=1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
	{
		/* >>chng 05 dec 14, zero Fe2LevelPop */
		Fe2LevelPop[ipHi] = 0.;
		for( ipLo=0; ipLo < ipHi; ipLo++ )
		{

			/* population of lower level with correction for stim emission */
			Fe2LevN[ipHi][ipLo].Emis->PopOpc = 0.;

			/* population of lower level */
			Fe2LevN[ipHi][ipLo].Lo->Pop = 0.;

			/* population of upper level */
			Fe2LevN[ipHi][ipLo].Hi->Pop = 0.;

			/* following two heat exchange excitation, deexcitation */
			Fe2LevN[ipHi][ipLo].Coll.cool = 0.;
			Fe2LevN[ipHi][ipLo].Coll.heat = 0.;

			/* intensity of line */
			Fe2LevN[ipHi][ipLo].Emis->xIntensity = 0.;

			Fe2LevN[ipHi][ipLo].Emis->phots = 0.;
			Fe2LevN[ipHi][ipLo].Emis->ots = 0.;
			Fe2LevN[ipHi][ipLo].Emis->ColOvTot = 0.;
		}
	}

	TauLines[ipTuv3].Lo->Pop = 0;
	TauLines[ipTuv3].Hi->Pop = 0;
	TauLines[ipTuv3].Emis->PopOpc = 0;
	TauLines[ipTuv3].Emis->phots = 0;
	TauLines[ipTuv3].Emis->xIntensity = 0;

	TauLines[ipTr48].Lo->Pop = 0;
	TauLines[ipTr48].Hi->Pop = 0;
	TauLines[ipTr48].Emis->PopOpc = 0;
	TauLines[ipTr48].Emis->phots = 0;
	TauLines[ipTr48].Emis->xIntensity = 0;

	FeII.for7 = 0;

	TauLines[ipTFe16].Lo->Pop = 0;
	TauLines[ipTFe16].Hi->Pop = 0;
	TauLines[ipTFe16].Emis->PopOpc = 0;
	TauLines[ipTFe16].Emis->phots = 0;
	TauLines[ipTFe16].Emis->xIntensity = 0;

	TauLines[ipTFe26].Lo->Pop = 0;
	TauLines[ipTFe26].Hi->Pop = 0;
	TauLines[ipTFe26].Emis->PopOpc = 0;
	TauLines[ipTFe26].Emis->phots = 0;
	TauLines[ipTFe26].Emis->xIntensity = 0;

	TauLines[ipTFe34].Lo->Pop = 0;
	TauLines[ipTFe34].Hi->Pop = 0;
	TauLines[ipTFe34].Emis->PopOpc = 0;
	TauLines[ipTFe34].Emis->phots = 0;
	TauLines[ipTFe34].Emis->xIntensity = 0;

	TauLines[ipTFe35].Lo->Pop = 0;
	TauLines[ipTFe35].Hi->Pop = 0;
	TauLines[ipTFe35].Emis->PopOpc = 0;
	TauLines[ipTFe35].Emis->phots = 0;
	TauLines[ipTFe35].Emis->xIntensity = 0;

	TauLines[ipTFe46].Lo->Pop = 0;
	TauLines[ipTFe46].Hi->Pop = 0;
	TauLines[ipTFe46].Emis->PopOpc = 0;
	TauLines[ipTFe46].Emis->phots = 0;
	TauLines[ipTFe46].Emis->xIntensity = 0;

	TauLines[ipTFe56].Lo->Pop = 0;
	TauLines[ipTFe56].Hi->Pop = 0;
	TauLines[ipTFe56].Emis->PopOpc = 0;
	TauLines[ipTFe56].Emis->phots = 0;
	TauLines[ipTFe56].Emis->xIntensity = 0;

	TauLines[ipT191].Lo->Pop = 0;
	TauLines[ipT191].Hi->Pop = 0;
	TauLines[ipT191].Emis->PopOpc = 0;
	TauLines[ipT191].Emis->phots = 0;
	TauLines[ipT191].Emis->xIntensity = 0;

	return;
}

/*
 *=====================================================================
 * save line data for FeII atom */
void FeIIPunData(
	/* io unit for save */
	FILE* ioPUN ,
	/* save all levels if true, only subset if false */
	bool lgDoAll )
{
	long int ipLo , ipHi;

	DEBUG_ENTRY( "FeIIPunData()" );

	if( lgDoAll )
	{
		fprintf( ioQQQ, 
			" FeIIPunData ALL option not implemented yet 1\n" );
		cdEXIT(EXIT_FAILURE);
	}
	else if( lgFeIIEverCalled )
	{
		long int nSkip=0, limit=MIN2(64, FeII.nFeIILevel_local);

		/* false, only save subset in above init */
		/* first 64 do all lines */
		for( ipHi=1; ipHi<limit; ++ipHi )
		{
			for( ipLo=0; ipLo<ipHi; ++ipLo )
			{
				fprintf(ioPUN,"%4li%4li ",ipLo,ipHi );
				Save1LineData( &Fe2LevN[ipHi][ipLo] , ioPUN, false);
			}
		}
		fprintf( ioPUN , "\n");

		if( limit==64 )
		{
			/* higher than 64 only do real transitions - the majority of those above
			 * n = 64 have no radiative transitions but still exist to hold collision,
			 * energy, and other line data - they will have Aul == 1e-5 */
			for( ipHi=limit; ipHi<FeII.nFeIILevel_local; ++ipHi )
			{
				/*>>chng 06 jun 23, ipLo had ranged from limit up to ipHi,
				 * so never printed lines from ipHi to ipLo<64 */
				for( ipLo=0; ipLo<ipHi; ++ipLo )
				{
					/*fprintf(ioPUN , "%li %li\n", ipHi , ipLo );
					if( ipHi==70 && ipLo==0 )
						Save1LineData( &Fe2LevN[ipHi][ipLo] , ioPUN , false); */
					/* ncs1 of 3 and Aul = 1e-5 indicate transition that is not optically
					 * allowed with fake cs */
					if( ncs1[ipHi][ipLo] == 3 && 
						(fabs(Fe2LevN[ipHi][ipLo].Emis->Aul-1e-5) < 1e-8 ) )
					{
						++nSkip;
					}
					else
					{
						/* add one to ipLo and ipHi so that printed number is on atomic, 
						 * not C, scale */
						fprintf(ioPUN,"%4li%4li ",ipLo+1,ipHi+1 );
						Save1LineData( &Fe2LevN[ipHi][ipLo] , ioPUN , false);
					}
				}
			}
			fprintf( ioPUN , " %li lines skiped\n",nSkip);
		}
	}

	return;
}

/*
 *=====================================================================
 */
void FeIIPunDepart(
	/* io unit for save */
	FILE* ioPUN ,
	/* save all levels if true, only subset if false */
	bool lgDoAll )
{
	/* numer of levels with dep coef that we will save out */
#	define NLEVDEP 11
	/* these are the levels on the physical, not c, scale (count from 1) */
	const int LevDep[NLEVDEP]={1,10,25,45,64,124,206,249,295,347,371};
	long int i;
	static bool lgFIRST=true;

	DEBUG_ENTRY( "FeIIPunDepart()" );

	/* on first call only, print levels that we will save later */
	if( lgFIRST && !lgDoAll )
	{
		/* but all the rest do */
		for( i=0; i<NLEVDEP; ++i )
		{
			fprintf( ioPUN , "%i\t", LevDep[i] );
		}
		fprintf( ioPUN , "\n");
		lgFIRST = false;
	}

	if( lgDoAll )
	{
		/* true, save all levels, one per line */
		for( i=1; i<=FeII.nFeIILevel_local; ++i )
		{
			FeIIPun1Depart( ioPUN , i );
			fprintf( ioPUN , "\n");
		}
	}
	else
	{
		/* false, only save subset in above init */
		for( i=0; i<NLEVDEP; ++i )
		{
			FeIIPun1Depart( ioPUN , LevDep[i] );
			fprintf( ioPUN , "\t");
		}
		fprintf( ioPUN , "\n");
	}

	return;
}

/*
 *=====================================================================
 */
void FeIIPun1Depart( 
	/* the io unit where the print should be directed */
	FILE * ioPUN , 
	/* the physical (not c) number of the level */
	long int nPUN )
{

	DEBUG_ENTRY( "FeIIPun1Depart()" );

	ASSERT( nPUN > 0 );
	/* need real assert to keep lint happy */
	assert( ioPUN != NULL );

	/* print the level departure coef on ioPUN if the level was computed,
	 * print a zero if it was not */
	if( nPUN <= FeII.nFeIILevel_local )
	{
		fprintf( ioPUN , "%e ",Fe2DepCoef[nPUN-1] );
	}
	else
	{
		fprintf( ioPUN , "%e ",0. );
	}

	return;
}

/*
 *=====================================================================
 */
void FeIIReset(void)
{

	DEBUG_ENTRY( "FeIIReset()" );

	/* space has been allocated, so reset number of levels to whatever it was */
	FeII.nFeIILevel_local = FeII.nFeIILevel_malloc;

	return;
}

/*
 *=====================================================================
 */
/*FeIIZero initialize some variables, called by zero one time before commands parsed*/
void FeIIZero(void)
{

	DEBUG_ENTRY( "FeIIZero()" );

	/* heating, cooling, and deriv wrt temperature */
	FeII.Fe2_large_cool = 0.;
	FeII.ddT_Fe2_large_cool = 0.;
	FeII.Fe2_large_heat = 0.;

	/* flag saying that lya pumping of FeII in large atom is turned on */
	FeII.lgLyaPumpOn = true;

	/*FeII. lgFeIION = false;*/
	/* >>chng 05 nov 29, test effects of always including FeII ground term with large atom
	 * but if ground term only is done, still also do simple UV approximation */
	/*FeII. lgFeIION = true;*/
	/* will not compute large FeII atom */
	FeII.lgFeIILargeOn = false;

	/* energy range of large FeII atom is zero to infinity */
	FeII.fe2ener[0] = 0.;
	FeII.fe2ener[1] = 1e8;

	/* pre mar 01, these had both been 1, ipPRD */
	/* resonance lines, ipCRD is -1 */
	FeII.ipRedisFcnResonance = ipCRD;
	/* subordinate lines, ipCRDW is 2 */
	FeII.ipRedisFcnSubordinate = ipCRDW;

	/* set zero for the threshold of weakest large FeII atom line to save */
	FeII.fe2thresh = 0.;

	/* normally do not constantly reevaluate the atom, set true with
	 * SLOW key on atom FeII command */
	FeII.lgSlow = false;

	/* option to print each call to FeIILevelPops, set with print option on atom FeII */
	FeII.lgPrint = false;

	/* option to only simulate calls to FeIILevelPops */
	FeII.lgSimulate = false;

	/* set number of levels for the atom
	 * changed with the atom FeII levels command */
	if( lgFeIIMalloc )
	{
		/* space has been allocated, so reset number of levels to whatever it was */
		FeII.nFeIILevel_local = FeII.nFeIILevel_malloc;
	}
	else
	{
		/* always include FeII ground term with large atom
		 * but if ground term only is done, still also do simple UV approximation 
		 * set this to only ground term, - will reset to NFE2LEVN when atom FeII parsed if levels not set */
		FeII.nFeIILevel_local = 16;
	}

	return;
}

/*FeIIPunPop - save level populations */
void FeIIPunPop(
	/* io unit for save */
	FILE* ioPUN ,
	/* save range of levels if true, only selected subset if false */
	bool lgPunchRange ,
	/* if ipPunchRange is true, this is lower bound of range on C scale */
	long int ipRangeLo ,
	/* if ipPunchRange is true, this is upper bound of range on C scale */
	long int ipRangeHi ,
	/* flag saying whether to save density (cm-3, true) or relative population (flase) */
	bool lgPunchDensity )
{
	/* numer of levels with dep coef that we will save out */
#	define NLEVPOP 11
	/* these are the levels on the physical, not c, scale (count from 1) */
	const int LevPop[NLEVPOP]={1,10,25,45,64,124,206,249,295,347,371};
	long int i;
	static bool lgFIRST=true;
	realnum denominator = 1.f;

	DEBUG_ENTRY( "FeIIPunPop()" );

	/* this implements the relative option on save FeII populations command */
	if( !lgPunchDensity )
		denominator = SDIV( dense.xIonDense[ipIRON][1] );

	/* on first call only, print levels that we will save later,
	 * but only if we will only save selected levels*/
	if( lgFIRST && !lgPunchRange )
	{
		/* but all the rest do */
		for( i=0; i<NLEVPOP; ++i )
		{
			/* indices for particular levels */
			fprintf( ioPUN , "%i\t", LevPop[i] );
		}
		fprintf( ioPUN , "\n");
		lgFIRST = false;
	}

	if( lgPunchRange )
	{
		/* confirm sane level indices */
		ASSERT( ipRangeLo>=0 && ipRangeLo<ipRangeHi  );

		/* true, save all levels across line,
		 * both call with physical level so that list is physical */
		long nHigh = MIN2(ipRangeHi, FeII.nFeIILevel_local );
		for( i=ipRangeLo; i<nHigh; ++i )
		{
			/* routine takes levels on physics scale */
			fprintf( ioPUN , "%.3e\t", Fe2LevelPop[i]/denominator );
		}
		fprintf( ioPUN , "\n");
	}
	else
	{
		/* false, only save subset in above init,
		 * both call with physical level so that list is physical  */
		for( i=0; i<NLEVPOP; ++i )
		{
			fprintf( ioPUN , "%.3e\t", Fe2LevelPop[LevPop[i]-1]/denominator );
		}
		fprintf( ioPUN , "\n");
	}

	return;
}

/*
 *=====================================================================
 */
#if 0
void FeIIPun1Pop( 
	/* the io unit where the print should be directed */
	FILE * ioPUN , 
	/* the physical (not c) number of the level */
	long int nPUN )
{
	DEBUG_ENTRY( "FeIIPun1Pop()" );

	ASSERT( nPUN > 0 );
	/* need real assert to keep lint happy */
	assert( ioPUN != NULL );

	/* print the level population on ioPUN if the level was computed,
	 * print a zero if it was not, 
	 * note that nPUN is on physical scale, so test is <=*/
	if( nPUN <= FeII.nFeIILevel_local )
	{
		fprintf( ioPUN , "%e ",Fe2LevelPop[nPUN-1] );
	}
	else
	{
		fprintf( ioPUN , "%e ",0. );
	}

	return;
}
#endif

/*
 *=====================================================================
 */
STATIC int FeIIBandsCreate(
	/* chFile is optional filename, if void then use default bands,
	 * if not void then use file specified,
	 * return value is 0 for success, 1 for failure */
	 const char chFile[] )
{

	char chLine[FILENAME_PATH_LENGTH_2];
	const char* chFile1;
	FILE *ioDATA;

	bool lgEOL;
	long int i,k;

	/* keep track of whether we have been called - want to be
	 * called a total of one time */
	static bool lgCalled=false;

	DEBUG_ENTRY( "FeIIBandsCreate()" );

	/* return previous number of bands if this is second or later call*/
	if( lgCalled )
	{
		/* success */
		return 0;
	}
	lgCalled = true;

	/* use default filename if void string, else use file specified */
	chFile1 = ( strlen(chFile) == 0 ) ? "FeII_bands.ini" : chFile;

	/* get FeII band data */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " FeIICreate opening %s:", chFile1 );
	}

	ioDATA = open_data( chFile1, "r" );

	/* now count how many bands are in the file */
	nFeIIBands = 0;

	/* first line is a version number and does not count */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " FeIICreate could not read first line of %s.\n", chFile1 );
		return 1;
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
			++nFeIIBands;
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " FeIICreate could not rewind %s.\n", chFile1 );
		return 1;
	}

	FeII_Bands = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)(nFeIIBands) );

	/* now make second dim, id wavelength, and lower and upper bounds */
	for( i=0; i<nFeIIBands; ++i )
	{
		FeII_Bands[i] = (realnum *)MALLOC(sizeof(realnum)*(unsigned)(3) );
	}

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " FeIICreate could not read first line of %s.\n", chFile1 );
		return 1;
	}
	{
		i = 1;
		const long int iyr = 9, imo=6 , idy = 11;
		long iyrread, imoread , idyread;
		iyrread = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		imoread = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		idyread = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

		if(( iyrread != iyr ) ||
		  (  imoread != imo ) ||
		  (  idyread != idy ) )
		{
			fprintf( ioQQQ, 
				" PROBLEM FeIIBandsCreate: the version of %s is not the current version.\n", chFile1 );
			fprintf( ioQQQ, 
				" FeIIBandsCreate: I expected the magic numbers %li %li %li but found %li %li %li.\n", 
				iyr, imo , idy ,iyrread, imoread , idyread  );
			return 1;
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
			i = 1;
			FeII_Bands[k][0] = (realnum)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this band line 1.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				return 1;
			}
			FeII_Bands[k][1] = (realnum)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this band line 2.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				return 1;
			}
			FeII_Bands[k][2] = (realnum)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this band line 3.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				return 1;
			}
			/*fprintf(ioQQQ,
			" band data %f %f %f \n", FeII_Bands[k][0],FeII_Bands[k][1],FeII_Bands[k][2]);*/
			++k;
		}
	}
	/* now validate this incoming data */
	for( i=0; i<nFeIIBands; ++i )
	{
		/* make sure all are positive */
		if( FeII_Bands[i][0] <=0. || FeII_Bands[i][1] <=0. || FeII_Bands[i][2] <=0. )
		{
			fprintf( ioQQQ, " FeII band %li has none positive entry.\n",i );
			return 1;
		}
		/* make sure bands bounds are in correct order, shorter - longer wavelength*/
		if( FeII_Bands[i][1] >= FeII_Bands[i][2] )
		{
			fprintf( ioQQQ, " FeII band %li has improper bounds.\n" ,i);
			return 1;
		}

	}

	fclose(ioDATA);

	/* return success */
	return 0;
}
/*
 *=====================================================================
 */
void AssertFeIIDep( double *pred , double *BigError , double *StdDev )
{
	long int n;
	double arg ,
		error ,
		sum2;

	DEBUG_ENTRY( "AssertFeIIDep()" );

	if( FeII.lgSimulate || !lgFeIIEverCalled )
	{

		*pred = 0.;
		*BigError = 0.;
		*StdDev = 0.;

		return;
	}

	/* sanity check */
	ASSERT( FeII.nFeIILevel_local > 0 );

	/* find sum of deviation of departure coef from unity */
	*BigError = 0.;
	*pred = 0.;
	sum2 = 0;
	for( n=0; n<FeII.nFeIILevel_local; ++n )
	{
		/* get mean */
		*pred += Fe2DepCoef[n];

		/* error is the largest deviation from unity for any single level*/
		error = fabs( Fe2DepCoef[n] -1. );
		/* remember biggest deviation */
		*BigError = MAX2( *BigError , error );

		/* get standard deviation */
		sum2 += POW2( Fe2DepCoef[n] );
	}

	/* now get standard deviation */
	arg = sum2 - POW2( *pred ) / (double)FeII.nFeIILevel_local;
	ASSERT( (arg >= 0.) );
	*StdDev = sqrt( arg / (double)(FeII.nFeIILevel_local - 1.) );

	/* this is average value, should be unity */
	*pred /= (double)(FeII.nFeIILevel_local);

	return;
}

/* do ots rates for FeII, called by RT_OTS */
void FeII_OTS( void )
{
	long int ipLo , 
		ipHi;

	DEBUG_ENTRY( "FeII_OTS()" );

	for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
	{
		for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_local; ipHi++ )
		{
			/* >>chng 02 feb 11, skip bogus transitions */
			if( Fe2LevN[ipHi][ipLo].ipCont < 1)
				continue;

			/* ots rates, the destp prob was set in hydropesc */
			Fe2LevN[ipHi][ipLo].Emis->ots = 
				Fe2LevN[ipHi][ipLo].Emis->Aul*
				Fe2LevN[ipHi][ipLo].Hi->Pop*
				Fe2LevN[ipHi][ipLo].Emis->Pdest;

			ASSERT( Fe2LevN[ipHi][ipLo].Emis->ots >= 0. );

			/* finally dump the ots rate into the stack */
			RT_OTS_AddLine(Fe2LevN[ipHi][ipLo].Emis->ots,
				Fe2LevN[ipHi][ipLo].ipCont );
		}
	}

	return;
}

/*
 *=====================================================================
 */
/*FeII_RT_Out - do outward rates for FeII, 
 * called by RT_diffuse, which is called by cloudy */
void FeII_RT_Out(void)
{
	long int ipLo , ipHi;

	DEBUG_ENTRY( "FeIIRTOut()" );

	/* only do this if Fe+ exists */
	if( dense.xIonDense[ipIRON][1] > 0. )
	{
		/* outward line photons */
		for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
		{
			for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_local; ipHi++ )
			{
				/* >>chng 02 feb 11, skip bogus transitions */
				if( Fe2LevN[ipHi][ipLo].ipCont < 1)
					continue;
				/* >>chng 03 sep 09, change to outline, ouutlin is
				 * not exactly the same in the two - tmn missing in outline */
				Fe2LevN[ipHi][ipLo].outline_resonance();

			}
		}
	}

	return;
}

/*
 *=====================================================================
 */
STATIC void FeIIContCreate(
	/* wavelength of low-lambda end */
	double xLamLow , 
	/* wavelength of high end */
	double xLamHigh , 
	/* number of cells to break this into */
	long int ncell )
{

	double step , wl1;
	long int i;

	/* keep track of whether we have been called - want to be
	 * called a total of one time */
	static bool lgCalled=false;

	DEBUG_ENTRY( "FeIIContCreate()" );

	/* return previous number of bands if this is second or later call*/
	if( lgCalled )
	{
		/* return value of number of bands, may be used by calling program*/
		return;
	}
	lgCalled = true;

	/* how many cells will be needed to go from xLamLow to xLamHigh in ncell steps */
	nFeIIConBins = ncell;

	FeII_Cont = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)(nFeIIConBins+1) );

	/* now make second dim, id wavelength, and lower and upper bounds */
	for( i=0; i<nFeIIConBins+1; ++i )
		FeII_Cont[i] = (realnum *)MALLOC(sizeof(realnum)*(unsigned)(3) );

	step = log10( xLamHigh/xLamLow)/ncell;
	wl1 = log10( xLamLow);
	for( i=0; i<nFeIIConBins+1; ++i)
		// [n][0] are bounds of cells in Angstroms 
		FeII_Cont[i][0] = (realnum)pow(10. , ( wl1 + i*step ) );

	return;
}

/*ParseAtomFeII parse the atom FeII command */
void ParseAtomFeII( Parser &p )
{
	DEBUG_ENTRY( "ParseAtomFeII()" );

	/* turn on the large verner atom */
	FeII.lgFeIILargeOn = true;
	/* >>chng 05 nov 29, reset number of levels when called, so increased above default of 16 */
	if( lgFeIIMalloc )
	{
		/* space has been allocated, so reset number of levels to whatever it was */
		FeII.nFeIILevel_local = FeII.nFeIILevel_malloc;
	}
	else
	{
		/* space not allocated yet, set to largest possible number of levels */
		FeII.nFeIILevel_local = NFE2LEVN;
	}

	/* levels keyword is to adjust number of levels.  But this only has effect
	 * BEFORE space is allocated for the FeII arrays */
	if( p.nMatch("LEVE") )
	{
		/* do option only if space not yet allocated */
		if( !lgFeIIMalloc )
		{
			/* number of levels for hydrogen at, 2s is this plus one */
			FeII.nFeIILevel_local = (long int)p.getNumberCheck("LEVEL");
			if( FeII.nFeIILevel_local <16 )
			{
				fprintf( ioQQQ, " This would be too few levels, must have at least 16.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			else if( FeII.nFeIILevel_local > NFE2LEVN )
			{
				fprintf( ioQQQ, " This would be too many levels.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* slow keyword means do not try to avoid evaluating atom */
	else if( p.nMatch("SLOW") )
	{
		FeII.lgSlow = true;
	}

	/* redistribution keyword changes form of redistribution function */
	else if( p.nMatch("REDI") )
	{
		int ipRedis=0;
		/* there are three functions, PRD_, CRD_, and CRDW,
		 * representing partial redistribution, 
		 * complete redistribution (doppler core only, no wings)
		 * and complete with wings */
		/* partial redistribution */
		if( p.nMatch(" PRD") )
		{
			ipRedis = ipPRD;
		}
		/* complete redistribution */
		else if( p.nMatch(" CRD") )
		{
			ipRedis = ipCRD;
		}
		/* complete redistribution with wings */
		else if( p.nMatch("CRDW") )
		{
			ipRedis = ipCRDW;
		}

		/* if not SHOW option (handled below) then we have a problem */
		else if( !p.nMatch("SHOW") )
		{
			fprintf(ioQQQ," There should have been a second keyword on this command.\n");
			fprintf(ioQQQ," Options are _PRD, _CRD, CRDW (_ is space).  Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* resonance lines */
		if( p.nMatch("RESO") )
		{
			FeII.ipRedisFcnResonance = ipRedis;
		}
		/* subordinate lines */
		else if( p.nMatch("SUBO") )
		{
			FeII.ipRedisFcnSubordinate = ipRedis;
		}
		/* the show option, say what we are assuming */
		else if( p.nMatch("SHOW") )
		{
			fprintf(ioQQQ," FeII resonance lines are ");
			if( FeII.ipRedisFcnResonance ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( FeII.ipRedisFcnResonance ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( FeII.ipRedisFcnResonance ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for ipRedisFcnResonance.\n");
				TotalInsanity();
			}

			fprintf(ioQQQ," FeII subordinate lines are ");
			if( FeII.ipRedisFcnSubordinate ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( FeII.ipRedisFcnSubordinate ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( FeII.ipRedisFcnSubordinate ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for ipRedisFcnSubordinate.\n");
				TotalInsanity();
			}
		}
		else
		{
			fprintf(ioQQQ," here should have been a second keyword on this command.\n");
			fprintf(ioQQQ," Options are RESONANCE, SUBORDINATE.  Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* trace keyword - print info for each call to FeIILevelPops */
	else if( p.nMatch("TRAC") )
	{
		FeII.lgPrint = true;
	}

	/* only simulate the FeII atom, do not actually call it */
	else if( p.nMatch("SIMU") )
	{
		/* option to only simulate calls to FeIILevelPops */
		FeII.lgSimulate = true;
	}

	/* only simulate the FeII atom, do not actually call it */
	else if( p.nMatch("CONT") )
	{
		/* the continuum output with the save FeII continuum command */

		/* number of levels for hydrogen at, 2s is this plus one */
		FeII.feconwlLo = (realnum)p.getNumberCheck("low wavelength");
		FeII.feconwlHi = (realnum)p.getNumberCheck("high wavelength");
		FeII.nfe2con = (long int)p.getNumberCheck("number of intervals");
		/* check that all are positive */
		if( FeII.feconwlLo<=0. || FeII.feconwlHi<=0. || FeII.nfe2con<= 0 )
		{
			fprintf(ioQQQ," there are three numbers on the FeII continuum command, start and end wavelengths, and number of intervals.\n");
			fprintf(ioQQQ," all three must be greater than zero, sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* make sure that wl1 is less than wl2 */
		if( FeII.feconwlLo>= FeII.feconwlHi )
		{
			fprintf(ioQQQ," there are three numbers on the FeII continuum command, start and end wavelengths, and number of intervals.\n");
			fprintf(ioQQQ," the second wavelength must be greater than the first, sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	/* note no fall-through error since routine can be called with no options,
	 * to turn on the large atom */

	return;
}

void PunFeII( FILE * io )
{
	long int n, ipHi;

	DEBUG_ENTRY( "PunFeII()" );

	for( n=0; n<FeII.nFeIILevel_local-1; ++n)
	{
		for( ipHi=n+1; ipHi<FeII.nFeIILevel_local; ++ipHi )
		{
			if( Fe2LevN[ipHi][n].ipCont > 0 ) 
				fprintf(io,"%li\t%li\t%.2e\n", n , ipHi , 
					Fe2LevN[ipHi][n].Emis->TauIn );
		}
	}

	return;
}

/* include FeII lines in punched optical depths, etc, called from SaveLineStuff */
void FeIIPunchLineStuff( FILE * io , realnum xLimit  , long index)
{
	long int n, ipHi;

	DEBUG_ENTRY( "FeIIPunchLineStuff()" );

	for( n=0; n<FeII.nFeIILevel_local-1; ++n)
	{
		for( ipHi=n+1; ipHi<FeII.nFeIILevel_local; ++ipHi )
		{
			Save1Line( &Fe2LevN[ipHi][n] , io , xLimit  , index, GetDopplerWidth(dense.AtomicWeight[ipIRON]) );
		}
	}

	return;
}

/* rad pre due to FeII lines called in PresTotCurrent*/
double FeIIRadPress(void)
{

	/* will be used to check on size of opacity, was capped at this value */
	realnum smallfloat=SMALLFLOAT*10.f;
	double press,
		RadPres1;
#	undef	DEBUGFE
#	ifdef DEBUGFE
	double RadMax;
	long ipLoMax , ipHiMax;
#	endif
	long int ipLo, ipHi;

	DEBUG_ENTRY( "FeIIRadPress()" );

	press = 0.;
	if( !lgFeIIEverCalled )
		return 0.;

#	ifdef DEBUGFE
	RadMax = 0;
	ipLoMax = -1;
	ipHiMax = -1;
#	endif
	/* loop over all levels to find radiation pressure */
	for( ipHi=1; ipHi<FeII.nFeIILevel_local; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo)
		{
			/* >>chng 04 aug 27, add test on ipCont for bogus line */
			if( Fe2LevN[ipHi][ipLo].ipCont > 0  &&
				Fe2LevN[ipHi][ipLo].Hi->Pop > 1e-30 )
			{
				if( Fe2LevN[ipHi][ipLo].Hi->Pop > smallfloat &&
					Fe2LevN[ipHi][ipLo].Emis->PopOpc > smallfloat )
				{
					RadPres1 =  PressureRadiationLine( &Fe2LevN[ipHi][ipLo], GetDopplerWidth(dense.AtomicWeight[ipIRON]) );

#					ifdef DEBUGFE
					if( RadPres1 > RadMax )
					{
						RadMax = RadPres1;
						ipLoMax = ipLo;
						ipHiMax = ipHi;
					}
#					endif
					press += RadPres1;
				}
			}
		}
	}

#	ifdef	DEBUGFE
	/* option to print radiation pressure */
	if( iteration > 1 || nzone > 1558 )
	{
		fprintf(ioQQQ,"DEBUG FeIIRadPress finds P(FeII) = %.2e %.2e %li %li width %.2e\n",
			press  ,
			RadMax ,
			ipLoMax ,
			ipHiMax ,
			RT_LineWidth(&Fe2LevN[9][0],GetDopplerWidth(dense.AtomicWeight[ipIRON])) );
		DumpLine( &Fe2LevN[9][0] );
	}
#	endif
#	undef	DEBUGFE

	return press;
}

#if 0
/* internal energy of FeII */
double FeII_InterEnergy(void)
{
	double energy;

	DEBUG_ENTRY( "FeII_InterEnergy()" );

	/* There is no stack of levels, so we access all levels 
	 * uniquely by via the line stack with Fe2LevN[ipHi][0].Hi */

	energy = 0.;
	for( long ipHi=1; ipHi<FeII.nFeIILevel_local; ++ipHi )
	{
		if( Fe2LevN[ipHi][0].Hi->Pop > 1e-30 )
		{
			energy += 
				Fe2LevN[ipHi][0].Hi->Pop * Fe2LevN[ipHi][0].EnergyErg;
		}
	}

	return energy;
}
#endif

/* HP cc cannot compile this routine with any optimization */
#if defined(__HP_aCC)
#	pragma OPT_LEVEL 1
#endif
/*FeIILyaPump find rate of Lya excitation of the FeII atom */
STATIC void FeIILyaPump(void)
{

	long int ipLo ,
		ipHi;
	double EnerLyaProf2, 
	  EnerLyaProf3, 
	  EnergyWN,
	  Gup_ov_Glo, 
	  PhotOccNum_at_nu, 
	  PumpRate, 
	  de, 
	  FeIILineWidthHz, 
	  taux;

	DEBUG_ENTRY( "FeIILyaPump()" );

	/* lgLyaPumpOn is false if no Lya pumping, with no FeII command */
	/* get rates FeII atom is pumped */

	/* elsewhere in this file the dest prob hydro.dstfe2lya is defined from
	 * quantites derived here, and the resulting populations */
	if( FeII.lgLyaPumpOn )
	{

		/*************trapeze form La profile:de,EnerLyaProf1,EnerLyaProf2,EnerLyaProf3,EnerLyaProf4*************************
		 * */
		/* width of Lya in cm^-1 */
		/* HLineWidth has units of cm/s, as was evaluated in PresTotCurrent */
		/* the factor is 1/2 of E(Lya, cm^-1_/c */
		de = 1.37194e-06*hydro.HLineWidth*2.0/3.0;
		/* 82259 is energy of Lya in wavenumbers, so these are the form of the trapezoid */
		EnerLyaProf1 = 82259.0 - de*2.0;
		EnerLyaProf2 = 82259.0 - de;
		EnerLyaProf3 = 82259.0 + de;
		EnerLyaProf4 = 82259.0 + de*2.0;

		/* find Lya photon occupation number */
		if( StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop > SMALLFLOAT )
		{
			/* This is the photon occupation number at the Lya line center */
			PhotOccNumLyaCenter = 
				MAX2(0.,1.0- Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Pelec_esc - 
				Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].Emis->Pesc)/
			  (StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH1s].Pop/StatesElemNEW[ipHYDROGEN][ipHYDROGEN-ipH_LIKE][ipH2p].Pop*3. - 1.0);
		}
		else
		{
			/* lya excitation temperature not available */
			PhotOccNumLyaCenter = 0.;
		}
	}
	else
	{
		PhotOccNumLyaCenter = 0.;
		de = 0.;
		EnerLyaProf1 = 0.;
		EnerLyaProf2 = 0.;
		EnerLyaProf3 = 0.;
		EnerLyaProf4 = 0.;
	}

	/* is Lya pumping enabled?  */
	if( FeII.lgLyaPumpOn )
	{
		for( ipLo=0; ipLo < (FeII.nFeIILevel_local - 1); ipLo++ )
		{
			for( ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_local; ipHi++ )
			{
				/* on first iteration optical depth in line is inward only, on later
				 * iterations is total optical depth */
				if( iteration == 1 )
				{
					taux = Fe2LevN[ipHi][ipLo].Emis->TauIn;
				}
				else
				{
					taux = Fe2LevN[ipHi][ipLo].Emis->TauTot;
				}

				/* Gup_ov_Glo is ratio of g values */
				Gup_ov_Glo = Fe2LevN[ipHi][ipLo].Hi->g/Fe2LevN[ipHi][ipLo].Lo->g;

				/* the energy of the FeII line */
				EnergyWN = Fe2LevN[ipHi][ipLo].EnergyWN;

				if( EnergyWN >= EnerLyaProf1 && EnergyWN <= EnerLyaProf4  &&  taux > 1 )
				{
					/* this branch, line is within the Lya profile */

					/*
					 * Lya source function, at peak is PhotOccNumLyaCenter,
					 *
					 *     Prof2    Prof3
					 *       ----------
					 *      /          \
					 *     /            \
					 *    /              \
					 *  ======================
					 * Prof1              Prof4
					 *
					 */

					if( EnergyWN < EnerLyaProf2 )
					{
						/* linear interpolation on edge of trapazoid */
						PhotOccNum_at_nu = PhotOccNumLyaCenter*(EnergyWN - EnerLyaProf1)/ de;
					}
					else if( EnergyWN < EnerLyaProf3 )
					{
						/* this is the central plateau */
						PhotOccNum_at_nu = PhotOccNumLyaCenter;
					}
					else
					{
						/* linear interpolation on edge of trapazoid */
						PhotOccNum_at_nu = PhotOccNumLyaCenter*(EnerLyaProf4 - EnergyWN)/de;
					}

					/* at this point Lya source function at FeII line energy is defined, but
					 * we need to multiply by line width in Hz,
					 * >>refer	fe2	pump rate	Netzer, H., Elitzur, M., & Ferland, G.J., 1985, ApJ, 299, 752-762*/

					/** \todo 2 change this number to speed of light. */
					/* width of FeII line in Hz  */
					FeIILineWidthHz = 1.e8/(EnergyWN*299792.5)*sqrt(log(taux))*GetDopplerWidth(dense.AtomicWeight[ipIRON]);

					/* final Lya pumping rate, s^-1*/
					PumpRate = FeIILineWidthHz * PhotOccNum_at_nu * Fe2LevN[ipHi][ipLo].Emis->Aul*
					  powi(82259.0f/EnergyWN,3);
					/* above must be bogus, use just occ num times A */
					PumpRate = Fe2LevN[ipHi][ipLo].Emis->Aul* PhotOccNum_at_nu;

					/* Lya pumping rate from ipHi to lower n */
					Fe2LPump[ipHi][ipLo] += PumpRate;

					/* Lya pumping rate from n to upper ipHi */
					Fe2LPump[ipLo][ipHi] += PumpRate*Gup_ov_Glo;
				}
			}
		}
	}

	return;
}

/* end work around bugs in HP compiler */
#if defined(__HP_aCC)
#pragma OPTIMIZE OFF
#pragma OPTIMIZE ON
#endif
