/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTable parse the table read command */
/*lines_table invoked by table lines command, check if we can find all lines in a given list */
/*read_hm05 read in the data files and interpolate the continuum to
 * the correct redshift */
#include "cddefines.h"
#include "cddrive.h"
#include "physconst.h"
#include "optimize.h"
#include "rfield.h"
#include "trace.h"
#include "radius.h"
#include "input.h"
#include "stars.h"
#include "lines.h"
#include "prt.h"
#include "parser.h"
#include "save.h"
#include "thirdparty.h"
#include "continuum.h"
/* HP cc cannot compile this routine with any optimization */
#if defined(__HP_aCC)
#	pragma OPT_LEVEL 1
#endif

/*ReadTable called by TABLE READ to read in continuum from PUNCH TRANSMITTED CONTINUUM */
STATIC void ReadTable(const char * fnam);

/* these will become the label and wl for a possible list of lines,
 * obtained when tables lines used */
static char **chLabel;
static realnum *wl;
static long int nLINE_TABLE = 0;
static char chLINE_LIST[FILENAME_PATH_LENGTH];

/* tables of various built in continue */
/* crab nebula */
const int NCRAB = 10;
static double tnucrb[NCRAB], 
  fnucrb[NCRAB];

/* Bob Rubin's corrected theta 1 Ori C continuum */
const int NRUBIN = 56;
double tnurbn[NRUBIN] = {1.05E-08,1.05E-07,1.05E-06,1.04E-05,1.00E-04,1.00E-03,1.00E-02,3.01E-02,1.00E-01,
	1.50E-01,2.50E-01,4.01E-01,6.01E-01,9.8E-01,9.96E-01,1.00E+00,1.02445,1.07266,1.12563,1.18411,1.23881,
	1.29328,1.35881,1.42463,1.48981,1.55326,1.6166,1.68845,1.76698,1.8019,1.808,1.84567,1.9317,2.04891,2.14533,
	2.19702,2.27941,2.37438,2.43137,2.49798,2.56113,2.59762,2.66597,2.80543,2.95069,3.02911,3.11182,3.22178,
	3.3155,3.42768,3.50678,3.56157,3.61811,3.75211,3.89643,4.}, 
	fnurbn[NRUBIN] = {1.56E-20,1.55E-17,1.54E-14,1.53E-11,1.35E-08,1.34E-05,1.35E-02,3.62E-01,1.27E+01,
	3.90E+01,1.48E+02,4.52E+02,1.02E+03,2.27E+03,8.69E+02,8.04E+02,6.58E+02,6.13E+02,6.49E+02,6.06E+02,
	6.30E+02,5.53E+02,5.55E+02,5.24E+02,4.86E+02,4.49E+02,4.42E+02,3.82E+02,3.91E+02,2.91E+02,2.61E+02,
	1.32E+02,1.20E+02,1.16E+02,9.56E+01,9.94E+01,9.10E+01,4.85E+01,7.53E+01,9.53E+01,8.52E+01,5.76E+01,
	6.72E+01,5.20E+01,8.09E+00,3.75E+00,5.57E+00,3.80E+00,2.73E+00,2.22E+00,3.16E-01,1.26E-01,1.39E-01,
	6.15E-02,3.22E-02,7.98E-03};

/* stores numbers for table cooling flow */
const int NCFL = 40;
static double cfle[NCFL], 
  cflf[NCFL];

/* Brad Peterson's AKN 120 continuum */
const int NKN120 = 11;
static double tnuakn[NKN120], 
  fnuakn[NKN120];

/* Black's ISM continuum, with He hole filled in */
const int NISM = 23;
static double tnuism[NISM], 
  fnuism[NISM];

/* z=2 background,
 * >>refer	continuum	background	Haardt, Francesco, & Madau, Piero, 1996, 
 * >>refercon	ApJ, 461, 20 */
const int NHM96 = 14;
/* log energy in Ryd */
static double tnuHM96[NHM96]={-8,-1.722735683,-0.351545683,-0.222905683,-0.133385683,
/* changeg these two energies to prevent degeneracy */
-0.127655683,-0.004575683,0.297544317,0.476753,0.476756,0.588704317,
0.661374317,1.500814317,2.245164317};
/*-0.127655683,-0.004575683,0.297544317,0.476754317,0.476754317,0.588704317,*/
/*log J in the units of (erg cm^{-2} s^{-1} Hz^{-1} sr^{-1})*/
static double fnuHM96[NHM96]={-32.53342863,-19.9789,-20.4204,-20.4443,-20.5756,-20.7546,
-21.2796,-21.6256,-21.8404,-21.4823,-22.2102,-22.9263,-23.32,-24.2865};

/* Mathews and Ferland generic AGN continuum */
const int NAGN = 8;
static double tnuagn[NAGN], 
  tslagn[NAGN];

/* table Draine ISM continuum */
const int NDRAINE = 15;
static double tnudrn[NDRAINE] , tsldrn[NDRAINE];

/* routine that stores values for above vectors */
STATIC void ZeroContin(void);

/*read_hm05 read in the data files and interpolate the Haardt & Madau continuum to
 * the correct redshift */
STATIC void read_hm05( const char chFile[] , double **tnuHM , double **fnuHM , 
	long int *nhm , double redshift )
{

	FILE *ioFILE;
	double *table_redshifts = NULL, 
		*table_wl = NULL , 
		**table_fn=NULL,
		frac_hi;
	char chLine[1000];
	long int nRedshift , i , n , nz , ipLo , ipHi;
	bool lgEOL;

	DEBUG_ENTRY( "read_hm05()" );

	ioFILE = open_data( chFile, "r", AS_LOCAL_DATA );

	/* this will be the number of continuum points in their table */
	*nhm = 0;
	/* the number of redshifts in their table */
	nRedshift = -1;
	/* first read past comments and get magic number */
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioFILE ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
		{
			++*nhm;
			/* check magic number if first valid line */
			if( *nhm==1 )
			{
				long mag_read;
				/* this is the magic number - read it in */
				i = 1;
				mag_read = (long)FFmtRead(chLine,&i,(int)sizeof(chLine),&lgEOL);
				if( mag_read != 50808 )
				{
					fprintf(ioQQQ,
						" Magic number in Haardt & Madau file is not correct, I expected 50808 and found %li\n",
						mag_read );
					cdEXIT(EXIT_FAILURE);
				}
			}
			/* second valid line count number of redshifts */
			else if( *nhm==2 )
			{
				double a;
				i = 2;
				lgEOL = false;
				nRedshift = 0;
				while( !lgEOL )
				{
					++nRedshift;
					a = FFmtRead(chLine,&i,(int)sizeof(chLine),&lgEOL);
					ASSERT( a >= 0. );
				}
				/* have over counted by one - back up one */
				--nRedshift;
				/* highest redshift data missing in file i received */
				--nRedshift;
				/*fprintf(ioQQQ," number of z is %li\n", nRedshift);*/
				/* malloc some space */
				table_redshifts = (double*)MALLOC( (size_t)nRedshift*sizeof(double) );
				/* now read in the redshifts */
				i = 2;
				for( n=0; n<nRedshift; ++n )
				{
					table_redshifts[n] = FFmtRead(chLine,&i,(int)sizeof(chLine),&lgEOL);
					/*fprintf(ioQQQ,"DEBUG %li z %.3f\n", n ,  table_redshifts[n] );*/
				}
				if( lgEOL )
					TotalInsanity();
			}
		}
	}

	/* malloc space for the arrays first wavelength array */
	table_wl = (double*)MALLOC( (size_t)*nhm*sizeof(double) );
	/* the spectrum array table_fn[z][nu] */
	table_fn = (double**)MALLOC( (size_t)nRedshift*sizeof(double*) );
	for(n=0; n<nRedshift; ++n )
	{
		table_fn[n] = (double*)MALLOC( (size_t)(*nhm)*sizeof(double) );
	}

	/* rewind the file so we can read it a second time*/
	if( fseek( ioFILE , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " read_hm05 could not rewind HM05 date file.\n");
		cdEXIT(EXIT_FAILURE);
	}
	n = 0;
	*nhm = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioFILE ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
		{
			/* this is number of non-comment lines - will skip magic number
			 * and line giving redshift */
			++n;
			if( n>2 )
			{
				i = 1;
				table_wl[*nhm] = FFmtRead(chLine,&i,(int)sizeof(chLine),&lgEOL);
				if( lgEOL )
					TotalInsanity();
				for( nz=0; nz<nRedshift; ++nz )
				{
					/* >>chng 07 feb 18, PvH change from last branck to first */
					table_fn[nz][*nhm] = FFmtRead(chLine,&i,(int)sizeof(chLine),&lgEOL);
				}
				++*nhm;
			}
		}
	}

	fclose( ioFILE );

	{
		/* change following to true to print their original table */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			fprintf(ioQQQ,"wavelength/z");
			for(nz=0; nz<nRedshift; ++nz )
				fprintf(ioQQQ,"\t%.3f", table_redshifts[nz] );
			fprintf(ioQQQ,"\n");
			for( i=0; i<*nhm; ++i )
			{
				fprintf(ioQQQ,"%.3e", table_wl[i] );
				for( nz=0; nz<nRedshift; ++nz )
					fprintf(ioQQQ,"\t%.3e", table_fn[nz][i] );
				fprintf(ioQQQ,"\n");
			}
		}
	}

	/* this is just to shut up the lint and must not be ASSERT */
	assert( table_redshifts!=NULL );

	/* first check that desired redshift is within range of table */
	if( redshift < table_redshifts[0] || 
		redshift > table_redshifts[nRedshift-1] )
	{
		fprintf(ioQQQ," The redshift requested on table HM05 is %g but is not within the range of the table, which goes from %g to %g.\n",
			redshift, table_redshifts[0] , table_redshifts[nRedshift-1] );
		fprintf(ioQQQ," Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	ipLo = -1;
	ipHi = -1;
	/* find which redshift bin we need */
	for( nz=0; nz<nRedshift-1; ++nz )
	{
		if( redshift >= table_redshifts[nz] &&
			redshift <= table_redshifts[nz+1] )
		{
			ipLo = nz;
			ipHi = nz+1;
			break;
		}
	}
	ASSERT( ipLo>=0 && ipHi >=0 );

	/* make sure that the wavelengths are in increasing order - they were not in the
	 * original data table, but had repeated values near important ionization edges */
	for( i=0; i<*nhm-1; ++i )
	{
		if( table_wl[i]>=table_wl[i+1] )
		{
			fprintf(ioQQQ," The wavelengths in the table HM05 data table are not in increasing order.  This is required.\n");
			fprintf(ioQQQ," Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* malloc the space for the returned arrays */
	*fnuHM = (double *)MALLOC( (size_t)(*nhm)*sizeof(double ) );
	*tnuHM = (double *)MALLOC( (size_t)(*nhm)*sizeof(double ) );

	/* now fill in the arrays with the interpolated table 
	 * we will interpolate linearly in redshift 
	 * in general the redshift will be between the tabulated redshift */
	frac_hi = (redshift-table_redshifts[ipLo]) / (table_redshifts[ipHi]-table_redshifts[ipLo]);
	for( i=0; i<*nhm; ++i )
	{
		/* we have wavelengths in angstroms and want log Ryd 
		 * original table was in decreasing energy order, must also
		 * flip it around since need increasing energy order */
		(*tnuHM)[*nhm-1-i] = RYDLAM / table_wl[i];
		/* original table in correct units so no conversion needed */
		(*fnuHM)[*nhm-1-i] = table_fn[ipLo][i]*(1.-frac_hi) + table_fn[ipHi][i]*frac_hi;
		/*fprintf(ioQQQ,"DEBUG1\t%.3e\t%.3e\n",(*tnuHM)[*nhm-1-i] , (*fnuHM)[*nhm-1-i] );*/
	}

	for( n=0; n<nRedshift; ++n )
		free( table_fn[n] );
	free( table_fn );
	free( table_wl );
	free( table_redshifts );
	return;
}

void ParseTable(Parser &p)
{
	char chFile[FILENAME_PATH_LENGTH_2];	/*file name for table read */

	IntMode imode = IM_ILLEGAL_MODE;
	bool lgHit, 
	  lgLogSet;
	static bool lgCalled=false;

	long int i, 
	  j, 
	  nstar,
	  ipNORM;

	double alpha, 
	  brakmm, 
	  brakxr, 
	  ConBreak, 
	  fac, 
	  scale, 
	  slopir, 
	  slopxr;

	bool lgNoContinuum = false,
		lgQuoteFound;

	DEBUG_ENTRY( "ParseTable()" );

	/* if first call then set up values for table */
	if( !lgCalled )
	{
		ZeroContin();
		lgCalled = true;
	}

	if( rfield.nShape >= LIMSPC )
	{
		fprintf( ioQQQ, " %ld is too many spectra entered.  Increase LIMSPC\n Sorry.\n", 
		  rfield.nShape );
		cdEXIT(EXIT_FAILURE);
	}

	/* three commands, tables line, read, and star, have quotes on the
	 * lines giving file names.  must get quotes first so that filename
	 * does not confuse parser */
	lgQuoteFound = true;
	if( p.GetQuote( chFile , false ) )
		lgQuoteFound = false;

	/* set flag telling interpolate */
	strcpy( rfield.chSpType[rfield.nShape], "INTER" );
	/* >>chng 06 jul 16, fix memory leak when optimizing, PvH */
	if( !rfield.lgContMalloc[rfield.nShape] )
	{
		rfield.tNu[rfield.nShape].resize( NCELL );
		rfield.tslop[rfield.nShape].resize( NCELL );
		rfield.tFluxLog[rfield.nShape].resize( NCELL );
		rfield.ncont[rfield.nShape] = 0;
		rfield.lgContMalloc[rfield.nShape] = true;
	}

	/* NB when adding more keys also change the comment at the end */
	if( p.nMatch(" AGN") )
	{
		/* do Mathews and Ferland generic AGN continuum */
		ASSERT( NAGN < NCELL);
		for( i=0; i < NAGN; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(tnuagn[i]);
			rfield.tslop[rfield.nShape][i] = (realnum)tslagn[i];
		}
		rfield.ncont[rfield.nShape] = NAGN;

		/* optional keyword break, to adjust IR cutoff */
		if( p.nMatch("BREA") )
		{
			ConBreak = p.FFmtRead();

			if( p.lgEOL() )
			{
				/* no break, set to low energy limit of code */
				if( p.nMatch(" NO ") )
				{
					ConBreak = rfield.emm*1.01f;
				}
				else
				{
					fprintf( ioQQQ, " There must be a number for the break.\n Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}

			if( ConBreak == 0. )
			{
				fprintf( ioQQQ, " The break must be greater than 0.2 Ryd.\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			if( p.nMatch("MICR") )
			{
				/*  optional keyword, ``microns'', convert to Rydbergs */
				ConBreak = 0.0912/ConBreak;
			}

			if( ConBreak < 0. )
			{
				/*  option to enter break as LOG10 */
				ConBreak = pow(10.,ConBreak);
			}

			else if( ConBreak == 0. )
			{
				fprintf( ioQQQ, " An energy of 0 is not allowed.\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			if( ConBreak >= rfield.tNu[rfield.nShape][2].Ryd() )
			{
				fprintf( ioQQQ, " The energy of the break cannot be greater than%10.2e Ryd.\n Sorry.\n", 
				  rfield.tNu[rfield.nShape][2].Ryd() );
				cdEXIT(EXIT_FAILURE);
			}

			else if( ConBreak <= rfield.tNu[rfield.nShape][0].Ryd() )
			{
				fprintf( ioQQQ, " The energy of the break cannot be less than%10.2e Ryd.\n Sorry.\n", 
				  rfield.tNu[rfield.nShape][0].Ryd() );
				cdEXIT(EXIT_FAILURE);
			}

			rfield.tNu[rfield.nShape][1].set(ConBreak);

			rfield.tslop[rfield.nShape][1] = 
				(realnum)(rfield.tslop[rfield.nShape][2] + 
			  log10(rfield.tNu[rfield.nShape][2].Ryd()/rfield.tNu[rfield.nShape][1].Ryd()));

			rfield.tslop[rfield.nShape][0] = 
				(realnum)(rfield.tslop[rfield.nShape][1] - 
			  2.5*log10(rfield.tNu[rfield.nShape][1].Ryd()/rfield.tNu[rfield.nShape][0].Ryd()));
		}
	}

	else if( p.nMatchErase("AKN1") )
	{
		/* AKN120 continuum from Brad Peterson */
		ASSERT( NKN120 < NCELL );
		for( i=0; i < NKN120; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(tnuakn[i]);
			rfield.tslop[rfield.nShape][i] = (realnum)log10(fnuakn[i]);
		}
		rfield.ncont[rfield.nShape] = NKN120;
	}

	else if( p.nMatch("CRAB") )
	{
		ASSERT( NCRAB < NCELL );
		/* Crab Nebula continuum from Davidson and Fesen 1985 Ann Rev article */
		for( i=0; i < NCRAB; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(tnucrb[i]);
			rfield.tslop[rfield.nShape][i] = (realnum)log10(fnucrb[i]);
		}
		rfield.ncont[rfield.nShape] = NCRAB;
	}

	else if( p.nMatch("COOL") )
	{
		ASSERT( NCFL < NCELL );
		/* cooling flow from Johnstone et al. 1992, MN 255, 431. */
		for( i=0; i < NCFL; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(cfle[i]);
			rfield.tslop[rfield.nShape][i] = (realnum)cflf[i];
		}
		rfield.ncont[rfield.nShape] = NCFL;
	}

	else if( p.nMatchErase("HM96") )
	{
		/* this is the old Haardt & Madau continuum, one set of points
		 * with only the quasars
		 * this command does not include the CMB - do that separately with the CMB command */
		/* set flag telling interpolate */
		strcpy( rfield.chSpType[rfield.nShape], "INTER" );

		ASSERT( NHM96 < NCELL );
		/* z=2 background,
		* >>refer	continuum	background	Haardt, Francesco, & Madau, Piero, 1996, ApJ, 461, 20 */
		for( j=0; j < NHM96; j++ )
		{
			/* frequency was stored as log of ryd */
			rfield.tNu[rfield.nShape][j].set(pow( 10. , tnuHM96[j] ));
			rfield.tslop[rfield.nShape][j] = (realnum)fnuHM96[j];
		}
		rfield.ncont[rfield.nShape] = NHM96;

		/* optional scale factor to change default intensity from their value
		 * assumed to be log if negative, and linear otherwise */
		scale = p.FFmtRead();
		if( scale > 0. )
			scale = log10(scale);

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::SYMMETRIC;

		/* this will be flux density at some frequency on the table.  the numbers
		 * are per Hz and sr so must multiply by 4 pi
		 * [2] is not special, could have been any within array*/
		rfield.range[p.m_nqh][0] = pow(10. , tnuHM96[2] )*1.0001;

		/* convert intensity HM96 give to current units of code */
		rfield.totpow[p.m_nqh] = (fnuHM96[2] + log10(PI4) + scale);

		/* set radius to very large value if not already set */
		/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
		if( !radius.lgRadiusKnown )
		{
			radius.Radius = pow(10.,radius.rdfalt);
		}
		++p.m_nqh;

		/* set radius to very large value if not already set */
		/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
		if( !radius.lgRadiusKnown )
		{
			radius.Radius = pow(10.,radius.rdfalt);
		}
	}

	else if( p.nMatchErase("HM05") )
	{
		double *tnuHM , *fnuHM;
		double redshift;
		long int nhm;
		/* the Haardt & Madau 2005 continuum, read in a table and interpolate
		 * for any redshift, background includes both quasars and galaxies
		 * >>refer	continuum	background	Haardt, Francesco, & Madau, Piero, 2005, in preparation
		 * this command does not include the CMB - do that separately with the CMB command
		 * set flag telling interpolate */
		strcpy( rfield.chSpType[rfield.nShape], "INTER" );
		if( p.nMatch("QUAS") )
		{
			/* quasar-only continuum */
			strcpy( chFile , "haardt_madau_quasar.dat" );
		}
		else
		{
			/* the default, quasar plus galaxy continuum */
			strcpy( chFile , "haardt_madau_galaxy.dat" );
		}

		/* find the redshift - must be within bounds of table, which we
		 * do not know yet */
		redshift = p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf(ioQQQ," The redshift MUST be specified on the table HM05 command.  I did not find one.\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* read in the data files and interpolate the continuum to
		 * the correct redshift */
		read_hm05( chFile , &tnuHM , &fnuHM , &nhm , redshift );
		/* now find a point near 1 Ryd where continuum may not be far too faint */
		ipNORM = -1;
		for( j=0; j < nhm-1; j++ )
		{
			/* looking for continuum frequency near 1 Ryd 
			 * does not need to be precise */
			if( tnuHM[j]<=1. && 1. <= tnuHM[j+1] )
			{
				ipNORM = j;
				break;
			}
		}
		ASSERT( ipNORM>=0 );

		/* optional scale factor to change default intensity from their value
		 * assumed to be log if negative, and linear otherwise 
		 * increment i to move past the 96 in the keyword */
		scale = p.FFmtRead();
		if( scale > 0. )
			scale = log10(scale);

		/* this will be flux density at some frequency on the table.  the numbers
		 * are per Hz and sr so must multiply by 4 pi */
		rfield.range[p.m_nqh][0] = tnuHM[ipNORM];

		/* convert intensity HM96 give to current units of code */
		rfield.totpow[p.m_nqh] = log10(fnuHM[ipNORM]) + log10(PI4) + scale;

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::SYMMETRIC;

		/* many of their values are outside the range of a 32 bit cpu and we don't
		 * want to fill array with zeros - so rescale to cont near 1 ryd */
		scale = SDIV( fnuHM[ipNORM] );
		/* their table does not extend to the low-energy limit of the code
		 * usually does not matter since CMB will take over, but there is
		 * an obvious gap at low, z = 0, redshift.  assume Rayleigh-Jeans tail
		 * and add a first continuum point */
		rfield.tNu[rfield.nShape][0].set(rfield.emm);
		rfield.tslop[rfield.nShape][0] = 
			(realnum)log10(fnuHM[0]*pow((double)(rfield.emm/tnuHM[0]),2.)/scale);
		/*fprintf(ioQQQ,"DEBUG2\t%.3e\t%.3e\n",
			rfield.tNuRyd[rfield.nShape][0]  , 
			rfield.tslop[rfield.nShape][0] );*/
		for( j=0; j < nhm; j++ )
		{
			/* frequency is already Ryd */
			rfield.tNu[rfield.nShape][j+1].set(tnuHM[j]);
			rfield.tslop[rfield.nShape][j+1] = (realnum)log10(fnuHM[j]/scale);
			/*fprintf(ioQQQ,"DEBUG2\t%.3e\t%.3e\n",
				rfield.tNuRyd[rfield.nShape][j+1]  , 
				rfield.tslop[rfield.nShape][j+1] );*/
		}
		++nhm;
		rfield.ncont[rfield.nShape] = nhm;

		/* now make sure they are in increasing order */
		for( j=0; j < nhm-1; j++ )
			ASSERT( rfield.tNu[rfield.nShape][j].Ryd() < rfield.tNu[rfield.nShape][j+1].Ryd() );

		/* set radius to very large value if not already set */
		/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
		if( !radius.lgRadiusKnown )
		{
			radius.Radius = pow(10.,radius.rdfalt);
		}
		++p.m_nqh;

		/* set radius to very large value if not already set */
		/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
		if( !radius.lgRadiusKnown )
		{
			radius.Radius = pow(10.,radius.rdfalt);
		}
		free( tnuHM );
		free( fnuHM );

	}
	else if( p.nMatch(" ISM") )
	{
		ASSERT( NISM < NCELL );
		/* local ISM radiation field from Black 1987, Interstellar Processes */
		/* >>chng 04 mar 16, rm CMB from field so that it can be used at
		 * any redshift */
		rfield.tNu[rfield.nShape][0].set(6.);
		rfield.tslop[rfield.nShape][0] = -21.21f - 6.f;
		for( i=6; i < NISM; i++ )
		{
			rfield.tNu[rfield.nShape][i-5].set(tnuism[i]);
			rfield.tslop[rfield.nShape][i-5] = (realnum)(fnuism[i] - 
			  tnuism[i]);
		}

		rfield.ncont[rfield.nShape] = NISM -5;

		/* optional scale factor to change default luminosity
		 * from observed value
		 * want final number to be log 
		 * assumed to be log if negative, and linear otherwise unless log option is present */
		scale = p.FFmtRead();
		if( scale > 0. && !p.nMatch(" LOG") )
			scale = log10(scale);

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %4ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			  p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::SYMMETRIC;

		/* this will be flux density at 1 Ryd
		 * >>chng 96 dec 18, from 1 Ryd to H mass Rydberg
		 * >>chng 97 jan 10, had HLevNIonRyd but not defined yet */
		rfield.range[p.m_nqh][0] = HIONPOT;

		/* interpolated from Black 1987 */
		rfield.totpow[p.m_nqh] = (-18.517 + scale);

		/* set radius to very large value if not already set */
		/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
		if( !radius.lgRadiusKnown )
		{
			radius.Radius = pow(10.,radius.rdfalt);
		}
		++p.m_nqh;

		if( optimize.lgVarOn )
		{
			optimize.nvarxt[optimize.nparm] = 1;
			strcpy( optimize.chVarFmt[optimize.nparm], "TABLE ISM LOG %f");
			/*  pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			/*  the scale factor */
			optimize.vparm[0][optimize.nparm] = (realnum)scale;
			optimize.vincr[optimize.nparm] = 0.2f;
			++optimize.nparm;
		}
	}
	else if( p.nMatch("DRAI") )
	{
		ASSERT( NDRAINE < NCELL );
		rfield.lgMustBlockHIon = true;
		/* local ISM radiation field from equation 23
		 *>>refer	ISM	continuum	Draine & Bertoldi 1996 */
		for( i=0; i < NDRAINE; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(tnudrn[i]);
			rfield.tslop[rfield.nShape][i] = (realnum)tsldrn[i];
		}

		rfield.ncont[rfield.nShape] = NDRAINE;

		/* optional scale factor to change default luminosity
		 * from observed value
		 * assumed to be log if negative, and linear otherwise unless log option is present */
		scale = p.FFmtRead();
		if( scale > 0. && !p.nMatch(" LOG") )
			scale = log10(scale);

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %4ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			  p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::SYMMETRIC;

		/* continuum normalization given by flux density at first point,
		 * must set energy a bit higher to make sure it is within energy bounds
		 * that results from float arithmetic */
		rfield.range[p.m_nqh][0] = tnudrn[0]*1.01;

		/* this is f_nu at this first point */
		rfield.totpow[p.m_nqh] = tsldrn[0] + scale;

		if( optimize.lgVarOn )
		{
			optimize.nvarxt[optimize.nparm] = 1;
			strcpy( optimize.chVarFmt[optimize.nparm], "TABLE DRAINE LOG %f");
			/*  pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			/*  the scale factor */
			optimize.vparm[0][optimize.nparm] = (realnum)scale;
			optimize.vincr[optimize.nparm] = 0.2f;
			++optimize.nparm;
		}

		/* set radius to very large value if not already set */
		/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
		if( !radius.lgRadiusKnown )
		{
			radius.Radius = pow(10.,radius.rdfalt);
		}
		++p.m_nqh;
	}

	else if( p.nMatch("LINE") )
	{
		/* table lines command - way to check that lines within a data
		 * file are still valid */
		static bool lgLines_Malloc = false;

		/* say that this is not a continuum command, so don't try to work with unallocated space */
		lgNoContinuum = true;

		/* this is not a continuum source - it is to read a table of lines */
		if( lgLines_Malloc )
		{
			fprintf(ioQQQ," sorry, only one table line per input stream\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* get file name within double quotes, if not present will use default
		 * return value of 1 indicates did not find double quotes on line */
		if( lgQuoteFound )
			strcpy( chLINE_LIST , chFile );
		else
			strcpy( chLINE_LIST , "LineList_BLR.dat" );

		lgLines_Malloc = true;
		/* this flag says this is not a continuum source - nearly all table XXX
		 * commands deal with continuum sources */
		rfield.ncont[rfield.nShape] = -9;

		/* actually get the lines, and malloc the space in the arrays */
		if( (nLINE_TABLE = cdGetLineList( chLINE_LIST , &chLabel , &wl) )<0 )
		{
			/* did not find file, abort */
			fprintf(ioQQQ,"\n DISASTER PROBLEM ParseTable could not find "
				"line list file %s\n", chLINE_LIST );
			fprintf(ioQQQ," Please check the spelling of the file name and that it "
				"is in either the local or data directory.\n\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("POWE") )
	{
		/* simple power law continuum between 10 micron and 50 keV
		 *  option to read in any slope for the intermediate continuum */
		alpha = p.FFmtRead();

		/* default (no number on line) is f_nu proportional nu^-1 */
		if( p.lgEOL() )
			alpha = -1.;

		/* this is low energy for code */
		rfield.tNu[rfield.nShape][0].set(rfield.emm);
		/* and the value of the flux at this point (f_nu units)*/
		rfield.tslop[rfield.nShape][0] = -5.f;

		lgLogSet = false;

		/* option to adjust sub-millimeter break */
		brakmm = p.FFmtRead();

		/* default is 10 microns */
		if( p.lgEOL() )
		{
			lgLogSet = true;
			brakmm = 9.115e-3;
		}

		else if( brakmm == 0. )
		{
			/* if second number on line is zero then set lower limit to
			 * low-energy limit of the code.  Also set linear mode,
			 * so that last number will also be linear. */
			lgLogSet = false;
			brakmm = rfield.tNu[rfield.nShape][0].Ryd()*1.001;
		}

		else if( brakmm < 0. )
		{
			/* if number is negative then this and next are logs */
			lgLogSet = true;
			brakmm = pow(10.,brakmm);
		}

		/* optional microns keyword - convert to Rydbergs */
		if( p.nMatch("MICR") )
			brakmm = RYDLAM / (1e4*brakmm);

		rfield.tNu[rfield.nShape][1].set(brakmm);

		/* check whether this is a reasonable mm break */
		if( brakmm > 1. )
			fprintf(ioQQQ,
			" Check the order of parameters on this table power law command - the low-energy break of %f Ryd seems high to me.\n",
			brakmm );

		/* this is spectral slope, in F_nu units, between the low energy limit 
		 * and the break that may have been adjusted above 
		 * this is the slope appropriate for self-absorbed synchrotron, see eq 6.54, p.190
		 *>>refer	continuum	synchrotron	Rybicki, G. B., & Lightman, A.P. 1979, 
		 *>>refercon	Radiative Processes in Astrophysics (New York: Wiley)*/
		slopir = 5./2.;

		/* now extrapolate a flux at this energy using the flux entered for
		 * the first point, and this slope */
		rfield.tslop[rfield.nShape][1] = 
		  (realnum)(rfield.tslop[rfield.nShape][0] + 
		  slopir*log10(rfield.tNu[rfield.nShape][1].Ryd()/rfield.tNu[rfield.nShape][0].Ryd()));

		/* this is energy of the high-energy limit to code */
		rfield.tNu[rfield.nShape][3].set(rfield.egamry);

		/* option to adjust hard X-ray break */
		brakxr = p.FFmtRead();

		/* default is 50 keV */
		if( p.lgEOL() )
		{
			brakxr = 3676.;
		}

		else if( brakxr == 0. )
		{
			brakxr = rfield.tNu[rfield.nShape][3].Ryd()/1.001;
		}

		else if( lgLogSet )
		{
			/* first number was negative this is a logs */
			brakxr = pow(10.,brakxr);
		}

		/* note that this second cutoff does not have the micron keyword */
		rfield.tNu[rfield.nShape][2].set(brakxr);

		/* >>chng 03 jul 19, check that upper energy is greater than lower energy,
		 * quit if this is not the case */
		if( brakmm >= brakxr )
		{
			fprintf( ioQQQ, " HELP!!  The lower energy for the power law, %f, is greater than the upper energy, %f. This is not possible.\n Sorry.\n",
				brakmm , brakxr );
			cdEXIT(EXIT_FAILURE);
		}

		/* alpha was first option on line, is slope of mid-range */
		rfield.tslop[rfield.nShape][2] = 
		  (realnum)(rfield.tslop[rfield.nShape][1] + 
		  alpha*log10(rfield.tNu[rfield.nShape][2].Ryd()/rfield.tNu[rfield.nShape][1].Ryd()));

		/* high energy range is nu^-2 */
		slopxr = -2.;

		rfield.tslop[rfield.nShape][3] = 
			(realnum)(rfield.tslop[rfield.nShape][2] + 
		  slopxr*log10(rfield.tNu[rfield.nShape][3].Ryd()/rfield.tNu[rfield.nShape][2].Ryd()));

		/* following is number of portions of continuum */
		rfield.ncont[rfield.nShape] = 4;
	}

	else if( p.nMatch("READ") )
	{
		/* set up eventual read of table of points previously punched by code 
		 * get file name within double quotes, return as null terminated string
		 * also blank out original, chCard  version of name so do not trigger */
		if( !lgQuoteFound )
		{
			fprintf( ioQQQ, " Name of file must appear on TABLE READ.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* set flag saying really just read in continuum exactly as punched */
		strcpy( rfield.chSpType[rfield.nShape], "READ " );

		ReadTable(chFile);

		/* this is flag saying continuum has been read in here */
		rfield.tslop[rfield.nShape][0] = 1.;
		rfield.tslop[rfield.nShape][1] = 0.;

		/* number of spectra shapes that have been specified */
		++rfield.nShape;
		return;
	}

	else if( p.nMatch("TLUSTY") && !p.nMatch("STAR") )
	{
		/* >>chng 04 nov 30, retired TABLE TLUSTY command, PvH */
		fprintf( ioQQQ, " The TABLE TLUSTY command is no longer supported.\n" );
		fprintf( ioQQQ, " Please use TABLE STAR TLUSTY instead. See Hazy for details.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( p.nMatch("RUBI") )
	{
		/* do Rubin modified theta 1 Ori c */
		for( i=0; i < NRUBIN; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(tnurbn[i]);
			rfield.tslop[rfield.nShape][i] = (realnum)log10(fnurbn[i] /tnurbn[i] );
		}
		rfield.ncont[rfield.nShape] = NRUBIN;
	}

	/* >>chng 06 jul 10, retired TABLE STARBURST command, PvH */

	else if( p.nMatch("STAR") )
	{
		char chMetalicity[5] = "", chODFNew[8], chVaryFlag[7] = "";
		bool lgHCa = false, lgHNi = false;
		long nval, ndim=0;
		double Tlow = -1., Thigh = -1.;
		double val[MDIM];

		/* >>chng 06 jun 22, add support for 3 and 4-dimensional grids, PvH */
		if( p.nMatchErase("1-DI") )
			ndim = 1;
		else if( p.nMatchErase("2-DI") )
			ndim = 2;
		else if( p.nMatchErase("3-DI") )
			ndim = 3;
		else if( p.nMatchErase("4-DI") )
			ndim = 4;

		if( ndim != 0 )
		{
			/* remember keyword for possible use in optimization command */
			sprintf(chVaryFlag,"%1ld-DIM",ndim);
		}

		/* time option to vary only some continua with time */
		rfield.lgTimeVary[p.m_nqh] = p.nMatch( "TIME" );
			
		static const char table[][2][10] = {
			{"Z+1.0 ", "p10"},
			{"Z+0.75", "p075"},
			{"Z+0.5 ", "p05"},
			{"Z+0.4 ", "p04"},
			{"Z+0.3 ", "p03"},
			{"Z+0.25", "p025"},
			{"Z+0.2 ", "p02"},
			{"Z+0.1 ", "p01"},
			{"Z+0.0 ", "p00"},
			{"Z-0.1 ", "m01"},
			{"Z-0.2 ", "m02"},
			{"Z-0.25", "m025"},
			{"Z-0.3 ", "m03"},
			{"Z-0.4 ", "m04"},
			{"Z-0.5 ", "m05"},
			{"Z-0.7 ", "m07"},
			{"Z-0.75", "m075"},
			{"Z-1.0 ", "m10"},
			{"Z-1.3 ", "m13"},
			{"Z-1.5 ", "m15"},
			{"Z-1.7 ", "m17"},
			{"Z-2.0 ", "m20"},
			{"Z-2.5 ", "m25"},
			{"Z-3.0 ", "m30"},
			{"Z-3.5 ", "m35"},
			{"Z-4.0 ", "m40"},
			{"Z-4.5 ", "m45"},
			{"Z-5.0 ", "m50"},
			{"Z-INF ", "m99"}
		};
	  
		strncpy( chMetalicity, "p00", sizeof(chMetalicity) );	// default
		for (i=0; i < (long)(sizeof(table)/sizeof(*table)); ++i)
		{
			if (p.nMatchErase(table[i][0]))
			{
				strncpy( chVaryFlag, table[i][0], sizeof(chVaryFlag) );
				strncpy( chMetalicity, table[i][1], sizeof(chMetalicity) );
				break;
			}
		}


		/* there are two abundance sets (solar and halo) for CoStar and Rauch H-Ca/H-Ni models.
		 * If halo keyword appears use halo, else use solar unless other metalicity was requested */
		bool lgHalo = p.nMatch("HALO");
		bool lgSolar = ( !lgHalo && strcmp( chMetalicity, "p00" ) == 0 );

		/* >>chng 05 oct 27, added support for PG1159 Rauch models, PvH */
		bool lgPG1159 = p.nMatchErase("PG1159");

		bool lgList = p.nMatch("LIST");

		if( p.nMatch("AVAI") )
		{
			AtmospheresAvail();
			cdEXIT(EXIT_SUCCESS);
		}

		/* now that all the keywords are out of the way, scan numbers from line image */
		for( nval=0; nval < MDIM; nval++ )
		{
			val[nval] = p.FFmtRead();
			if( p.lgEOL() )
				break;
		}

		if( nval == 0 && !lgList )
		{
			fprintf( ioQQQ, " The stellar temperature MUST be entered.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* option for log if less than or equal to 10 */
		/* Caution: For CoStar models this implicitly assumes that
		 * the minimum ZAMS mass and youngest age is greater than 10.
		 * As of June 1999 this is the case. However, this is not
		 * guaranteed for possible future upgrades. */
		/* match on "LOG " rather than " LOG" to avoid confusion with "LOG(G)" */
		if( ( val[0] <= 10. && !p.nMatch("LINE") ) || p.nMatch("LOG ") )
		{
			if( val[0] < log10(BIGFLOAT) )
				val[0] = pow(10.,val[0]);
			else
			{
				fprintf(ioQQQ," DISASTER the log of the temperature was specified, "
						"but the numerical value is too large.\n Sorry.\n\n");
				cdEXIT(EXIT_FAILURE );
			}
		}

		if( lgQuoteFound )
		{
			nstar = GridInterpolate(val,&nval,&ndim,chFile,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("ATLA") )
		{
			/* this sub-branch added by Kevin Volk, to read in large
			 * grid of Kurucz atmospheres */
			/* >>chng 05 nov 19, option TRACE is no longer supported, PvH */

			if( p.nMatch("ODFN") )
				strncpy( chODFNew, "_odfnew", sizeof(chODFNew) );
			else
				strncpy( chODFNew, "", sizeof(chODFNew) );

			/* >>chng 05 nov 19, add support for non-solar metalicities, PvH */
			nstar = AtlasInterpolate(val,&nval,&ndim,chMetalicity,chODFNew,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("COST") )
		{
			/* >>chng 99 apr 30, added CoStar stellar atmospheres */
			/* this subbranch reads in CoStar atmospheres, no log(g),
			 * but second parameter is age sequence, a number between 1 and 7,
			 * default is 1 */
			if( p.nMatch("INDE") )
			{
				imode = IM_COSTAR_TEFF_MODID;
				if( nval == 1 )
				{
					val[1] = 1.;
					nval++;
				}

				/* now make sure that second parameter is within allowed range -
				 * we do not have enough information at this time to verify temperature */
				if( val[1] < 1. || val[1] > 7. )
				{
					fprintf( ioQQQ, " The second number must be the id sequence number, 1 to 7.\n" );
					fprintf( ioQQQ, " reset to 1.\n" );
					val[1] = 1.;
				}
			}
			else if( p.nMatch("ZAMS") )
			{
				imode = IM_COSTAR_MZAMS_AGE;
				if( nval == 1 )
				{
					fprintf( ioQQQ, " A second number, the age of the star, must be present.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else if( p.nMatch(" AGE") )
			{
				imode = IM_COSTAR_AGE_MZAMS;
				if( nval == 1 )
				{
					fprintf( ioQQQ, " A second number, the ZAMS mass of the star, must be present.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else
			{
				if( nval == 1 )
				{
					imode = IM_COSTAR_TEFF_MODID;
					/* default is to use ZAMS models, i.e. use index 1 */
					val[1] = 1.;
					nval++;
				}
				else
				{
					/* Caution: the code in CoStarInterpolate implicitly
					 * assumes that the dependence of log(g) with age is
					 * strictly monotonic. As of June 1999 this is the case. */
					imode = IM_COSTAR_TEFF_LOGG;
				}
			}

			if( ! ( lgSolar || lgHalo ) )
			{
				fprintf( ioQQQ, " Please choose SOLAR or HALO abundances.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			nstar = CoStarInterpolate(val,&nval,&ndim,imode,lgHalo,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("KURU") )
		{
			nstar = Kurucz79Interpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("MIHA") )
		{
			nstar = MihalasInterpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("RAUC") )
		{
			if( ! ( lgSolar || lgHalo ) )
			{
				fprintf( ioQQQ, " Please choose SOLAR or HALO abundances.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* >>chng 97 aug 23, K Volk added Rauch stellar atmospheres */
			/* >>chng 05 oct 27, added support for PG1159 Rauch models, PvH */
			/* >>chng 06 jun 26, added support for H, He, H+He Rauch models, PvH */
			if( p.nMatch("H-CA") || p.nMatch(" OLD") )
			{
				lgHCa = true;
				nstar = RauchInterpolateHCa(val,&nval,&ndim,lgHalo,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("HYDR") )
			{
				nstar = RauchInterpolateHydr(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("HELI") )
			{
				nstar = RauchInterpolateHelium(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("H+HE") )
			{
				nstar = RauchInterpolateHpHe(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( lgPG1159 ) /* the keyword PG1159 was already matched and erased above */
			{
				nstar = RauchInterpolatePG1159(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else
			{
				lgHNi = true;
				nstar = RauchInterpolateHNi(val,&nval,&ndim,lgHalo,lgList,&Tlow,&Thigh);
			}
		}

		else if( p.nMatch("TLUS") )
		{
			if( p.nMatch("OBST") )
			{
				/* >>chng 09 dec 24, this sub-branch added to read in
				 * merged Tlusty O-star & B-star atmospheres, PvH */
				nstar = TlustyInterpolate(val,&nval,&ndim,TL_OBSTAR,chMetalicity,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("BSTA") )
			{
				/* >>chng 06 oct 19, this sub-branch added to read in
				 * large 2006 grid of Tlusty B-star atmospheres, PvH */
				nstar = TlustyInterpolate(val,&nval,&ndim,TL_BSTAR,chMetalicity,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("OSTA") )
			{
				/* >>chng 05 nov 21, this sub-branch added to read in
				 * large 2002 grid of Tlusty O-star atmospheres, PvH */
				nstar = TlustyInterpolate(val,&nval,&ndim,TL_OSTAR,chMetalicity,lgList,&Tlow,&Thigh);
			}
			else
			{
				fprintf( ioQQQ, " There must be a third key on TABLE STAR TLUSTY;" );
				fprintf( ioQQQ, "  the options are OBSTAR, BSTAR, OSTAR.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("WERN") )
		{
			/* this subbranch added by Kevin Volk, to read in large
			 * grid of hot PN atmospheres computed by Klaus Werner */
			nstar = WernerInterpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("WMBA") )
		{
			/* >>chng 06 jun 27, this subbranch added to read in
			 * grid of hot atmospheres computed by Pauldrach */
			nstar = WMBASICInterpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else
		{
			fprintf( ioQQQ, " There must be a second key on TABLE STAR;" );
			fprintf( ioQQQ, "  the options are ATLAS, KURUCZ, MIHALAS, RAUCH, WERNER, and WMBASIC.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* set flag saying really just read in continuum exactly as punched */
		strcpy( rfield.chSpType[rfield.nShape], "VOLK " );

		rfield.ncont[rfield.nShape] = nstar;

		/* vary option */
		if( optimize.lgVarOn )
		{
			optimize.vincr[optimize.nparm] = (realnum)0.1;

			if( lgQuoteFound )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = nval;
				sprintf( optimize.chVarFmt[optimize.nparm], "TABLE STAR \"%s\"", chFile );
				strcat( optimize.chVarFmt[optimize.nparm], " %f LOG" );
				for( i=1; i < nval; i++ )
					strcat( optimize.chVarFmt[optimize.nparm], " %f" );
			}

			if( p.nMatch("ATLA") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = ndim;
				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR ATLAS " );
				strcat( optimize.chVarFmt[optimize.nparm], chVaryFlag );
				if( p.nMatch("ODFN") )
					strcat( optimize.chVarFmt[optimize.nparm], " ODFNEW" );

				strcat( optimize.chVarFmt[optimize.nparm], " %f LOG %f" );

				if( ndim == 3 )
					strcat( optimize.chVarFmt[optimize.nparm], " %f" );

			}

			else if( p.nMatch("COST") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 2;

				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR COSTAR " );
				if( lgHalo )
					strcat( optimize.chVarFmt[optimize.nparm], "HALO " );

				if( imode == IM_COSTAR_TEFF_MODID )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "%f LOG , INDEX %f" );
				}
				else if( imode == IM_COSTAR_TEFF_LOGG )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "%f LOG %f" );
				}
				else if( imode == IM_COSTAR_MZAMS_AGE )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "ZAMS %f LOG %f" );
				}
				else if( imode == IM_COSTAR_AGE_MZAMS )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "AGE %f LOG %f" );
					optimize.vincr[optimize.nparm] = (realnum)0.5;
				}
			}

			else if( p.nMatch("KURU") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR KURUCZ %f LOG" );
			}

			else if( p.nMatch("MIHA") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR MIHALAS %f LOG" );
			}

			else if( p.nMatch("RAUC") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = ndim;

				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR RAUCH " );

				if( p.nMatch("HYDR") )
					strcat( optimize.chVarFmt[optimize.nparm], "HYDR " );
				else if( p.nMatch("HELI") )
					strcat( optimize.chVarFmt[optimize.nparm], "HELIUM " );
				else if( p.nMatch("H+HE") )
					strcat( optimize.chVarFmt[optimize.nparm], "H+HE " );
				else if( lgPG1159 )
					strcat( optimize.chVarFmt[optimize.nparm], "PG1159 " );
				else if( lgHCa )
					strcat( optimize.chVarFmt[optimize.nparm], "H-CA " );

				if( ( lgHCa || lgHNi ) && ndim == 2 )
				{
					if( lgHalo )
						strcat( optimize.chVarFmt[optimize.nparm], "HALO " );
					else
						strcat( optimize.chVarFmt[optimize.nparm], "SOLAR " );
				}

				strcat( optimize.chVarFmt[optimize.nparm], "%f LOG %f" );

				if( ndim == 3 )
				{
					if( p.nMatch("H+HE") )
						strcat( optimize.chVarFmt[optimize.nparm], " %f" );
					else
						strcat( optimize.chVarFmt[optimize.nparm], " %f 3-DIM" );
				}
			}

			if( p.nMatch("TLUS") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = ndim;
				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR TLUSTY " );
				if( p.nMatch("OBST") )
					strcat( optimize.chVarFmt[optimize.nparm], "OBSTAR " );
				else if( p.nMatch("BSTA") )
					strcat( optimize.chVarFmt[optimize.nparm], "BSTAR " );
				else if( p.nMatch("OSTA") )
					strcat( optimize.chVarFmt[optimize.nparm], "OSTAR " );
				else
					TotalInsanity();
				strcat( optimize.chVarFmt[optimize.nparm], chVaryFlag );
				strcat( optimize.chVarFmt[optimize.nparm], " %f LOG %f" );

				if( ndim == 3 )
					strcat( optimize.chVarFmt[optimize.nparm], " %f" );
			}

			else if( p.nMatch("WERN") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 2;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR WERNER %f LOG %f" );
			}

			else if( p.nMatch("WMBA") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 3;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR WMBASIC %f LOG %f %f" );
			}

			if( rfield.lgTimeVary[p.m_nqh] )
				strcat( optimize.chVarFmt[optimize.nparm], " TIME" );

			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;

			ASSERT( nval <= LIMEXT );

			/* log of temp will be pointer */
			optimize.vparm[0][optimize.nparm] = (realnum)log10(val[0]);
			for( i=1; i < nval; i++ )
				optimize.vparm[i][optimize.nparm] = (realnum)val[i];

			/* following are upper and lower limits to temperature range */
			optimize.varang[optimize.nparm][0] = (realnum)log10(Tlow);
			optimize.varang[optimize.nparm][1] = (realnum)log10(Thigh);

			/* finally increment this counter */
			++optimize.nparm;
		}
	}

	else
	{
		fprintf( ioQQQ, " There MUST be a keyword on this line. The keys are:"
			"AGN, AKN120, CRAB, COOL, DRAINE, HM96, HM05, _ISM, LINE, POWERlaw, "
			"READ, RUBIN, and STAR.\n  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* table star Werner and table star atlas are special
	 * cases put in by Kevin Volk - they are not really tables
	 * at all, so return if chSpType is Volk */
	if( strcmp(rfield.chSpType[rfield.nShape],"VOLK ") == 0 )
	{ 
		++rfield.nShape;
		return;
	}

	/* this flag set true if we did not parse a continuum source, thus creating
	 * the arrays that are needed - return at this point */
	if( lgNoContinuum )
	{
		return;
	}

	/* ncont only positive when real continuum entered 
	 * convert from log(Hz) to Ryd if first nu>5 */
	if( rfield.tNu[rfield.nShape][0].Ryd() >= 5. )
	{
		for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
		{
			rfield.tNu[rfield.nShape][i].set( 
				pow(10.,rfield.tNu[rfield.nShape][i].Ryd() - 15.51718));
		}
	}

	else if( rfield.tNu[rfield.nShape][0].Ryd() < 0. )
	{
		/* energies entered as logs */
		for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
		{
			rfield.tNu[rfield.nShape][i].set( 
				(realnum)pow(10.,rfield.tNu[rfield.nShape][i].Ryd()));
		}
	}

	/* tFluxLog will be log(fnu) at that spot, read into tslop */
	for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
	{
		rfield.tFluxLog[rfield.nShape][i] = rfield.tslop[rfield.nShape][i];
	}

	for( i=0; i < rfield.ncont[rfield.nShape]-1; i++ )
	{
		rfield.tslop[rfield.nShape][i] = 
			(realnum)((rfield.tslop[rfield.nShape][i+1] - 
		  rfield.tslop[rfield.nShape][i])/
		  log10(rfield.tNu[rfield.nShape][i+1].Ryd()/
		  rfield.tNu[rfield.nShape][i].Ryd()));
	}

	if( rfield.ncont[rfield.nShape] > 1 && (rfield.ncont[rfield.nShape] + 1 < rfield.nupper) )
	{
		/* zero out remainder of array */
		for( i=rfield.ncont[rfield.nShape]; i < rfield.nupper; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(0.);
		}
	}

	if( trace.lgConBug && trace.lgTrace )
	{
		fprintf( ioQQQ, " Table for this continuum; TNU, TFAC, TSLOP\n" );
		for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
		{
			fprintf( ioQQQ, "%12.4e %12.4e %12.4e\n", rfield.tNu[rfield.nShape][i].Ryd(), 
			  rfield.tFluxLog[rfield.nShape][i], rfield.tslop[rfield.nShape][i] );
		}
	}

	/* renormalize continuum so that quantity passes through unity at 1 Ryd
	 * lgHit will be set false when we get a hit in following loop,
	 * then test made to make sure it happened*/
	lgHit = false;
	/*following will be reset when hit occurs*/
	fac = -DBL_MAX;
	/* >>chng 04 mar 16, chng loop so breaks when hit, previously had cycled
	 * until last point reached, so last point always used */
	for( i=0; i < rfield.ncont[rfield.nShape] && !lgHit; i++ )
	{
		if( rfield.tNu[rfield.nShape][i].Ryd() > 1. )
		{
			fac = rfield.tFluxLog[rfield.nShape][i];
			lgHit = true;
		}
	}

	if( rfield.ncont[rfield.nShape] > 1 && lgHit )
	{
		/* do the renormalization */
		for( i=0; i < rfield.ncont[rfield.nShape] ; i++ )
		{
			rfield.tFluxLog[rfield.nShape][i] -= (realnum)fac;
		}
	}

	if( rfield.ncont[rfield.nShape] > 1 )
		++rfield.nShape;
	return;
}


/* this allows the low energy point of any built in array to be reset to the
 * current low energy point in the code - nothing need be done if this is reset 
 * tnu is array of energies, [0] is first, and we want it to be lower
 * fluxlog is flux at tnu, and may or may not be log 
 * lgLog says whether it is */
STATIC void resetBltin( double *tnu , double *fluxlog , bool lgLog )
{
	/* this will multiply low-energy bounds of code and go into element[0] 
	 * ensures that energy range is fully covered */
	const double RESETFACTOR = 0.98;
	double power;
	/* this makes sure we are called after emm is defined */
	ASSERT( rfield.emm  > 0. );

	if( lgLog )
	{
		/* continuum comes in as log of flux */
		/* this is current power-law slope of low-energy continuum */
		power = (fluxlog[1] - fluxlog[0] ) / log10( tnu[1]/tnu[0] );
		/* this will be new low energy bounds to this continuum */
		tnu[0] = rfield.emm*RESETFACTOR;
		fluxlog[0] = fluxlog[1] + power * log10( tnu[0] /tnu[1] );
	}
	else
	{
		/* continuum comes in as linear flux */
		/* this is current power-law slope of low-energy continuum */
		power = log10( fluxlog[1]/fluxlog[0]) / log10( tnu[1]/tnu[0] );
		/* this will be new low energy bounds to this continuum */
		tnu[0] = rfield.emm*RESETFACTOR;
		fluxlog[0] = log10(fluxlog[1]) + power * log10( tnu[0] /tnu[1] );
		/* flux is not really log, we want linear */
		fluxlog[0] = pow(10. , fluxlog[0]);
	}
	/*fprintf(ioQQQ," power is %f lgLog is %i\n", power, lgLog );*/
	return;
}

/*ZeroContin store sets of built in continua */
STATIC void ZeroContin(void)
{

	DEBUG_ENTRY( "ZeroContin()" );

	/* Draine 1978 ISM continuum */
	/* freq in ryd */
	tnudrn[0] = 0.3676;
	tnudrn[1] = 0.4144;
	tnudrn[2] = 0.4558;
	tnudrn[3] = 0.5064;
	tnudrn[4] = 0.5698;
	tnudrn[5] = 0.6511;
	tnudrn[6] = 0.7012;
	tnudrn[7] = 0.7597;
	tnudrn[8] = 0.8220;
	tnudrn[9] = 0.9116;
	tnudrn[10] = 0.9120;
	tnudrn[11] = 0.9306;
	tnudrn[12] = 0.9600;
	tnudrn[13] = 0.9806;
	/* >>chng 05 aug 03, this energy is so close to 1 ryd that it spills over
	 * into the H-ionizing continuum - move down to 0.99 */
	/* >>chng 05 aug 08, this destabilized pdr_leiden_hack_v4 (!) so put back to
	 * original value and include extinguish command */
	tnudrn[14] = 0.9999;
	/*tnudrn[14] = 0.99;*/

	/* f_nu from equation 23 of 
	 * >>refer	ism	field	Draine, B.T. & Bertoldi, F. 1996, ApJ, 468, 269 */
	int i;
	i= 0;
	tsldrn[i] = -17.8063;
	++i;tsldrn[i] = -17.7575;
	++i;tsldrn[i] = -17.7268;
	++i;tsldrn[i] = -17.7036;
	++i;tsldrn[i] = -17.6953;
	++i;tsldrn[i] = -17.7182;
	++i;tsldrn[i] = -17.7524;
	++i;tsldrn[i] = -17.8154;
	++i;tsldrn[i] = -17.9176;
	++i;tsldrn[i] = -18.1675;
	++i;tsldrn[i] = -18.1690;
	++i;tsldrn[i] = -18.2477;
	++i;tsldrn[i] = -18.4075;
	++i;tsldrn[i] = -18.5624;
	++i;tsldrn[i] = -18.7722;

	/* generic AGN continuum taken from 
	 * >>refer	AGN	cont	Mathews and Ferland ApJ Dec15 '87
	 * except self absorption at 10 microns, nu**-2.5 below that */
	tnuagn[0] = 1e-5;
	tnuagn[1] = 9.12e-3;
	tnuagn[2] = .206;
	tnuagn[3] = 1.743;
	tnuagn[4] = 4.13;
	tnuagn[5] = 26.84;
	tnuagn[6] = 7353.;
	tnuagn[7] = 7.4e6;

	tslagn[0] = -3.388;
	tslagn[1] = 4.0115;
	tslagn[2] = 2.6576;
	tslagn[3] = 2.194;
	tslagn[4] = 1.819;
	tslagn[5] = -.6192;
	tslagn[6] = -2.326;
	tslagn[7] = -7.34;
	resetBltin( tnuagn , tslagn , true );

	/* Crab Nebula continuum from 
	 *>>refer	continuum	CrabNebula	Davidson, K., & Fesen, 1985, ARAA,  
	 * second vector is F_nu as seen from Earth */
	tnucrb[0] = 1.0e-5;
	tnucrb[1] = 5.2e-4;
	tnucrb[2] = 1.5e-3;
	tnucrb[3] = 0.11;
	tnucrb[4] = 0.73;
	tnucrb[5] = 7.3;
	tnucrb[6] = 73.;
	tnucrb[7] = 7300.;
	tnucrb[8] = 1.5e6;
	tnucrb[9] = 7.4e6;

	fnucrb[0] = 3.77e-21;
	fnucrb[1] = 1.38e-21;
	fnucrb[2] = 2.10e-21;
	fnucrb[3] = 4.92e-23;
	fnucrb[4] = 1.90e-23;
	fnucrb[5] = 2.24e-24;
	fnucrb[6] = 6.42e-26;
	fnucrb[7] = 4.02e-28;
	fnucrb[8] = 2.08e-31;
	fnucrb[9] = 1.66e-32;
	resetBltin( tnucrb , fnucrb , false );

	/* Bob Rubin's theta 1 Ori C continuum, modified from Kurucz to
	 * get NeIII right */
	/* >>chng 02 may 02, revise continuum while working with Bob Rubin on NII revisit */
	resetBltin( tnurbn , fnurbn , false );

	/* cooling flow continuum from Johnstone et al. MNRAS 255, 431. */
	cfle[0] = 0.0000100;
	cflf[0] = -0.8046910;
	cfle[1] = 0.7354023;
	cflf[1] = -0.8046910;
	cfle[2] = 1.4708046;
	cflf[2] = -0.7436830;
	cfle[3] = 2.2062068;
	cflf[3] = -0.6818757;
	cfle[4] = 2.9416091;
	cflf[4] = -0.7168990;
	cfle[5] = 3.6770115;
	cflf[5] = -0.8068384;
	cfle[6] = 4.4124136;
	cflf[6] = -0.6722584;
	cfle[7] = 5.1478162;
	cflf[7] = -0.7626385;
	cfle[8] = 5.8832183;
	cflf[8] = -1.0396487;
	cfle[9] = 6.6186204;
	cflf[9] = -0.7972314;
	cfle[10] = 7.3540230;
	cflf[10] = -0.9883416;
	cfle[11] = 14.7080460;
	cflf[11] = -1.1675659;
	cfle[12] = 22.0620689;
	cflf[12] = -1.1985949;
	cfle[13] = 29.4160919;
	cflf[13] = -1.2263466;
	cfle[14] = 36.7701149;
	cflf[14] = -1.2918345;
	cfle[15] = 44.1241379;
	cflf[15] = -1.3510833;
	cfle[16] = 51.4781609;
	cflf[16] = -1.2715496;
	cfle[17] = 58.8321838;
	cflf[17] = -1.1098027;
	cfle[18] = 66.1862030;
	cflf[18] = -1.4315782;
	cfle[19] = 73.5402298;
	cflf[19] = -1.1327956;
	cfle[20] = 147.080459;
	cflf[20] = -1.6869649;
	cfle[21] = 220.620681;
	cflf[21] = -2.0239367;
	cfle[22] = 294.160919;
	cflf[22] = -2.2130392;
	cfle[23] = 367.701141;
	cflf[23] = -2.3773901;
	cfle[24] = 441.241363;
	cflf[24] = -2.5326197;
	cfle[25] = 514.7816162;
	cflf[25] = -2.5292389;
	cfle[26] = 588.3218384;
	cflf[26] = -2.8230250;
	cfle[27] = 661.8620605;
	cflf[27] = -2.9502323;
	cfle[28] = 735.4022827;
	cflf[28] = -3.0774822;
	cfle[29] = 1470.8045654;
	cflf[29] = -4.2239799;
	cfle[30] = 2206.2067871;
	cflf[30] = -5.2547927;
	cfle[31] = 2941.6091309;
	cflf[31] = -6.2353640;
	cfle[32] = 3677.0114746;
	cflf[32] = -7.1898708;
	cfle[33] = 4412.4135742;
	cflf[33] = -8.1292381;
	cfle[34] = 5147.8159180;
	cflf[34] = -9.0594845;
	cfle[35] = 5883.2182617;
	cflf[35] = -9.9830370;
	cfle[36] = 6618.6206055;
	cflf[36] = -10.9028034;
	cfle[37] = 7354.0229492;
	cflf[37] = -11.8188877;
	cfle[38] = 7400.0000000;
	cflf[38] = -30.0000000;
	cfle[39] = 10000000.0000000;
	cflf[39] = -30.0000000;
	resetBltin( cfle , cflf , true );

	/* AKN120 continuum from Brad Peterson's plot
	 * second vector is F_nu*1E10 as seen from Earth */
	tnuakn[0] = 1e-5;
	tnuakn[1] = 1.9e-5;
	tnuakn[2] = 3.0e-4;
	tnuakn[3] = 2.4e-2;
	tnuakn[4] = 0.15;
	tnuakn[5] = 0.30;
	tnuakn[6] = 0.76;
	tnuakn[7] = 2.0;
	tnuakn[8] = 76.0;
	tnuakn[9] = 760.;
	tnuakn[10] = 7.4e6;
	fnuakn[0] = 1.5e-16;
	fnuakn[1] = 1.6e-16;
	fnuakn[2] = 1.4e-13;
	fnuakn[3] = 8.0e-15;
	fnuakn[4] = 1.6e-15;
	fnuakn[5] = 1.8e-15;
	fnuakn[6] = 7.1e-16;
	fnuakn[7] = 7.9e-17;
	fnuakn[8] = 1.1e-18;
	fnuakn[9] = 7.1e-20;
	fnuakn[10] = 1.3e-24;
	resetBltin( fnuakn , fnuakn , false );

	/* interstellar radiation field from Black 1987, "Interstellar Processes"
	 * table of log nu, log nu*fnu taken from his figure 2 */
	/* >>chng 99 jun 14 energy range lowered to 1e-8 ryd */
	tnuism[0] = 6.00;
	/*tnuism[0] = 9.00;*/
	tnuism[1] = 10.72;
	tnuism[2] = 11.00;
	tnuism[3] = 11.23;
	tnuism[4] = 11.47;
	tnuism[5] = 11.55;
	tnuism[6] = 11.85;
	tnuism[7] = 12.26;
	tnuism[8] = 12.54;
	tnuism[9] = 12.71;
	tnuism[10] = 13.10;
	tnuism[11] = 13.64;
	tnuism[12] = 14.14;
	tnuism[13] = 14.38;
	tnuism[14] = 14.63;
	tnuism[15] = 14.93;
	tnuism[16] = 15.08;
	tnuism[17] = 15.36;
	tnuism[18] = 15.43;
	tnuism[19] = 16.25;
	tnuism[20] = 17.09;
	tnuism[21] = 18.00;
	tnuism[22] = 23.00;
	/* these are log nu*Fnu */
	fnuism[0] = -16.708;
	/*fnuism[0] = -7.97;*/
	fnuism[1] = -2.96;
	fnuism[2] = -2.47;
	fnuism[3] = -2.09;
	fnuism[4] = -2.11;
	fnuism[5] = -2.34;
	fnuism[6] = -3.66;
	fnuism[7] = -2.72;
	fnuism[8] = -2.45;
	fnuism[9] = -2.57;
	fnuism[10] = -3.85;
	fnuism[11] = -3.34;
	fnuism[12] = -2.30;
	fnuism[13] = -1.79;
	fnuism[14] = -1.79;
	fnuism[15] = -2.34;
	fnuism[16] = -2.72;
	fnuism[17] = -2.55;
	fnuism[18] = -2.62;
	fnuism[19] = -5.68;
	fnuism[20] = -6.45;
	fnuism[21] = -6.30;
	fnuism[22] = -11.3;

	return;
}

/*lines_table invoked by table lines command, check if we can find all lines in a given list */
int lines_table(void)
{
	long int n,
		miss;

	if( !nLINE_TABLE )
		return 0;

	DEBUG_ENTRY( "lines_table()" );

	fprintf( ioQQQ , "lines_table checking lines within data table %s\n", chLINE_LIST );
		miss = 0;

	/* \todo	2 DOS carriage return on linux screws this up.  Can we overlook the CR?  Skip in cdgetlinelist? */
	for( n=0; n<nLINE_TABLE; ++n )
	{
		double relative , absolute;
		if( (cdLine( chLabel[n], wl[n] , &relative , &absolute ))<=0 )
		{
			++miss;
			fprintf(ioQQQ,"lines_table in parse_table.cpp did not find line %4s ",chLabel[n]);
			prt_wl(ioQQQ,wl[n]);
			fprintf(ioQQQ,"\n");
		}
	}
	if( miss )
	{
		/* this is so that we pick up problem automatically */
		fprintf(  ioQQQ , "  BOTCHED MONITORS!!!   Botched Monitors!!! lines_table could not find a total of %li lines\n\n", miss );
	}
	else
	{
		fprintf(  ioQQQ , "lines_table found all lines\n\n" );
	}
	return miss;

}

/*ReadTable called by TABLE READ to read in continuum from SAVE TRANSMITTED CONTINUUM */
STATIC void ReadTable(const char *fnam)
{
	char chLine[INPUT_LINE_LENGTH];
	long int i;
	FILE *io;

	DEBUG_ENTRY( "ReadTable()" );

	/* make unused array has valid zeros */
	for( i=0; i < NCELL; i++ )
	{
		rfield.tFluxLog[rfield.nShape][i] = -70.;
	}

	io = open_data( fnam, "r", AS_LOCAL_ONLY );

	string unit;
	char *last;

	/* read in first line of header and parse for units, if present */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )==NULL )
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}
	
	unit = "Ryd"; // default
	last = strchr_s(chLine,'\t');
	if (last)
	{
		*last = '\0';
		char *first = strrchr(chLine,'/');
		if (first) 
		{
			unit = string(first+1);
		}
		*last = '\t';
	};

	/* line 2: empty comment line */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )==NULL )
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}

	/* line 3: the version number */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )!=NULL )
	{
		long version;
		sscanf( chLine, "%ld", &version );
		if( version != VERSION_TRNCON )
		{
			fprintf( ioQQQ, " the input continuum file has the wrong version number.\n" );
			goto error;
		}
	}
	else
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}

	char md5sum[NMD5];
	long nflux;
	double mesh_lo, mesh_hi;
	union {
		double x;
		uint32 i[2];
	} u;

	/* line 4: the MD5 checksum */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )!=NULL )
	{
		strncpy( md5sum, chLine, NMD5 );
	}
	else
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}
	
	/* line 5: the lower limit of the energy grid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )!=NULL )
	{
		if( cpu.big_endian() )
			sscanf( chLine, "%x %x", &u.i[0], &u.i[1] );
		else
			sscanf( chLine, "%x %x", &u.i[1], &u.i[0] );
		mesh_lo = u.x;
	}
	else
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}
	
	/* line 6: the upper limit of the energy grid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )!=NULL )
	{
		if( cpu.big_endian() )
			sscanf( chLine, "%x %x", &u.i[0], &u.i[1] );
		else
			sscanf( chLine, "%x %x", &u.i[1], &u.i[0] );
		mesh_hi = u.x;
	}
	else
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}
	
	if( strncmp( md5sum, continuum.mesh_md5sum.c_str(), NMD5 ) != 0 ||
	    !fp_equal( mesh_lo, double(rfield.emm) ) ||
	    !fp_equal( mesh_hi, double(rfield.egamry) ) )
	{
		fprintf( ioQQQ, " the input continuum file has an incompatible energy grid.\n" );
		goto error;
	}

	/* line 7: the energy grid resolution scale factor */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )!=NULL )
	{
		if( cpu.big_endian() )
			sscanf( chLine, "%x %x", &u.i[0], &u.i[1] );
		else
			sscanf( chLine, "%x %x", &u.i[1], &u.i[0] );
		rfield.RSFCheck[rfield.nShape] = u.x;
	}
	else
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}
	
	/* line 8: the number of frequency grid points contained in the file */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )!=NULL )
	{
		sscanf( chLine, "%ld", &nflux );
	}
	else
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}
	
	/* line 9: empty comment line */
	if( read_whole_line( chLine , (int)sizeof(chLine) , io )==NULL )
	{
		fprintf( ioQQQ, " the input continuum file has been truncated.\n" );
		goto error;
	}

	/* now read in the file of numbers */
	i = 0;
	/* keep reading until we hit eol or run out of room in the array */
	while( (read_whole_line(chLine, (int)sizeof(chLine),io)!=NULL) && (i<NCELL) )
	{
		double help[2];
		sscanf( chLine, "%lf%lf ", &help[0], &help[1] );
		rfield.tNu[rfield.nShape][i].set(help[0],unit.c_str());
		// Convert to standard FluxLog units
		if (help[1] > 0.0)
		{
			rfield.tFluxLog[rfield.nShape][i] = (realnum) log10(help[1]/
				rfield.tNu[rfield.nShape][i].Ryd());
		}
		++i;
	}
	/* put pointer at last good value */
	rfield.ncont[rfield.nShape] = i;

	/* check that sane number of values entered */
	if( rfield.ncont[rfield.nShape] != nflux )
	{
		fprintf( ioQQQ, " the input continuum file has the wrong number of points: %ld\n", 
		  rfield.ncont[rfield.nShape] );
		goto error;
	}

	fclose(io);
	return;

error:
	fprintf( ioQQQ, " please recreate this file using the SAVE TRANSMITTED CONTINUUM command.\n" );
	cdEXIT(EXIT_FAILURE);
}

#if defined(__HP_aCC)
#pragma OPTIMIZE OFF
#pragma OPTIMIZE ON
#endif
