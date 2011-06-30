/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseDrive parse the drive command - drive calls to various subs */
/*DrvCaseBHS allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
/*DrvHyas allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
/*dgaunt drive gaunt factor routines by letting user query values */
#include "cddefines.h"
#include "trace.h"
#include "hydro_bauman.h"
#include "hydroeinsta.h"
#include "thirdparty.h"
#include "phycon.h"
#include "atmdat.h"  
#include "taulines.h"
#include "thermal.h"
#include "abund.h"
#include "rt.h"
#include "continuum.h"
#include "parser.h"
#include "physconst.h"

/*dgaunt drive gaunt factor routines by letting user query values */
STATIC void dgaunt(void);

/*DrvHyas allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
STATIC void DrvHyas(void);

/* drive escape probability routines */
STATIC void DrvEscP( void );

/*DrvCaseBHS allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
STATIC void DrvCaseBHS(void);

void ParseDrive(Parser &p )
{
	bool lgEOL;
	long int n, 
	  i;
	double fac, 
	  zed;

	DEBUG_ENTRY( "ParseDrive()" );

	/* NB evolve all following names to style DrvSomething */

	/* option to drive cloudy, which one? */
	if( p.nMatch("FFMT") )
	{
		/* free format parser */
		char chInput[INPUT_LINE_LENGTH];
		fprintf( ioQQQ, " FFmtRead ParseDrive entered.  Enter number.\n" );
		lgEOL = false;
		while( !lgEOL )
		{
			if( read_whole_line( chInput , (int)sizeof(chInput) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " ParseDrive.dat error getting magic number\n" );
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			fac = FFmtRead(chInput,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " FFmtRead hit the EOL with no value, return=%10.2e\n", 
				  fac );
				break;
			}
			else if( fac == 0. )
			{
				break;
			}
			else
			{
				fprintf( ioQQQ, " FFmtRead returned with value%11.4e\n", 
				  fac );
			}
			fprintf( ioQQQ, " Enter 0 to stop, or another value.\n" );
		}
		fprintf( ioQQQ, " FFmtRead ParseDrive exits.\n" );
	}

	else if( p.nMatch("CASE") )
	{
		/* option to interpolate on Hummer and Storey case b hydrogenic emission routines */
		DrvCaseBHS( );
	}

	else if( p.nMatch("CDLI") )
	{
		/* drive cdLine to check that it finds all the right lines, routine is in lines.c */
		trace.lgDrv_cdLine = true;
	}

	else if( p.nMatch(" E1 ") )
	{
		// option to drive exponential integral routines
		// first, special case given in Abramovitc & Stegan
		double tau = 1.275;
		for( i=0; i<50; ++i )
		{
			fprintf(ioQQQ,"tau\t%.3e\t exp-tau\t%.5e\t e1 tau\t%.5e  \t e2 "
				"\t%.5e \te2n %.5e \t e3\t%.5e \t e4\t%.5e \n",
				tau , sexp(tau) , ee1(tau) , e2(tau ), expn(2, tau) ,
				expn(3 , tau ), expn(4 , tau ) );
			tau = pow( 10. , ((double)i/4. - 9.) );
		}
		cdEXIT(EXIT_SUCCESS);
	}

	else if( p.nMatch("ESCA") )
	{
		/* option to drive escape probability routines */
		DrvEscP( );
	}

	else if( p.nMatch("HYAS") )
	{
		/* option to drive Jason's hydrogen transition probabilities */
		DrvHyas();
	}

	else if( p.nMatch("GAUN") )
	{
		/* drive gaunt factor routine */
		dgaunt();
	}

	else if( p.nMatch("POIN") )
	{
		/* later on, check cell pointers, centers, widths */
		fprintf( ioQQQ, " Define entire model first, later I will ask for energies.\n" );
		trace.lgPtrace = true;
	}

	else if( p.nMatch("PUMP") )
	{
		char chInput[INPUT_LINE_LENGTH];
		lgEOL = false;
		fprintf( ioQQQ, " Continuum pump ParseDrive entered - Enter log tau\n" );
		while( !lgEOL )
		{
			if( read_whole_line( chInput , (int)sizeof(chInput) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " ParseDrive.dat error getting magic number\n" );
				cdEXIT(EXIT_FAILURE);
			}
			/* get tau */
			i = 1;
			fac = FFmtRead(chInput,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL )
				break;
			fac = pow(10.,fac);
			fprintf( ioQQQ, " Tau =%11.4e\n", fac );
			TauDummy.Emis->TauIn = (realnum)fac;
			// not sure what this is supposed to do, but TauDummy.Hi->nelem and therefore the doppler width are ill-defined here!
			// just put 1 for hydrogen
			fac = DrvContPump(&TauDummy,1.f); 
			fprintf( ioQQQ, " ContPump =%11.4e\n", fac );
			fprintf( ioQQQ, " Enter null to stop, or another value.\n" );
		}
		fprintf( ioQQQ, " ContPump ParseDrive exits.\n" );
	}

	else if( p.nMatch("STAR") )
	{
		char chInput[INPUT_LINE_LENGTH];
		/* get starburst abundances */
		for( i=0; i < 40; i++ )
		{
			zed = ((double)i+1.)/4. + 0.01;
			sprintf( chInput, "starburst, zed=%10.4f", zed );
			p.setline(chInput);
			abund_starburst(p);
			fprintf( ioQQQ, "%8.1e", zed );
			for(n=0; n < LIMELM; n++)
				fprintf( ioQQQ, "%8.1e", abund.solar[n] );
			fprintf( ioQQQ, "\n" );
		}
	}

	else
	{
		fprintf( ioQQQ, 
			" Unrecognized key; keys are FFmtRead, CASE, GAUNT, HYAS, POINT, MOLE, STAR, "
			"PUMP, and ESCApe.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/*DrvEscP user queries escape probability routines, which return values */
STATIC void DrvEscP( void )
{
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL;
	long i;
	double tau;

	DEBUG_ENTRY( "DrvEscP()" );

	/* this routine is enterd with the command escape probability, and
	 * drives the escape probability routine to compare answers */
	fprintf( ioQQQ, " Enter the log of the one-sided optical depth; line with no number to stop.\n" );

	lgEOL = false;
	while( !lgEOL )
	{
		if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
		{
			break;
		}

		i = 1;
		tau = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			break;
		}

		tau = pow(10.,tau);
		fprintf( ioQQQ, "tau was %e\n", tau );
		fprintf( ioQQQ, " ESCINC=%11.3e\n", esc_PRD_1side(tau,1e-4) );
		fprintf( ioQQQ, " ESCCOM=%11.3e\n", esc_CRDwing_1side(tau,1e-4 ) );
		fprintf( ioQQQ, " ESCA0K2=%11.3e\n", esca0k2(tau) );

	}
	return;
}

/*DrvCaseBHS allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
STATIC void DrvCaseBHS(void)
{
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL,
		lgHit;
	long int i, 
	  n1, 
	  nelem ,
	  n2;
	double Temperature,
		Density;

	DEBUG_ENTRY( "DrvCaseBHS()" );

	/* this routine is entered with the command DRIVE CASEB, and
	 * interpolates on the Hummer & Storey case b data set */

	/* read in some external data files, but only if this is first call */
	fprintf(ioQQQ," I will get needed H data files. This will take a second.\n");
	atmdat_readin();

	{
		/* following should be set true to print input lines */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			double xLyman , alpha;
			long int ipHi;
			nelem = 2;
			Temperature = 2e4;
			Density = 1e2;
			for( ipHi=3; ipHi<25; ++ipHi )
			{
				double photons = (1./POW2(ipHi-1.)-1./POW2((double)ipHi) ) /(1.-1./ipHi/ipHi );
				xLyman = atmdat_HS_caseB(1,ipHi, nelem,Temperature , Density , 'A' );
				alpha = atmdat_HS_caseB(ipHi-1,ipHi, nelem,Temperature , Density , 'A' );
				fprintf(ioQQQ,"%li\t%.3e\t%.3e\n", ipHi, xLyman/alpha*photons, photons );
			}
			cdEXIT(EXIT_SUCCESS);
		}
	}

	/* first get the charge, only H and He at present */
	lgHit = false;
	nelem = 0;
	while( !lgHit )
	{
		fprintf( ioQQQ, " Enter atomic number of species, either 1(H) or 2(He).\n" );
		if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
		{
			fprintf( ioQQQ, " error getting species \n" );
			cdEXIT(EXIT_FAILURE);
		}

		i = 1;
		nelem = (long int)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL || nelem< 1 || nelem > 2 )
		{
			fprintf( ioQQQ, " must be either 1 or 2!\n" );
		}
		else
		{
			lgHit = true;
		}
	}

	fprintf(ioQQQ," In the following temperatures <10 are log, >=10 linear.\n");
	fprintf(ioQQQ," The density is always a log.\n");
	fprintf(ioQQQ," The order of the quantum numbers do not matter.\n");
	fprintf(ioQQQ," The smallest must not be smaller than 2,\n");
	fprintf(ioQQQ," and the largest must not be larger than 25.\n");
	fprintf(ioQQQ," Units of emissivity are erg cm^3 s^-1\n\n");
	fprintf(ioQQQ," The limits of the HS tables are 2 <= n <= 25.\n");

	lgHit = true;
	/* this is always true */
	while( lgHit )
	{
		fprintf( ioQQQ, " Enter 4 numbers, temperature, density, 2 quantum numbers, null line stop.\n" );
		if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
		{
			fprintf( ioQQQ, " Thanks for interpolating on the Hummer & Storey data set!\n" );
			cdEXIT(EXIT_FAILURE);
		}

		i = 1;
		Temperature = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " error getting temperature!\n" );
			break;
		}

		/* log if less than 10 */
		if( Temperature < 10. )
		{
			Temperature = pow(10., Temperature );
		}
		fprintf(ioQQQ," Temperature is %g\n", Temperature );

		Density = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " error getting density!\n" );
			break;
		}
		Density = pow(10., Density );
		fprintf(ioQQQ," Density is %g\n", Density );

		/* these quantum numbers can be in any order */
		n1 = (long)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " error getting quantum number!\n" );
			break;
		}

		n2 = (long)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " error getting quantum number!\n" );
			break;
		}

		if( MAX2( n1 , n2 ) > 25 )
		{
			fprintf( ioQQQ," The limits of the HS tables are 2 <= n <= 25.  Sorry.\n");
			break;
		}

		fprintf( ioQQQ, 
			" 4pJ(%ld,%ld)/n_e n_p=%11.3e\n", 
			 n1, n2, 
			 atmdat_HS_caseB(n1,n2, nelem,Temperature , Density , 'B' ) );

		/* this is check that we were in bounds */
#if	0
		{
			long j;
			double tempTable[33] = {
					11870.,12490.,12820.,
					11060.,17740.,12560.,
					16390.,16700.,11360.,
					10240.,20740.,12030.,
					14450.,19510.,12550.,
					16470.,16560.,12220.,
					15820.,12960.,10190.,
					12960.,14060.,12560.,
					11030.,10770.,13360.,
					10780.,11410.,11690.,
					12500.,13190.,21120. };
			double edenTable[33] = {
					10.,270.,80.,10.,70.,
					110.,200.,10.,40.,90.,
					340.,80.,60.,340.,30.,
					120.,10.,50.,450.,30.,
					180.,20.,170.,60.,20.,
					40.,30.,20.,100.,130.,
					10.,10.,110. };


			for( j=0; j<33; j++ )
			{
				double halpha, hbeta, hgamma;

				halpha = atmdat_HS_caseB(2,3, 1,tempTable[j] , edenTable[j] , 'B' );
				hbeta =  atmdat_HS_caseB(2,4, 1,tempTable[j] , edenTable[j] , 'B' );
				hgamma = atmdat_HS_caseB(2,5, 1,tempTable[j] , edenTable[j] , 'B' );

				fprintf( ioQQQ, "%e\t%e\t%e\t%e\n",
					tempTable[j], 
					edenTable[j],
					halpha/hbeta,
					hgamma/hbeta );
			}
		}
#endif
	}

	fprintf( ioQQQ, " Thanks for interpolating on the Hummer & Storey data set!\n" );
	cdEXIT(EXIT_FAILURE);

}

/*DrvHyas allow user to query hydrogen A's, asks for up, low level, gives A, drive hyas */
STATIC void DrvHyas(void)
{
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL;
	long int i, nHi, lHi, nLo, lLo;

	DEBUG_ENTRY( "DrvHyas()" );

	/* this routine is entered with the command DRIVE HYAS, and
	 * drives Jason's hydrogen einstein A routines */

	nHi = 1;
	/* nHi never lt 1 */
	while( nHi != 0 )
	{
		fprintf( ioQQQ, " Enter four quantum numbers (n, l, n', l'), null line to stop.\n" );
		if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
		{
			fprintf( ioQQQ, " error getting drvhyas \n" );
			cdEXIT(EXIT_FAILURE);
		}

		i = 1;
		nHi = (long int)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
			break;

		lHi = (long int)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " must be four numbers!\n" );
			break;
		}

		if( lHi >= nHi )
		{
			fprintf( ioQQQ, " l must be less than n!\n" );
			break;
		}

		nLo = (long int)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " must be four numbers!\n" );
			break;
		}

		lLo = (long int)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
		if( lgEOL )
		{
			fprintf( ioQQQ, " must be four numbers!\n" );
			break;
		}

		if( lLo >= nLo )
		{
			fprintf( ioQQQ, " l must be less than n!\n" );
			break;
		}

		if( nLo > nHi )
		{
			long nTemp, lTemp;

			/* swap hi and lo */
			nTemp = nLo;
			lTemp = lLo;
			nLo = nHi;
			lLo = lHi;
			nHi = nTemp;
			lHi = lTemp;
		}

		fprintf( ioQQQ, " A(%3ld,%3ld->%3ld,%3ld)=%11.3e\n", 
			nHi, lHi, nLo, lLo,
			H_Einstein_A( nHi, lHi, nLo, lLo, 1 ) );

	}
	fprintf( ioQQQ, " Driver exits, enter next line.\n" );

	return;
}

/*dgaunt drive gaunt factor routines by letting user query values */
STATIC void dgaunt(void)
{
	char chCard[INPUT_LINE_LENGTH];
	bool lgEOL;
	int inputflag;
	long int i, 
	  ierror;
	realnum enerlin[1];
	double SaveTemp;
	double z,mygaunt=0.;
	double loggamma2, logu;

	DEBUG_ENTRY( "dgaunt()" );

	SaveTemp = phycon.te;

	/* this routine is entered with the command DRIVE GAUNT, and
	 * drives the gaunt factor routine to check range
	 * */
	fprintf( ioQQQ, " Enter 0 to input temp, energy, and net charge, or 1 for u and gamma**2.\n" );
	if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
	{
		fprintf( ioQQQ, " dgaunt error getting magic number\n" );
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	inputflag = (int)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);

	if( inputflag == 0 )
	{
		fprintf( ioQQQ, " Enter the temperature (log if <=10), energy (Ryd), and net charge. Null line to stop.\n" );
		/* >>chng 96 july 07, got rid of statement labels replacing with do while
		 * */
		ierror = 0;
		while( ierror == 0 )
		{
			if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " dgaunt error getting magic number\n" );
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			phycon.alogte = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
			/* the line may be trash but ierror will pick it up  */
			if( lgEOL  )
			{
				fprintf( ioQQQ, " Gaunt driver exits, enter next line.\n" );
				break;
			}
			/* numbers less than or equal to 10 are the log of the temperature */
			double TeNew;
			if( phycon.alogte > 10. )
			{
				TeNew = phycon.alogte;
			}
			else
			{
				TeNew = pow(10.,phycon.alogte);
			}
			TempChange(TeNew , false);

			enerlin[0] = (realnum)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL || enerlin[0] == 0. )
			{
				fprintf( ioQQQ, " Sorry, but there should be two more numbers, energy and charge.\n" );
			}

			z = (double)FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL || z == 0. )
			{
				fprintf( ioQQQ, " Sorry, but there should be a third number, charge.\n" );
			}

			/* This is non-thermally averaged gaunt factors.	*/
			mygaunt = cont_gaunt_calc( (double)phycon.te, z, enerlin[0] );

			fprintf( ioQQQ, " Using my  routine, Gff= \t" );
			fprintf( ioQQQ, "%11.3e\n", mygaunt );

		}
	}
	else
	{
		/* this routine is entered with the command DRIVE GAUNT, and
		 * drives the gaunt factor routine to check range
		 * */
		fprintf( ioQQQ, " Enter log u and log gamma2. Null line to stop.\n" );
		ierror = 0;
		while( ierror == 0 )
		{
			if( read_whole_line( chCard , (int)sizeof(chCard) , ioStdin ) == NULL )
			{
				fprintf( ioQQQ, " dgaunt error getting magic number\n" );
				cdEXIT(EXIT_FAILURE);
			}
			i = 1;
			logu = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
			/* the line may be trash but ierror will pick it up  */
			if( lgEOL  )
			{
				fprintf( ioQQQ, " Gaunt driver exits, enter next line.\n" );
				break;
			}

			loggamma2 = FFmtRead(chCard,&i,INPUT_LINE_LENGTH,&lgEOL);
			if( lgEOL  )
			{
				fprintf( ioQQQ, " Sorry, but there should be another numbers, log gamma2.\n" );
			}

			/* This is my attempt to calculate non-thermally averaged gaunt factors.	*/
			mygaunt = cont_gaunt_calc( TE1RYD/pow(10.,loggamma2), 1., pow(10.,logu-loggamma2) );

			TempChange(TE1RYD/pow(10.,loggamma2) , false);

			fprintf( ioQQQ, " Using my  routine, Gaunt factor is" );
			fprintf( ioQQQ, "%11.3e\n", mygaunt );
		}
	}

	TempChange(SaveTemp , false);
	return;
}
