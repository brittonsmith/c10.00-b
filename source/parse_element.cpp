/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseElement parse options on element command */
#include "cddefines.h"
#include "optimize.h"
#include "abund.h"
#include "dense.h"
#include "input.h"
#include "iso.h"
#include "hmi.h"
#include "mole.h"
#include "elementnames.h"
#include "parser.h"

void ParseElement( Parser &p)
{
	/* this will be used to remember how many levels are active in any element we turn off,
	 * so that we can retain state if turned back on  */
	static bool lgFirst = true;
	static long levels[NISO][LIMELM];
	bool lgEnd, 
	  lgHIT;

	bool lgForceLog=false, lgForceLinear=false;

	long int nelem, 
	  j, 
	  k;
	double param=0.0;

	DEBUG_ENTRY( "ParseElement()" );

	/* zero out this array if first call */
	if( lgFirst )
	{
		lgFirst = false;
		for( long i=0; i<NISO; ++i )
		{
			for( nelem=0; nelem<LIMELM; ++nelem )
			{
				levels[i][nelem] = 0;
			}
		}
	}
	/* say that abundances have been changed */
	abund.lgAbnSolar = false;

	/* read in list of elements for abundances command */
	if( p.nMatch("READ") )
	{
		abund.npSolar = 0;
		for( long i=1; i < LIMELM; i++ )
		{
			p.getline();

			if( p.m_lgEOF )
			{
				fprintf( ioQQQ, " Hit EOF while reading element list; use END to end list.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* this would be a line starting with END to say end on list */
			if( p.strcmp("END") == 0 )
			{
				return;
			}

			j = 1;
			lgHIT = false;
			while( j < LIMELM && !lgHIT )
			{
				j += 1;

				if( p.strcmp( elementnames.chElementNameShort[j-1] ) == 0 )
				{
					abund.npSolar += 1;
					abund.ipSolar[abund.npSolar-1] = j;
					lgHIT = true;
				}
			}

			if( !lgHIT )
			{
				p.PrintLine( ioQQQ);
				fprintf( ioQQQ, " Sorry, but I did not recognize element name on this line.\n" );
				fprintf( ioQQQ, " Here is the list of names I recognize.\n" );
				fprintf( ioQQQ, " " );

				for( k=2; k <= LIMELM; k++ )
				{
					fprintf( ioQQQ, "%4.4s\n", elementnames.chElementNameShort[k-1] );
				}

				cdEXIT(EXIT_FAILURE);
			}
		}

		/* fell through, make sure one more line, the end line */
		p.getline();

		if( p.strcmp("END") == 0 )
		{
			return;
		}

		else
		{
			fprintf( ioQQQ, " Too many elements were entered.\n" );
			fprintf( ioQQQ, " I only know about%3d elements.\n", 
			  LIMELM );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* find which element - will be used in remainder of routine to 
	 * adjust aspects of this element */
	nelem = p.GetElem( );

	bool lgElementSet = false;
	if( nelem >= ipHYDROGEN )
	{
		lgElementSet = true;
		/* nelem is now the atomic number of the element, and must be correct */
		ASSERT( nelem>=0 && nelem < LIMELM );
	}

	/* look for log or linear keywords */
	if( p.nMatch(" LOG") )
		lgForceLog = true;
	else if( p.nMatch("LINE") )
		lgForceLinear = true;

	if( p.nMatch("SCAL") )
	{
		param = p.FFmtRead();

		/* enter abundance as scale factor, relative to what is in now */
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be a number on this line.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			/* interpret as log unless forced linear */
			if( (lgForceLog || param <= 0.) && !lgForceLinear )
			{
				/* option for log of scale factor */
				param = pow(10.,param);
			}
			abund.ScaleElement[nelem] = (realnum)param;
		}
	}

	else if( p.nMatch("ABUN") )
	{
		param = p.FFmtRead();

		/* log of actual abundance */
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be a number on this line.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			/* interpret as log unless forced linear */
			if( !lgForceLinear )
			{
				/* option for log of scale factor */
				param = pow(10.,param);
			}
			abund.solar[nelem] = (realnum)param;

			if( abund.solar[nelem] > 1. )
			{
				fprintf( ioQQQ, 
					" Please check the abundance of this element.  It seems high to me.\n" );
			}
		}
	}

	else if( p.nMatch(" OFF") )
	{
		param = p.FFmtRead();
	
		/* keyword LIMIT sets lower limit to abundance of element
		 * that will be included */
		/* option to turn off this element, set abundance to zero */
		if( p.nMatch( "LIMI") )
		{
			if( p.lgEOL() )
			{
				fprintf( ioQQQ, " There must be an abundances on the ELEMENT OFF LIMIT command.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			dense.AbundanceLimit = (realnum)pow(10., param);
		}
		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine ( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			/* specify element name */
			dense.lgElmtOn[nelem] = false;
			/* another flag that element is off */
			dense.IonHigh[nelem] = -1;
			dense.lgSetIoniz[nelem] = false;
			/* set no levels for all elements that are turned off, except for helium itself, which always exists */
			if( nelem > ipHELIUM )
			{
				levels[ipHYDROGEN][nelem] = iso.numLevels_max[ipH_LIKE][nelem];
				levels[ipHELIUM][nelem] = iso.numLevels_max[ipHE_LIKE][nelem];
				iso.numLevels_max[ipH_LIKE][nelem] = 0;
				iso.numLevels_max[ipHE_LIKE][nelem] = 0;
			}

			if( nelem == ipHYDROGEN )
			{
				fprintf( ioQQQ, " It is not possible to turn hydrogen off.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* specify an ionization distribution */
	else if( p.nMatch("IONI") )
	{
		bool lgLogSet = false;
		int ion;
		int ihi , low;


		if( !lgElementSet )
		{
			fprintf( ioQQQ, 
				" ParseElement did not find an element on the following line:\n" );
			p.PrintLine( ioQQQ );
			cdEXIT(EXIT_FAILURE);
		}
		if( !dense.lgElmtOn[nelem] )
		{
			fprintf(ioQQQ,"Sorry, you cannot set the ionization of %s since it has been turned off.\n" ,
				elementnames.chElementName[nelem] );
			cdEXIT(EXIT_FAILURE);
		}

		/* option to specify ionization distribution for this element */
		dense.lgSetIoniz[nelem] = true;
		ion = 0;
		while( ion<nelem+2 )
		{
			/* the ionization fractions that are set when above set true,
			 * gas phase abundance is this times total abundance
			 * Ionization fraction for [nelem][ion] */
			dense.SetIoniz[nelem][ion] = (realnum)p.FFmtRead();

			if( p.lgEOL() ) 
				break;

			/* all are log if any are negative */
			if( dense.SetIoniz[nelem][ion] < 0. )
				lgLogSet = true;
			++ion;
		}
		param = dense.SetIoniz[nelem][0];

		/* zero out ones not specified */
		for( long i=ion; i<nelem+2; ++i )
			dense.SetIoniz[nelem][i] = 0.;

		/* convert rest to linear if any were negative */
		if( lgLogSet )
		{
			for( long i=0; i<ion; ++i )
				dense.SetIoniz[nelem][i] = (realnum)pow((realnum)10.f, dense.SetIoniz[nelem][i]);
		}

		/* now check that no zero abundances exist between lowest and highest non-zero values */
		low = 0;
		while( dense.SetIoniz[nelem][low]==0 && low < nelem+1 )
			++low;
		ihi = nelem+1;
		while( dense.SetIoniz[nelem][ihi]==0 && ihi > low )
			--ihi;

		if( ihi==low && dense.SetIoniz[nelem][ihi]==0 )
		{
			fprintf(ioQQQ," element ionization command has all zero ionization fractions.  This is not possible.\n Sorry\n");
			cdEXIT(EXIT_FAILURE);
		}
		realnum AbundSum=0.;
		for(long i=low; i<=ihi; ++i )
		{
			if( dense.SetIoniz[nelem][i]==0 )
			{
				fprintf(ioQQQ," element abundance command has zero abundance between positive values.  This is not possible.\n Sorry\n");
				cdEXIT(EXIT_FAILURE);
			}
			AbundSum += dense.SetIoniz[nelem][i];
		}

		// renormalize so that sum is unity - needed for various sum rules
		for(long i=low; i<=ihi; ++i )
			dense.SetIoniz[nelem][i] /= AbundSum;

		/* >>chng 05 mar 13, add this
		 * finally turn off any part of the chemical networks that might use this element,
		 * since chemistry must come into equilibrium with the ionization, and vice vese -
		 * ion distribution cannot come into equilibrium with chemistry when ions are set,
		 * so much turn off chemistry */
		if( nelem==ipHYDROGEN && !hmi.lgNoH2Mole )
		{
			fprintf(ioQQQ," Hydrogen ionization has been set, H chemistry is disabled.\n");
			/* turn off only H2 */
			hmi.lgNoH2Mole = true;
		}
		/* loop over all atoms in co chem net - these occur as the last NUM_ELEMENTS 
		 * elements in the full set of mole.num_heavy_molec */
		/* >>chng 06 mar 19, do not print if network already turned off */
		if( mole.lgElem_in_chemistry[nelem] && !co.lgNoCOMole )
		{
			fprintf(ioQQQ," The ionization of an element in the CO chemistry network has been set, CO chemistry is disabled.\n");
			/* turn off CO network */
			co.lgNoCOMole = true;
		}
	}

	/* option to turn element on - some ini files turn off elements,
	 * this can turn them back on */
	else if( p.nMatch(" ON ") )
	{

		if( !lgElementSet )
		{
			fprintf( ioQQQ, 
				" ParseElement did not find an element on the following line:\n" );
			p.PrintLine( ioQQQ );
			cdEXIT(EXIT_FAILURE);
		}
		/* option to turn off this element, set abundance to zero */
		dense.lgElmtOn[nelem] = true;
		/* reset levels to default if they were ever turned off with element off command */
		if( levels[ipHYDROGEN][nelem] )
		{
			iso.numLevels_max[ipH_LIKE][nelem] = levels[ipHYDROGEN][nelem];
			iso.numLevels_max[ipHE_LIKE][nelem] = levels[ipHELIUM][nelem];
		}
	}

	else if( p.nMatch("TABL") )
	{

		if( !lgElementSet )
		{
			fprintf( ioQQQ, 
				" ParseElement did not find an element on the following line:\n" );
			p.PrintLine( ioQQQ );
			cdEXIT(EXIT_FAILURE);
		}
		/* read in table of position-dependent abundances for a particular element. */
		abund.lgAbunTabl[nelem] = true;

		/* general flag saying this option turned on */
		abund.lgAbTaON = true;

		/* does the table give depth or radius?  keyword DEPTH
		 * says it is depth, default is radius */
		if( p.nMatch("DEPT") )
			abund.lgAbTaDepth[nelem] = true;
		else
			abund.lgAbTaDepth[nelem] = false;

		/* make sure not trying to change hydrogen */
		if( nelem == ipHYDROGEN )
		{
			fprintf( ioQQQ, " cannot change abundance of hydrogen.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* read pair giving radius and abundance */
		p.getline();	  
		abund.AbTabRad[0][nelem] = (realnum)p.FFmtRead();
		abund.AbTabFac[0][nelem] = (realnum)p.FFmtRead();

		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " no pairs entered - cannot interpolate\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* number of points in the table */
		abund.nAbunTabl = 2;

		lgEnd = false;
		/* LIMTAB is limit to number of points we can store and is
		 * set to 500 in abundances */
		while( !lgEnd && abund.nAbunTabl < LIMTABD )
		{
			p.getline();
			lgEnd = p.m_lgEOF;
			if( !lgEnd )
			{
				/* convert first 4 char to caps, into chCap */
				if( p.strcmp("END") == 0 )
					lgEnd = true;
			}

			/* lgEnd may have been set within above if, if end line encountered*/
			if( !lgEnd )
			{
				abund.AbTabRad[abund.nAbunTabl-1][nelem] = 
					(realnum)p.FFmtRead();

				abund.AbTabFac[abund.nAbunTabl-1][nelem] = 
					(realnum)p.FFmtRead();
				abund.nAbunTabl += 1;
			}
		}

		abund.nAbunTabl -= 1;

		/* now check that radii are in increasing order */
		for( long i=1; i < abund.nAbunTabl; i++ )
		{
			/* the radius values are assumed to be strictly increasing */
			if( abund.AbTabRad[i][nelem] <= abund.AbTabRad[i-1][nelem] )
			{
				fprintf( ioQQQ, "ParseElement: TABLE ELEMENT TABLE radii "
					"must be in increasing order\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else
	{
		fprintf( ioQQQ, "ParseElement: ELEMENT command - there must be "
			"a keyword on this line.\n" );
		fprintf( ioQQQ, " The keys I know about are TABLE, SCALE, _OFF, "
			"_ON_, IONIZATION, and ABUNDANCE.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		if( p.nMatch("SCAL") )
		{
			/* vary scale factor */
			sprintf( optimize.chVarFmt[optimize.nparm], "ELEMENT %4.4s SCALE %%f LOG", 
			  elementnames.chElementNameShort[nelem] );

		}
		else if( p.nMatch("ABUN") )
		{
			/* vary absolute abundance */
			sprintf( optimize.chVarFmt[optimize.nparm], "ELEMENT %4.4s ABUNDANCE %%f LOG", 
			  elementnames.chElementNameShort[nelem] );
		}
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		if( param >=0. )
			optimize.vparm[0][optimize.nparm] = (realnum)log10(param);
		else
			optimize.vparm[0][optimize.nparm] = (realnum)param;
		optimize.vincr[optimize.nparm] = 0.2f;
		optimize.nvarxt[optimize.nparm] = 1;
		++optimize.nparm;
	}
	return;
}
