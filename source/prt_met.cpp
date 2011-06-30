/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*prtmet print all line optical depths at end of iteration */
#include "cddefines.h"
#include "taulines.h"
#include "h2.h"
#include "iso.h"
#include "lines_service.h"
#include "dense.h"
#include "prt.h"
#include "mole.h"

void prtmet(void)
{
	long int i,
		nelem , 
		ipHi , 
		ipLo , 
		ipISO;

	DEBUG_ENTRY( "prtmet()" );

	/* default is to not print optical depths, turn on with
	 * print optical depths on command */
	if( prt.lgPrtTau )
	{
		fprintf( ioQQQ, "                                                    Line Optical Depths\n");

		/* "IN" - initialize */
		prme("IN",&TauLines[0]);

		/* iso sequences */
		for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] )
				{
					/* print Lyman, Balmer, Paschen, etc sequence optical depths */	
					for( ipLo=0; ipLo < iso.numLevels_local[ipISO][nelem]-1; ipLo++ )
					{
						for( ipHi=ipLo+1; ipHi < iso.numLevels_local[ipISO][nelem]; ipHi++ )
						{
							prme(" c",&Transitions[ipISO][nelem][ipHi][ipLo]);
						}
					}
				}
			}
		}

		/* print main lines optical depths */
		for( i=1; i <= nLevel1; i++ )
		{
			prme(" c",&TauLines[i]);
		}

		for( i=0; i < nWindLine; i++ )
		{
			if( TauLine2[i].Hi->IonStg < TauLine2[i].Hi->nelem+1-NISO )
			{
				prme(" c",&TauLine2[i]);
			}
		}

		for( i=0; i < nUTA; i++ )
		{
			prme(" c",&UTALines[i]);
		}

		/* print H2 line optical depths */
		H2_Prt_line_tau();

		for( i=0; i < nHFLines; i++ )
		{
			prme(" c",&HFLines[i]);
		}

		/* data base lines */
		for( i=0; i < linesAdded2; i++)
		{
			prme("DB",dBaseLines[i].tran);
		}

		fprintf( ioQQQ, "\n");
	}
	return;
}

/* prme - print optical depth */
void prme(
  /* flag saying what to do 
   * "IN" initialize
   * " c" add to list of old style lines
   * "DB" add to list of database lines
   */
  const char *chDoIt, 
  transition * t)
{
	char chAtMolWL[20],chAtMol[35];
	static long int n ;

	DEBUG_ENTRY( "prme()" );

	if( t->ipCont <= 0 )
	{
		/* line is not transferred */
		return;
	}

	/* "In" is to initialize for new printout */
	if( strncmp(chDoIt,"IN",2) == 0 )
	{
		n = 0;
	}

	else if( strncmp(chDoIt,"DB",2) == 0)
	{
		/* database lines, - cannot now address species labels for atoms
		 * and molecules in one simple way so separate from most lines
		 * print optical depth if greater than lower limit, or significantly negative */
		if( t->Emis->TauIn > prt.PrtTauFnt || t->Emis->TauIn < -1e-5 )
		{
			sprt_wl(chAtMolWL,t->WLAng);
			strcpy(chAtMol,t->Hi->chLabel);
			strcat(chAtMol," ");
			strcat(chAtMol,chAtMolWL);
			fprintf( ioQQQ, "  %10.15s",chAtMol);
			fprintf( ioQQQ,PrintEfmt("%9.2e", t->Emis->TauIn));
			fprintf( ioQQQ, " ");
			// throw CR after printing 6 numbers
			++n;
			if(n == 6)
			{
				n = 0;
				fprintf( ioQQQ, " \n");
			}
		}
	}

	else if( strncmp(chDoIt," c",2) == 0)
	{
		/* print optical depth if greater than lower limit, or significantly negative */
		if( t->Emis->TauIn > prt.PrtTauFnt || t->Emis->TauIn < -1e-5 )
		{
			/* PrtTauFnt is threshold for printing it */
			fprintf( ioQQQ, "  %10.10s",chLineLbl(t));
			fprintf( ioQQQ, PrintEfmt("%9.2e", t->Emis->TauIn ));

			// throw CR after printing 6 numbers
			++n;
			if(n == 6)
			{
				n = 0;
				fprintf( ioQQQ, " \n");
			}
		}
	}
	else
		TotalInsanity();
	return;
}
