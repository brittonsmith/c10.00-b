/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtAllTau master routine controlling printout of optical depths at end of calculation */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "opacity.h"
#include "dense.h"
#include "colden.h"
#include "elementnames.h"
#include "geometry.h"
#include "prt.h"

void PrtAllTau(void)
{
	long int i, 
		n,
	  nelem;
	realnum fcon, 
	  flin;

	DEBUG_ENTRY( "PrtAllTau()" );

	/* optical depths used by code are total through model,
	 * when sphere is set, this is twice optical depth through
	 * computed structure */
	if( geometry.lgSphere )
	{
		fcon = 2.;
		if( geometry.lgStatic )
		{
			flin = 2.;
		}
		else
		{
			flin = 1.;
		}
	}
	else
	{
		fcon = 1.;
		flin = 1.;
	}

	/* print out optical depths and column densities */

	/* R(1300) is Rayleigh scattering */
	fprintf( ioQQQ, "\n Contin Optical Depths: COMP:");

	fprintf( ioQQQ,PrintEfmt("%9.2e", opac.telec));
	fprintf( ioQQQ, "    H-:");
	fprintf( ioQQQ,PrintEfmt("%9.2e",opac.thmin ));

	fprintf( ioQQQ, " R(1300):");
	fprintf( ioQQQ,PrintEfmt("%9.2e", colden.colden[ipCOL_H0]*6.71e-24));

	fprintf( ioQQQ, "  H2+:");
	fprintf( ioQQQ,PrintEfmt("%9.2e", colden.colden[ipCOL_H2p]*7e-18));

	fprintf( ioQQQ, "  Pfa:");
	if( iso.n_HighestResolved_local[ipH_LIKE][ipHYDROGEN] >= 5 )
	{
		long ip5p = iso.QuantumNumbers2Index[ipH_LIKE][ipHYDROGEN][5][1][2];
		ASSERT( Transitions[ipH_LIKE][ipHYDROGEN][ip5p][ipH4s].ipCont > 0 );
		PrintE82( ioQQQ , opac.TauTotalGeo[0][Transitions[ipH_LIKE][ipHYDROGEN][ip5p][ipH4s].ipCont-1]/fcon);
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}
	fprintf( ioQQQ, "\n" );

	fprintf( ioQQQ, "                          Pa:");
	/* 06 aug 28, from numLevels_max to _local. */
	if( iso.numLevels_local[ipH_LIKE][ipHYDROGEN] > ipH4p )
	{
		fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH3s].ipCont-1]/fcon));
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}

	fprintf( ioQQQ, "    Ba:");
	/* 06 aug 28, from numLevels_max to _local. */
	if( iso.numLevels_local[ipH_LIKE][ipHYDROGEN] > 3 )
	{
		fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][Transitions[ipH_LIKE][ipHYDROGEN][ipH3p][ipH2s].ipCont-1]/fcon));
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}

	fprintf( ioQQQ, "      Hb:");
	/* 06 aug 28, from numLevels_max to _local. */
	if( iso.numLevels_local[ipH_LIKE][ipHYDROGEN] > 4 )
	{
		fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][Transitions[ipH_LIKE][ipHYDROGEN][ipH4p][ipH2s].ipCont-1]/fcon));
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}

	fprintf( ioQQQ, "   La:");
	fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][Transitions[ipH_LIKE][ipHYDROGEN][ipH2p][ipH1s].ipCont-1]/fcon));

	fprintf( ioQQQ, "     1r:");
	PrintE93( ioQQQ , opac.TauTotalGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][ipHYDROGEN][ipH1s]-1]/fcon);

	fprintf( ioQQQ, "  1.8:");
	PrintE82( ioQQQ , opac.TauTotalGeo[0][iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][0]-1]/fcon);

	fprintf( ioQQQ, " 4.:");
	PrintE93( ioQQQ , opac.TauTotalGeo[0][iso.ipIsoLevNIonCon[ipH_LIKE][1][ipH1s]-1]/fcon);
	fprintf( ioQQQ, "\n");

	/* print optical depths of some metal lines */
	prtmet();

	/* print H-like H, He+ optical depths */
	for( nelem=ipHYDROGEN; nelem<=ipHELIUM; ++nelem )
	{
		/* helium may be turned off */
		if( dense.lgElmtOn[nelem] )
		{
#			define NUMB_PER_LINE	8
			fprintf( ioQQQ, "\n Old, new %s%2li continuum optical depths:\n",
				elementnames.chElementSym[nelem] ,
				nelem+1);
			/* absorption continuum optical depths are energies of the h-like ion continua
			 * loop over old, then new optical depths */
			for( i=1; i>=0; --i )
			{
				/* print ground, skip t2, then do 2p */
				for( n=ipH1s; n < iso.numLevels_max[ipH_LIKE][nelem] - iso.nCollapsed_max[ipH_LIKE][nelem]; n++ )
				{
					if( n==ipH2s )
						continue;
					if( n%NUMB_PER_LINE ==1)
						fprintf(ioQQQ,"\n");
					/* this, combined with "continue" above, ensures that we print
					 * 1 (1s), 2(tot 2), then 3 */
					fprintf( ioQQQ , "%6ld",MAX2(1,n));
					fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauAbsGeo[i][iso.ipIsoLevNIonCon[ipH_LIKE][nelem][n]-1]/fcon));
				}
				fprintf( ioQQQ, "\n" );
			}

			/* now do h-like line optical depths */
			fprintf( ioQQQ, "\n Old, new %s%2li line optical depths:\n",
				elementnames.chElementSym[nelem] ,
				nelem+1);
			/* Lya is a special case due to 2s-2p resolution - explictly print it first */
			fprintf( ioQQQ, "%3i-%2i",2, 1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauTot/flin ));
			/* total optical depth in 3-2s and 3-2p, is total of both so 2-1 is correct for 3-2*/
			/* 06 aug 28, from numLevels_max to _local. */
			for( n=3; n <= iso.n_HighestResolved_local[ipH_LIKE][nelem]; n++ )
			{
				if( n%NUMB_PER_LINE ==1)
					fprintf(ioQQQ,"\n");
				fprintf( ioQQQ, "%3ld-%2ld",n, n-1 );
				fprintf( ioQQQ,PrintEfmt("%9.2e", 
					/* just do nP - n'S */
					Transitions[ipH_LIKE][nelem][ iso.QuantumNumbers2Index[ipH_LIKE][nelem][n][1][2] ][ iso.QuantumNumbers2Index[ipH_LIKE][nelem][n-1][0][2] ].Emis->TauTot/flin ));
			}
			if( prt.lgPrnIsoCollapsed )
			{
				/* above flag set with print line iso collapsed command
				 * not done by default */
				for( n=iso.numLevels_local[ipH_LIKE][nelem] - iso.nCollapsed_local[ipH_LIKE][nelem]; n < iso.numLevels_local[ipH_LIKE][nelem]; n++ )
				{
					if( StatesElemNEW[nelem][nelem-ipH_LIKE][n].n % NUMB_PER_LINE ==1)
						fprintf(ioQQQ,"\n");
					fprintf( ioQQQ, "%3ld-%2ld", StatesElemNEW[nelem][nelem-ipH_LIKE][n].n, StatesElemNEW[nelem][nelem-ipH_LIKE][n-1].n );
					fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][n][n-1].Emis->TauTot/flin ));
				}
			}

			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "%3i-%2i",2, 1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][ipH2p][ipH1s].Emis->TauIn/flin ));
			/* 06 aug 28, from numLevels_max to _local. */
			for( n=3; n <= iso.n_HighestResolved_local[ipH_LIKE][nelem]; n++ )
			{
				if( n%NUMB_PER_LINE ==1)
					fprintf(ioQQQ,"\n");
				fprintf( ioQQQ, "%3ld-%2ld",n, n-1 );
				fprintf( ioQQQ,PrintEfmt("%9.2e", 
					/* just do nP - n'S */
					Transitions[ipH_LIKE][nelem][ iso.QuantumNumbers2Index[ipH_LIKE][nelem][n][1][2] ][ iso.QuantumNumbers2Index[ipH_LIKE][nelem][n-1][0][2] ].Emis->TauIn/flin ));
			}
			if( prt.lgPrnIsoCollapsed )
			{
				for( n=iso.numLevels_local[ipH_LIKE][nelem] - iso.nCollapsed_local[ipH_LIKE][nelem]; n < iso.numLevels_local[ipH_LIKE][nelem]; n++ )
				{
					if( StatesElemNEW[nelem][nelem-ipH_LIKE][n].n % NUMB_PER_LINE ==1)
						fprintf(ioQQQ,"\n");
					fprintf( ioQQQ, "%3ld-%2ld", StatesElemNEW[nelem][nelem-ipH_LIKE][n].n, StatesElemNEW[nelem][nelem-ipH_LIKE][n-1].n );
					fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipH_LIKE][nelem][n][n-1].Emis->TauIn/flin ));
				}
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* ================================================================================ */

	/* print helium lines if helium exists */
	if( dense.lgElmtOn[ipHELIUM] )
	{
		fprintf( ioQQQ, "\n Old He Is optical depths:" );
		for( i=0; i < 5; i++ )
		{
			fprintf( ioQQQ, "%5ld", i+1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauAbsGeo[1][iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][i]-1]/fcon) );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " New He Is optical depths:" );
		for( i=0; i < 5; i++ )
		{
			fprintf( ioQQQ, "%5ld", i+1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauAbsGeo[0][iso.ipIsoLevNIonCon[ipHE_LIKE][ipHELIUM][i]-1]/fcon ));
		}
		fprintf( ioQQQ, "\n" );

		/* ================================================================================*/

		fprintf( ioQQQ, "          Old He Is Lines:" );
		fprintf( ioQQQ, " %4d",584 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe2p1P][ipHe1s1S].Emis->TauTot/flin ));
		fprintf( ioQQQ, " %4d",3889 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe3p3P][ipHe2s3S].Emis->TauTot/flin ));
		fprintf( ioQQQ, " %4d",5016 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe3p1P][ipHe2s1S].Emis->TauTot/flin ));
		fprintf( ioQQQ, " %4d",5876 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P2].Emis->TauTot/flin ));
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "          New He Is Lines:" );
		fprintf( ioQQQ, " %4d",584 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe2p1P][ipHe1s1S].Emis->TauIn/flin ));
		fprintf( ioQQQ, " %4d",3889 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe3p3P][ipHe2s3S].Emis->TauIn/flin ));
		fprintf( ioQQQ, " %4d",5016 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe3p1P][ipHe2s1S].Emis->TauIn/flin ));
		fprintf( ioQQQ, " %4d",5876 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", Transitions[ipHE_LIKE][ipHELIUM][ipHe3d3D][ipHe2p3P2].Emis->TauIn/flin ));
		fprintf( ioQQQ, "\n" );

		/* ================================================================================*/
	}
	return;
}
