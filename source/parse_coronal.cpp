/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCoronal parse parameters off coronal equilibrium command */
#include "cddefines.h"
#include "rfield.h"
#include "thermal.h"
#include "input.h"
#include "optimize.h"
#include "phycon.h"
#include "radius.h"
#include "dynamics.h"
#include "parser.h"

/*ParseCoronal parse parameters off coronal equilibrium command */
void ParseCoronal(Parser &p)
{
	double a;

	DEBUG_ENTRY( "ParseCoronal()" );

	/* use coronal command to establish initial conditions in a cooling
	 * time-varying cloud */
	if( p.nMatch( "INIT"  ) && p.nMatch( "TIME"  ) )
	{
		dynamics.lg_coronal_time_init = true;
		if( p.nMatch( "TRAC"  ) )
			dynamics.lgTracePrint = true;
	}

	/* coronal equilibrium; set constant temperature to number on line */
	thermal.lgTemperatureConstant = true;
	thermal.lgTemperatureConstantCommandParsed = true;
	a = p.FFmtRead();
	if( p.lgEOL() )
	{
		fprintf( ioQQQ, " There should be a temperature on this line.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* numbers less than or equal to 10 are the log of the temperature */
	if( (a <= 10. && !p.nMatch("LINE")) || p.nMatch(" LOG") )
	{
		thermal.ConstTemp = (realnum)pow(10.,a);
	}
	else
	{
		thermal.ConstTemp = (realnum)a;
	}

	/* check temperature bounds */
	if( thermal.ConstTemp < phycon.TEMP_LIMIT_LOW )
	{
		thermal.ConstTemp = (realnum)(1.0001*phycon.TEMP_LIMIT_LOW);
		fprintf( ioQQQ, " PROBLEM Te too low, reset to %g K.\n",
			 thermal.ConstTemp );
	}
	if( thermal.ConstTemp > phycon.TEMP_LIMIT_HIGH )
	{
		thermal.ConstTemp = (realnum)(0.9999*phycon.TEMP_LIMIT_HIGH);
		fprintf( ioQQQ, " PROBLEM Te too high, reset to %g K.\n",
			 thermal.ConstTemp );
	}

	/* now simulate a LASER line */
	strcpy( rfield.chSpType[rfield.nShape], "LASER" );
	/* scan off the laser's energy ion Rydbergs */
	rfield.slope[rfield.nShape] = rfield.egamry*0.99;
	/* default width is 0.05 */
	rfield.cutoff[rfield.nShape][0] = 0.05;

	/* simulate an ionization parameter line */
	strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
	strcpy( rfield.chSpNorm[p.m_nqh], "IONI" );

	/* >>chng 96 jun 17, to stop mole network from crashing */
	/* >>chng 05 aug 15, this sets ionization parameter, in test case ism_hot_brems the
	 * value of 1e-10 was enough to dominate the ionization of he-like N - it's ionization
	 * then jumped due to large optical depth in the continuum - change U from -10 to -15 */
	/* >>chng 05 aug 16, this very strongly affected the coll_t4 sim - apparently there
	 * was a significant photoionization contribution from the -10 continuum,
	 * this was close to a 'no photoionization' case, but lower further to insure
	 * no photo contribution 
	 * chang from -15 to -20 */
	/* 
	   8-29-2012 - Britton changed the heating value to so runs with no 
	   external radiation don't crash at low temperatures.
	rfield.totpow[p.m_nqh] = -20.f;
	*/
	rfield.totpow[p.m_nqh] = -13.f;

	++rfield.nShape;
	if( rfield.nShape >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	++p.m_nqh;
	if( p.m_nqh >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* set R to large value if U specified but R is not */
	/* set radius to very large value if not already set */
	/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
	if( !radius.lgRadiusKnown )
	{
		radius.Radius = pow(10.,radius.rdfalt);
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/*  no luminosity options on vary */
		optimize.nvarxt[optimize.nparm] = 1;
		strcpy( optimize.chVarFmt[optimize.nparm], "COROnal equilibrium %f LOG" );
		if( dynamics.lg_coronal_time_init )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME INIT" );

		/*  pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;

		/*  log of temp will be pointer */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(thermal.ConstTemp);
		optimize.vincr[optimize.nparm] = 0.1f;
		optimize.varang[optimize.nparm][0] = (realnum)log10(1.00001*phycon.TEMP_LIMIT_LOW);
		optimize.varang[optimize.nparm][1] = (realnum)log10(0.99999*phycon.TEMP_LIMIT_HIGH);
		++optimize.nparm;
	}
	return;
}
