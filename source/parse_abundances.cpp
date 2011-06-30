/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAbundances parse and read in composition as set by abundances command */
#include "cddefines.h"
#include "grains.h"
#include "grainvar.h"
#include "abund.h"
#include "phycon.h"
#include "called.h"
#include "elementnames.h"
#include "input.h"
#include "parser.h"

void ParseAbundances(Parser &p)
			/* following set true by grains command,
			* so this will not set more grains if grains already on. */
{
	bool lgLog;
	long int i;
	double absav[LIMELM], 
	  chk;
	GrainPar gp;

	DEBUG_ENTRY( "ParseAbundances()" );

	if( p.nMatch("STAR") )
	{
		/* Fred Hamann's star burst galaxy mixture -- includes a number which isn't an abundance */
		abund_starburst(p);
		return;
	}

	absav[0] = p.FFmtRead();
	/* this branch at least one number on the line - an abundance */
	if( !p.lgEOL() )
	{
		absav[1] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* this branch, we found one, but not two, numbers -
			 * must be one of a few special cases */
			if( p.nMatch(" ALL") )
			{
				/* special option, all abundances will be this number */
				if( absav[0] <= 0. )
				{
					absav[0] = pow(10.,absav[0]);
				}
				for( i=1; i < LIMELM; i++ )
				{
					abund.solar[i] = (realnum)absav[0];
				}

			}
			else if( p.nMatch("OLD ") && p.nMatch("SOLA") )
			{
				i = (int)absav[0];
				/* old solar - better be the number "84" on the line */
				if( i!=84 )
				{
					fprintf( ioQQQ, 
						" The only old abundance set I have is for version 84 - %3ld was requested.  Sorry.\n", 
						i );
					cdEXIT(EXIT_FAILURE);
				}
				for( i=1; i < LIMELM; i++ )
				{
					/* these are the old abundances */
					abund.SolarSave[i] = abund.OldSolar84[i];
					abund.solar[i] = abund.OldSolar84[i];
				}
			}

			else
			{
				fprintf( ioQQQ, 
					" I did not recognize a sub-keyword - options are ALL and OLD SOLAR 84.  Sorry.\n");
				cdEXIT(EXIT_FAILURE);
			}

			/* normal return */
			return;
		}

		/* we get here if there is a second number - read in all abundances */
		for( i=2; i < abund.npSolar; i++ )
		{
			absav[i] = p.FFmtRead();
			if( p.lgEOL() )
			{
				/*   read CONTINUE line if J scanned off end before reading all abund */
				do 
				{
					p.getline();
					if( p.m_lgEOF )
					{
						fprintf( ioQQQ, " There MUST be%3ld abundances entered, there were only%3ld.  Sorry.\n", 
									abund.npSolar, i );
						cdEXIT(EXIT_FAILURE);
					}
				} while( p.isComment() ); 

				p.echo();

				if( p.strcmp("CONT") != 0 )
				{
					fprintf( ioQQQ, " There MUST be%3ld abundances entered, there were only%3ld.  Sorry.\n", 
						abund.npSolar, i );
					cdEXIT(EXIT_FAILURE);
				}
				else
				{
					absav[i] = p.FFmtRead();
					if( p.lgEOL() )
					{
						fprintf( ioQQQ, " There MUST be%3ld abundances entered, there were only%3ld.  Sorry.\n", 
							abund.npSolar, i);
						cdEXIT(EXIT_FAILURE);
					}
				}
			}
		}

		/* fell through to here after reading in N abundances for N elements
		 * check that there are no more abundances on the line - that would
		 * be an error - a typo */
		chk = p.FFmtRead();
		if( !p.lgEOL() || (chk!=0.) )
		{
			/* got another number, not lgEOL, so too many numbers */
			fprintf( ioQQQ, " There were more than %ld abundances entered\n", 
				abund.npSolar );
			fprintf( ioQQQ, " Could there have been a typo somewhere?\n" );
		}

		/* are numbers scale factors, or log of abund rel to H?? */
		lgLog = false;
		for( i=0; i < abund.npSolar; i++ )
			if( absav[i] < 0. )
				lgLog = true;

		if( lgLog )
		{
			/* entered as log of number rel to hydrogen */
			for( i=0; i < abund.npSolar; i++ )
				abund.solar[abund.ipSolar[i]-1] = (realnum)pow(10.,absav[i]);
		}
		else
		{
			/* scale factors relative to solar */
			for( i=0; i < abund.npSolar; i++ )
				abund.solar[abund.ipSolar[i]-1] *= (realnum)absav[i];
		}

		/* check whether the abundances are reasonable */
		for( i=1; i < LIMELM; i++ )
		{
			if( abund.solar[i] > 0.2 )
			{
				fprintf( ioQQQ, " Is an abundance of %.3e relative to H reasonable for %2.2s?\n", 
					abund.solar[i], elementnames.chElementSym[i] );
			}
		}
		return;
	}

	/* following set of clauses -  no numbers on the line - 
	 * use one of several stored abundance sets */
	if( p.nMatch(" AGB") || p.nMatch("AGB ") || p.nMatch("PLAN") )
	{
		/* AGB/planetary nebula abundances */
		/* only turn on grains if "no grains" is not present */
		if( !p.nMatch("NO GR") )
		{
			/* turn on grains if not already done */
			if( !p.m_lgDSet )
			{
				/* return bins allocated by previous abundances ... commands */
				gv.clear();
				/* now allocate new grain bins */
				gp.dep = 1.;
				// standard dust to gas ratio
				gp.nDustFunc = DF_STANDARD;
				gp.lgRequestQHeating = true;
				gp.lgForbidQHeating = false;
				gp.lgGreyGrain = false;

				/* NO QHEAT option to turn off quantum heating for grains */
				if( p.nMatch("NO QH") ) 
				{
					gp.lgForbidQHeating = true;
					gp.lgRequestQHeating = false;
					phycon.lgPhysOK = false;
				}

				/* actually set up the grains */
				mie_read_opc("graphite_ism_10.opc",gp);
				mie_read_opc("silicate_ism_10.opc",gp);
			}
		}

		for( i=0; i < LIMELM; i++ )
		{
			abund.solar[i] = abund.apn[i];
			if( !abund.lgElmONapn[i] )
			{
				/* turn off elements - do this way to make sure,
				* that iso sequence limits are handled properly */
				char chDUMMY[INPUT_LINE_LENGTH];
				sprintf(chDUMMY,"element %s off ", elementnames.chElementName[i] );
				p.setline(chDUMMY);				
				ParseElement( p );
			}
		}
	}

	else if( p.nMatch("CAME") )
	{
		/* mix from Cameron 1982, "Essays on Nuclear Astrophysics" */
		/* now turn off the heavy elements */
		for( i=0; i < LIMELM; i++ )
			abund.solar[i] = abund.camern[i];
	}

	else if( p.nMatch("CRAB") )
	{
		// Crab Nebula "typical" abundances
		for( i=0; i < LIMELM; i++ )
		{
			abund.solar[i] = abund.aCrab[i];
			if( !abund.lgElmONaCrab[i] )
			{
				/* turn off elements - do this way to make sure,
				* that iso sequence limits are handled properly */
				char chDUMMY[INPUT_LINE_LENGTH];
				sprintf(chDUMMY,"element %s off ", elementnames.chElementName[i] );
				p.setline(chDUMMY);
				ParseElement( p );
			}
		}
	}
	else if( p.nMatch("HII ") || p.nMatch("H II") || p.nMatch("ORIO") )
	{
		/* H II region abundances - turn on grains by default
		 * "no grains" to not do this */
		if( !p.nMatch("NO GR") )
		{
			/* option to turn on grains */
			if( !p.m_lgDSet )
			{
				/* return bins allocated by previous abundances ... commands */
				gv.clear();
				/* now allocate new grain bins */
				gp.dep = 1.;
				// standard dust to gas ratio
				gp.nDustFunc = DF_STANDARD;
				gp.lgRequestQHeating = true;
				gp.lgForbidQHeating = false;
				gp.lgGreyGrain = false;

				/* NO QHEAT option to turn off quantum heating for grains */
				if( p.nMatch("NO QH") ) 
				{
					gp.lgForbidQHeating = true;
					gp.lgRequestQHeating = false;
					phycon.lgPhysOK = false;
				}
				/* This scales the Orion grain abundance so that the observed 
				 * dust to gas ratio that Cloudy predicts is in agreement with
				 * that observed in the Veil,
				 *>>refer	grain	Abel, N., Brogan, C., Ferland, G., O'Dell, C.R., 
				 *>>refercon	Shaw, G., Troland, T., 2004, ApJ, submitted */
				gp.dep *= 0.85;

				mie_read_opc("graphite_orion_10.opc",gp);
				mie_read_opc("silicate_orion_10.opc",gp);
			}
		}

		for( i=0; i < LIMELM; i++ )
		{
			abund.solar[i] = abund.ahii[i];
			if( !abund.lgElmONahii[i] )
			{
				/* turn off elements - do this way to make sure,
				* that iso sequence limits are handled properly */
				char chDUMMY[INPUT_LINE_LENGTH];
				sprintf(chDUMMY,"element %s off ", elementnames.chElementName[i] );
				p.setline(chDUMMY);
				ParseElement( p );
			}
		}
	}

	else if( p.nMatch("ISM ") || p.nMatch(" ISM") )
	{
		/* ISM abundances from Cowie and Songaila Ann Rev '86 */
		/* only turn on grains if "no grains" is not present */
		if( !p.nMatch("NO GR") )
		{
			if( !p.m_lgDSet )
			{
				/* return bins allocated by previous abundances ... commands */
				gv.clear();
				/* now allocate new grain bins */
				gp.dep = 1.;
				// standard dust to gas ratio
				gp.nDustFunc = DF_STANDARD;
				gp.lgRequestQHeating = true;
				gp.lgForbidQHeating = false;
				gp.lgGreyGrain = false;

				/* NO QHEAT option to turn off quantum heating */
				if( p.nMatch("NO QH") ) 
				{
					gp.lgForbidQHeating = true;
					gp.lgRequestQHeating = false;
					phycon.lgPhysOK = false;
				}

				/* actually set up the grains */
				mie_read_opc("graphite_ism_10.opc",gp);
				mie_read_opc("silicate_ism_10.opc",gp);
			}
		}

		for( i=0; i < LIMELM; i++ )
		{
			abund.solar[i] = abund.aism[i];
			if( !abund.lgElmONaism[i] )
			{
				/* turn off elements - do this way to make sure,
				 * that iso sequence limits are handled properly */
				char chDUMMY[INPUT_LINE_LENGTH];
				sprintf(chDUMMY,"element %s off ", elementnames.chElementName[i] );
				p.setline(chDUMMY);
				ParseElement( p );
			}
		}
	}

	else if( p.nMatch("NOVA") )
	{
		/* Nova Cyg abundances */
		for( i=0; i < LIMELM; i++ )
			abund.solar[i] = abund.anova[i];
	}

	else if( p.nMatch("PRIM") )
	{
		/* roughly primordial abundances: He/H=.072 */
		for( i=0; i < 4; i++ )
		{
			abund.solar[i] = abund.aprim[i];
		}

		/* now turn off the heavy elements */
		for( i=4; i < LIMELM; i++ )
		{
			/* turn off heavy elements - do this way to make sure,
			 * that H-like and He-like level limits are handled properly */
			char chDUMMY[INPUT_LINE_LENGTH];
			sprintf(chDUMMY,"element %s off ", elementnames.chElementName[i] );
			p.setline(chDUMMY);
			ParseElement( p );
		}
	}

	else
	{
		fprintf( ioQQQ, 
			" ABUNDances must have PLAN, H II, CAMERON, CRAB, NOVA, ALL, STARBURST, OLD SOLAR 84 or PRIMORDIAL.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
