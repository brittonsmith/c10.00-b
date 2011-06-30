/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAtomISO parse information from the atom XX-like command line */
#include "cddefines.h"
#include "elementnames.h"
#include "optimize.h"
#include "hydrogenic.h"
#include "input.h"
#include "iso.h"
#include "parser.h"
#include "phycon.h"
#include "rfield.h"
#include "taulines.h"
#include "thirdparty.h"

/*ParseAtomISO parse parameters off the XX-like command */
void ParseAtomISO(long ipISO, Parser &p )
{
	long int 
		numLevels,
		nelem;

	DEBUG_ENTRY( "ParseAtomISO()" );

	/* look for the name of an element - if we don't find one do the entire
	 * iso sequence - returns negative number if element not found */
	nelem = p.GetElem( );

	/* H-like Helium is not possible */
	if( ipISO==ipHE_LIKE && nelem==ipHYDROGEN )
	{
		fprintf(ioQQQ," Sorry, H-like He does not exist.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* collisions - don't pick up the levels collapsed command */
	if( p.nMatch("COLL") && !p.nMatch("LEVE"  ) )
	{
		/* option to turn collisions off, all are
		 * set to 1 in zero.  command can accept only one option at a time */
		if( p.nMatch("EXCI") )
		{
			/* turn off collisional excitation */
			iso.lgColl_excite[ipISO] = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("IONI") )
		{
			/* turn off collisional ionization */
			iso.lgColl_ionize[ipISO] = false;
			phycon.lgPhysOK = false;
		}

		else if( p.nMatch("2S2P") || ( p.nMatch("2P2S") && ipISO == ipH_LIKE ) )
		{
			/* >>chng 02 feb 07, change from 2s2p to l-mixing */
			/* this is the atom h-like collision l-mixing command */
			fprintf(ioQQQ,"This command changed to ATOM H-LIKE COLLISIONS L-MIXING\n");
			fprintf(ioQQQ,"I will parse it for now, but may not in the future.\n");
			/* turn off 2s - 2p collisions */
			iso.lgColl_l_mixing[ipISO] = false;
			phycon.lgPhysOK = false;
		}

		else if( p.nMatch("L-MI") )
		{
			if( ipISO == ipH_LIKE )
			{
				/* >>chng 02 feb 07, change from 2s2p to l-mixing */
				/* this is the atom h-like collision l-mixing command */
				iso.lgColl_l_mixing[ipISO] = false;
				phycon.lgPhysOK = false;
			}
			else if( p.nMatch("THER") )
			{
				/* use l-mix from 
				 *>>refer	l-mix	all Vrinceanu, D. & Flannery, M. R. 2001, PhysRevA 63, 032701 */
				if( p.nMatch("NO T") )			
				{
					/* This is the "NO Thermal average" command.  It
					 * causes collisions strengths to be evaluated at kT rather than
					 * integrated over a Maxwellian peaked at kT. */
					iso.lgCS_therm_ave[ipISO] = false;
				}
				else
				{
					iso.lgCS_therm_ave[ipISO] = true;
				}
			}
			else if( p.nMatch("PENG") )
			{
				iso.lgCS_Vrinceanu[ipISO] = false;
			}
			else if( p.nMatch(" OFF"  ) )
			{
				/* this is the atom xx-like collision l-mixing command */
				/* turn off same-n collisions */
				iso.lgColl_l_mixing[ipISO] = false;
				phycon.lgPhysOK = false;
				iso.lgCS_Vrinceanu[ipISO] = false;
				iso.lgCS_therm_ave[ipISO] = false;
			}
			else
			{
				fprintf( ioQQQ, " needs parameter\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( p.nMatch(" OFF"  ) )
		{
			/* turn everything off, since no keyword given */
			iso.lgColl_excite[ipISO] = false;
			iso.lgColl_ionize[ipISO] = false;
			iso.lgColl_l_mixing[ipISO] = false;
			phycon.lgPhysOK = false;
		}
		else
		{
			fprintf( ioQQQ, " needs parameter\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("CONT") && p.nMatch("LOWE") )
	{
		/* disable continuum lowering for this isoelectronic sequence */
		if( p.nMatch("OFF") )
			iso.lgContinuumLoweringEnabled[ipISO] = false;
		else
			iso.lgContinuumLoweringEnabled[ipISO] = true;
	}

	else if( p.nMatch("DAMP") )
	{
		if( ipISO == ipHE_LIKE )
		{
			fprintf(ioQQQ," Sorry, the DAMPING option is not implemented for the he-like sequence.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* turn off absorption due to Ly alpha damping wings */
		hydro.DampOnFac = 0.;
	}

	else if( p.nMatch("DIEL") )
	{
		if( ipISO == ipH_LIKE )
		{
			fprintf(ioQQQ," Sorry, but dielectronic recombination onto the h-like sequence is not possible.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* This sets which set of data to use for dielectronic recombination.	*/
		if( p.nMatch(" OFF") )
		{
			iso.lgDielRecom[ipISO] = false;
		}
		else 
			iso.lgDielRecom[ipISO] = true;
	}

	else if( p.nMatch("LEVE") )
	{
		/* the number of levels read in is n, the principal quantum number
		 * only lines with upper levels less than n will be printed */

		/* number of levels for iso-sequence */
		/* there are two options here,
		 * when keyword ELEMENT appears, scan off element name and change levels only for
		 * that one.
		 * when there is no ELEMENT then set all in iso to same number */

		/* lgHydroMalloc is false at start of calculation, set true when space 
		 * allocated for the hydrogen and helium lines.  Once done we must ignore all 
		 * future changes in the number of levels */
		if( p.nMatch("PRIN") )
		{
			/* only print - do not change levels */
			iso.lgPrintNumberOfLevels = true;
		}
		else if( !lgHydroMalloc )
		{
			numLevels = (long int)p.FFmtRead();

			if( p.lgEOL() )
			{
				/* must be a number, or a key, either large (or limit) or small
				 * if no number but LARGER or SMALL then set to preset number */
				if( p.nMatch("LARG") )
				{
					/* there is no real limit, but this includes all levels with tabulated rec coefficient */
					numLevels = RREC_MAXN;
				}

				/* this is small or compact keyword */
				else if( p.nMatch("SMAL") || p.nMatch("COMP") )
				{
					if( ipISO == ipH_LIKE )
						numLevels = 5;
					else if( ipISO == ipHE_LIKE )
						numLevels = 3;
				}

				else
					/* punch out if no number */
					p.NoNumb("levels");
			}
			else if( ( numLevels < 3 ) && !p.nMatch("COLL") )
			{
				/* only enforce this limit if not working on collapsed levels */
				fprintf( ioQQQ, " cannot have fewer than 3 levels, the requested number was %li\n" ,
					numLevels  );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			if( ipISO == ipH_LIKE && numLevels > NHYDRO_MAX_LEVEL-2 )
			{
				fprintf( ioQQQ, " Not possible to set nhlvl to >NHYDRO_MAX_LEVEL-2= %i\n", 
				  NHYDRO_MAX_LEVEL-2 );
				fprintf( ioQQQ, " change NHYDRO_MAX_LEVEL\n");
				cdEXIT(EXIT_FAILURE);
			}

			/* this is check that alpha transition of highest level is within energy
			 * bounds of continuum */
			if( ipISO == ipH_LIKE &&
				( 2. / POW3((double)numLevels) < rfield.emm ) )
			{
				fprintf( ioQQQ, " Not possible to set iso.numLevels_max[ipH_LIKE][ipHYDROGEN] to such a high value, since "
					"alpha transition not within energy bounds of code\n");

				fprintf( ioQQQ, " lowest energy is %e and corresponding highest level is %li\n" ,
					rfield.emm, (long)pow(2./rfield.emm, 0.3333) );
				cdEXIT(EXIT_FAILURE);
			}

			if( p.nMatch("COLL") )
			{
				if( numLevels < 1 )
				{
					fprintf( ioQQQ, "There must be at least one collapsed level.\n");
					cdEXIT(EXIT_FAILURE);
				}

				/* set collapsed levels for entire sequence or just one element */
				if( p.nMatch("ELEM") )
				{
					/* check which element is on the line, nelem is atomic number on C scale */
					nelem = p.GetElem();
					iso.nCollapsed_max[ipISO][nelem] = numLevels;
					iso_update_num_levels( ipISO, nelem );
				}
				else
				{
					/* keyword element did not appear, so set all iso species to this number */
					for( nelem=ipISO; nelem<LIMELM; ++nelem )
					{
						iso.nCollapsed_max[ipISO][nelem] = numLevels;
						iso_update_num_levels( ipISO, nelem );
					}
				}
			}

			/* we now have the desired number of levels, and it has been fully qualified. 
			 * now check whether ELEMENT keyword is on line, if not then set all to
			 * the number, if so then only element in use */
			else if( p.nMatch("ELEM") )
			{
				/* check which element is on the line, nelem is atomic number on C scale */
				nelem = p.GetElem();
				iso.n_HighestResolved_max[ipISO][nelem] = numLevels;
				iso_update_num_levels( ipISO, nelem );
			}
			else
			{
				/* keyword element did not appear, so set all species to this number */
				for( nelem=ipISO; nelem<LIMELM; ++nelem )
				{
					iso.n_HighestResolved_max[ipISO][nelem] = numLevels;
					iso_update_num_levels( ipISO, nelem );
				}
			}
		}
	}

	else if( p.nMatch("ERRO") && p.nMatch("GENE"  ) )
	{
		/* Rates will be modified by a randomly generated error that falls within
		 * the range specifically set for each rate (or set of rates).	*/
		iso.lgRandErrGen[ipISO] = true;
		iso.modelRank[ipISO] = (int)p.FFmtRead();

		iso.modelRank[ipISO] = MAX2( 0, iso.modelRank[ipISO] );
		if( p.lgEOL() )
			/* Changed to avoid lint complaint */
			/* iso.modelRank[ipISO] = (unsigned)time(NULL); */
			iso.modelRank[ipISO] = abs((int)time(NULL));

		/* this allows a seed that is dependent upon the processor rank
		 * in a parallel run. */
		/* We add 2 so that the seed is always greater than 1, which would reset the generator. */
		init_genrand( (unsigned)iso.modelRank[ipISO] + 2);

		if( p.nMatch("PESS") )
			iso.lgPessimisticErrors = true;
		else
			iso.lgPessimisticErrors = false;
	}

	else if( p.nMatch(" FSM") )
	{
		if( ipISO == ipH_LIKE )
		{
			fprintf(ioQQQ," Sorry, but fine-structure mixing can only be implemented for the He-like sequence.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* turn on fine structure mixing of spontaneous decays.  
		 * >>refer	Helike	FSM	Bauman, R., MacAdams, K., and Ferland, G. (2003).	*/
		if( p.nMatch(" OFF") )
			iso.lgFSM[ipISO] = false;
		else
			iso.lgFSM[ipISO] = true;
	}

	else if( p.nMatch("GBAR") )
	{
		if( ipISO == ipH_LIKE )
		{
			fprintf(ioQQQ," Sorry, the GBAR option is only implemented for the He-like sequence.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* the HEGBAR command - to change cs of higher levels */
		/* first turn all off, one will be turned back on */
		iso.lgCS_Vriens[ipISO] = false;
		iso.lgCS_None[ipISO] = false;
		iso.nCS_new[ipISO] = false;

		/* now turn one on */
		if( p.nMatch("VRIE") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso.lgCS_Vriens[ipISO] = true;
		}
		else if( p.nMatch(" NEW") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso.nCS_new[ipISO] = (int)p.FFmtRead();

			/* there are two options, 1 and 2, for which fit - default (no number)
			 * will be 1, the broken power law fit */
			if( p.lgEOL() )
				iso.nCS_new[ipISO] = 1;

			ASSERT( iso.nCS_new[ipISO] );
		}
		else if( p.nMatch(" OFF") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso.lgCS_None[ipISO] = true;
		}
		else
		{
			fprintf( ioQQQ, " needs parameter\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}


	else if( p.nMatch("LYMA") )
	{
		if( ipISO == ipH_LIKE && p.nMatch("PUMP") )
		{
			/* >>chng 05 jul 08, separate out Lyman pump commands */
			if( p.nMatch(" OFF") )
			{
				/* option to turn off all continuum pumping of Lyman lines */
				hydro.lgLymanPumping = false;
			}
			else if( p.nMatch("SCALE") )
			{
				/* multiplicative factor for all continuum pumping of H I Lyman lines,
				 * account for possible emission in the line - only affects H I
				 * not entire H-like iso sequence */
				hydro.xLymanPumpingScaleFactor = 
					(realnum)p.FFmtRead();
				/* scale factor is log if <=0, 
				 * represents line in emission if >1
				 * LOG keyword forces interpretation as a log */
				if( hydro.xLymanPumpingScaleFactor <= 0. ||
					p.nMatch(" LOG") )
				{
					hydro.xLymanPumpingScaleFactor = 
						(realnum)pow((realnum)10.f , hydro.xLymanPumpingScaleFactor );
				}

				/* vary option */
				if( optimize.lgVarOn )
				{
					optimize.nvarxt[optimize.nparm] = 1;
					strcpy( optimize.chVarFmt[optimize.nparm], "ATOM H-LIKE LYMAN PUMPING SCALE %f LOG" );

					/*  pointer to where to write */
					optimize.nvfpnt[optimize.nparm] = input.nRead;

					/*  current parameters - always log so steps are log  */
					optimize.vincr[optimize.nparm] = 0.1f;
					optimize.vparm[0][optimize.nparm] = (realnum)log10(hydro.xLymanPumpingScaleFactor);
					++optimize.nparm;
				}
			}
			else
			{
				fprintf(ioQQQ," Sorry, I didn\'t recognize an option on this ATOM H-LIKE LYMAN PUMP command.\n");
				fprintf(ioQQQ," The options are \" OFF\", and \"SCALE\".\n");  
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			/* option to set number of "extra" Lyman lines, used for optical depths only */
			iso.nLyman[ipISO] = (long int)p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("'extra' Lyman lines");
		}
	}

	/* don't interpolate on look-up tables but compute recombination on the fly instead */
	else if( p.nMatch("RECO") &&
		p.nMatch(" NO ") && p.nMatch("INTE") )
	{
		/* flag set by atom xx-like no recombination interp command, 
		 * says to generate recombination coefficients
		 * on the fly */
		iso.lgNoRecombInterp[ipISO] = true;
	}

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
		/* complete redistribution no wings */
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

		/* resonance lines - not Lya*/
		if( p.nMatch("ALPH") )
		{
			iso.ipLyaRedist[ipISO] = ipRedis;
		}
		/* Lya itself */
		else if( p.nMatch("RESO") )
		{
			iso.ipResoRedist[ipISO] = ipRedis;
		}
		/* subordinate lines */
		else if( p.nMatch("SUBO") )
		{
			iso.ipSubRedist[ipISO] = ipRedis;
		}
		/* the show option, say what we are assuming */
		else if( p.nMatch("SHOW") )
		{
			fprintf(ioQQQ," Ly a is ");
			if( iso.ipLyaRedist[ipISO] ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( iso.ipLyaRedist[ipISO] ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( iso.ipLyaRedist[ipISO] ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else if( iso.ipLyaRedist[ipISO] ==ipLY_A )
			{
				fprintf(ioQQQ,"special Lya.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for iso.ipLyaRedist.\n");
				TotalInsanity();
			}

			fprintf(ioQQQ," Other %s resonance lines are ",
				elementnames.chElementSym[ipISO] );

			if( iso.ipResoRedist[ipISO] ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( iso.ipResoRedist[ipISO] ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( iso.ipResoRedist[ipISO] ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for iso.ipResoRedist.\n");
				TotalInsanity();
			}

			fprintf(ioQQQ," %s subordinate lines are ",
				elementnames.chElementSym[ipISO] );

			if( iso.ipSubRedist[ipISO] ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( iso.ipSubRedist[ipISO] ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( iso.ipSubRedist[ipISO] ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for iso.ipSubRedist.\n");
				TotalInsanity();
			}
		}
		else
		{
			fprintf(ioQQQ," here should have been another keyword on this command.\n");
			fprintf(ioQQQ," Options are ALPHA, RESONANCE, SUBORDINATE.  Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* which type of matrix solution, populations or lowt */
	else if( p.nMatch("MATR") )
	{
		fprintf( ioQQQ, " This command has been deprecated.  It now has no effect.\n" );
	}

	else if( p.nMatch("TOPO") )
	{
		if( ipISO == ipH_LIKE )
		{
			fprintf( ioQQQ, " The behavior of the TOPOFF option for the H-like sequence has changed since version 07.02.\n");
			fprintf( ioQQQ, " Originally, it accepted a number and distributed \"topoff\" equally into that many levels at the top of the model.\n");
			fprintf( ioQQQ, " Now, the option works the same as in the He-like sequence, which is to turn topoff on or (with keyword _OFF) off.\n" );
		}

		if( p.nMatch(" OFF") )
		{
			iso.lgTopoff[ipISO] = false;
			fprintf( ioQQQ, "ISO %li TOPOFF is OFF\n", ipISO );
		}
		else
			iso.lgTopoff[ipISO] = true;
	}


	else
	{
		fprintf( ioQQQ, " There should have been a keyword - STOP\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
