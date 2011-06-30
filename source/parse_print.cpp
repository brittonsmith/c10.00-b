/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParsePrint parse the print command  */
/*prt_constants print physical constants */
/* PrtMacros - print macros in cddefines.h and their status -  *print macros* */
#include "cddefines.h"
#include "physconst.h"
#include "rfield.h"
#include "iso.h"
#include "iterations.h"
#include "lines.h"
#include "called.h"
#include "elementnames.h"
#include "prt.h"
#include "h2.h"
#include "parser.h"
#include "version.h"
/*prt_constants print physical constants */
STATIC void prt_constants(void);


// PrtMacros - print macros in cddefines.h and their status -  *print macros*
STATIC void PrtMacros( void )
{

	DEBUG_ENTRY( "PrtMacros()" );

	fprintf( ioQQQ," PrtMacros:\n FLT_IS_DBL\t");
#	ifdef FLT_IS_DBL
	fprintf( ioQQQ,"defined.\n");
#	endif

	fprintf( ioQQQ , "\n DMALLOC\t");
#	ifdef DMALLOC
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n SYS_CONFIG\t");
#	ifdef SYS_CONFIG
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n MPI_GRID_RUN\t");
#	ifdef MPI_GRID_RUN
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n USE_GPROF\t");
#	ifdef USE_GPROF
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_FUNC\t");
#	ifdef HAVE_FUNC
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n _MSC_VER\t");
#	ifdef _MSC_VER
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n __INTEL_COMPILER\t");
#	ifdef __INTEL_COMPILER
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n OLD_ASSERT\t");
#	ifdef OLD_ASSERT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POWI\t");
#	ifdef HAVE_POWI
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_DOUBLE_INT\t");
#	ifdef HAVE_POW_DOUBLE_INT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_DOUBLE_LONG\t");
#	ifdef HAVE_POW_DOUBLE_LONG
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_FLOAT_INT\t");
#	ifdef HAVE_POW_FLOAT_INT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_FLOAT_LONG\t");
#	ifdef HAVE_POW_FLOAT_LONG
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_FLOAT_DOUBLE\t");
#	ifdef HAVE_POW_FLOAT_DOUBLE
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_DOUBLE_FLOAT\t");
#	ifdef HAVE_POW_DOUBLE_FLOAT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n");

}

void ParsePrint(
	/* input line was converted to caps by calling routine */
	Parser &p )
{
	int ipISO;
	long int
	  j, 
	  nelem,
	  num1;
	double a;

	DEBUG_ENTRY( "ParsePrint()" );

	/* >>chng 01 aug 91, had been series of if branches, and could hit more than
	 * one key - dangerous!  changed to else if so only one hit per line possible */
	if( p.nMatch("AGES") )
	{
		/* print all estimates of cloud timescales */
		prt.lgPrnAges = true;
	}

	else if( p.nMatch("ARRA") )
	{
		/* print arrays for ionization balance of heavy elements */
		if( p.nMatch( "ONLY"  ) )
		{
			/* returns element number on C scale */
			if( (nelem = p.GetElem())<0 )
			{
				fprintf(ioQQQ,"An element name must appear on this PRINT ARRAYS ONLY xx command.\n");
				cdEXIT(EXIT_FAILURE);
			}
			/* have the element number, turn on its print */
			prt.lgPrtArry[nelem] = true;
		}
		else
		{
			/* this flag, print arrays for all elements */
			for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
			{
				prt.lgPrtArry[nelem] = true;
			}
		}
	}

	else if( p.nMatch("CITA") )
	{
		prt.lgPrtCitations = true;
	}

	else if( p.nMatch("COLU") && p.nMatch("DENS") )
	{
		if( p.nMatch(" ON ") )
		{
			/* print column densities of elements - this is default */
			prt.lgPrintColumns = true;
		}
		else if( p.nMatch(" OFF") )
		{
			/* print column densities of elements */
			prt.lgPrintColumns = false;
		}
	}

	/* option to print departure coefficients in addition to level pops 
	 * keywords He-like to do He-like sequence element, else do h-like
	 * element name, if not recognized, does hydrogen
	 * so with no options prints hydrogen itself */
	else if( p.nMatch("DEPA") )
	{
		if( p.nMatch("HE-L") )
		{
			ipISO = ipHE_LIKE;
		}
		else
		{
			ipISO = ipH_LIKE;
		}

		/* now check for element name */
		nelem = p.GetElem( );

		/* return value is < 0 if no element recognized - in this case use root of sequence */
		nelem = MAX2( nelem, ipISO );

		/* print departure coefficients instead of hydrogen level populations */
		iso.lgPrtDepartCoef[ipISO][nelem] = true;
	}

	else if( p.nMatch("CONS") )
	{
		/* print physical constants, routine is below */
		prt_constants();
	}

	else if( p.nMatch("ERRO") )
	{
		/* print errors to special window */
		lgPrnErr = true;
	}

	else if( p.nMatch("HEAT") )
	{
		/* print heat arrays */
		prt.lgPrintHeating = true;
	}

	else if( p.nMatch("PATH") )
	{
		/* print the path */
		cpu.printDataPath();
	}

	/*else if( p.nMatch("H-LI"))*/
	else if( p.nMatch("POPU"))
	{
		if( p.nMatch("HE-L") )
		{
			ipISO = ipHE_LIKE;
		}
		else
		{
			ipISO = ipH_LIKE;
		}

		/* now check for element name */
		nelem = p.GetElem( );
		/* return value is < 0 if no element recognized - in this case use H */
		nelem = MAX2(0,nelem);

		/* if no element specified but he-like iso, then use helium */
		if( nelem==0 && ipISO==ipHE_LIKE )
			nelem = ipHELIUM;

		if( nelem < ipISO )
		{
			fprintf(ioQQQ,"This iso-sequence (%s) and element (%s) are impossible.\n",
				elementnames.chElementName[ipISO],
				elementnames.chElementName[nelem]);
			cdEXIT(EXIT_FAILURE);
		}

		/* print hydrogenic H-like level populations */
		iso.lgPrtLevelPops[ipISO][nelem] = true;
	}

	/* option to only print last iteration */
	else if( p.nMatch("LAST") )
	{
		prt.lgPrtLastIt = true;
	}

	/* the print line command as several options */
	else if( p.nMatch("LINE") )
	{
		if( p.nMatch(" ALL") )
		{
			/* turn on all printed components */
			prt.lgPrnPump = true;
			prt.lgPrnColl = true;
			prt.lgPrnHeat = true;
		}

		else if( p.nMatch("CELL") )
		{
			/* print line cell on physics scale, first cell in continuum is 1
			 * give all lines in this cell */
			prt.lgPrnLineCell = true;
			prt.nPrnLineCell = (long)p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("cell for line print" );
			if( prt.nPrnLineCell < 1 )
			{
				/* non-positive cells are not allowed */
				fprintf(ioQQQ , "The cell number on the PRINT LINE CELL command must be positive.\n");
				fprintf(ioQQQ , "The cell number was %li.\n" , prt.nPrnLineCell);
			}
		}

		else if( p.nMatch("COLL") )
		{
			/* either print line collisions or print line iso collapsed */
			/* also print collisional contributions */
			if( p.nMatch(" ISO") )
			{
				/* print predictions from collapsed levels of iso sequences */
				prt.lgPrnIsoCollapsed = true;
			}
			else
			{
				/* print line collisions */
				prt.lgPrnColl = true;
			}
		}

		else if( p.nMatch("COLU") )
		{
			/* option to print main line array as a single column */
			prt.lgPrtLineArray = false;
			/* this also has an option - liNEAR - to print linear quantity 
			 * in exponential format */
			if( p.nMatch("NEAR") )
				prt.lgPrtLineLog = false;
		}

		else if( p.nMatch("FAIN") && !(p.nMatch("OPTI")&&p.nMatch("DEPT")) )
		{
			/* print line faint - above do not trigger on optical depth 
			 * option to adjust intensity of faintest line to print */
			/* >> 01 feb 01, move print faint into print line faint */
			/* faintest line, rel to norm line, to print; either linear of log */
			a = p.FFmtRead();

			/* option for, if no number, keyword=" OFF", to print all lines */
			if( p.lgEOL() )
			{
				if( p.nMatch(" OFF") )
				{
					prt.lgFaintOn = false;
				}
				else
				{
					fprintf( ioQQQ, 
						" There faintest line to print must be on this line, sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}

			prt.lgFntSet = true;
			if( a <= 0. || p.nMatch(" LOG") )
			{
				prt.TooFaint = (realnum)pow(10.,a);
			}
			else
			{
				prt.TooFaint = (realnum)a;
			}
		}

		else if( p.nMatch("FLUX") && p.nMatch("EART"))
		{
			/* print line flux seen at earth */
			prt.lgPrintFluxEarth = true;
		}

		else if( p.nMatch(" H2") && p.nMatch("ELEC") )
		{
			/* print H2 electronic lines too - -1 since number of electronic
			 * levels is not yet known, will set when H2 actually called */
			h2.nElecLevelOutput = -1;
		}

		else if( p.nMatch("HEAT") )
		{
			/* also print heating contributions */
			prt.lgPrnHeat = true;
		}

		else if( p.nMatch("INWA") )
		{
			/* also print inward contributions */
			prt.lgPrnInwd = true;
		}

		else if( p.nMatch("OPTI") && p.nMatch("DEPT") )
		{
			/* print line optical depths, with option for smallest to print */
			if( p.nMatch(" OFF") )
			{
				/* turn off or on printing of optical depths - default off */
				prt.lgPrtTau = false;
			}
			else
			{
				prt.lgPrtTau = true;
			}
			if( p.nMatch("FAIN") )
			{
				/* log of faintest optical depth, default is linear value of 0.1 */
				prt.PrtTauFnt = (realnum)pow(10.,p.FFmtRead());
				if( p.lgEOL() )
				{
					fprintf( ioQQQ, " There must be a number for the FAINT option.  They are HEAD and ZONE.  Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
		}

		else if( p.nMatch("PUMP") )
		{
			/* also print pump contributions */
			prt.lgPrnPump = true;
		}

		else if( p.nMatch("SORT") )
		{
			/* >>chng 01 aug 18, print sort command works after all these years,
			 * sort by wavelength or intensity */
			/* turn on sorting with respect to wavelength */
			prt.lgSortLines = true;
			if( p.nMatch("WAVE") )
			{
				/* sort by wavelength */
				/* remember which one to do */
				prt.lgSortLineIntensity = false;
				prt.lgSortLineWavelength = true;

				/* wavelength has range option */
				/* option to only certain print range of lines */
				if( p.nMatch("RANG") )
				{
					prt.wlSort1 = (realnum)p.getWaveOpt();

					prt.wlSort2 = (realnum)p.getWaveOpt();

					if( p.lgEOL() )
					{
						fprintf( ioQQQ, " There must be two numbers for the RANGE option, the lower and upper wavelength.  Sorry.\n" );
						cdEXIT(EXIT_FAILURE);
					}
					if( prt.wlSort1 <0. || prt.wlSort2 <0. || 
						prt.wlSort1 >= prt.wlSort2 )
					{
						fprintf( ioQQQ, " The lower and upper wavelength must be positive and in the correct order.  Sorry.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
				else
				{
					prt.wlSort1 = -1;
					prt.wlSort2 = 1e30f;
				}
			}
			else if( p.nMatch("INTE") )
			{
				/* sort by intensity/luminosity */
				/* remember which one to do */
				prt.lgSortLineIntensity = true;
				prt.lgSortLineWavelength = false;
			}
			else
			{
				fprintf( ioQQQ, "I can sort by wavelength or intensity - one must be specified.\nSorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch(" SUM") )
		{
			/* option to read in set of lines to sum over */
			ParsePrtLineSum( p );
		}

		else if( p.nMatch("SURF") && p.nMatch("BRIG") )
		{
			/* print surface brightness rather than 4pi J */
			prt.lgSurfaceBrightness = true;
			/* default is per sr, arcsec option changes to sq arcsec */
			if( p.nMatch("ARCS" ) )
			{
				/* use sr */
				prt.lgSurfaceBrightness_SR = false;
			}
			else
			{
				/* use sq arcsec */
				prt.lgSurfaceBrightness_SR = true;
			}
		}

		else if( p.nMatch("CUMU") )
		{
			/* print lines cumulative - integral of line emission over time */
			prt.lgPrintLineCumulative = true;
		}

		else
		{
			fprintf( ioQQQ, "One of the keys should have appeared.  \nPlease consult Hazy.\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* print maser lines when TAV is called */
	else if( p.nMatch("MASE") )
	{
		prt.lgPrtMaser = true;
	}

	else if( p.nMatch("ONLY") )
	{
		if( p.nMatch("ZONE") )
			prt.lgOnlyZone = true;

		else if( p.nMatch("HEAD") )
			prt.lgOnlyHead = true;

		else
		{
			fprintf( ioQQQ, " There must be a keyword for the ONLY option.  They are HEAD and ZONE.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("STAR") )
	{
		/* start printout at specified zone */
		called.lgTalk = false;
		prt.lgPrtStart = true;
		prt.nstart = (long int)p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, 
				" The zone on which the print is to start MUST be entered on this line.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* print continuum command */
	else if( p.nMatch("CONT") )
	{
		/* >>chng 01 jun 30, logic of keys swapped around, so no keyword gets nFnu info */
		/* print continuum diffuse turns on nFnu continuum at many energies */
		if( p.nMatch("BLOC") )
		{
			/* option to print emergent continuum at end of calculation*/
			prt.lgPrtCont = true;
		}
		else if( p.nMatch("INDI"  ))
		{
			/* option to print lines and continuum that go into each continuum 
			 * index the continuum index is the cell within the continuum 
			 * array - this identifies lines that occur within each 
			 * continuum cell */
			prt.lgPrtContIndices = true;
			/* these are lower and upper limits to the energy range in Rydbergs.
			* they are the first and second number on the command line, lower and
			* upper bounds of the code are used if not specified */
			/* if no number on line then zero is returned, this is fine, since
			 * we want the lower energy bound of the code */
			prt.lgPrtContIndices_lo_E = (realnum)p.FFmtRead();
			prt.lgPrtContIndices_hi_E = (realnum)p.FFmtRead();
			/* if we hit end of line then use high-energy limit of code - that is,
			 * include all energies */
			if( p.lgEOL() )
				prt.lgPrtContIndices_hi_E = (realnum)rfield.egamry;
		}
		else
		{
			/* option to print continuum points within emission lines block */
			prt.lgPrnDiff = true;
		}
	}

	else if( p.nMatch("COOL") )
	{
		/* print cooling array for a specified one */
		prt.nzdump = (long int)p.FFmtRead();

		/* dump all zones if argument is zero or not present  */
		if( p.lgEOL() )
		{
			prt.nzdump = 0;
		}
	}

	else if( p.nMatch("QUIE") || (p.nMatch(" OFF") && 
		!p.nMatch("FAIN" )) )
	{
		/* in above, there is also a 'print faint off' command
		 * QUIET or OFF means turn off printout */
		called.lgTalk = false;
	}

	else if( p.nMatch("MACR") )
	{
		// print status of macros in cddefines.ht */
		PrtMacros();
	}

	else if( p.nMatch(" ON ") )
	{
		/* on means turn on printout, lgTalkIsOK is set false in grid_do.cpp.
		 * this keeps printout quiet during optimize, even when init files are parsed */
		/* called.lgTalkForcedOff was set true with cdTalk(false), if this was
		 * set then do not respect this command.  this is to prevent print on at end
		 * of init file from turning on print in grids when print is turned off */
		if( called.lgTalkIsOK && !called.lgTalkForcedOff )
		{
			called.lgTalk = cpu.lgMPI_talk();
		}
	}

	else if( p.nMatch("SHOR") )
	{
		/* make short printout, don't print last */
		prt.lgPrtShort = true;
		if( !prt.lgFntSet )
			prt.TooFaint = 0.001f;
	}

	else if( p.nMatch("VERS") )
	{
		/* print compiler and code version information */
		fprintf( ioQQQ, "\nThis is Cloudy %s\n%s\n\n" ,
			t_version::Inst().chVersion,
			t_version::Inst().chInfo );
	}

	else if( p.nMatch("ZONE") || p.nMatch("EVER") )
	{
		/* print every nth zone - command was originally print every but
		 * is being changed to print zone */
		num1 = (long int)p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " The number of zones to print MUST be entered on this line.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		iterations.IterPrnt[0] = MAX2(num1,1);

		for( j=1; j < iterations.iter_malloc; j++ )
		{
			iterations.IterPrnt[j] = (long int)p.FFmtRead();
			if( p.lgEOL() )
			{
				iterations.IterPrnt[j] = iterations.IterPrnt[j-1];
			}
		}
	}

	/* check if no keywords were recognized. */
	else
	{
		fprintf( ioQQQ, " There MUST be a keyword on the following line.  Sorry.\n" );
		fprintf( ioQQQ, " The PRINT FAINT command is now the PRINT LINE FAINT command.\n" );
		p.PrintLine(ioQQQ);
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/*prt_constants print physical and machine constants */
STATIC void prt_constants(void)
{

	DEBUG_ENTRY( "prt_constants()" );

	fprintf(ioQQQ,"\n\nPhysical constants used by Cloudy, taken from physconst.h\n");

	fprintf(ioQQQ,"EE\t%.15g\n",EE);
	fprintf(ioQQQ,"EULER\t%.15g\n",EULER);
	fprintf(ioQQQ,"PI\t%.15g\n",PI);
	fprintf(ioQQQ,"PI2\t%.15g\n",PI2);
	fprintf(ioQQQ,"PI4\t%.15g\n",PI4);
	fprintf(ioQQQ,"PI8\t%.15g\n",PI8);
	fprintf(ioQQQ,"SQRT2\t%.15g\n",SQRT2);
	fprintf(ioQQQ,"SQRTPI\t%.15g\n",SQRTPI);
	fprintf(ioQQQ,"SQRTPIBY2\t%.15g\n",SQRTPIBY2);
	fprintf(ioQQQ,"LN_TWO\t%.15g\n",LN_TWO);
	fprintf(ioQQQ,"LN_TEN\t%.15g\n",LN_TEN);
	fprintf(ioQQQ,"LOG10_E\t%.15g\n",LOG10_E);
	fprintf(ioQQQ,"OPTDEP2EXTIN\t%.15g\n",OPTDEP2EXTIN);
	fprintf(ioQQQ,"RADIAN\t%.15g\n",RADIAN);
	fprintf(ioQQQ,"SOLAR_MASS\t%.15g\n",SOLAR_MASS);
	fprintf(ioQQQ,"SOLAR_LUMINOSITY\t%.15g\n",SOLAR_LUMINOSITY);
	fprintf(ioQQQ,"AU\t%.15g\n",AU);
	fprintf(ioQQQ,"ATOMIC_MASS_UNIT\t%.15g\n",ATOMIC_MASS_UNIT);
	fprintf(ioQQQ,"ELECTRON_MASS\t%.15g\n",ELECTRON_MASS);
	fprintf(ioQQQ,"PROTON_MASS\t%.15g\n",PROTON_MASS);
	fprintf(ioQQQ,"BOLTZMANN\t%.15g\n",BOLTZMANN);
	fprintf(ioQQQ,"SPEEDLIGHT\t%.15g\n",SPEEDLIGHT);
	fprintf(ioQQQ,"HPLANCK\t%.15g\n",HPLANCK);
	fprintf(ioQQQ,"GRAV_CONST\t%.15g\n",GRAV_CONST);
	fprintf(ioQQQ,"ELEM_CHARGE\t%.15g\n",ELEM_CHARGE);
	fprintf(ioQQQ,"RYD_INF\t%.15g\n",RYD_INF);
	fprintf(ioQQQ,"HIONPOT\t%.15g\n",HIONPOT);
	fprintf(ioQQQ,"AS1RAD\t%.15g\n",AS1RAD);
	fprintf(ioQQQ,"SQAS1SR\t%.15g\n",SQAS1SR);
	fprintf(ioQQQ,"SQAS_SKY\t%.15g\n",SQAS_SKY);
	fprintf(ioQQQ,"PARSEC\t%.15g\n",PARSEC);
	fprintf(ioQQQ,"H_BAR \t%.15g\n",H_BAR );
	fprintf(ioQQQ,"ELEM_CHARGE_ESU \t%.15g\n",ELEM_CHARGE_ESU );
	fprintf(ioQQQ,"ELECTRIC_CONST\t%.15g\n",ELECTRIC_CONST);
	fprintf(ioQQQ,"HION_LTE_POP\t%.15g\n",HION_LTE_POP);
	fprintf(ioQQQ,"SAHA\t%.15g\n",SAHA);
	fprintf(ioQQQ,"ERG1CM\t%.15g\n",ERG1CM);
	fprintf(ioQQQ,"T1CM\t%.15g\n",T1CM);
	fprintf(ioQQQ,"WAVNRYD\t%.15g\n",WAVNRYD);
	fprintf(ioQQQ,"RYDLAM\t%.15g\n",RYDLAM);
	fprintf(ioQQQ,"EN1RYD\t%.15g\n",EN1RYD);
	fprintf(ioQQQ,"TE1RYD\t%.15g\n",TE1RYD);
	fprintf(ioQQQ,"EVDEGK\t%.15g\n",EVDEGK);
	fprintf(ioQQQ,"EVRYD\t%.15g\n",EVRYD);
	fprintf(ioQQQ,"EN1EV\t%.15g\n",EN1EV);
	fprintf(ioQQQ,"FR1RYD\t%.15g\n",FR1RYD);
	fprintf(ioQQQ,"HNU3C2\t%.15g\n",HNU3C2);
	fprintf(ioQQQ,"FR1RYDHYD\t%.15g\n",FR1RYDHYD );
	fprintf(ioQQQ,"HBAReV\t%.15g\n",HBAReV );
	fprintf(ioQQQ,"RYDLAMHYD\t%.15g\n",RYDLAMHYD );
	fprintf(ioQQQ,"STEFAN_BOLTZ\t%.15g\n",STEFAN_BOLTZ);
	fprintf(ioQQQ,"FREQ_1EV\t%.15g\n",FREQ_1EV);
	fprintf(ioQQQ,"FINE_STRUCTURE\t%.15g\n",FINE_STRUCTURE);
	fprintf(ioQQQ,"BOHR_RADIUS_CM\t%.15g\n",BOHR_RADIUS_CM);
	fprintf(ioQQQ,"TWO_PHOT_CONST\t%.15g\n",TWO_PHOT_CONST);
	fprintf(ioQQQ,"COLL_CONST\t%.15g\n",COLL_CONST);
	fprintf(ioQQQ,"MILNE_CONST\t%.15g\n",MILNE_CONST);
	fprintf(ioQQQ,"TRANS_PROB_CONST\t%.15g\n",TRANS_PROB_CONST);
	fprintf(ioQQQ,"\n");

	fprintf(ioQQQ,"Some other interesting sizes:\n");
	fprintf(ioQQQ,"bool\t%lu\n",(unsigned long)sizeof(bool));
	fprintf(ioQQQ,"char\t%lu\n",(unsigned long)sizeof(char));
	fprintf(ioQQQ,"int\t%lu\n",(unsigned long)sizeof(int));
	fprintf(ioQQQ,"long int\t%lu\n",(unsigned long)sizeof(long int));
	fprintf(ioQQQ,"unsigned int\t%lu\n",(unsigned long)sizeof(unsigned int));
	fprintf(ioQQQ,"float\t%lu\n",(unsigned long)sizeof(sys_float));
	fprintf(ioQQQ,"realnum\t%lu\n",(unsigned long)sizeof(realnum));
	fprintf(ioQQQ,"double\t%lu\n",(unsigned long)sizeof(double));
	fprintf(ioQQQ,"double*\t%lu\n",(unsigned long)sizeof(double*));
	fprintf(ioQQQ,"\n");

	fprintf(ioQQQ,"Some constants from float.h.\n");
	/* some constants from float.h */
	fprintf(ioQQQ,"DBL_DIG \t%i\n", DBL_DIG);         /* # of decimal digits of precision */
	fprintf(ioQQQ,"DBL_EPSILON \t%.15g\n",DBL_EPSILON);   /* smallest such that 1.0+DBL_EPSILON != 1.0 */
	fprintf(ioQQQ,"DBL_MANT_DIG\t%i\n",DBL_MANT_DIG); /* # of bits in mantissa */
	fprintf(ioQQQ,"DBL_MAX\t%.15g\n", DBL_MAX);           /* max value */
	fprintf(ioQQQ,"DBL_MAX_10_EXP\t%i\n", DBL_MAX_10_EXP); /* max decimal exponent */
	fprintf(ioQQQ,"DBL_MAX_EXP\t%i\n", DBL_MAX_EXP);  /* max binary exponent */
	fprintf(ioQQQ,"DBL_MIN\t%.15g\n", DBL_MIN);           /* min positive value */

	fprintf(ioQQQ,"FLT_DIG\t%i\n", FLT_DIG);          /* # of decimal digits of precision */
	fprintf(ioQQQ,"FLT_EPSILON\t%.15g\n", FLT_EPSILON);   /* smallest such that 1.0+FLT_EPSILON != 1.0 */
	fprintf(ioQQQ,"FLT_MANT_DIG\t%i\n", FLT_MANT_DIG); /* # of bits in mantissa */
	fprintf(ioQQQ,"FLT_MAX\t%.15g\n", FLT_MAX);            /* max value */
	fprintf(ioQQQ,"FLT_MAX_10_EXP\t%i\n", FLT_MAX_10_EXP);/* max decimal exponent */
	fprintf(ioQQQ,"FLT_MAX_EXP\t%i\n", FLT_MAX_EXP);   /* max binary exponent */
	fprintf(ioQQQ,"FLT_MIN\t%.15g\n", FLT_MIN);            /* min positive value */

	fprintf(ioQQQ,"BIGFLOAT\t%.15g\n", BIGFLOAT);
	fprintf(ioQQQ,"SMALLFLOAT\t%.15g\n", SMALLFLOAT);

	return;
}
