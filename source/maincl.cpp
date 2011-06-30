/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that reads input and calls cloudy to compute a single model, or
 * try to optimize an observed model.  Routine returns 0 if model is ok, 
 * and 1 if problems occurred. */
#include "cddefines.h"
#include "cddrive.h"
#include "input.h"
#include "prt.h"
#include "save.h"
#include "called.h"
#include "monitor_results.h"
#ifdef MPI_ENABLED
#include "mpi_utilities.h"
#include "grid.h"
#endif

int cdMain( int argc, const char* argv[] );

/** main: this is an MPI wrapper around cdMain. It should ONLY do MPI stuff!
 * for non-MPI runs, this should do nothing more than call cdMain and exit. */
int main( int argc, char *argv[] )
{
	int exit_status = 0;

	DEBUG_ENTRY( "main()" );

#	ifdef MPI_ENABLED

	MPI::Init( argc, argv );

	cpu.set_MPI();
	cpu.set_nCPU( MPI::COMM_WORLD.Get_size() );
	cpu.set_nRANK( MPI::COMM_WORLD.Get_rank() );

#	ifdef MPI_GRID_RUN

	// this will generate input files for each grid point,
	// or execute the Phymir run, whichever is appropriate
	exit_status = cdMain( argc, (const char**)argv );

	// wait for writing of input files to finish
	MPI::COMM_WORLD.Barrier();

	// process the individual grid points
	if( grid.lgGrid )
	{
		// this was set to true after we wrote the last input script
		grid.lgGridDone = false;

		// from now on each rank will run its own model
		cpu.set_MPISingleRankMode( true );

		MPI::load_balance lb( grid.totNumModels );

		// Each MPI rank will get jobs assigned by lb and execute them.
		// If there are no jobs left, lb.next_job() will return -1.
		while( ( optimize.nOptimiz = lb.next_job() ) >= 0 )
		{
			const char** mpi_argv = new const char*[argc+2];

			string jobName = GridPointPrefix( optimize.nOptimiz );
			for( int i=0; i < argc; ++i )
				mpi_argv[i] = argv[i];
			mpi_argv[argc] = "-g";
			mpi_argv[argc+1] = jobName.c_str();

			int retval = cdMain( argc+2, mpi_argv );

			exit_status = max( retval, exit_status );
			delete[] mpi_argv;
		}

		lb.finalize();

		// gather the spectral information from all ranks for saving FITS files
		grid.lgGridDone = true;
		optimize.nOptimiz = grid.totNumModels;
		GridGatherInCloudy();

		// and concatenate the output
		if( cpu.lgMaster() )
			process_output();
	}

	// remove empty output files from slave ranks
	if( cpu.lgMaster() )
	{
		for( long n=1; n < cpu.nCPU(); ++n )
		{
			ostringstream oss;
			oss << ".err" << setfill('0') << setw(2) << n;
			string slave_output = save.chRedirectPrefix + oss.str();
			FILE *io = open_data( slave_output.c_str(), "a", AS_LOCAL_ONLY );
			bool lgEmpty = ( ftell(io) == 0 );
			fclose( io );
			if( lgEmpty )
				remove( slave_output.c_str() );
		}
	}

#	endif /* MPI_GRID_RUN */

	MPI::Finalize();

#	else /* MPI_ENABLED */

	// do non-MPI run
	exit_status = cdMain( argc, (const char**)argv );

#	endif /* MPI_ENABLED */

	return exit_status;
}

/** cdMain: this is the main entry point for Cloudy */
int cdMain( int argc, const char* argv[] )
{
	/* these will be used to count number of various problems */
	long int NumberWarnings, 
	  NumberCautions, 
	  NumberNotes, 
	  NumberSurprises, 
	  NumberTempFailures, 
	  NumberPresFailures,
	  NumberIonFailures, 
	  NumberNeFailures;

	bool lgAbort_exit,
	  lgFileIO,
	  lgBadExit;
	/* number of lines we can still read in */
	int nread=0;

	int i;
	const char *s, 
	  *prefix = "",
	  *gprefix = "", // grid prefix
	  *pprefix = "", // save prefix
	  *rprefix = ""; // redirect prefix

	/* the length of the following vector will be the longest line image
	 * the code will be able to read here.  Cloudy itself will ignore anything 
	 * beyond INPUT_LINE_LENGTH, and checks that no information exists beyond it. 
	 * The code will stop if the input line is longer than INPUT_LINE_LENGTH
	 * since extra characters would become a new command line due to buffer overrun */
	char chLine[INPUT_LINE_LENGTH];

	/* indicates that a command line flag to redirect I/O has been used */
	lgFileIO = false;

	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "cdMain()" );

	try {
		/* 
		 * Handle argument input -- written generally, but presently handles
		 * only `-p prefix' or `-pprefix' to set input file as `prefix.in',
		 * output file as `prefix.out' and the save prefix.
		 */
		for( i=1; i < argc; i++ ) 
		{
			s = argv[i];
			if( *s != '-' || s[1] == '\0' ) 
			{
				fprintf( ioQQQ, "%s: argument %d `%s' not understood\n", argv[0], i, s );
				cdEXIT( EXIT_FAILURE );
			}
			else
			{
				while( s != NULL && *(++s) )
				{
					switch( *s ) 
					{
					case 'a':
						cpu.setAssertAbort( true );
						break;
					case 'g':
					case 'p':
					case 'r':
						if( s[1] != '\0' )
						{
							prefix = s+1;
						}					
						else
						{
							if( ++i == argc || argv[i][0] == '-' )
							{
								fprintf( ioQQQ, "%s: no argument given for -%c flag\n",
									 argv[0], *s );
								cdEXIT( EXIT_FAILURE );
							}
							prefix = argv[i];
						}
						if( *s == 'g' )
							gprefix = prefix;
						else if( *s == 'p' )
						{
							pprefix = prefix;
							rprefix = prefix;
						}
						else if( *s == 'r' )
							rprefix = prefix;
						else
							TotalInsanity();
						s = NULL;
						lgFileIO = true;
						break;
					default:
						fprintf( ioQQQ, "%s: argument %d, `%s': flag -%c not understood\n",
							 argv[0], i, argv[i], *s );
						cdEXIT( EXIT_FAILURE );
					}
				}
			}
		}

		/* initialize the code for this run */
		cdInit();

		save.chGridPrefix = gprefix;
		save.chFilenamePrefix = pprefix;
		save.chRedirectPrefix = rprefix;

		if( cpu.lgMPI() && save.chRedirectPrefix.empty() )
		{
			// Only the master rank should give this error message to prevent ranks from clobbering
			// each others output. After the next statements each rank will have its own output file.
			if( cpu.lgMPI_talk() )
			{
				fprintf( ioQQQ, " The -p or -r flag to redirect I/O is mandatory in MPI runs.\n" );
				fprintf( ioQQQ, " Usage: mpirun -np <n> /path/to/cloudy.exe [ -p | -r ] <prefix>\n" );
			}
			cdEXIT(EXIT_FAILURE);
		}

		/* following should be set true to print to file instead of std output */
		if( lgFileIO )
		{
			string Base = save.chGridPrefix + save.chRedirectPrefix;
			string InName = Base + ".in";
			string OutName;
			if( cpu.lgMPI_talk() )
				OutName = Base + ".out";
			else
			{
				ostringstream oss;
				oss << ".err" << setfill('0') << setw(2) << cpu.nRANK();
				OutName = Base + oss.str();
			}
			cdInput( InName.c_str(), "r" );
			cdOutput( OutName.c_str(), "w" );
		}

		nread = 1;
		/* keep reading input lines until end of file */
		while( read_whole_line(chLine, (int)sizeof(chLine), ioStdin)!= NULL )
		{
			char *chChar;
			/* when running from command prompt, last line can be \n or 
			 * \r (this can happen with gcc under cygwin in windows) 
			 * or space, so check for each here, and break if present,
			 * check on chLine[23] is for case where we are reading cloudy output */
			/*fprintf(ioQQQ,"DEBUG char0 %i\n", (int)chLine[0] );*/
			if( chLine[0] == '\n' || chLine[0] == '\r' || 
			    ( chLine[0] == ' ' && chLine[23] != '*' ) )
				break;
			
			/* read entire line of input - lgAbort set true if lines too long */
			if( !lgInputComment(chLine) )
			{

				if( (chChar = strchr_s(chLine , '\"' ) ) ==NULL )
				{
					/* check for underscore, probably meant as a space */
					while( (chChar = strchr_s(chLine , '_' ) ) !=NULL )
					{
						*chChar = ' ';
						input.lgUnderscoreFound = true;
					}
				}

				/* change _, [, and ] to space if no filename occurs on line 
				 * >>chng 06 sep 04 use routine to check for comments 
				 * do not remove _, [, and ] in comments */
				/* check for left or right bracket, probably meant as a space */
				while( (chChar = strchr_s(chLine , '[' ) ) !=NULL )
				{
					*chChar = ' ';
					input.lgBracketFound = true;
				}

				while( (chChar = strchr_s(chLine , ']' ) ) !=NULL )
				{
					*chChar = ' ';
					input.lgBracketFound = true;
				}
			}

			/* this is trick so that cloudy input can read cloudy output */
			/* are first 25 char of input string equal to output? */
			if( strncmp(chLine,"                       * ",25) == 0 )
			{
				/* reading cloudy output, send in shifted input */
				nread = cdRead( chLine+25 );
			}
			else
			{
				/* stuff the command line into the internal stack */
				nread = cdRead( chLine  );
			}
		}

		if( lgAbort )
		{
			/* input parser hit something REALLY bad */
			lgBadExit = true;
			return(lgBadExit);
		}

		if( nread <= 0 )
			fprintf(ioQQQ," Warning: limit to number of lines exceeded, %i\n", nread);

		/* actually call the code.  This routine figures out whether the code will do
		 * a single model or be used to optimize on a spectrum, by looking for the
		 * keyword VARY on command lines.  It will call routine cloudy if no vary commands
		 * occur, and lgOptimize_do if VARY does occur.  
		 * cdDrive returns 0 if calculation is ok, 1 if problems happened */
		if( cdDrive() )
			lgBadExit = true;
		else
			lgBadExit = false;

		/* the last line of output will contain some interesting information about the model*/
		cdNwcns(
			/* abort status, this better be false, 0 */
			&lgAbort_exit,
			/* the number of warnings, cautions, notes, and surprises */
			&NumberWarnings, 
			&NumberCautions, 
			&NumberNotes, 
			&NumberSurprises, 
			/* the number of temperature convergence failures */
			&NumberTempFailures, 
			/* the number of pressure convergence failures */
			&NumberPresFailures,
			/* the number of ionization convergence failures */
			&NumberIonFailures, 
			/* the number of electron density convergence failures */
			&NumberNeFailures );

		ostringstream finalMsg;

		finalMsg << " Cloudy ends: " << nzone << " zone";
		if( nzone > 1 )
			finalMsg << "s";

		finalMsg << ", " << iteration << " iteration";
		if( iteration > 1 )
			finalMsg << "s";

		if( lgAbort_exit )
			finalMsg << ", ABORT DISASTER PROBLEM";

		if( NumberWarnings > 0 )
		{
			finalMsg << ", " << NumberWarnings << " warning";
			if( NumberWarnings > 1 )
				finalMsg << "s";
			/* this indicates error */
			lgBadExit = 1;
		}

		if( NumberCautions > 0 )
		{
			finalMsg << ", " << NumberCautions << " caution";
			if( NumberCautions > 1 )
				finalMsg << "s";
		}

		/* this flag was set in lgCheckMonitors*/
		if( !lgMonitorsOK )
		{
			finalMsg << ", ";
			/* some botches were three sigma */
			if( lgBigBotch  )
				finalMsg << "BIG ";
			finalMsg << "BOTCHED MONITORS!!!";
			/* this indicates error */
			lgBadExit = 1;
		}

		if( NumberTempFailures+NumberPresFailures+NumberIonFailures+NumberNeFailures > 0 )
		{
			finalMsg << ". Failures: " << NumberTempFailures << " thermal, ";
			finalMsg << NumberPresFailures << " pressure, ";
			finalMsg << NumberIonFailures << " ionization, ";
			finalMsg << NumberNeFailures << " electron density";
		}

		/* NB DO NOT CHANGE ANY ASPECT OF THE FOLLOWING STRINGS - THEY ARE USED TO RECORD
		 * EXEC TIME BY A PERL SCRIPT */
		if( prt.lgPrintTime )
		{
			/* print execution time [s] by default,
			 * need spaces around number so that logging perl script picks up correct number 
			 * ir_extime.pl script will delete through "ExecTime(s)" and remainder of line must be number */
			finalMsg << ".  ExecTime(s) " << fixed << setprecision(2) << cdExecTime();
		}

		if( called.lgTalk )
			fprintf( ioQQQ, "%s\n", finalMsg.str().c_str() );

		/* cdDrive returned 1 if something bad happened, and 0 if everything is ok.  We will
		 * return 0 if everything is ok, and 1 if something bad happened.*/
		cdEXIT(lgBadExit);
	}
	catch( bad_alloc )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Bailing out...\n" );
		cdPrepareExit();
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		cdPrepareExit();
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		cdPrepareExit();
	}
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT )
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
			else if( e.sig() == SIGTERM )
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
			else if( e.sig() == SIGILL )
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
			else if( e.sig() == SIGFPE )
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
			else if( e.sig() == SIGSEGV )
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
#			endif
			else
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );

		}
		cdPrepareExit();
	}
#endif
	catch( cloudy_exit& e )
	{
		if( called.lgTalk )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		cdPrepareExit();
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		cdPrepareExit();
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		cdPrepareExit();
	}

	return exit_status;
}
