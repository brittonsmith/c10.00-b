/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* Template main program for calling Cloudy as a subroutine
 * Routine returns 0 if model is ok, and 1 if problems occurred. */
#include "cddefines.h"
#include "cddrive.h"

int main( int argc, char *argv[] )
{
	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "main()" ); // do not remove this!

	try {
		// ============================================================================
		// START ENTERING YOUR CODE AFTER THIS LINE
		// inside this try block you can enter your code to call Cloudy as a subroutine
		// replace the code between here and the next ===== line with your own program

		bool lgBad;
		long nleft;

		// you can call Cloudy in a loop if you want...
		for( int i=0; i < 1; ++i )
		{
			// the code always needs to be initialized first
			cdInit();
			// replace this with a series of calls to cdRead to define the input script
			// this particular command line exercises the smoke test
			nleft = cdRead( "test" );
			// this calls Cloudy to execute the input script you defined above
			lgBad = cdDrive();
		}

		cdEXIT(lgBad); // always exit with cdEXIT, this assures files are properly closed.
		// ============================================================================
		// THIS IS THE END OF THE PROGRAM YOU WROTE
		// this ends the try block - your code calling Cloudy ends above this line
	}
	// here we catch all the possible exceptions that the code can throw, do not remove this!
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
		MyAssert( e.file(), e.line(), e.comment() );
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
		if( ioQQQ != NULL )
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
