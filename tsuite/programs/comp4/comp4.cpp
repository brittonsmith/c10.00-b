/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that calls cloudy with 2 different models, called twice 
to confirm that code is properly initialized */
#include "cddefines.h"
#include "cddrive.h"
/*
This model runs two sets of models that have identical boundary conditions.  
The output and results should also be identical.  This tests that all 
variables have been properly reset.  If the results of the second 
calculation differ from the first then the code has some memory of 
previous calculations, a bad feature if it is to be used to compute 
a grid of calculations.

There are two sets of output.  The files comp4a.out and comp4b.out have the main
output from the code.  The files file1.txt and file2.txt have the punched results
files.  These files should be identical if all is well.
*/

int main( void )
{
	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "main()" );

	try {
		int lgOK ;
		long int i;
		char chCard[200];

		/* this is limit to the number of command chLines we can still put in */
		long int nleft;

		for( i=0; i<2; ++i )
		{
			if( i==0 )
				cdOutput( "comp4a.out" );
			else
				cdOutput( "comp4b.out" );

			/* initialize the code for this run */
			cdInit();
			/*this is a very simple constant temp model */
			nleft = cdRead("test " );
			/* must have some grains to malloc their space in this grid */
			nleft = cdRead("grains ism abundance -10 " );
			nleft = cdRead("no times " );
			/* write results to either file1.txt or file2.txt */
			nleft = cdRead("print lines column  " );
			sprintf(chCard , "punch results column \"file%li.txt\" hide ",i+1);
			nleft = cdRead( chCard );
			lgOK = cdDrive();
			/* end of the first model */

			/* start of the second model, fully molecular */
			cdInit();
			nleft = cdRead("blackbody 5000 " );
			nleft = cdRead("luminosity total solar linear 2 " );
			nleft = cdRead("brems 6 " );
			nleft = cdRead("luminosity total solar log -2.7 " );
			nleft = cdRead("hden 10 " );
			nleft = cdRead("abundances ism " );
			nleft = cdRead("stop temperature 10K linear " );
			nleft = cdRead("radius 15.8  " );
			nleft = cdRead("stop zone 1 " );
			nleft = cdRead("no times " );
			nleft = cdRead("print lines column  " );
			/* write results to either file1.txt or file2.txt */
			sprintf(chCard , "punch results column \"file%li.txt\" hide ",i+1);
			lgOK = cdDrive();
			/* end of the second model */

			/* start of the third model */
			cdInit();
			/* inputs */
			nleft = cdRead( "ioniz -1 "  );
			nleft = cdRead( "sphere "  );
			nleft = cdRead( "abundances ism "  );
			nleft = cdRead( "table agn "  );
			nleft = cdRead( "hden 11 " );
			nleft = cdRead( "stop column density 19 " );
			nleft = cdRead( "stop zone 2 " );
			nleft = cdRead( "no times " );
			nleft = cdRead( "print lines column  " );
			/* write results to either file1.txt or file2.txt */
			sprintf(chCard , "punch results column \"file%li.txt\" hide ",i+1);
			/* actually call the code */
			lgOK = cdDrive();
		}

		// redirect to stdout
		cdOutput();

		printf("calculations complete, now compare comp4a.out and comp4b.out, and file1.txt and file2.txt\n");
		printf("vsdiff comp4a.out comp4b.out\n");
		printf("vsdiff file1.txt file2.txt\n");

		cdEXIT(lgOK);
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

