/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* test case to vary density at constant temperature, to show how line ratios change */
#include "cddefines.h"
#include "cddrive.h"
/*int main( int argc, char *argv[] )*/
int main( void )
{
	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "main()" );

	try {
		int lgBAD ;
		double telog ;
		double hden , hdeninc;

		FILE *ioRES ;
		char chLine[100];

		/* this will be limit to the number of command chLines we can still put in */
		long int nleft;

		/* calculation's results */
		if( (ioRES = fopen("varyn.txt","w")) == NULL )
		{
			printf(" could not open collion.txt for writing.\n");
			cdEXIT(EXIT_FAILURE);
		}
		fprintf(ioRES,"density\t[OII] 3726\t3729\t2471\t7323\t7332\n");

		/* the first temperature */
		telog = 4.;
		/* the increment in the temperature */
		hdeninc = 0.5;
		/* the hydrogen density that will be used */
		hden = -1.;

		/* this is limit on 32 bit double */
		while( hden < 8. )
		{
			double relint , absint ;
			/* initialize the code for this run */
			cdInit();
			cdTalk(false);
			/*cdNoExec( );*/
			printf("hden %g\n",hden);
			fprintf(ioRES,"%g",hden);

			/* inputs */
			nleft = cdRead( "title vary density at constant temperature"  );

			sprintf(chLine,"hden %f ",hden);
			nleft = cdRead( chLine  );

			nleft = cdRead( "table agn "  );
			nleft = cdRead( "stop zone 1 "  );
			nleft = cdRead( "ionization parameter -3 "  );
			nleft = cdRead( "constant temperature 4 "  );
			nleft = cdRead( "normalize to \"o ii\" 3726 "  );

			/* actually call the code */
			lgBAD = cdDrive();

			/* now get the lines */
			cdLine( "o ii" , 3726 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "o ii" , 3729 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "o ii" , 2471 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "o ii" , 7323 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "o ii" , 7332 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			fprintf(ioRES, "\n");

			hden += hdeninc;
		}

		cdEXIT(lgBAD);
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

