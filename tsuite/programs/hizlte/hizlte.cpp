/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* test case to show behavior of very high metals environmet going to lte
 * cloud irradiated by STE blackbody, should equilibriate at continuum temperature
 * output hizlte.txt gives input continuum temperature, deduced electron temperature, error */
#include "cddefines.h"
#include "cddrive.h"
/*int main( int argc, char *argv[] )*/
int main( void )
{
	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "main()" );

	try {
		int lgBAD ;
		double telog , teinc;
		double hden , te, error,BigError;

		FILE *ioRES ;
		char chLine[100];

		/* this will be limit to the number of command chLines we can still put in */
		long int nleft;

		/* calculation's results */
		ioRES = fopen("hizlte.txt","w");
		if( ioRES == NULL )
		{
			printf(" could not open collion.txt for writing.\n");
			cdEXIT(EXIT_FAILURE);
		}
		fprintf(ioRES,"T(BB)\tT(equil)\terror\n");

		/* the first temperature */
		telog = 4.;
		/* the increment in the temperature */
		teinc = 0.5;
		/* the hydrogen density that will be used */
		hden = 10.;

		BigError = 0.;
		/* this is limit on 32 bit double */
		while( telog < 7.99 )
		{
			/* initialize the code for this run */
			cdInit();
			cdTalk(false);
			/*cdNoExec( );*/
			printf("te %g ",telog);

			/* inputs */
			nleft = cdRead( "title high metal clouds into lte"  );

			nleft = cdRead( "iterate"  );
			nleft = cdRead( "set dr 0"  );
			nleft = cdRead( "abundances starburst 10 "  );

			sprintf(chLine,"blackbody lte %f ",telog);
			nleft = cdRead( chLine  );

			sprintf(chLine,"hden %f ",hden);
			nleft = cdRead( chLine  );
			nleft = cdRead( "stop zone 1 "  );
			sprintf(chLine,"hden %f ",hden);
			nleft = cdRead( chLine  );

			/* actually call the code */
			lgBAD = cdDrive();

			te = cdTemp_last();
			error = (pow(10.,telog)-te)/te;
			BigError = MAX2( BigError , fabs(error) );
			/*remember largest error */
			fprintf(ioRES,"%e\t%e\t%e\n", telog , log10(te) ,error);
			printf(" %g %g\n",log10(te),error );
			telog += teinc;
		}
		printf("biggest error was %g\n", BigError );
		fprintf(ioRES,"biggest error was %g\n", BigError );
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
