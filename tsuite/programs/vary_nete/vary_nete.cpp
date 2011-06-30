/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that calls cloudy when used as a stand-alone program */
#include "cddefines.h"
#include "cddrive.h"

/*int main( int argc, char *argv[] )*/
int main( void )
{
	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "main()" );

	try {
		int lgOK ;
		double hdenLimit , hdenInit , hden, TeInc ,
			r5007,r4363,r1665,ro3_88,ro3_52, r1661 ,
			absolute , temp, TeLimit, TeInit , hdenInc ;

		FILE *ioDATA ;
		char chLine[100];

		/* this is limit to the number of command chLines we can still put in */
		long int nleft;

		cdOutput( "vary_nete.out" );

		/* calculation's results will go to this file*/
		if( (ioDATA = fopen("vary_nete.txt","w")) == NULL )
		{
			printf(" could not open varyTF.txt for writing.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* the range of density, and the increment in density, for this grid */
		hdenInit = 0.;
		hdenLimit = 7.;
		hdenInc = 1.;

		/* the range of temperatures for the grid */
		TeInit = 5000.;
		TeLimit = 30000.;
		/* multiplicative inc */
		TeInc = 1.5;

		/* set the density and temperature to the initial values */
		hden = hdenInit;
		temp = TeInit;

		/* print the header for the data file */
		fprintf(ioDATA,
			"density\ttemp\tOIII 5007\t4363\t1665\t88\t52\t1661\n");

		while( hden < 1.01 * hdenLimit )
		{
			while( temp < 1.01*TeLimit )
			{

				/* initialize the code for this run */
				cdInit();

				/* option to not execute the code, uncomment when debugging setup */
				/*cdNoexec( );*/

				/* gas temperature for this calculation */
				sprintf(chLine,"constant temperature %f", temp);
				nleft = cdRead( chLine  );

				/* gas density for this calculation */
				sprintf(chLine,"hden %f", hden);
				nleft = cdRead( chLine  );

				/* only want the first zone */
				nleft = cdRead( "stop zone 1 "  );

				/* speed up calculation by not included some elements and lines */
				nleft = cdRead( "init file \"fast.ini\" "  );

				nleft = cdRead( "normalize to \"o  3\" 5007 "  );

				/* an incident continuum must be specified to get the code
				 * to run at all - not very important since we will set
				 * the gas temperature */
				nleft = cdRead( "blackbody 40000  "  );
				nleft = cdRead( "ionization parameter -3  "  );

				/* actually call the code */
				lgOK = cdDrive();
				/* flush the output so we see it on the screen */
				fflush(ioQQQ);

				fprintf(ioDATA,"%.3e\t%.3e",	hden, temp );
				fprintf(stderr,"%.3e\t%.3e",	hden, temp );

				/************************ O III ******************/
				lgOK = cdLine( "O  3", 5007 , &r5007 , &absolute );
				if( lgOK==false )
					printf("did not find O  3 5007\n");
				/* this should be unity */
				fprintf(ioDATA,"\t%.3e",	r5007 );
				fprintf(stderr,"\t%.3e",	r5007 );

				lgOK = cdLine( "totl", 4363 , &r4363 , &absolute );
				if( lgOK==false )
					printf("did not find totl 4363\n");
				fprintf(ioDATA,"\t%.3e",	r4363 );
				fprintf(stderr,"\t%.3e",	r4363 );

				lgOK = cdLine( "totl", 1665 , &r1665 , &absolute );
				if( lgOK==false )
					printf("did not find totl 1665\n");
				fprintf(ioDATA,"\t%.3e",	r1665 );
				fprintf(stderr,"\t%.3e",	r1665 );

				lgOK = cdLine( "O  3", 883300 , &ro3_88 , &absolute );
				if( lgOK==false )
					printf("did not find O  3 88\n");
				fprintf(ioDATA,"\t%.3e",	ro3_88 );
				fprintf(stderr,"\t%.3e",	ro3_88 );

				lgOK = cdLine( "O  3", 518000 , &ro3_52 , &absolute );
				if( lgOK==false )
					printf("did not find O  3 52\n");
				fprintf(ioDATA,"\t%.3e",	ro3_52 );
				fprintf(stderr,"\t%.3e",	ro3_52 );

				lgOK = cdLine( "O  3", 1661 , &r1661 , &absolute );
				if( lgOK==false )
					printf("did not find O  3 1661\n");
				fprintf(ioDATA,"\t%.3e",	r1661 );
				fprintf(stderr,"\t%.3e",	r1661 );


				/************************* end lines with lf and flush it ***************/
				fprintf(ioDATA,"\n" );
				fflush(ioDATA );
				fprintf(stderr,"\n" );

				temp *= TeInc;
			}
			temp = TeInit;
			hden += hdenInc;
		}

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

