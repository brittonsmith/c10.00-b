/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "mpi_utilities.h"
#include "save.h"
#include "grid.h"

#ifdef MPI_ENABLED

namespace MPI
{
	void load_balance::init(int nJobs)
	{
		ASSERT( nJobs > 0 );
		p_jobs.resize( nJobs );
		p_ptr = COMM_WORLD.Get_rank();
		// the master rank will now set up a random sequence for the jobs
		// this way we hope to get statistical load balancing of the ranks
		if( p_ptr == 0 )
		{
			for( int i=0; i < nJobs; ++i )
				p_jobs[i] = i;
			// This may or may not seed the random number generator used by
			// random_shuffle. There is no portable C++ interface to do this :-(
			srand( unsigned( time(NULL) ) );
			random_shuffle( p_jobs.begin(), p_jobs.end() );
		}
		// now broadcast the random sequence to the other ranks...
		COMM_WORLD.Bcast( &p_jobs[0], nJobs, type(p_jobs[0]), 0 );
	}
}

#endif /* MPI_ENABLED */

/** process_output: concatenate output files produced in MPI grid run */
void process_output()
{
	DEBUG_ENTRY( "process_output()" );

	// NOTE: when this routine is called all file handles have already been closed

	string main_input = save.chRedirectPrefix + ".in";
	string main_output = save.chRedirectPrefix + ".out";

	// first process main output files
	FILE* main_output_handle = open_data( main_output.c_str(), "a", AS_LOCAL_ONLY );
	for( int j=0; j < grid.totNumModels; ++j )
	{
		string Base = GridPointPrefix(j) + save.chRedirectPrefix;
		string in_name = Base + ".in";
		remove( in_name.c_str() );
		string out_name = Base + ".out";
		append_file( main_output_handle, out_name.c_str() );
		remove( out_name.c_str() );
	}
	fclose( main_output_handle );

	fstream main_input_handle;
	open_data( main_input_handle, main_input.c_str(), mode_r, AS_LOCAL_ONLY );
	string line;

	int ipPun = 0;

	while( getline( main_input_handle, line ) )
	{
		string caps_line;
		// create all caps version
		string::const_iterator p = line.begin();
		while( p != line.end() )
			caps_line.push_back( toupper(*p++) );
		if( caps_line.compare( 0, 4, "SAVE" ) == 0 || caps_line.compare( 0, 4, "PUNC" ) == 0 )
		{
			ASSERT( ipPun < save.nsave );
			string fnam = save.chFilenamePrefix;
			string::size_type p = line.find( '"' );
			fnam += line.substr( ++p );
			fnam.erase( fnam.find( '"' ) );
			// open in binary mode in case we are writing a FITS file
			FILE *dest = open_data( fnam.c_str(), "ab", AS_LOCAL_ONLY_TRY );
			if( dest != NULL )
			{
				if( save.lgSaveToSeparateFiles[ipPun] )
				{
					// keep the save files for each grid point separate
					// the main save file contains the save header
					// salvage it by prepending it to the first save file
					// this gives the same behavior as in non-MPI runs
					string gridnam = GridPointPrefix(0) + fnam;
					append_file( dest, gridnam.c_str() );
					fclose( dest );
					dest = NULL;
					// this will overwrite the old file gridnam
					rename( fnam.c_str(), gridnam.c_str() );
				}
				else
				{
					// concatenate the save files for each grid point
					for( int j=0; j < grid.totNumModels; ++j )
					{
						string gridnam = GridPointPrefix(j) + fnam;
						append_file( dest, gridnam.c_str() );
						remove( gridnam.c_str() );
					}
				}
				if( caps_line.find( "XSPE", 4 ) != string::npos )
				{
					// dest points to an empty file, so generate the complete FITS file now
					ASSERT( save.FITStype[ipPun] >= 0 &&
						save.FITStype[ipPun] < NUM_OUTPUT_TYPES );
					saveFITSfile( dest, save.FITStype[ipPun] );
					fseek( dest, 0, SEEK_END );
					ASSERT( ftell(dest)%2880 == 0 );
				}
				if( dest != NULL )
				{
					fclose( dest );
				}
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM - could not open file %s\n", fnam.c_str() );
			}
			++ipPun;
		}
	}
}

/** append_file: append output produced on file <source> to open file descriptor <dest> */
void append_file( FILE *dest, const char *source )
{
	DEBUG_ENTRY( "append_file()" );

	FILE *src = open_data( source, "rb", AS_LOCAL_ONLY_TRY );
	if( src == NULL )
		return;
	int chr = fgetc(src);
	while( ! feof(src) )
	{
		fputc( chr, dest );
		chr = fgetc(src);
	}
	fclose(src);
	return;
}
