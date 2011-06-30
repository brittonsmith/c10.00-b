/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"

/** dump the array to a file in binary format; the file must
 *  already have been opened prior to calling this method */
void do_dump_state(const void* buf, size_t nelem, size_t size, FILE* out, int32 magic)
{
	DEBUG_ENTRY( "do_dump_state()" );

	bool lgErr = ( fwrite( &magic, sizeof(int32), 1, out ) != 1 );
	int32 help = (int32)sizeof(size_t);
	lgErr = lgErr || ( fwrite( &help, sizeof(int32), 1, out ) != 1 );
	lgErr = lgErr || ( fwrite( &size, sizeof(size_t), 1, out ) != 1 );
	lgErr = lgErr || ( fwrite( buf, size, nelem, out ) != nelem );
	if( lgErr )
	{
		fprintf( ioQQQ, " I/O error while dumping state!\n" );
		cdEXIT(EXIT_FAILURE);
	}
}

/** restore the array from a file in binary format; the file must already
 *  have been opened prior to calling this method and the array must already
 *  have been allocated in exactly the same way as when it was dumped; some
 *  checks are performed, but not every error is excluded */
void do_restore_state(void* buf, size_t nelem, size_t size, FILE *in, int32 magic)
{
	DEBUG_ENTRY( "do_restore_state()" );

	int32 help = 0;
	size_t help2 = 0;
	bool lgErr = ( fread( &help, sizeof(int32), 1, in ) != 1 );
	// this checks for correct version and prevents mixing up old style and new style data
	// it also prevents mixing up data from big-endian and little-endian machines.
	lgErr = lgErr || ( help != magic );
	lgErr = lgErr || ( fread( &help, sizeof(int32), 1, in ) != 1 );
	// this prevents mixing up data from 32-bit and 64-bit systems
	lgErr = lgErr || ( help != (int32)sizeof(size_t) );
	lgErr = lgErr || ( fread( &help2, sizeof(size_t), 1, in ) != 1 );
	// this may guard against reading an older, incompatible version of the array
	lgErr = lgErr || ( help2 != size );
	lgErr = lgErr || ( fread( buf, size, nelem, in ) != nelem );
	if( lgErr )
	{
		fprintf( ioQQQ, " Error while restoring state!\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
