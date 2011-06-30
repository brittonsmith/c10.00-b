/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* input_readarray read input commands from array where images are stored *
 * returns chCard, which will have <=80 characters before eol    *
 * line image is up and low case                                 */
/*input_init initial input_readarray array for storing line images at start of calculation */
/*lgInputComment - parse comment - check if argument is comment string */
#include "cddefines.h"
#include "trace.h"
#include "input.h"

/*lgInputComment - parse comment - check if argument is comment string, 
 * either upper or lower case -
 * returns true if line is a comment, false if not 
 * a comment is any line starting with "C ", *, %, //, or # */
bool lgInputComment( const char *chLine )
{
	bool lgReturn;

	DEBUG_ENTRY( "lgInputComment()" );

	/* should not call this routine with null line */
	if( chLine[0] == 0 )
		TotalInsanity();

	/* first case - the special characters that can start a line */
	if( chLine[0] == '#' || chLine[0] == '*' || chLine[0] == '%' || chLine[0] == ' ' )
	{
		lgReturn = true;
	}
	else if( strncmp(chLine,"//", 2 ) == 0 )
	{
		lgReturn = true;
	}
	/* second case is line that starts with c */
	else if( chLine[0] == 'C' || chLine[0] == 'c' )
	{
		/* line starts with C, could be a command or a comment,
		 * if a comment then line is "C ", but there could be a newline '\n')
		 * or carriage return '\r' in the [1] position 
		 * '\r' is carriage return, happens on cygwin gcc */
		if( chLine[1] == '\n' || chLine[1] == ' ' || chLine[1] == '\r' ) 
		{
			lgReturn = true;
		}
		else
		{
			lgReturn = false;
		}
	}
	else
	{
		lgReturn = false;
	}
	/*fprintf(ioQQQ,"DEBUG %c \n", TorF(lgReturn ) );*/

	return lgReturn;
}

/*input_init initial input_readarray array for storing line images at start of calculation */
void t_input::init(void)
{

	DEBUG_ENTRY( "t_input::init()" );

	/* this sub must be called before calling READAR to get line images
	 * it simply sets the pointer to set up reading the images
	 * */
	if( iReadWay > 0 )
	{
		/* this is usual case, read from the start of array, the commands */
		nRead = -1;
	}
	else if( iReadWay < 0 )
	{
		/* this is special case where we read from end of array, the ini file */
		/* save the current counter so we can reset it when done */
		nReadSv = nRead;

		/* and set current counter to the bottom of the stack */
		nRead = NKRD;
	}

	return;
}

void t_input::echo( FILE *ipOUT)
{
	char chCard[INPUT_LINE_LENGTH];

	/* start the file with the input commands */
	init();

	bool lgEOF = false;
	while( !lgEOF )
	{
		readarray(chCard,&lgEOF);
		if( !lgEOF )
		{
			char chCAPS[INPUT_LINE_LENGTH];
			strcpy( chCAPS , chCard );
			caps( chCAPS );
			/* keyword HIDE means to hide the command - do not print it */
			if( !nMatch( "HIDE" , chCAPS ) )
				fprintf( ipOUT, "%s\n", chCard );
		}
	}

}
/*input_readarray read input commands from array where images are stored *
 * returns chCard, which will have <=80 characters before eol    */
void t_input::readarray(char *chCard, 
  bool *lgEOF)
{
	long int last;

	DEBUG_ENTRY( "t_input::readarray()" );

	if( iReadWay > 0 )
	{
		/* usual case, reading commands from start of array
		 * nRead points to one plus the array element with the next line, it is
		 * one on the first call, which references line[0] */
		++nRead;

		/* nSave points to the last line array element that was saved,
		 * so it is one less than the number of lines read.  the last element
		 * containing a line image is [input.nSave].  There is a -1 for
		 * nRead to bring it onto the same c counting scale as nSave */
		if( nRead > nSave )
		{
			*lgEOF = true;
		}
		else
		{
			/* get the line image */
			strcpy( chCard, chCardSav[nRead] );

			*lgEOF = false;
		}
	}
	else
	{
		/* this is special case of reading cloudy.ini file, 
		 * nRead was set to 1+last image in input_init, so first time
		 * we get here it is very large.  decrement counter from end of file */
		nRead -= 1;

		/* last one with real data is NKRD+1-nSaveIni */
		last = NKRD - nSaveIni;

		/* this read is eof eof */
		if( nRead < last )
		{
			/* reset counter so we read in the proper direction */
			iReadWay = 1;
			/* pointer to next line to read.  this is on the scale where nRead-1
			 * is the actual array element */
			nRead = nReadSv+1;
		}

		/* check if we hit eof while reading in forward direction */
		if( iReadWay == 1 && nRead > nSave )
		{
			*lgEOF = true;
		}
		else
		{
			strcpy( chCard, chCardSav[nRead] );

			/* did not hit eof */
			*lgEOF = false;
		}
	}

	/* if any "trace" appeared on a command line, then this flag was set
	 * so print the input command before it is parsed */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "t_input::readarray returns=%s=\n",chCard );
	}

	return;
}

/** input_readvector: read n numbers from the file chFile and store them in vector[] */
void input_readvector(const char* chFile, /**< file name to read from */
		      double vector[],    /**< vector[n] - the numbers that were read from the input line(s) */
		      long n,             /**< number of elements in vector[] that we need to read */
		      bool* lgEOF)        /**< was EOF reached before enough numbers were read? */
{
	DEBUG_ENTRY( "input_readvector()" );

	fstream ioDATA;
	open_data( ioDATA, chFile, mode_r, AS_LOCAL_ONLY );

	for( long i=0; i < n; ++i )
		ioDATA >> vector[i];

	*lgEOF = !ioDATA.good();
	return;
}
