/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* 
* a set of  routines that are widely used across the code for various
* housekeeping chores.  These do not do any physics and are unlikely to
* change over time.  The prototypes are in cddefines.h and so are 
* automatically picked up by all routines 
*/
/*FFmtRead scan input line for free format number */
/*e2 second exponential integral */
/*caps convert input command line (through eol) to ALL CAPS */
/*ShowMe produce request to send information to GJF after a crash */
/*AnuUnit produce continuum energy in arbitrary units */
/*cap4 convert first 4 char of input line chLab into chCAP all in caps, null termination */
/*insane set flag saying that insanity has occurred */
/*nMatch determine whether match to a keyword occurs on command line,
 * return value is 0 if no match, and position of match within string if hit */
/*fudge enter fudge factors, or some arbitrary number, with fudge command*/
/*GetQuote get any name between double quotes off command line
 * return string as chLabel, is null terminated */
/*qip compute pow(x,n) for positive integer n through repeated squares */
/*dsexp safe exponential function for doubles */
/*sexp safe exponential function */
/*TestCode set flag saying that test code is in place */
/*CodeReview - placed next to code that needs to be checked */
/*fixit - say that code needs to be fixed */
/*broken set flag saying that the code is broken, */
/*dbg_printf is a debug print routine that was provided by Peter Teuben,
 * as a component from his NEMO package.  It offers run-time specification
 * of the level of debugging */
/*qg32 32 point Gaussian quadrature, original Fortran given to Gary F by Jim Lattimer */
/*TotalInsanity general error handler for something that cannot happen */
/*BadRead general error handler for trying to read data, but failing */
/*MyMalloc wrapper for malloc().  Returns a good pointer or dies. */
/*MyCalloc wrapper for calloc().  Returns a good pointer or dies. */
/*spsort netlib routine to sort array returning sorted indices */
/*chLineLbl use information in line transfer arrays to generate a line label *
 * this label is null terminated */
/*chIonLbl use information in line array to generate a null terminated ion label in "Fe 2" */
/*csphot returns photoionization cross section from opacity stage using std pointers */
/*MyAssert a version of assert that fails gracefully */
/*RandGauss normal random variate generator */
/*MyGaussRand a wrapper for RandGauss, see below */
#include <cstdarg>	/* ANSI variable arg macros */
#include "cddefines.h"
#include "physconst.h"
#include "cddrive.h"
#include "called.h"
#include "opacity.h"
#include "rfield.h"
#include "hextra.h"
#include "struc.h"
#include "hmi.h"
#include "fudgec.h"
#include "broke.h"
#include "trace.h"
#include "input.h"
#include "save.h"
#include "version.h"
#include "warnings.h"
#include "conv.h"
#include "thirdparty.h"

/*read_whole_line - safe version of fgets - read a line, 
 * return null if cannot read line or if input line is too long */
char *read_whole_line( char *chLine , int nChar , FILE *ioIN )
{
	char *chRet;

	DEBUG_ENTRY( "read_whole_line()" );

	/* wipe the buffer to prevent the code from accidentally picking up on previous input */
	memset( chLine, 0, (size_t)nChar );

	/* this always writes a '\0' character, even if line does not fit in buffer
	 * the terminating newline is copied only if the line does fit in the buffer */
	if( (chRet = fgets( chLine, nChar, ioIN )) == NULL )
		return NULL;

	/* check where the '\n' character was written, if it is not there, the line was truncated */
	char *chEOL = (char*)memchr( chLine, '\n', nChar );

	/* return null if input string is longer than nChar-1 (including terminating newline),
	 * the longest we can read. Print and return null but chLine still has as much of
	 * the line as could be placed in the buffer */
	if( chEOL == NULL )
	{
		if( called.lgTalk )
			fprintf( ioQQQ, "DISASTER PROBLEM read_whole_line found input"
			" with a line too long to be read, limit is %i char.  "
			"Start of line follows:\n%s\n",
			nChar , chLine );

		lgAbort = true;
		return NULL;
	}
	return chRet;
}

/** Split: split a string into substrings using "sep" as separator */
void Split(const string& str,   // input string
	   const string& sep,   // separator, may be multiple characters
	   vector<string>& lst, // the separated items will be appended here
	   split_mode mode)     // SPM_RELAX, SPM_KEEP_EMPTY, or SPM_STRICT; see cddefines.h
{
	DEBUG_ENTRY( "Split()" );

	bool lgStrict = ( mode == SPM_STRICT );
	bool lgKeep = ( mode == SPM_KEEP_EMPTY );
	bool lgFail = false;
	string::size_type ptr1 = 0;
	string::size_type ptr2 = str.find( sep );
	string sstr = str.substr( ptr1, ptr2-ptr1 );
	if( sstr.length() > 0 )
		lst.push_back( sstr );
	else {
		if( lgStrict ) lgFail = true;
		if( lgKeep ) lst.push_back( sstr );
	}
	while( ptr2 != string::npos ) {
		// the separator is skipped
		ptr1 = ptr2 + sep.length();
		if( ptr1 < str.length() ) {
			ptr2 = str.find( sep, ptr1 );
			sstr = str.substr( ptr1, ptr2-ptr1 );
			if( sstr.length() > 0 )
				lst.push_back( sstr );
			else {
				if( lgStrict ) lgFail = true;
				if( lgKeep ) lst.push_back( sstr );
			}
		}
		else {
			ptr2 = string::npos;
			if( lgStrict ) lgFail = true;
			if( lgKeep ) lst.push_back( "" );
		}
	}
	if( lgFail )
	{
		fprintf( ioQQQ, " A syntax error occurred while splitting the string: \"%s\"\n", str.c_str() );
		fprintf( ioQQQ, " The separator is \"%s\". Empty substrings are not allowed.\n", sep.c_str() );
		cdEXIT(EXIT_FAILURE);
	}
}

/* a version of assert that fails gracefully */
void MyAssert(const char *file, int line, const char *comment)
{
	DEBUG_ENTRY( "MyAssert()" );

	fprintf(ioQQQ,"\n\n\n PROBLEM DISASTER\n An assert has been thrown, this is bad.\n");
	fprintf(ioQQQ," %s\n",comment);
	fprintf(ioQQQ," It happened in the file %s at line number %i\n", file, line );
	fprintf(ioQQQ," This is iteration %li, nzone %li, fzone %.2f, lgSearch=%c.\n", 
		iteration , 
		nzone ,
		fnzone ,
		TorF(conv.lgSearch) );

	ShowMe();
#	ifdef OLD_ASSERT
	cdEXIT(EXIT_FAILURE);
#	endif
}

/*AnuUnit produce continuum energy in arbitrary units */
double AnuUnit(realnum energy_ryd)
/*static double AnuUnit(long int ip)*/
{
	DEBUG_ENTRY( "AnuUnit()" );

	return Energy((double)energy_ryd).get(save.chConPunEnr[save.ipConPun]);
}

/*ShowMe produce request to send information to GJF after a crash */
void ShowMe(void)
{

	DEBUG_ENTRY( "ShowMe()" );

	/* print info if output unit is defined */
	if( ioQQQ != NULL )
	{
		/* >>chng 06 mar 02 - check if molecular but cosmic rays are ignored */
		if( (hextra.cryden == 0.) && hmi.BiggestH2 > 0.1 )
		{
			fprintf( ioQQQ, " >>> \n >>> \n >>> Cosmic rays are not included and the gas is molecular.  "
				"THIS IS KNOWN TO BE UNSTABLE.  Add cosmic rays and try again.\n >>> \n >>>\n\n");
		}
		else
		{
			fprintf( ioQQQ, "\n\n\n" );
			fprintf( ioQQQ, "           vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv \n" );
			fprintf( ioQQQ, "          > PROBLEM DISASTER PROBLEM DISASTER.      <\n" );
			fprintf( ioQQQ, "          > Sorry, something bad has happened.      <\n" );
			fprintf( ioQQQ, "          > Please post this on the Cloudy web site <\n" );
			fprintf( ioQQQ, "          > discussion board at www.nublado.org     <\n" );
			fprintf( ioQQQ, "          > Please send all following information:  <\n" );
			fprintf( ioQQQ, "           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" );
			fprintf( ioQQQ, "\n\n" );


			fprintf( ioQQQ, " Cloudy version number is %s\n", 
				t_version::Inst().chVersion );
			fprintf( ioQQQ, " %s\n\n", t_version::Inst().chInfo );

			fprintf( ioQQQ, "%5ld warnings,%3ld cautions,%3ld temperature failures.  Messages follow.\n", 
			  warnings.nwarn, warnings.ncaun, conv.nTeFail );

			/* print the warnings first */
			cdWarnings(ioQQQ);

			/* now print the cautions */
			cdCautions(ioQQQ);

			/* now output the commands */
			cdPrintCommands(ioQQQ);

			/* if init command was present, this is the number of lines in it -
			 * if no init then still set to zero as done in cdInit */
			if( input.nSaveIni )
			{
				fprintf(ioQQQ," This input stream included an init file.\n");
				fprintf(ioQQQ," If this init file is not part of the standard Cloudy distribution\n"); 
				fprintf(ioQQQ," then I will need a copy of it too.\n");
			}
		}
	}
	return;
}

/*cap4 convert first 4 char of input line chLab into chCAP all in caps, null termination */
void cap4(
		char *chCAP ,	/* output string, cap'd first 4 char of chLab, */
						/* with null terminating */
		const char *chLab)	/* input string ending with eol*/
{
	long int /*chr,*/ 
	  i;

	DEBUG_ENTRY( "cap4()" );

	/* convert 4 character string in chLab to ALL CAPS in chCAP */
	for( i=0; i < 4; i++ )
	{
		/* toupper is function in ctype that converts to upper case */
		chCAP[i] = toupper( chLab[i] );
	}

	/* now end string with eol */
	chCAP[4] = '\0';
	return;
}

/*uncaps convert input command line (through eol) to all lowercase */
void uncaps(char *chCard )
{
	long int i;

	DEBUG_ENTRY( "caps()" );

	/* convert full character string in chCard to ALL CAPS */
	i = 0;
	while( chCard[i]!= '\0' )
	{
		chCard[i] = tolower( chCard[i] );
		++i;
	}
	return;
}

/*caps convert input command line (through eol) to ALL CAPS */
void caps(char *chCard )
{
	long int i;

	DEBUG_ENTRY( "caps()" );

	/* convert full character string in chCard to ALL CAPS */
	i = 0;
	while( chCard[i]!= '\0' )
	{
		chCard[i] = toupper( chCard[i] );
		++i;
	}
	return;
}

/*e2 second exponential integral */
/*>>chng 07 jan 17, PvH discover that exp-t is not really
 * exp-t - this changed results in several tests */
double e2(
	/* the argument to E2 */
	double t )
{
	/* use recurrence relation */
	/* ignore exp_mt, it is *very* unreliable */
	double hold = sexp(t) - t*ee1(t);
	DEBUG_ENTRY( "e2()" );
	/* guard against negative results, this can happen for very large t */
	return max( hold, 0. );
}

/*ee1 first exponential integral */
double ee1(double x)
{
	double ans, 
	  bot, 
	  top;
	static double a[6]={-.57721566,.99999193,-.24991055,.05519968,-.00976004,
	  .00107857};
	static double b[4]={8.5733287401,18.0590169730,8.6347608925,.2677737343};
	static double c[4]={9.5733223454,25.6329561486,21.0996530827,3.9584969228};

	DEBUG_ENTRY( "ee1()" );

	/* computes the exponential integral E1(x),
	 * from Abramowitz and Stegun
	 * stops with error condition for negative argument,
	 * returns zero in large x limit 
	 * */

	/* error - does not accept negative arguments */
	if( x <= 0 )
	{
		fprintf( ioQQQ, " DISASTER negative argument in function ee1, x<0\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* branch for x less than 1 */
	else if( x < 1. )
	{
		/* abs. accuracy better than 2e-7 */
		ans = ((((a[5]*x + a[4])*x + a[3])*x + a[2])*x + a[1])*x + a[0] - log(x);
	}

	/* branch for x greater than, or equal to, one */
	else
	{
		/* abs. accuracy better than 2e-8 */
		top = (((x + b[0])*x + b[1])*x + b[2])*x + b[3];
		bot = (((x + c[0])*x + c[1])*x + c[2])*x + c[3];
		ans = top/bot/x*exp(-x);
	}
	return ans;
}

/* same as ee1, except without factor of exp(x) in returned value	*/
double ee1_safe(double x)
{
	double ans, 
	  bot, 
	  top;
	/*static double a[6]={-.57721566,.99999193,-.24991055,.05519968,-.00976004,
	  .00107857};*/
	static double b[4]={8.5733287401,18.0590169730,8.6347608925,.2677737343};
	static double c[4]={9.5733223454,25.6329561486,21.0996530827,3.9584969228};

	DEBUG_ENTRY( "ee1_safe()" );

	ASSERT( x > 1. );

	/* abs. accuracy better than 2e-8 */
	/*	top = powi(x,4) + b[0]*powi(x,3) + b[1]*x*x + b[2]*x + b[3]; */
	top = (((x + b[0])*x + b[1])*x + b[2])*x + b[3];
	/*	bot = powi(x,4) + c[0]*powi(x,3) + c[1]*x*x + c[2]*x + c[3]; */
	bot = (((x + c[0])*x + c[1])*x + c[2])*x + c[3];

	ans = top/bot/x;
	return ans;
}

/*FFmtRead scan input line for free format number */
double FFmtRead(const char *chCard, 
		long int *ipnt, 
		/* the contents of this array element is the last that will be read */
		long int last, 
		bool *lgEOL)
{
	DEBUG_ENTRY( "FFmtRead()" );

	char chr = '\0';
	const char *eol_ptr = &chCard[last]; // eol_ptr points one beyond last valid char
	const char *ptr = min(&chCard[*ipnt-1],eol_ptr); // ipnt is on fortran scale

	ASSERT( *ipnt > 0 && last - *ipnt < INPUT_LINE_LENGTH );

	while( ptr < eol_ptr && ( chr = *ptr++ ) != '\0' )
	{
		const char *lptr = ptr;
		char lchr = chr;
		if( lchr == '-' || lchr == '+' )
			lchr = *lptr++;
		if( lchr == '.' )
			lchr = *lptr;
		if( isdigit(lchr) )
			break;
	}

	if( ptr == eol_ptr || chr == '\0' )
	{
		*ipnt = last+1;
		*lgEOL = true;
		return 0.;
	}

	string chNumber;
	bool lgCommaFound = false, lgLastComma = false;
	do
	{
		lgCommaFound = lgLastComma;
		if( chr != ',' )
		{
			chNumber += chr;
		}
		else
		{
			/* don't complain about comma if it appears after number,
				as determined by exiting loop before this sets lgCommaFound */
			lgLastComma = true;

		}
		if( ptr == eol_ptr )
			break;
		chr = *ptr++;
	}
	while( isdigit(chr) || chr == '.' || chr == '-' || chr == '+' || chr == ',' || chr == 'e' || chr == 'E' );

	if( lgCommaFound )
	{
		fprintf( ioQQQ, " PROBLEM - a comma was found embedded in a number, this is deprecated.\n" );
		fprintf(ioQQQ, "== %-80s ==\n",chCard);
	}

	double value = atof(chNumber.c_str());

	*ipnt = (long)(ptr - chCard); // ptr already points 1 beyond where next read should start
	*lgEOL = false;
	return value;
}

/*nMatch determine whether match to a keyword occurs on command line,
 * return value is 0 if no match, and position of match within string if hit */
long nMatch(const char *chKey, 
	    const char *chCard)
{
	const char *ptr;
	long Match_v;

	DEBUG_ENTRY( "nMatch()" );

	ASSERT( strlen(chKey) > 0 );

	if( ( ptr = strstr_s( chCard, chKey ) ) == NULL )
	{
		/* did not find match, return 0 */
		Match_v = 0L;
	}
	else
	{
		/* return position within chCard (fortran scale) */
		Match_v = (long)(ptr-chCard+1);
	}
	return Match_v;
}

/* fudge enter fudge factors, or some arbitrary number, with fudge command
 * other sections of the code access these numbers by calling fudge
 * fudge(0) returns the first number that was entered
 * prototype for this function is in cddefines.h so that it can be used without
 * declarations 
 * fudge(-1) queries the routine for the number of fudge parameters that were entered,
 * zero returned if none */
double fudge(long int ipnt)
{
	double fudge_v;

	DEBUG_ENTRY( "fudge()" );

	if( ipnt < 0 )
	{
		/* this is special case, return number of arguments */
		fudge_v = fudgec.nfudge;
		fudgec.lgFudgeUsed = true;
	}
	else if( ipnt >= fudgec.nfudge )
	{
		fprintf( ioQQQ, " FUDGE factor not entered for array number %3ld\n", 
		  ipnt );
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		fudge_v = fudgec.fudgea[ipnt];
		fudgec.lgFudgeUsed = true;
	}
	return fudge_v;
}


/* GetQuote get any name between double quotes off command line
 * sets chLabel - null terminated string as chLabel
 * sets string between quotes to spaces
 * returns zero for success, 1 for did not find double quotes 
 * lgAbort - true then above, false only set flag*/
int GetQuote(
		/* we will generate a label and stuff it here */
		char *chStringOut,	
		/* line image, we will blank out anything between quotes */
		char *chCard,
		char *chCardRaw,
		/* if true then abort if no double quotes found, 
		 * if false then return null string in this case */
		bool lgAbort )
{
	char *i0,       /* pointer to first " */
	  *i1,          /* pointer to second ", name is in between */
	  *iLoc;        /* pointer to first " in local version of card in calling routine */
	size_t len;

	DEBUG_ENTRY( "GetQuote()" );

	/* find first quote start of string, string begins and ends with quotes */
	i0 = strchr_s( chCardRaw,'\"' );

	if( i0 != NULL ) 
	{
		/* get pointer to next quote */
		i1 = strchr_s( i0+1 ,'\"' );
	}
	else 
	{
		i1 = NULL;
	}

	/* check that pointers were not NULL */
	/* >>chng 00 apr 27, check for i0 and i1 equal not needed anymore, by PvH */
	if( i0 == NULL || i1 == NULL )
	{
		if( lgAbort )
		{
			/* lgAbort true, must abort if no string found */
			fprintf( ioQQQ, 
				" A filename or label must be specified within double quotes, but no quotes were encountered on this command.\n" );
			fprintf( ioQQQ, " Name must be surrounded by exactly two double quotes, like \"name.txt\". \n" );
			fprintf( ioQQQ, " The line image follows:\n" );
			fprintf( ioQQQ, " %s\n", chCardRaw);
			fprintf( ioQQQ, " Sorry\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			/* this branch, ok if not present, return null string in that case */
			chStringOut[0] = '\0';
			/* return value of 1 indicates did not find double quotes */
			return 1;
		}
	}

	/* now copy the text in between quotes */
	len = (size_t)(i1-i0-1);
	strncpy(chStringOut,i0+1,len);
	/* strncpy doesn't terminate the label */
	chStringOut[len] = '\0';

	/* get pointer to first quote in local copy of line image in calling routine */
	iLoc = strchr_s( chCard, '\"' );
	if( iLoc == NULL )
		TotalInsanity();

	// blank out label once finished, to not be picked up later
	// erase quotes as well, so that we can find second label, by PvH
	while( i0 <= i1 )
	{
		*i0 = ' ';
		*iLoc = ' ';
		++i0;
		++iLoc;
	}
	/* return condition of 0 indicates success */
	return 0;
}

/* want to define this only if no native os support exists */
#ifndef HAVE_POWI

/* powi.c - calc x^n, where n is an integer! */

/* Very slightly modified version of power() from Computer Language, Sept. 86,
	pg 87, by Jon Snader (who extracted the binary algorithm from Donald Knuth,
	"The Art of Computer Programming", vol 2, 1969).
	powi() will only be called when an exponentiation with an integer
	exponent is performed, thus tests & code for fractional exponents were 
	removed.
 */

double powi( double x, long int n )	/* returns:  x^n */
/* x;	 base */
/* n;	 exponent */
{
	double p;	/* holds partial product */

	DEBUG_ENTRY( "powi()" );

	if( x == 0 )
		return 0.;

	/* test for negative exponent */
	if( n < 0 )
	{	
		n = -n;
		x = 1/x;
	}

	p = is_odd(n) ? x : 1;	/* test & set zero power */

	/*lint -e704 shift right of signed quantity */
	/*lint -e720 Boolean test of assignment */
	while( n >>= 1 )
	{	/* now do the other powers */
		x *= x;			/* sq previous power of x */
		if( is_odd(n) )	/* if low order bit set */
			p *= x;		/*	then, multiply partial product by latest power of x */
	}
	/*lint +e704 shift right of signed quantity */
	/*lint +e720 Boolean test of assignment */
	return p;
}

#endif /* HAVE_POWI */

long ipow( long m, long n )	/* returns:  m^n */
/* m;		 base */
/* n;		 exponent */
{
	long p;	/* holds partial product */

	DEBUG_ENTRY( "ipow()" );

	if( m == 0 || (n < 0 && m > 1) )
		return 0L;
	/* NOTE: negative exponent always results in 0 for integers!
	 * (except for the case when m==1 or -1) */

	if( n < 0 )
	{	/* test for negative exponent */
		n = -n;
		m = 1/m;
	}

	p = is_odd(n) ? m : 1;	/* test & set zero power */

	/*lint -e704 shift right of signed quantity */
	/*lint -e720 Boolean test of assignment */
	while( n >>= 1 )
	{	/* now do the other powers */
		m *= m;			/* sq previous power of m */
		if( is_odd(n) )	/* if low order bit set */
			p *= m;		/*	then, multiply partial product by latest power of m */
	}
	/*lint +e704 shift right of signed quantity */
	/*lint +e720 Boolean test of assignment */
	return p;
}

/*PrintE82 - series of routines to mimic 1p, e8.2 fortran output */
/***********************************************************
 * contains the following sets of routines to get around   *
 * the MS C++ compilers unusual exponential output.        *
 * PrintEfmt <= much faster, no overhead with unix         *
 * PrintE93                                                *
 * PrintE82                                                *
 * PrintE71                                                *
 **********************************************************/

#ifdef _MSC_VER
/**********************************************************/
/*
 * Instead of printf'ing with %e or %.5e or whatever, call
 * efmt("%e", val) and print the result with %s.  This lets
 * us work around bugs in VS C 6.0.
 */
char *PrintEfmt(const char *fmt, double val /* , char *buf */) 
{
	static char buf[30]; /* or pass it in */

	DEBUG_ENTRY( "PrintEfmt()" );

	/* create the string */
	sprintf(buf, fmt, val);

	/* code to fix incorrect ms v e format.  works only for case where there is
	 * a leading space in the format - for formats like 7.1, 8.2, 9.3, 10.4, etc
	 * result will have 1 too many characters */
	char *ep , buf2[30];

	/* msvc behaves badly in different ways for positive vs negative sign vals,
	 * if val is positive must create a leading space */
	if( val >= 0.)
	{
		strcpy(buf2 , " " );
		strcat(buf2 , buf);
		strcpy( buf , buf2);
	}

	/* allow for both e and E formats */
	if((ep = strchr_s(buf, 'e')) == NULL)
	{
		ep = strchr_s(buf, 'E');
	}

	/* ep can still be NULL if val is Inf or NaN */
	if(ep != NULL) 
	{
		/* move pointer two char past the e, to pick up the e and sign */
		ep += 2;

		/* terminate buf where the e is, *ep points to this location */
		*ep = '\0';

		/* skip next char, */
		++ep;

		/* copy resulting string to return string */
		strcat( buf, ep );
	}
	return buf;
}
#endif

/**********************************************************/
void PrintE82( FILE* ioOUT, double value )
{
	double frac , xlog , xfloor , tvalue;
	int iExp;

	DEBUG_ENTRY( "PrintE82()" );

	if( value < 0 )
	{
		fprintf(ioOUT,"********");
	}
	else if( value == 0 )
	{
		fprintf(ioOUT,"0.00E+00");
	}
	else
	{
		/* round number off for 8.2 format (not needed since can't be negative) */
		tvalue = value;
		xlog = log10( tvalue );
		xfloor = floor(xlog);
		/* this is now the fractional part */
		frac = tvalue*pow(10.,-xfloor);
		/*this is the possibly signed exponential part */
		iExp = (int)xfloor;
		if( frac>9.9945 )
		{
			frac /= 10.;
			iExp += 1;
		}
		/* print the fractional part*/
		fprintf(ioOUT,"%.2f",frac);
		/* E for exponent */
		fprintf(ioOUT,"E");
		/* if positive throw a + sign*/
		if(iExp>=0 )
		{
			fprintf(ioOUT,"+");
		}
		fprintf(ioOUT,"%.2d",iExp);
	}
	return;
}
/*
 *==============================================================================
 */
void PrintE71( FILE* ioOUT, double value )
{
	double frac , xlog , xfloor , tvalue;
	int iExp;

	DEBUG_ENTRY( "PrintE71()" );

	if( value < 0 )
	{
		fprintf(ioOUT,"*******");
	}
	else if( value == 0 )
	{
		fprintf(ioOUT,"0.0E+00");
	}
	else
	{
		/* round number off for 8.2 format (not needed since can't be negative) */
		tvalue = value;
		xlog = log10( tvalue );
		xfloor = floor(xlog);
		/* this is now the fractional part */
		frac = tvalue*pow(10.,-xfloor);
		/*this is the possibly signed exponential part */
		iExp = (int)xfloor;
		if( frac>9.9945 )
		{
			frac /= 10.;
			iExp += 1;
		}
		/* print the fractional part*/
		fprintf(ioOUT,"%.1f",frac);
		/* E for exponent */
		fprintf(ioOUT,"E");
		/* if positive throw a + sign*/
		if(iExp>=0 )
		{
			fprintf(ioOUT,"+");
		}
		fprintf(ioOUT,"%.2d",iExp);
	}
	return;
}

/*
 *==============================================================================
 */
void PrintE93( FILE* ioOUT, double value )
{
	double frac , xlog , xfloor, tvalue;
	int iExp;

	DEBUG_ENTRY( "PrintE93()" );

	if( value < 0 )
	{
		fprintf(ioOUT,"*********");
	}
	else if( value == 0 )
	{
		fprintf(ioOUT,"0.000E+00");
	}
	else
	{
		/* round number off for 9.3 format, neg numb not possible */
		tvalue = value;
		xlog = log10( tvalue );
		xfloor = floor(xlog);
		/* this is now the fractional part */
		frac = tvalue*pow(10.,-xfloor);
		/*this is the possibly signed exponential part */
		iExp = (int)xfloor;
		if( frac>9.99949 )
		{
			frac /= 10.;
			iExp += 1;
		}
		/* print the fractional part*/
		fprintf(ioOUT,"%5.3f",frac);
		/* E for exponent */
		fprintf(ioOUT,"E");
		/* if positive throw a + sign*/
		if(iExp>=0 )
		{
			fprintf(ioOUT,"+");
		}
		fprintf(ioOUT,"%.2d",iExp);
	}
	return;
}

/*TotalInsanity general error handler for something that cannot happen */
NORETURN void TotalInsanity(void)
{
	DEBUG_ENTRY( "TotalInsanity()" );

	/* something that cannot happen, happened,
	 * if this message is triggered, simply place a breakpoint here
	 * and debug the error */
	fprintf( ioQQQ, " Something that cannot happen, has happened.\n" );
	fprintf( ioQQQ, " This is TotalInsanity, I live in %s.\n", __FILE__ );
	ShowMe();

	cdEXIT(EXIT_FAILURE);
}

/*BadRead general error handler for trying to read data, but failing */
NORETURN void BadRead(void)
{
	DEBUG_ENTRY( "BadRead()" );

	/* read failed */
	fprintf( ioQQQ, " A read of internal input data has failed.\n" );
	fprintf( ioQQQ, " This is BadRead, I live in %s.\n", __FILE__ );
	ShowMe();

	cdEXIT(EXIT_FAILURE);
}

/*sexp safe exponential function */
sys_float sexp(sys_float x)
{
	sys_float sexp_v;

	DEBUG_ENTRY( "sexp()" );

	/* SEXP_LIMIT is 84 in cddefines.h */
	if( x < SEXP_LIMIT )
	{
		sexp_v = exp(-x);
	}
	else
	{
		sexp_v = 0.f;
	}
	return sexp_v;
}

/*sexp safe exponential function */
double sexp(double x)
{
	double sexp_v;

	DEBUG_ENTRY( "sexp()" );

	/* SEXP_LIMIT is 84 in cddefines.h */
	if( x < SEXP_LIMIT )
	{
		sexp_v = exp(-x);
	}
	else
	{
		sexp_v = 0.;
	}
	return sexp_v;
}


/*dsexp safe exponential function for doubles */
double dsexp(double x)
{
	double dsexp_v;

	DEBUG_ENTRY( "dsexp()" );

	if( x < DSEXP_LIMIT )
	{
		dsexp_v = exp(-x);
	}
	else
	{
		dsexp_v = 0.;
	}
	return dsexp_v;
}

/*TestCode set flag saying that test code is in place 
 * prototype in cddefines.h */
void TestCode(void)
{
	DEBUG_ENTRY( "TestCode( )" );

	/* called if test code is in place */
	lgTestCodeCalled = true;
	return;
}

/*broken set flag saying that the code is broken, */
void broken(void)
{
	DEBUG_ENTRY( "broken( )" );

	broke.lgBroke = true;
	return;
}

/*fixit say that code needs to be fixed */
void fixit(void)
{
	DEBUG_ENTRY( "fixit( )" );

	broke.lgFixit = true;
	return;
}

/*CodeReview placed next to code that needs to be checked */
void CodeReview(void)
{
	DEBUG_ENTRY( "CodeReview( )" );

	broke.lgCheckit = true;
	return;
}

/** dprintf -- version of fprintf which prepends DEBUG */
int dprintf(FILE *fp, const char *format, ...)
{
	va_list ap;
	int i1, i2;

	DEBUG_ENTRY( "dprintf()" );
	va_start(ap,format);
	i1 = fprintf(fp,"DEBUG ");
	if (i1 >= 0)
		i2 = vfprintf(fp,format,ap);
	else
		i2 = 0;
	if (i2 < 0)
		i1 = 0;
	va_end(ap);

	return i1+i2;
}

/* dbg_printf is a debug print routine that was provided by Peter Teuben,
 * as a component from his NEMO package.  It offers run-time specification
 * of the level of debugging */
int dbg_printf(int debug, const char *fmt, ...)
{
	va_list ap;
	int i=0;

	DEBUG_ENTRY( "dbg_printf()" );

	/* print this debug message? (debug_level not currently used)*/
	if(debug <= trace.debug_level) 
	{		
		va_start(ap, fmt);	

		i = vfprintf(ioQQQ, fmt, ap);
		/* drain ioQQQ */
		fflush(ioQQQ);
		va_end(ap);
	}
	return i;
}


/*qg32 32 point Gaussian quadrature, originally given to Gary F by Jim Lattimer */
double qg32(
	double xl, /*lower limit to integration range*/
	double xu, /*upper limit to integration range*/
	/*following is the pointer to the function that will be evaluated*/
	double (*fct)(double) )
{
	double a = 0.5*(xu+xl), 
	  b = xu-xl, 
	  y = 0.;

	DEBUG_ENTRY( "qg32()" );

	/********************************************************************************
	 *                                                                              *
	 *  32-point Gaussian quadrature                                                *
	 *  xl  : the lower limit of integration                                        *
	 *  xu  : the upper limit                                                       *
	 *  fct : the (external) function                                               *
	 *  returns the value of the integral                                           *
	 *                                                                              *
	 * simple call to integrate sine from 0 to pi                                   *
	 * double agn = qg32( 0., 3.141592654 ,  sin );                                 *
	 *                                                                              *
	 *******************************************************************************/

	double weights[16] = {
		.35093050047350483e-2, .81371973654528350e-2, .12696032654631030e-1, .17136931456510717e-1,
		.21417949011113340e-1, .25499029631188088e-1, .29342046739267774e-1, .32911111388180923e-1,
		.36172897054424253e-1, .39096947893535153e-1, .41655962113473378e-1, .43826046502201906e-1,
		.45586939347881942e-1, .46922199540402283e-1, .47819360039637430e-1, .48270044257363900e-1};

	double c[16] = {
		.498631930924740780, .49280575577263417, .4823811277937532200, .46745303796886984000,
		.448160577883026060, .42468380686628499, .3972418979839712000, .36609105937014484000,
		.331522133465107600, .29385787862038116, .2534499544661147000, .21067563806531767000,
		.165934301141063820, .11964368112606854, .7223598079139825e-1, .24153832843869158e-1};

	for( int i=0; i<16; i++)
	{
		y += b * weights[i] * ((*fct)(a+b*c[i]) + (*fct)(a-b*c[i]));
	}

	/* the answer */
	return y;
}

/*spsort netlib routine to sort array returning sorted indices */
void spsort(
	  /* input array to be sorted */
	  realnum x[], 
	  /* number of values in x */
	  long int n, 
	  /* permutation output array */
	  long int iperm[], 
	  /* flag saying what to do - 1 sorts into increasing order, not changing
	   * the original vector, -1 sorts into decreasing order. 2, -2 change vector */
	  int kflag, 
	  /* error condition, should be 0 */
	  int *ier)
{
	/*
	 ****BEGIN PROLOGUE  SPSORT
	 ****PURPOSE  Return the permutation vector generated by sorting a given
	 *            array and, optionally, rearrange the elements of the array.
	 *            The array may be sorted in increasing or decreasing order.
	 *            A slightly modified quicksort algorithm is used.
	 ****LIBRARY   SLATEC
	 ****CATEGORY  N6A1B, N6A2B
	 ****TYPE      SINGLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
	 ****KEY WORDS NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
	 ****AUTHOR  Jones, R. E., (SNLA)
	 *           Rhoads, G. S., (NBS)
	 *           Wisniewski, J. A., (SNLA)
	 ****DESCRIPTION
	 *
	 *   SPSORT returns the permutation vector IPERM generated by sorting
	 *   the array X and, optionally, rearranges the values in X.  X may
	 *   be sorted in increasing or decreasing order.  A slightly modified
	 *   quicksort algorithm is used.
	 *
	 *   IPERM is such that X(IPERM(I)) is the Ith value in the rearrangement
	 *   of X.  IPERM may be applied to another array by calling IPPERM,
	 *   SPPERM, DPPERM or HPPERM.
	 *
	 *   The main difference between SPSORT and its active sorting equivalent
	 *   SSORT is that the data are referenced indirectly rather than
	 *   directly.  Therefore, SPSORT should require approximately twice as
	 *   long to execute as SSORT.  However, SPSORT is more general.
	 *
	 *   Description of Parameters
	 *      X - input/output -- real array of values to be sorted.
	 *          If ABS(KFLAG) = 2, then the values in X will be
	 *          rearranged on output; otherwise, they are unchanged.
	 *      N - input -- number of values in array X to be sorted.
	 *      IPERM - output -- permutation array such that IPERM(I) is the
	 *              index of the value in the original order of the
	 *              X array that is in the Ith location in the sorted
	 *              order.
	 *      KFLAG - input -- control parameter:
	 *            =  2  means return the permutation vector resulting from
	 *                  sorting X in increasing order and sort X also.
	 *            =  1  means return the permutation vector resulting from
	 *                  sorting X in increasing order and do not sort X.
	 *            = -1  means return the permutation vector resulting from
	 *                  sorting X in decreasing order and do not sort X.
	 *            = -2  means return the permutation vector resulting from
	 *                  sorting X in decreasing order and sort X also.
	 *      IER - output -- error indicator:
	 *          =  0  if no error,
	 *          =  1  if N is zero or negative,
	 *          =  2  if KFLAG is not 2, 1, -1, or -2.
	 ****REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
	 *                 for sorting with minimal storage, Communications of
	 *                 the ACM, 12, 3 (1969), pp. 185-187.
	 ****ROUTINES CALLED  XERMSG
	 ****REVISION HISTORY  (YYMMDD)
	 *   761101  DATE WRITTEN
	 *   761118  Modified by John A. Wisniewski to use the Singleton
	 *           quicksort algorithm.
	 *   870423  Modified by Gregory S. Rhoads for passive sorting with the
	 *           option for the rearrangement of the original data.
	 *   890620  Algorithm for rearranging the data vector corrected by R.
	 *           Boisvert.
	 *   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
	 *   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
	 *   920507  Modified by M. McClain to revise prologue text.
	 *   920818  Declarations section rebuilt and code restructured to use
	 *           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
	 ****END PROLOGUE  SPSORT
	 *     .. Scalar Arguments ..
	 */
	long int i, 
	  ij, 
	  il[21], 
	  indx, 
	  indx0, 
	  istrt, 
	  istrt_, 
	  iu[21], 
	  j, 
	  k, 
	  kk, 
	  l, 
	  lm, 
	  lmt, 
	  m, 
	  nn;
	realnum r, 
	  ttemp;

	DEBUG_ENTRY( "spsort()" );

	/*     .. Array Arguments .. */
	/*     .. Local Scalars .. */
	/*     .. Local Arrays .. */
	/*     .. External Subroutines .. */
	/*     .. Intrinsic Functions .. */
	/****FIRST EXECUTABLE STATEMENT  SPSORT */
	*ier = 0;
	nn = n;
	if( nn < 1 )
	{
		*ier = 1;
		return;
	}
	else
	{
		kk = labs(kflag);
		if( kk != 1 && kk != 2 )
		{
			*ier = 2;
			return;
		}
		else
		{

			/* Initialize permutation vector to index on f scale
			 * */
			for( i=0; i < nn; i++ )
			{
				iperm[i] = i+1;
			}

			/* Return if only one value is to be sorted */
			if( nn == 1 )
			{ 
				--iperm[0];
				return;
			}

			/* Alter array X to get decreasing order if needed */
			if( kflag <= -1 )
			{
				for( i=0; i < nn; i++ )
				{
					x[i] = -x[i];
				}
			}

			/* Sort X only */
			m = 1;
			i = 1;
			j = nn;
			r = .375e0;
		}
	}

	while( true )
	{
		if( i == j )
			goto L_80;
		if( r <= 0.5898437e0 )
		{
			r += 3.90625e-2;
		}
		else
		{
			r -= 0.21875e0;
		}

L_40:
		k = i;

		/*     Select a central element of the array and save it in location L
		 * */
		ij = i + (long)((j-i)*r);
		lm = iperm[ij-1];

		/*     If first element of array is greater than LM, interchange with LM
		 * */
		if( x[iperm[i-1]-1] > x[lm-1] )
		{
			iperm[ij-1] = iperm[i-1];
			iperm[i-1] = lm;
			lm = iperm[ij-1];
		}
		l = j;

		/*     If last element of array is less than LM, interchange with LM
		 * */
		if( x[iperm[j-1]-1] < x[lm-1] )
		{
			iperm[ij-1] = iperm[j-1];
			iperm[j-1] = lm;
			lm = iperm[ij-1];

			/*        If first element of array is greater than LM, interchange
			 *        with LM
			 * */
			if( x[iperm[i-1]-1] > x[lm-1] )
			{
				iperm[ij-1] = iperm[i-1];
				iperm[i-1] = lm;
				lm = iperm[ij-1];
			}
		}

		/* Find an element in the second half of the array which is smaller
		 * than LM */
		while( true )
		{
			l -= 1;
			if( x[iperm[l-1]-1] <= x[lm-1] )
			{

				/* Find an element in the first half of the array which is greater
				 * than LM */
				while( true )
				{
					k += 1;
					if( x[iperm[k-1]-1] >= x[lm-1] )
						break;
				}

				/* Interchange these elements */
				if( k > l )
					break;
				lmt = iperm[l-1];
				iperm[l-1] = iperm[k-1];
				iperm[k-1] = lmt;
			}
		}

		/* Save upper and lower subscripts of the array yet to be sorted */
		if( l - i > j - k )
		{
			il[m-1] = i;
			iu[m-1] = l;
			i = k;
			m += 1;
		}
		else
		{
			il[m-1] = k;
			iu[m-1] = j;
			j = l;
			m += 1;
		}

L_90:
		if( j - i >= 1 )
			goto L_40;
		if( i == 1 )
			continue;
		i -= 1;

		while( true )
		{
			i += 1;
			if( i == j )
				break;
			lm = iperm[i];
			if( x[iperm[i-1]-1] > x[lm-1] )
			{
				k = i;

				while( true )
				{
					iperm[k] = iperm[k-1];
					k -= 1;

					if( x[lm-1] >= x[iperm[k-1]-1] )
						break;
				}
				iperm[k] = lm;
			}
		}

		/* Begin again on another portion of the unsorted array */
L_80:
		m -= 1;
		if( m == 0 )
			break;
		/*lint -e644 not explicitly initialized */
		i = il[m-1];
		j = iu[m-1];
		/*lint +e644 not explicitly initialized */
		goto L_90;
	}

	/* Clean up */
	if( kflag <= -1 )
	{
		for( i=0; i < nn; i++ )
		{
			x[i] = -x[i];
		}
	}

	/* Rearrange the values of X if desired */
	if( kk == 2 )
	{

		/* Use the IPERM vector as a flag.
		 * If IPERM(I) < 0, then the I-th value is in correct location */
		for( istrt=1; istrt <= nn; istrt++ )
		{
			istrt_ = istrt - 1;
			if( iperm[istrt_] >= 0 )
			{
				indx = istrt;
				indx0 = indx;
				ttemp = x[istrt_];
				while( iperm[indx-1] > 0 )
				{
					x[indx-1] = x[iperm[indx-1]-1];
					indx0 = indx;
					iperm[indx-1] = -iperm[indx-1];
					indx = labs(iperm[indx-1]);
				}
				x[indx0-1] = ttemp;
			}
		}

		/* Revert the signs of the IPERM values */
		for( i=0; i < nn; i++ )
		{
			iperm[i] = -iperm[i];
		}
	}

	for( i=0; i < nn; i++ )
	{
		--iperm[i];
	}
	return;
}

/*MyMalloc wrapper for malloc().  Returns a good pointer or dies. 
 * memory is filled with NaN
 * >>chng 05 dec 14, do not set to NaN since tricks debugger 
 * routines within code do not call this or malloc, but rather MALLOC
 * which is resolved into MyMalloc or malloc depending on whether 
 * NDEBUG is set by the compiler to indicate "not debugging",
 * in typical negative C style */
void *MyMalloc( 
	/*use same type as library function MALLOC*/ 
	size_t size ,
	const char *chFile, 
	int line
	)
{
	void *ptr;

	DEBUG_ENTRY( "MyMalloc()" );

	ASSERT( size > 0 );

	/* debug branch for printing malloc args */
	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			static long int kount=0, nTot=0;
			nTot += (long)size;
			fprintf(ioQQQ,"%li\t%li\t%li\n", 
				kount , 
				(long)size , 
				nTot );
			++kount;
		}
	}

	if( ( ptr = malloc( size ) ) == NULL )
	{
		fprintf(ioQQQ,"DISASTER MyMalloc could not allocate %lu bytes.  Exit in MyMalloc.",
			(unsigned long)size );
		fprintf(ioQQQ,"MyMalloc called from file %s at line %i.\n",
			chFile , line );

		if( struc.nzlim>2000 )
			fprintf(ioQQQ,"This may have been caused by the large number of zones."
			" %li zones were requested.  Is this many zones really necessary?\n",
			struc.nzlim );

		cdEXIT(EXIT_FAILURE);
	}

	/* flag -DNOINIT will turn off this initialization which can fool valgrind/purify */
#	if !defined(NDEBUG) && !defined(NOINIT)

	size_t nFloat = size/4;
	size_t nDouble = size/8;
	sys_float *fptr = static_cast<sys_float*>(ptr);
	double *dptr = static_cast<double*>(ptr);

	/* >>chng 04 feb 03, fill memory with invalid numbers, PvH */
	/* on IA32/AMD64 processors this will generate NaN's for both float and double;
	 * on most other (modern) architectures it is likely to do the same... */
	/* >>chng 05 dec 14, change code to generate signaling NaN's for most cases (but not all!) */
	if( size == nDouble*8 )
	{
		/* this could be an array of doubles as well as floats -> we will hedge our bets
		 * we will fill the array with a pattern that is interpreted as all signaling
		 * NaN's for doubles, and alternating signaling and quiet NaN's for floats:
		 * byte offset:  0      4      8     12     16
		 * double        | SNaN        | SNan        |
		 * float         | SNaN | QNaN | SNan | QNaN |  (little-endian, e.g. Intel, AMD, alpha)
		 * float         | QNaN | SNaN | QNan | SNaN |  (big-endian, e.g. Sparc, PowerPC, MIPS) */
		set_NaN( dptr, (long)nDouble );
	}
	else if( size == nFloat*4 )
	{
		/* this could be an arrays of floats, but not doubles -> init to all float SNaN */
		set_NaN( fptr, (long)nFloat );
	}
	else
	{
		memset( ptr, 0xff, size );
	}

#	endif /* !defined(NDEBUG) && !defined(NOINIT) */
	return ptr;
}


/* wrapper for calloc().  Returns a good pointer or dies. 
 * routines within code do not call this or malloc, but rather CALLOC
 * which is resolved into MyCalloc or calloc depending on whether 
 * NDEBUG is set in cddefines. \h */
void *MyCalloc( 
	/*use same type as library function CALLOC*/ 
	size_t num ,
	size_t size )
{
	void *ptr;

	DEBUG_ENTRY( "MyCalloc()" );

	ASSERT( size > 0 );

	/* debug branch for printing malloc args */
	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			static long int kount=0;
			fprintf(ioQQQ,"%li\tcall\t%li\tbytes\n", kount, 
				(long)size );
			++kount;
		}
	}

	if( ( ptr = calloc( num , size ) ) == NULL )
	{
		fprintf(ioQQQ,"MyCalloc could not allocate %lu bytes.  Exit in MyCalloc.",
			(unsigned long)size );
		cdEXIT(EXIT_FAILURE);
	}
	return ptr;
}

/* wrapper for realloc().  Returns a good pointer or dies. 
 * routines within code do not call this or malloc, but rather REALLOC
 * which is resolved into MyRealloc or realloc depending on whether 
 * NDEBUG is set in cddefines.h */
void *MyRealloc( 
	/*use same type as library function realloc */ 
	void *p ,
	size_t size )
{
	void *ptr;

	DEBUG_ENTRY( "MyRealloc()" );

	ASSERT( size > 0 );

	/* debug branch for printing malloc args */
	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
			static long int kount=0;
			fprintf(ioQQQ,"%li\tcall\t%li\tbytes\n", kount, 
				(long)size );
			++kount;
		}
	}

	if( ( ptr = realloc( p , size ) ) == NULL )
	{
		fprintf(ioQQQ,"MyRealloc could not allocate %lu bytes.  Exit in MyRealloc.",
			(unsigned long)size );
		cdEXIT(EXIT_FAILURE);
	}
	return ptr;
}

/* function to facilitate addressing opacity array */
double csphot(
	/* INU is array index pointing to frequency where opacity is to be evaluated
	 * on f not c scale */
	long int inu, 
	/* ITHR is pointer to threshold*/
	long int ithr, 
	/* IOFSET is offset as defined in opac0*/
	long int iofset)
{
	double csphot_v;

	DEBUG_ENTRY( "csphot()" );

	csphot_v = opac.OpacStack[inu-ithr+iofset-1];
	return csphot_v;
}

/*RandGauss normal Gaussian random number generator
 * the user must call srand to set the seed before using this routine.

 * the random numbers will have a mean of xMean and a standard deviation of s

 * The convention is for srand to be called when the command setting 
 * the noise is parsed 

 * for very small dispersion there are no issues, but when the dispersion becomes
 * large the routine will find negative values - so most often used in this case
 * to find dispersion in log space 
 * this routine will return a normal Gaussian - must be careful in how this is
 * used when adding noise to physical quantity */
/*
NB - following from Ryan Porter:
I discovered that I unintentionally created an antisymmetric skew in my 
Monte Carlo.  RandGauss is symmetric in log space, which means it is not 
symmetric in linear space.  But to get the right standard deviation you 
have to take 10^x, where x is the return from RandGauss.  The problem is 
10^x will happen less frequently than 10^-x, so without realizing it, the 
average "tweak" to every piece of atomic data in my Monte Carlo run was 
not 1.0 but something greater than 1.0, causing every single line to have 
an average Monte Carlo emissivity greater than the regular value.  Any place
that takes 10^RandGauss() needs to be adjusted if what is intended is +/- x. */
double RandGauss(
	/* mean value */
	double xMean, 
	/*standard deviation s */
	double s )
{
	double  x1, x2, w, yy1;
	static double yy2=-BIGDOUBLE;
	static int use_last = false;

	DEBUG_ENTRY( "RandGauss()" );

	if( use_last )
	{
		yy1 = yy2;
		use_last = false;
	}
	else
	{
		do {
			x1 = 2.*genrand_real3() - 1.;
			x2 = 2.*genrand_real3() - 1.;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt((-2.0*log(w))/w);
		yy1 = x1 * w;
		yy2 = x2 * w;
		use_last = true;
	}
	return xMean + yy1 * s;
}

/* MyGaussRand takes as input a percent uncertainty less than 50%
 * (expressed as 0.5). The routine then assumes this input variable represents one
 * standard deviation about a mean of unity, and returns a random number within
 * that range.  A hard cutoff is imposed at two standard deviations, which 
 * eliminates roughly 5% of the normal distribution.  In other words, the routine
 * returns a number in a normal distribution with standard deviation equal to
 * the input.  The number will be between 1-3*stdev and 1+3*stdev. */
double MyGaussRand( double PctUncertainty )
{
	double StdDev;
	double result;

	DEBUG_ENTRY( "MyGaussRand()" );

	ASSERT( PctUncertainty < 0.5 );
	/* We want this "percent uncertainty" to represent one standard deviation */
	StdDev = PctUncertainty;

	do
	{
		/*result = pow( 10., RandGauss( 0., logStdDev ) );*/
		result = 1.+RandGauss( 0., StdDev );
	}
	/* only allow values that are within 3 standard deviations */
	while( (result < 1.-3.*PctUncertainty) || (result > 1.+3.*PctUncertainty) );

	ASSERT( result>0. && result<2. );
	return result;
}

/*plankf evaluate Planck function for any cell at current electron temperature */
double plankf(long int ip)
{
	double plankf_v;

	DEBUG_ENTRY( "plankf()" );

	/* evaluate Planck function
	 * argument is pointer to cell energy in ANU
	 * return photon flux for cell IP */
	if( rfield.ContBoltz[ip] <= 0. )
	{
		plankf_v = 1e-35;
	}
	else
	{
		plankf_v = 6.991e-21*POW2(FR1RYD*rfield.anu[ip])/
			(1./rfield.ContBoltz[ip] - 1.)*FR1RYD*4.;
	}
	return plankf_v;
}
