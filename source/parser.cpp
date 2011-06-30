/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "parser.h"
#include "called.h"
#include "energy.h"
#include "flux.h"
#include "input.h"
#include "elementnames.h"

#include <deque>
class Token
{
public:
	enum symType { symNull, symNumber, symOp };
	string s;
	symType t;
	explicit Token(enum symType type) : s(""), t(type) {}
	explicit Token() : s(""), t(symNull) {}
};
STATIC bool ParseExpr(deque<Token> &chTokens, vector<double> &valstack);

const char *Parser::nWord(const char *chKey) const
{
	return ::nWord(chKey, m_card);
}

/*nWord determine whether match to a keyword occurs on command line,
 * return value is 0 if no match, and position of match within string if hit */
const char *nWord(const char *chKey, 
	    const char *chCard)
{
	DEBUG_ENTRY( "nWord()" );

	// Ignore leading space in chKey -- logic below is designed
	// to avoid the need to include this in the first place
	while (isspace(*chKey))
	{
		++chKey;
	}

	const long lenkey = strlen(chKey);
	ASSERT( lenkey > 0 );

	bool atBoundary = true;
	for (const char *ptr = chCard; *ptr; ++ptr)
	{
		if ( atBoundary && strncmp( ptr, chKey, lenkey) == 0 )
		{
			return ptr;
		}

		atBoundary = isBoundaryChar(*ptr);
	}

	return NULL;
}

bool isBoundaryChar(char c)
{
	const bool lgAnyWhitespacePrecedesWord = false;
	
	if (lgAnyWhitespacePrecedesWord)
		return isspace(c) ? true : false ;
	else 	// Words are strings starting with A-Z, a-z or _
		return (! isalpha(c) ) && c != '_';
}

bool Parser::isComment(void) const
{
	return lgInputComment(m_card);
}
bool Parser::isCommandComment(void) const
{
	return ( m_card[0]=='C' && (m_card[1]==' ' || m_card[1]== '\0')) ||
		isComment();
}

void Parser::echo(void) const
{
	/* >>chng 04 jan 21, add HIDE option, mostly for print quiet command */
	if( called.lgTalk && !::nMatch("HIDE",m_card) )
		fprintf( ioQQQ, "%23c* %-80s*\n", ' ', m_card_raw );
}

NORETURN void Parser::CommandError(void) const
{
	DEBUG_ENTRY(" Parser::CommandError()");
	fprintf( ioQQQ, "  Unrecognized command. Key=\"%4.4s\".  This is routine ParseCommands.\n", 
				m_card );
	fprintf( ioQQQ, " The line image was\n");
	PrintLine(ioQQQ);
	fprintf( ioQQQ, " Sorry.\n" );
	cdEXIT(EXIT_FAILURE);
}
bool Parser::getline(void)
{
	input.readarray(m_card_raw,&m_lgEOF);
	newlineProcess();
	if (m_lgEOF)
		return false;
	else
		return true;
}
	
const char *Parser::StandardEnergyUnit(void) const
{
	return ::StandardEnergyUnit(m_card);
}
string Parser::StandardFluxUnit(void) const
{
	return ::StandardFluxUnit(m_card);
}
void Parser::help(FILE *fp) const
{
	DEBUG_ENTRY("Parser::help()");
	fprintf(fp,"Available commands are:\n\n");
	long int i=0, l=0, len;
	while (1)
	{
		len = strlen(m_Commands[i].name);
		if (l+len+2 > 80)
		{
			fprintf(fp,"\n");
			l = 0;
		}
		l += len+2;
		fprintf(fp,"%s",m_Commands[i].name);
		++i;
		if (m_Commands[i].name == NULL)
			break;
		fprintf(fp,", ");
	}
				  
	fprintf(fp,"\n\nSorry, no further help available yet -- try Hazy.\n\n");
	cdEXIT(EXIT_SUCCESS);
}

/*GetElem scans line image, finds element. returns atomic number j, 
 * on C scale, -1 if no hit.  chCARD_CAPS must be in CAPS to hit element */
long int Parser::GetElem(void ) const
{
	int i;

	DEBUG_ENTRY( "GetElem()" );

	/* find which element */

	/* >>>chng 99 apr 17, lower limit to loop had been 1, so search started with helium,
	 * change to 0 so we can pick up hydrogen.  needed for parseasserts command */
	/* find match with element name, start with helium */
	for( i=0; i<(int)LIMELM; ++i )
	{
		if( nMatch( elementnames.chElementNameShort[i] ) )
		{
			/* return value is in C counting, hydrogen would be 0*/
			return i;
		}
	}
	/* fall through, did not hit, return -1 as error condition */
	return (-1 );
}

/*NoNumb general error handler for no numbers on input line */
NORETURN void Parser::NoNumb(const char * chDesc) const
{
	DEBUG_ENTRY( "NoNumb()" );

	/* general catch-all for no number when there should have been */
	fprintf( ioQQQ, " There is a problem on the following command line:\n" );
	fprintf( ioQQQ, " %s\n", m_card_raw );
	fprintf( ioQQQ, " A value for %s should have been on this line.\n   Sorry.\n",chDesc );
	cdEXIT(EXIT_FAILURE);
 }

double Parser::getWaveOpt()
{
	double val = FFmtRead();
	/* check for optional micron or cm units, else interpret as Angstroms */
	if( chPoint() == 'M' )
	{
		/* microns */
		val *= 1e4;
	}
	else if( chPoint() == 'C' )
	{
		/* centimeters */
		val *= 1e8;
	}
	return val;
}
double Parser::getWave()
{
	double val = getWaveOpt();
	if( lgEOL() )
	{
		NoNumb("wavelength");
	}
	return val;
}
double Parser::getNumberPlain( const char * )
{
	return FFmtRead();
}
double Parser::getNumberCheck( const char *chDesc )
{
	double val = FFmtRead();
	if( lgEOL() )
	{
		NoNumb(chDesc);
	}
	return val;
}
double Parser::getNumberDefault( const char *, double fdef )
{
	double val = FFmtRead();
	if( lgEOL() )
	{
		val = fdef;
	}
	return val;
}
double Parser::getNumberCheckLogLinNegImplLog( const char *chDesc )
{
	double val = getNumberCheck(chDesc);
	if( nMatch(" LOG") )
	{
		val = pow(10.,val);
	}		
	else if(! nMatch("LINE") )
	{
		/* log, linear not specified, neg so log */
		if( val <= 0. )
		{
			val = pow(10.,val);
		}
	}
	return val;
}
double Parser::getNumberCheckAlwaysLog( const char *chDesc )
{
	double val = getNumberCheck(chDesc);
	val = pow(10., val);
	return val;
}
double Parser::getNumberCheckAlwaysLogLim( const char *chDesc, double flim )
{
	double val = getNumberCheck(chDesc);
	if ( val > flim )
	{
		fprintf(ioQQQ,"WARNING - the log of %s is too "
				  "large, I shall probably crash.  The value was %.2e\n",
				  chDesc, val );
		fflush(ioQQQ);		
	}
	val = pow(10., val);
	return val;
}
double Parser::getNumberDefaultAlwaysLog( const char *, double fdef )
{
	double val = pow(10.,FFmtRead());
	if ( lgEOL() )
	{
		val = fdef;
	}
	return val;
}
double Parser::getNumberDefaultNegImplLog( const char *, double fdef )
{
	double val = FFmtRead();
	if ( lgEOL() )
	{
		val = fdef;
	}
	if (val < 0.0)
	{
		val = pow(10.,val);
	}
	return val;
}

/*FFmtRead scan input line for free format number */


double Parser::FFmtRead(void)
{

	DEBUG_ENTRY( "Parser::FFmtRead()" );

	char chr = '\0';
	// eol_ptr points one beyond last valid char
	const char * const eol_ptr = m_card+m_len; 

	// Look for start of next expression
	while( m_ptr < eol_ptr && ( chr = *m_ptr++ ) != '\0' )
	{
		const char *lptr = m_ptr;
		char lchr = chr;
		if( lchr == '-' || lchr == '+' )
			lchr = *lptr++;
		if( lchr == '.' )
			lchr = *lptr;
		if( isdigit(lchr) )
			break;
	}

	if( m_ptr == eol_ptr || chr == '\0' )
	{
		m_lgEOL = true;
		return 0.;
	}

	// Lexer for expression
	deque<Token> chTokens(0);
	bool lgCommaFound = false, lgLastComma = false;
	do
	{
		lgCommaFound = lgLastComma;
		if( chr != ',' )
		{
			if (chr != '^' && chr != '*') 
			{
				if (chTokens.size() == 0 || chTokens.back().t != Token::symNumber)
					chTokens.push_back(Token(Token::symNumber));
				chTokens.back().s += chr;
			}
			else
			{
				chTokens.push_back(Token(Token::symOp));
				chTokens.back().s += chr;
			}
		}
		else
		{
			/* don't complain about comma if it appears after number,
				as determined by exiting loop before this sets lgCommaFound */
			lgLastComma = true;

		}
		if( m_ptr == eol_ptr )
			break;
		chr = *m_ptr++;
	}
	while( isdigit(chr) || chr == '.' || chr == '-' || chr == '+' || chr == ',' || chr == 'e' || chr == 'E' || chr == '^' || chr == '*' );

	if( lgCommaFound )
	{
		fprintf( ioQQQ, " PROBLEM - a comma was found embedded in a number, this is deprecated.\n" );
		fprintf(ioQQQ, "== %-80s ==\n",m_card);
	}

	// Parse tokens
	vector<double> valstack;
	const bool lgParseOK = ParseExpr(chTokens, valstack);
	if (!lgParseOK || 1 != valstack.size())
	{
		fprintf(ioQQQ," PROBLEM - syntax error in number\n");
		fprintf(ioQQQ, "== %-80s ==\n",m_card);
	}

	double value = valstack[0];

	m_lgEOL = false;
	m_ptr--;	 // m_ptr already points 1 beyond where next read should start

	return value;
}

STATIC bool ParseNumber(deque<Token> &chTokens, vector<double> &valstack)
{
	if ( chTokens.size() < 1 || 
		  Token::symNumber != chTokens[0].t)
		return false;

	valstack.push_back(atof(chTokens[0].s.c_str()));
	chTokens.pop_front();
	return true;
}

STATIC bool doop(vector<double> &valstack, const string &op)
{
	const double v2 = valstack.back();
	valstack.pop_back();
	const double v1 = valstack.back();
	valstack.pop_back();
	double result;
	if (op == "^")
	{
		result = pow(v1,v2);
	}
	else if (op == "*")
	{
		result = v1*v2;
	}
	else
	{
		fprintf(ioQQQ,"Unknown operator '%s'\n",op.c_str());
		return false;
	}
	valstack.push_back(result);
	return true;
}

STATIC bool ParseExp(deque<Token> &chTokens, vector<double> &valstack)
{
	if (!ParseNumber(chTokens, valstack))
		return false;

	if ( chTokens.size() == 0 )
		return true;

	if ( chTokens.size() < 2 )
		return false;
	if ( Token::symOp != chTokens[0].t || "^" != chTokens[0].s )
		return true;

	string op =  chTokens[0].s;
	chTokens.pop_front();
	
	if (!ParseExp(chTokens, valstack))
		return false;
	
	return doop(valstack, op);
}

STATIC bool ParseProduct(deque<Token> &chTokens, vector<double> &valstack)
{
	if (!ParseExp(chTokens, valstack))
		return false;

	if ( chTokens.size() == 0 )
		return true;

	if ( chTokens.size() < 2 )
		return false;

	if ( Token::symOp != chTokens[0].t || "*" != chTokens[0].s )
		return true;

	string op =  chTokens[0].s;
	chTokens.pop_front();
	
	if (!ParseProduct(chTokens, valstack))
		return false;
	
	return doop(valstack, op);
}

STATIC bool ParseExpr(deque<Token> &chTokens, vector<double> &valstack)
{
	if (ParseProduct(chTokens, valstack))
		return true;
	return false;
}
